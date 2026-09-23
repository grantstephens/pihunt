//! Exact π digit extraction via Machin's formula (`π = 16·arctan(1/5) − 4·arctan(1/239)`),
//! evaluated as an integer fixed-point sum on [`Big`]. This is the exact fallback used by
//! [`crate::nthdigit`] and [`crate::nthdigit2`] wherever their fixed-point certification loop
//! can't be trusted directly — see `nthdigit::digits_fallback` for when and why.
//!
//! # Method
//!
//! `arctan(1/x) = Σ_{k=0}^∞ (-1)^k / ((2k+1)·x^(2k+1))`. Scaling by `S = 10^P` and truncating
//! every division to an integer (floor — every term is positive before its sign is applied)
//! gives an integer approximation to `arctan(1/x)·S` with a per-term error bounded by a small
//! constant (see [`error_bound`]), so the total error after `terms` loop iterations is
//! `O(terms)` — utterly negligible next to `S` for any `P` this is used at.
//!
//! `16·atan_scaled(5) − 4·atan_scaled(239)` (exact `Big` arithmetic — `mul_u64`/`sub` don't
//! round, so this combination step adds no further error) approximates `floor(π·S)` to within
//! that same `O(terms)` bound. [`digits`] computes this at `P = n + count + guard` decimal
//! digits, then checks that the requested digit window doesn't change when the error bound is
//! added or subtracted (the same guard-band-and-retry certification `nthdigit::digits` uses)
//! before trusting it, widening the guard and retrying otherwise.

use crate::bignum::Big;

/// `Σ_{k=0}^{K} (-1)^k · S / ((2k+1)·x^(2k+1))`, `K` chosen so the running term underflows to
/// 0 at this scale (`S = scale`). Kept as a running term `t = S/x^(2k+1)`, divided by `x²`
/// each step, per Gourdon-style fixed-point accumulation — each loop iteration costs two
/// truncating [`Big::div_u64`] calls on a shrinking `Big`.
///
/// Returns the sum — always non-negative, since partial sums of this alternating,
/// term-magnitude-decreasing series stay within `[0, first term]` — and the number of terms
/// added, which [`error_bound`] uses to bound the total truncation error.
fn atan_scaled(x: u64, scale: &Big) -> (Big, u64) {
    let x2 = x * x;
    let mut term = scale.div_u64(x); // k=0: S / x^1
    let mut sum = Big::zero();
    let mut add = true;
    let mut k: u64 = 0;
    let mut terms: u64 = 0;
    while !term.is_zero() {
        let contribution = term.div_u64(2 * k + 1);
        sum = if add {
            sum.add(&contribution)
        } else {
            sum.sub(&contribution)
        };
        terms += 1;
        add = !add;
        k += 1;
        term = term.div_u64(x2);
    }
    (sum, terms)
}

/// Bound on the total truncation error of an [`atan_scaled`] sum with `terms` terms, in
/// units of the scale `S` it was computed at.
///
/// Each loop iteration performs two truncating divisions. Let `e_k` be the absolute error
/// (in units of `S`) between the true `S/x^(2k+1)` and the computed running term at step
/// `k`: `e_0 < 1` (one division from `S`), and `e_{k+1} <= e_k/x² + 1` (dividing an
/// under-estimate by `x²` shrinks the inherited error by `x²`, and that one new division
/// adds less than 1 more unit). Since Machin only ever calls this with `x ∈ {5, 239}`, so
/// `x² >= 25`, this recursion is bounded by `e_k < x²/(x²-1) <= 25/24` for every `k`. The
/// contribution `term/(2k+1)` then adds less than `e_k + 1 < 2.05` more absolute error on
/// top. Summed over `terms` terms, worst case additive (no cancellation assumed even though
/// the sum itself alternates sign): `< 2.05 * terms`. Rounded up generously — this costs
/// nothing, `S` has orders of magnitude more digits than this bound ever will — to
/// `4 * terms + 16`.
fn error_bound(terms: u64) -> u64 {
    4 * terms + 16
}

/// `x`'s decimal string, left-padded with zeros to `width` digits, then sliced to
/// `[start, end)`. `x` is `floor(π·S)` or a nearby perturbation by the error bound: always
/// positive, and in practice always exactly `width` digits (`π·S ∈ [3·S, 4·S)`, and the
/// error bound is astronomically smaller than `S`); the padding only guards the
/// vanishingly-unlikely edge case where it isn't.
fn extract(x: &Big, width: usize, start: usize, end: usize) -> String {
    let s = x.to_decimal_string();
    let s = if s.len() < width {
        format!("{}{s}", "0".repeat(width - s.len()))
    } else {
        s
    };
    s[start..end].to_string()
}

/// `count` decimal digits of π at positions `n+1..=n+count` (position 1 is the '1' in
/// 3.14159...), computed exactly via Machin's formula: guard digits plus a boundary check
/// against the error bound, retrying with more guard digits whenever the requested window is
/// too close to a rounding boundary (e.g. a run of repeated digits) to certify.
pub fn digits(n: u64, count: usize) -> String {
    if count == 0 {
        return String::new();
    }
    let mut guard: u64 = 25;
    loop {
        let total = n + count as u64 + guard;
        let scale = Big::pow_u64(10, total as u32);
        let (atan5, terms5) = atan_scaled(5, &scale);
        let (atan239, terms239) = atan_scaled(239, &scale);
        let pi_scaled = atan5.mul_u64(16).sub(&atan239.mul_u64(4));
        let err = Big::from_u64(16 * error_bound(terms5) + 4 * error_bound(terms239));

        let start = n as usize + 1;
        let end = start + count;
        let width = total as usize + 1;
        let base = extract(&pi_scaled, width, start, end);
        let plus = extract(&pi_scaled.add(&err), width, start, end);
        let minus = extract(&pi_scaled.sub(&err), width, start, end);
        if base == plus && base == minus {
            return base;
        }
        guard += 20;
    }
}
