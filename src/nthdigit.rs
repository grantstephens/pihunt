//! Low-memory n-th decimal digit extraction for π.
//!
//! Implements Xavier Gourdon, "Computation of the n-th decimal digit of π with low
//! memory" (2003), Theorem 1 / Algorithm 1 / Algorithm 2. See `docs/nthdigit.md` for a
//! summary and `docs/findings-bbp-hunt.md` for how we got here (this is unrelated to the
//! PSLQ-based BBP hunt elsewhere in this crate).
//!
//! # The formula
//!
//! Accelerating `π/4 = Σ (-1)^k/(2k+1)` with the Cohen-Villegas-Zagier process and
//! `P(x) = x^M (1-x)^N` gives (Gourdon eq. 7)
//!
//! ```text
//! S = Σ_{k=0}^{(M+1)N-1} (-1)^k 4/(2k+1)  -  Σ_{k=0}^{N-1} (-1)^k 4 s_k / (2^N (2MN+2k+1))
//! s_k = Σ_{j=0}^{k} C(N,j)
//! |S - π| <= π/(2eM)^N
//! ```
//!
//! With `M`, `N` chosen as in [`Params::new`], `frac(10^n π)` is approximated (Gourdon
//! Proposition 1, error < `10^-n0`) by `frac(B - C)`, where every numerator below is an
//! *integer* (that's the whole trick: `N <= n + 2` makes `4*10^n` and
//! `5^(N-2) 10^(n-N+2)` integers, so the sums only need modular arithmetic on
//! machine-word-sized numbers, not big floats):
//!
//! ```text
//! B = Σ_{k=0}^{(M+1)N-1} (-1)^k (4*10^n mod (2k+1)) / (2k+1)
//! C = Σ_{k=0}^{N-1}      (-1)^k (5^(N-2) 10^(n-N+2) s_k mod m_k) / m_k,   m_k = 2MN+2k+1
//! ```
//!
//! # Memory bound
//!
//! Every term of `B` and `C` is computed independently from `n`, `M`, `N` and `k` alone,
//! using only fixed-size (`u64`/`u128`) scratch state: no array is ever sized by `N` or
//! `n`. Algorithm 2 (below) needs `O(log m)` extra `u64`s for the prime factors of `m`
//! that are `<= k`. So per-thread memory is `O(log^2 n)` (Gourdon's theorem 1), and total
//! memory is that times the (small, fixed) thread count.
//!
//! # Fixed-point accumulation and its error budget
//!
//! Each term `(x mod m)/m` is stored as a `u128` fixed-point fraction `f` with
//! `f / 2^128` approximating the true rational value, computed by [`frac_fixed_point`] via
//! exact 128-by-64 long division (`floor(x * 2^128 / m)`, so the representation error is
//! strictly less than one ulp = `2^-128`). Terms are accumulated with `wrapping_add` /
//! `wrapping_sub`: since we only ever want the *fractional part* of `B - C`, arithmetic
//! modulo 1 is exactly arithmetic modulo `2^128` on this fixed-point encoding, so wrapping
//! add/sub is exact (order-independent) arithmetic on `frac(sum) mod 1` — no error is
//! introduced by accumulation itself, only by the one `floor` per term.
//!
//! Total representation error after summing `(M+1)N + N` terms is therefore
//! `< ((M+1)N + N) * 2^-128`, astronomically below `2^-128` relative to any error bound
//! `10^-n0` we ever ask for in practice (`(M+1)N` stays well under `2^64` for any n this
//! implementation is remotely capable of running in a human lifetime). The error that
//! actually matters is Gourdon's own truncation bound `π/(2eM)^N < 10^-(n+n0)`, which is
//! what `n0` is chosen to satisfy.
//!
//! # Digit-boundary safety
//!
//! [`digits`] never returns a digit it can't certify: it computes `frac_10n_pi` with a
//! guard band of extra digits (`n0 = count + guard`), and re-extracts the requested
//! `count` digits from `x`, `x + err` and `x - err` (`err` = the `10^-n0` error bound
//! converted to fixed-point units). If those three don't agree digit-for-digit, the true
//! value is too close to a rounding boundary (e.g. a run of `9`s or `0`s) for this guard to
//! resolve, so the guard is doubled and the whole computation retried.

use rug::{Float, Integer, float::Constant, ops::Pow};

// ---------------------------------------------------------------------------------
// Parameters (Gourdon Algorithm 1, step 1 / eq. 8)
// ---------------------------------------------------------------------------------

/// `M` and `N` for one Algorithm 1 run: extracting digits certified to `10^-n0` at decimal
/// position `n` (i.e. computing `frac(10^n π)`).
#[derive(Debug, Clone, Copy)]
pub struct Params {
    /// Degree parameter `M` of `P(x) = x^M (1-x)^N`. Always even, `>= 4`.
    pub m: u64,
    /// Degree parameter `N`. Always even. Gourdon requires `N <= n + 2`, checked here.
    pub big_n: u64,
}

impl Params {
    /// Computes `M`, `N` per Gourdon eq. (8). Panics if `N > n + 2` (Gourdon: this holds
    /// once `n >= 4*n0` or so; too-small `n` for the requested `n0` should use the direct
    /// MPFR fallback instead, see [`digits`]).
    pub fn new(n: u64, n0: u32) -> Self {
        assert!(n >= 1, "n must be >= 1");
        let nf = n as f64;
        let ln_n = nf.ln();
        let raw_m = (nf / ln_n.powi(3)).ceil().max(1.0) as u64;
        let mut m = 2 * raw_m;
        if m < 4 {
            m = 4;
        }

        let log_2em = (2.0 * std::f64::consts::E * m as f64).ln();
        let raw_n = (((n + n0 as u64 + 1) as f64) * 10f64.ln() / log_2em).ceil();
        let mut big_n = raw_n as u64;
        if !big_n.is_multiple_of(2) {
            big_n += 1;
        }
        if big_n < 2 {
            big_n = 2;
        }

        assert!(
            big_n <= n + 2,
            "N={big_n} exceeds n+2={} for n={n}, n0={n0}: n is too small for the requested \
             precision (Gourdon requires roughly n >= 4*n0); use the MPFR fallback instead",
            n + 2
        );

        Params { m, big_n }
    }

    /// `log10` of Gourdon's bound `π/(2eM)^N` on `|S - π|` (always negative).
    pub fn error_bound_log10(&self) -> f64 {
        let log_2em = (2.0 * std::f64::consts::E * self.m as f64).ln();
        (std::f64::consts::PI.ln() - self.big_n as f64 * log_2em) / 10f64.ln()
    }
}

// ---------------------------------------------------------------------------------
// Modular arithmetic on u64
// ---------------------------------------------------------------------------------

/// `a * b mod m`, exact, via a `u128` intermediate.
#[inline]
pub fn mulmod(a: u64, b: u64, m: u64) -> u64 {
    ((a as u128 * b as u128) % m as u128) as u64
}

/// `base^exp mod m`, by binary exponentiation.
pub fn powmod(mut base: u64, mut exp: u64, m: u64) -> u64 {
    if m == 1 {
        return 0;
    }
    let mut result = 1u64 % m;
    base %= m;
    while exp > 0 {
        if exp & 1 == 1 {
            result = mulmod(result, base, m);
        }
        base = mulmod(base, base, m);
        exp >>= 1;
    }
    result
}

/// Modular inverse of `a` mod `m` (which must be coprime to `m`), via the extended
/// Euclidean algorithm on `i128` (large enough to hold every intermediate for `u64` inputs
/// without overflow).
fn mod_inverse(a: u64, m: u64) -> u64 {
    fn ext_gcd(a: i128, b: i128) -> (i128, i128, i128) {
        if b == 0 {
            (a, 1, 0)
        } else {
            let (g, x1, y1) = ext_gcd(b, a % b);
            (g, y1, x1 - (a / b) * y1)
        }
    }
    let (g, x, _) = ext_gcd(a as i128, m as i128);
    debug_assert_eq!(g, 1, "{a} is not invertible mod {m}");
    (((x % m as i128) + m as i128) % m as i128) as u64
}

/// `floor((x mod m) * 2^128 / m)` as an exact `u128`: `x/m` encoded as a 128-bit
/// fixed-point fraction in `[0, 1)`. Requires `x < m`.
///
/// Computed by explicit long division of the 192-bit numerator `x * 2^128` (three 64-bit
/// limbs: `x`, `0`, `0`, most significant first) by the 64-bit divisor `m`, one limb at a
/// time, entirely in exact `u128` arithmetic (each per-limb quotient and remainder fits
/// comfortably in `u128` since the running remainder is always `< m <= 2^64`). This is
/// exact: the only "error" is the final `floor`, i.e. strictly less than one ulp = `2^-128`.
pub fn frac_fixed_point(x: u64, m: u64) -> u128 {
    debug_assert!(m > 0);
    debug_assert!(x < m);
    let m128 = m as u128;
    // First limb is x itself (remainder after it is just x, since x < m means the
    // quotient contribution of this limb is 0).
    let r = x as u128;
    // Second limb (value 0).
    let cur = r << 64;
    let q_hi = cur / m128;
    let r = cur % m128;
    // Third limb (value 0).
    let cur = r << 64;
    let q_lo = cur / m128;
    (q_hi << 64) | q_lo
}

// ---------------------------------------------------------------------------------
// Algorithm 2: Σ_{j=0}^{k} C(N,j) mod m, with low memory
// ---------------------------------------------------------------------------------

/// Primes `<= k` dividing `m`, found by trial division of `m` against `2..=k` (stopping
/// early once `m`'s cofactor is fully spent). `O(k)` time, `O(log m)` memory — this is the
/// "acceptable in v1" trial division the spec calls out; a sieve would be faster but would
/// need `O(k)` memory, which we're avoiding on purpose.
fn primes_le_k_dividing(m: u64, k: u64) -> Vec<u64> {
    let mut temp = m;
    let mut primes = Vec::new();
    let mut d = 2u64;
    while d <= k && d <= temp {
        if temp.is_multiple_of(d) {
            primes.push(d);
            while temp.is_multiple_of(d) {
                temp /= d;
            }
        }
        d += 1;
    }
    primes
}

/// Strips every factor of `p` out of `*x`, returning how many there were.
fn strip_power(x: &mut u64, p: u64) -> u32 {
    let mut count = 0u32;
    while (*x).is_multiple_of(p) {
        *x /= p;
        count += 1;
    }
    count
}

/// `Σ_{j=0}^{k} C(N,j) mod m` via Gourdon's Algorithm 2, for `k <= N/2` (see
/// [`sum_binomials_mod`] for the complement trick that handles `k > N/2`).
///
/// Maintains `A = Π a* mod m`, `B = Π b* mod m` (numerator/denominator of the binomial
/// product with every prime factor `p_i <= k` of `m` stripped out), and for each such
/// prime `p_i` a running *exact* integer `R_i = p_i^{e_i}` where `e_i` is the cumulative
/// sum of `v_{p_i}(a) - v_{p_i}(b)` so far. Gourdon proves `e_i` is always `>= 0` and
/// `R_i <= N`, so `R_i` fits in a `u64` and never needs to be reduced mod `m` on its own
/// (only the product `Π R_i` is, when it's folded into `C`). `C` accumulates
/// `C := C * b* + A * (Π R_i) mod m`, so that after step `j`, `C/B mod m = Σ_{i=0}^{j}
/// C(N,i) mod m`. At the end, `s_k mod m = C * B^-1 mod m` (`B` is coprime to `m` because
/// every prime factor of any `b* <= k` that also divides `m` has already been stripped).
fn s_k_mod_direct(big_n: u64, k: u64, m: u64) -> u64 {
    if m == 1 {
        return 0;
    }
    let primes = primes_le_k_dividing(m, k);
    let mut r_vals = vec![1u64; primes.len()];
    let (mut a_acc, mut b_acc, mut c_acc) = (1u64 % m, 1u64 % m, 1u64 % m);

    for j in 1..=k {
        let a = big_n - j + 1;
        let mut a_star = a;
        let mut b_star = j;
        for (idx, &p) in primes.iter().enumerate() {
            let alpha = strip_power(&mut a_star, p);
            let beta = strip_power(&mut b_star, p);
            for _ in 0..alpha {
                r_vals[idx] *= p;
            }
            for _ in 0..beta {
                r_vals[idx] /= p;
            }
        }
        debug_assert!(
            primes.iter().enumerate().all(|(i, _)| r_vals[i] <= big_n),
            "R_i must stay <= N (Gourdon lemma 1)"
        );

        let a_star_m = a_star % m;
        let b_star_m = b_star % m;
        a_acc = mulmod(a_acc, a_star_m, m);
        b_acc = mulmod(b_acc, b_star_m, m);

        let mut r_prod = 1u64 % m;
        for &r in &r_vals {
            r_prod = mulmod(r_prod, r % m, m);
        }
        c_acc = (mulmod(c_acc, b_star_m, m) + mulmod(a_acc, r_prod, m)) % m;
    }

    let b_inv = mod_inverse(b_acc, m);
    mulmod(c_acc, b_inv, m)
}

/// `Σ_{j=0}^{k} C(N,j) mod m` (Gourdon Algorithm 2), for `0 <= k < N`.
///
/// Uses the complement identity `Σ_{j<=k} C(N,j) = 2^N - Σ_{j<=N-k-1} C(N,j)` when
/// `k > N/2`, halving the work in that regime.
pub fn sum_binomials_mod(big_n: u64, k: u64, m: u64) -> u64 {
    if m == 1 {
        return 0;
    }
    if 2 * k > big_n {
        let comp_k = big_n - k - 1;
        let s_comp = s_k_mod_direct(big_n, comp_k, m);
        let two_pow_n = powmod(2, big_n, m);
        (two_pow_n + m - s_comp % m) % m
    } else {
        s_k_mod_direct(big_n, k, m)
    }
}

// ---------------------------------------------------------------------------------
// The two sums, and frac(10^n π)
// ---------------------------------------------------------------------------------

#[inline]
fn signed(term: u128, k: u64) -> u128 {
    if k.is_multiple_of(2) {
        term
    } else {
        term.wrapping_neg()
    }
}

fn b_sum(n: u64, terms: u64) -> u128 {
    use rayon::prelude::*;
    (0..terms)
        .into_par_iter()
        .map(|k| {
            let m = 2 * k + 1;
            let pow10 = powmod(10, n, m);
            let x = mulmod(4 % m, pow10, m);
            signed(frac_fixed_point(x, m), k)
        })
        .reduce(|| 0u128, u128::wrapping_add)
}

fn c_sum(n: u64, p: Params) -> u128 {
    use rayon::prelude::*;
    let (big_n, big_m) = (p.big_n, p.m);
    (0..big_n)
        .into_par_iter()
        .map(|k| {
            let m = 2 * big_m * big_n + 2 * k + 1;
            let s = sum_binomials_mod(big_n, k, m);
            let e1 = big_n - 2;
            let e2 = n - big_n + 2;
            let pow5 = powmod(5, e1, m);
            let pow10 = powmod(10, e2, m);
            let y = mulmod(mulmod(pow5, pow10, m), s, m);
            signed(frac_fixed_point(y, m), k)
        })
        .reduce(|| 0u128, u128::wrapping_add)
}

/// `frac(10^n π)` to within `10^-n0`, as a 128-bit fixed-point fraction (value = `x /
/// 2^128`). Requires `n` large enough that [`Params::new`] doesn't panic (roughly `n >=
/// 4*n0`); see [`digits`] for a version with an MPFR fallback for small `n`.
pub fn frac_10n_pi(n: u64, n0: u32) -> u128 {
    let p = Params::new(n, n0);
    let terms = (p.m + 1) * p.big_n;
    let (b, c) = rayon::join(|| b_sum(n, terms), || c_sum(n, p));
    b.wrapping_sub(c)
}

// ---------------------------------------------------------------------------------
// Decimal digit extraction
// ---------------------------------------------------------------------------------

/// Multiplies the fixed-point fraction `*x` (value `= *x / 2^128`) by 10 in place,
/// returning the digit that carries out (`0..=9`) — i.e. one step of the standard
/// "multiply by the base, take the integer part" digit-extraction loop, done exactly in
/// `u128` by splitting into 64-bit halves so the intermediate `x * 10` (up to 132 bits)
/// never overflows.
fn next_digit(x: &mut u128) -> u8 {
    let x_hi = (*x >> 64) as u64;
    let x_lo = *x as u64;
    let lo_prod = (x_lo as u128) * 10;
    let carry = lo_prod >> 64;
    let new_lo = lo_prod as u64;
    let hi_prod = (x_hi as u128) * 10 + carry;
    let digit = (hi_prod >> 64) as u8;
    let new_hi = hi_prod as u64;
    *x = ((new_hi as u128) << 64) | (new_lo as u128);
    digit
}

fn extract_digits(mut x: u128, count: usize) -> String {
    let mut s = String::with_capacity(count);
    for _ in 0..count {
        s.push((b'0' + next_digit(&mut x)) as char);
    }
    s
}

/// A conservative bound, in fixed-point units (`err / 2^128`), on the total error of
/// `frac_10n_pi`: the series truncation (`< 10^-n0`, Gourdon's Proposition 1, rounded up)
/// plus one ulp of fixed-point rounding for each of the `terms` accumulated fractions
/// (each is floored once, then added or subtracted exactly).
fn error_units(n0: u32, terms: u64) -> u128 {
    let num = Integer::from(1) << 128u32;
    let den = Integer::from(10).pow(n0);
    let mut q = Integer::from(&num / &den);
    let rem = Integer::from(&num - &q * &den);
    if rem > 0 {
        q += 1;
    }
    q += terms;
    q.to_u128().unwrap_or(u128::MAX)
}

/// Largest `n0` worth attempting: rounding costs up to one ulp (2^-128) per term, and there
/// are ~2^45 terms by n ~ 10^9 (~10^-25), so a truncation bound much below 10^-24 can't be
/// certified anyway. `error_units` still accounts for the rounding exactly; this only stops
/// the guard-doubling loop from chasing precision the accumulator doesn't have.
const MAX_N0: u32 = 24;

/// Most digits one evaluation hands out; longer requests are split (see [`digits`]).
const CHUNK: usize = 16;

/// Below this position, just ask MPFR for π directly — it's cheap there, and it sidesteps
/// [`Params::new`]'s `N <= n + 2` precondition for small `n`.
const SMALL_N_THRESHOLD: u64 = 2000;

/// `count` decimal digits of π at positions `n+1 ..= n+count` after the decimal point
/// (position 1 is the '1' in 3.14159...).
///
/// For small `n` (below [`SMALL_N_THRESHOLD`], or too small relative to `count` for
/// Gourdon's method to apply), falls back to computing π directly with MPFR. Otherwise
/// uses [`frac_10n_pi`] with a guard band of extra digits, doubling the guard and retrying
/// whenever the requested digits are too close to a rounding boundary to certify (see the
/// module docs).
pub fn digits(n: u64, count: usize) -> String {
    if count == 0 {
        return String::new();
    }
    if count > CHUNK {
        // One fixed-point evaluation can only certify ~MAX_N0 digits, so long requests are
        // stitched together from independent chunks.
        return (0..count)
            .step_by(CHUNK)
            .map(|i| digits(n + i as u64, CHUNK.min(count - i)))
            .collect();
    }
    let mut guard: u64 = 4;
    loop {
        let n0 = count as u64 + guard;
        if n < SMALL_N_THRESHOLD || n < 4 * n0 || n0 > MAX_N0 as u64 {
            // Past MAX_N0 the u128 accumulator can't certify the digits (a pathological run of
            // ~20 equal digits right after the block). MPFR is exact but needs O(n) memory.
            return digits_via_mpfr(n, count);
        }
        let n0 = n0 as u32;
        let x = frac_10n_pi(n, n0);
        let p = Params::new(n, n0);
        let err = error_units(n0, (p.m + 1) * p.big_n + p.big_n);
        let base = extract_digits(x, count);
        let plus = extract_digits(x.wrapping_add(err), count);
        let minus = extract_digits(x.wrapping_sub(err), count);
        if base == plus && base == minus {
            return base;
        }
        guard *= 2;
    }
}

/// Computes `count` digits of π at position `n+1..=n+count` directly with MPFR (cheap for
/// small `n`, used as the fallback in [`digits`]).
fn digits_via_mpfr(n: u64, count: usize) -> String {
    let guard: u64 = 20;
    let total = n + count as u64 + guard;
    let bits = crate::pslq::digits_to_bits(total as u32) + 8;
    let pi = Float::with_val(bits, Constant::Pi);
    let scale = Float::with_val(bits, Integer::from(10).pow(total as u32));
    let scaled = Float::with_val(bits, &pi * &scale);
    let int_part = scaled.to_integer().expect("pi * 10^total is finite");
    let s = int_part.to_string(); // "3" followed by `total` decimal digits
    let start = (n + 1) as usize;
    s[start..start + count].to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn error_bound_includes_rounding_per_term() {
        let base = error_units(20, 0);
        assert_eq!(error_units(20, 1_000), base + 1_000);
        // At n0 = 30 the truncation bound is below what 10^10 rounded terms can promise.
        assert!(error_units(30, 10_000_000_000) > 2 * error_units(30, 0));
    }
}
