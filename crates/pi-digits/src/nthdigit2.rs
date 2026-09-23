//! Xavier Gourdon's unpublished "Theorem 2": decimal digits of π with `O(m)`-bit memory
//! (`m` a tunable budget), trading memory for time against [`crate::nthdigit`]'s Theorem 1
//! (`O(log² n)` memory, `O(n²)` time). See `docs/nthdigit-theorem2.md` for the full
//! reconstruction, proofs and prototype measurements this module ports to Rust (the
//! reference implementation is `research/thm2/thm2.py`, validated against MPFR to
//! `n = 1.28e6`). Section numbers cited below (`§4.2` etc.) refer to that doc.
//!
//! # The formula (unchanged from Theorem 1, doc §2)
//!
//! `frac(10^n π) ≈ frac(B - C)` (Gourdon Proposition 1, error `< 10^-n0`):
//!
//! ```text
//! B = Σ_{k<(M+1)N} (-1)^k (4·10^n mod (2k+1)) / (2k+1)
//! C = Σ_{k<N}      (-1)^k (X·s_k mod m_k) / m_k,   X = 5^(N-2)·10^(n-N+2), m_k = 2MN+2k+1
//! s_k = Σ_{j<=k} C(N,j)
//! ```
//!
//! `B` is computed exactly as in Theorem 1 ([`crate::nthdigit::b_sum`], reused verbatim: the
//! formula, the fixed-point encoding and the rayon parallelisation are all independent of how
//! `C` gets computed). Theorem 2 replaces Theorem 1's `O(k)`-per-term Algorithm 2 for `s_k mod
//! m_k` with a **chunked accumulating remainder tree** (doc §4.4): binary splitting of the
//! recurrence `P_t = N!/(N-t)!, D_t = t!, T_t = t!·s_t` on numbers of `~mem_bits` bits, so the
//! per-target cost drops from `Θ(t)` to `Θ(polylog(m)/m)` amortised. `M` is re-balanced against
//! `mem_bits` (doc §5.4) instead of Theorem 1's `M ≈ n/log³n`, since a smaller `M` (hence
//! larger `N`) is cheaper once `C` itself got cheaper.
//!
//! # Why every "big number" in this module fits in a `u64` except the ART state
//!
//! `m_k = 2MN + 2k + 1` is bounded by `O(n^1.5/log n)` for `mem_bits ∝ √n` (the headline case;
//! see doc §5.4), which stays far below `2^64` for every `n` this implementation can run in
//! practice (as with Theorem 1, this is an implicit domain limit shared with
//! [`crate::nthdigit`]). So every individual modulus (`m_k`, or one of its prime-power parts
//! `q`) fits in a `u64`, and all of Theorem 1's `u64` modular-arithmetic helpers
//! ([`crate::nthdigit::mulmod`], `powmod`, `mod_inverse`, `frac_fixed_point`) are reused
//! unchanged. The *only* place that needs `rug::Integer` bignums is the ART: the recurrence
//! state `(P, T, D)` and the chunk modulus `Q = Π q` genuinely need `~mem_bits` bits, because
//! that's the whole point (binary-splitting many small factors into one big exact product
//! before reducing, doc §3 "obstacle 1").
//!
//! # Precision and certification
//!
//! Exactly as Theorem 1 (module docs on [`crate::nthdigit`]): each of the `(M+1)N` B-terms and
//! every C-part partial-fraction contribution is stored as a `u128` fixed-point fraction via
//! [`crate::nthdigit::frac_fixed_point`] (exact `floor(x·2^128/m)`, `< 1` ulp error), and
//! accumulated with `wrapping_add`/`wrapping_sub` (exact arithmetic mod 1 on this encoding).
//! Unlike Theorem 1, a single `k` in the C-part can contribute *more than one* rounded term
//! (one per coprime prime-power part of `m_k`: the main ART part, plus one per small prime via
//! Lucas or the p-adic recursion) — see [`c_part`], which counts every actual call to
//! [`add_contribution`] and returns that exact count so [`error_units`][crate::nthdigit] gets
//! charged the true number of roundings, not an estimate.
//!
//! # Memory: what's `O(mem_bits)` and what isn't (be honest, doc §7)
//!
//! The C-part item generation is **streamed** (doc §4.4's last paragraph, implemented here):
//! there is no global `O(N)` factorisation table and no global sorted item list. Instead:
//!
//! * [`stream_needs_and_chunks`] makes one sequential pass over target `t = 0..N/2` (both
//!   `k = t` and `k = N-1-t` per target, doc §4.4: "Main items come out naturally in target
//!   order"), factoring `O(sqrt(N))`-sized windows at a time via [`factor_window`] (freed once
//!   processed) and deciding Main-item **chunk boundaries by running bit-count alone** — the
//!   items themselves aren't kept, only their bit lengths, until a chunk's worth (`~mem_bits`)
//!   has been seen. [`art_chunk_by_range`] then *regenerates* each chunk's items from its target
//!   window (re-factoring that one small window) instead of slicing a pre-built array, so at
//!   most one chunk's items (`O(mem_bits / log n)`) are ever live at once. The ART's own bignum
//!   working set — product tree, `(P,T,D)` state, one group's binary-splitting stack, per chunk
//!   — stays `O(mem_bits log(mem_bits))` bits (doc §5.1), run in parallel across chunks.
//! * Small-prime Lucas/p-adic needs (`p <= sqrt(max m_k)`, doc §4.3) are collected into
//!   `HashMap<u64, SideRanges>`s (see [`SideRanges`]'s docs) bounded by `O(pi(sqrt(max m_k)))`
//!   entries and resolved by their own small dedicated ART sub-pass ([`resolve_lucas_items`]) —
//!   doc §4.4's "process small primes in their own pass with their own small ART". **This used
//!   to be the dominant term** (see "What was actually dominant" below) until the `SideRanges`
//!   fix: earlier, each key held a `Vec<(k, t[, e])>` of every matching `k`, which sums to `O(N
//!   loglog sqrt(max m_k))` words over all keys (Mertens' third theorem: `Σ_{p<=P} 1/p ~ loglog
//!   P`), not the `O(pi(sqrt(max m_k)))` this section always *claimed*. `SideRanges` fixes that:
//!   for a fixed prime and side (low `k=t` / high `k=N-1-t`), the matching `k` form one
//!   arithmetic progression, so only its first/last `t` needs storing (two `u64`s), and
//!   consumers ([`lucas_consumers`], [`padic_consumers`]) regenerate `(k, t)` pairs on the fly by
//!   walking that range. p-adic tables (`O(p)` words per prime, doc §5.2) are still built and
//!   dropped one prime at a time, so their peak is `O(max p)` over primes actually used, not
//!   `O(Σp)`; here `p <= sqrt(max m_k)`, so that peak is `O(sqrt(max m_k))`, not `O(N)`.
//! * The small-primes sieve itself ([`primes_upto`], up to `sqrt(max m_k)`) is
//!   `O(pi(sqrt(max m_k)))` words — explicitly within budget (it doesn't grow with `N` the way
//!   the old per-`k` factor table did; for `mem_bits ∝ sqrt(n)` it's `O(n^{0.75}/log n)`-ish, far
//!   below the other terms in practice — see the measured table in `docs/nthdigit.md`).
//!
//! ## What was actually dominant (found by profiling, not by re-reading the complexity argument)
//!
//! Measured peak RSS at `n = 1e7` was 812 MiB before streaming, 447 MiB after streaming (the
//! caveat this section used to end on: the cofactor `Vec` below). Streaming's own complexity
//! argument said the *only* remaining `O(N)` piece was `cofactor` — but `cofactor` measures only
//! tens of MiB at these `n` (see below), nowhere near 447 MiB. Profiling with `--features
//! mem-profile` (a counting global allocator plus an explicit per-structure breakdown, both in
//! [`crate::mem_profile`] and `log_needs_sizes`, feature-gated so not always in scope) at
//! `n = 1e6`/`3e6` found two real culprits,
//! **neither of them `cofactor`**:
//!
//! 1. **The actual dominant term**, by a wide margin: [`PadicBinom`]'s `memo_s`/`memo_c` were
//!    memoising across an entire prime's worth of queries instead of within one query's
//!    recursion tree (a bug, not a documented tradeoff — nothing in doc §5.2's "O(e² p log_p N)
//!    per query" cost model called for this). For small `p` (where a prime's `klist` is large —
//!    e.g. `p=3` at `n=1e6` had 34 557 queries), this grew the memo maps to roughly *10x the
//!    query count* instead of the intended `O(e² log_p N)` per query, measured at ~40 MiB of the
//!    ~58 MiB peak at `n=1e6` from `p=3` alone. Fixed by [`PadicBinom::clear_query_memo`],
//!    called before every top-level query in [`padic_consumers`]: cheap (`HashMap::clear` keeps
//!    the allocation) and correctness-neutral (`s`/`c` are pure functions of their arguments, so
//!    discarding memo entries only forces recomputation, never changes the answer).
//! 2. `lucas_small`/`padic_small` themselves, per the `SideRanges` fix described above: measured
//!    12.9 MiB / 2.4 MiB at `n=1e6`, growing to 36.7 MiB / 6.1 MiB at `n=3e6` before the fix, now
//!    a fraction of a MiB at both.
//!
//! **Result:** peak RSS 67 MiB -> 14 MiB at `n=1e6`; 62 MiB -> 23 MiB at `n=3e6`; 447 MiB (the
//! streaming-only figure) -> 69 MiB at `n=1e7` (see `docs/nthdigit.md`'s updated table for the
//! full before/after/after-after numbers). The remaining ~69 MiB at `n=1e7` is now genuinely
//! dominated by the one piece below, `cofactor`, plus a roughly-comparable amount of transient
//! `rug::Integer`/rayon/glibc-arena overhead this section doesn't itemise further.
//!
//! **What's still not `O(mem_bits)` (doc §7's one remaining honest caveat):** a prime `p` can
//! divide `m_k` *without* being `<= sqrt(max m_k)` — it's `m_k`'s single larger leftover
//! cofactor (there's at most one per `k`, `m_k`'s factorisation leaves at most one prime factor
//! above its square root). When that cofactor is itself `<= t_k`, it still needs Lucas
//! treatment, but unlike a small sieve prime it's essentially unique to one or two `k` (not
//! shared by a residue-class arithmetic progression the way a `p <= sqrt(max m_k)` prime is, so
//! most of its `SideRanges` entries — yes, it's routed through the same machinery as
//! `lucas_small` now, see [`c_part`] — have exactly one member and don't compact). It can't be
//! resolved by the same bounded per-prime iteration. [`stream_needs_and_chunks`] collects the
//! raw needs into `cofactor`, a `Vec` whose size is **not** bounded independent of `N` — measured
//! at 45 K entries at `n=1e6` (~1.5 MiB) and 133 K entries at `n=3e6` (~6 MiB), 24 bytes each.
//! It's still `O(N)` words. Eliminating it for real would need either an external (disk-backed)
//! sort of the cofactor needs by target, or a smarter per-prime classification that doesn't
//! require discovering the cofactor's value before knowing whether it's "small" — both are
//! future work, not implemented here. In practice, at the `n` this implementation is run at, it
//! no longer dominates: see the measured table in `docs/nthdigit.md`.
//!
//! **Overall bound achieved:** peak memory is `O(mem_bits log(mem_bits) * threads + pi(sqrt(max
//! m_k)) + chunks + cofactor_count)` where `chunks = O(N log n / mem_bits)` (tiny: a few thousand
//! `(u64,u64)` pairs even at `n = 1e7`) and `cofactor_count` is the one term above that scales
//! with `N` — now, empirically, the dominant *named* term, but itself only tens of MiB through
//! `n = 1e7` (measured; see `docs/nthdigit.md`).

use crate::bignum::Big;
use crate::nthdigit::{
    self, MAX_N0, extract_digits, frac_fixed_point, mod_inverse, mulmod, powmod, signed,
};
use crate::par::{maybe_join, maybe_par_iter, maybe_reduce};
#[cfg(feature = "parallel")]
use rayon::iter::ParallelIterator as _;
use std::collections::HashMap;

/// Cache of Lucas' theorem's per-digit binomial row + prefix sums, keyed `(prime, digit
/// position)` (valid as long as `N` doesn't change between calls — see [`lucas_s`]).
type LucasRowCache = HashMap<(u64, u32), (Vec<u64>, Vec<u64>)>;
/// `(s_r, C(N,r)) mod p`, keyed `(p, r)`: what a Lucas ART item resolves to.
type LucasVal = (u64, u64);
/// One resolved Lucas leaf, ready to feed into [`lucas_consumers`].
type LucasLeaf = ((u64, u64), LucasVal);

// ---------------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------------

/// `M`, `N` for one Theorem-2 run at memory budget `mem_bits` (doc §4, §5.4).
///
/// Unlike Theorem 1's `M ≈ n/log³n` (chosen to make the `O(k)`-per-term `B`/`C` balance when
/// `C` costs `Θ(k)`), Theorem 2 makes `C` cost roughly `Θ(N log n log²(mem_bits)/mem_bits)`
/// (doc §5.1), so a *smaller* `M` (hence larger `N`, since `N ∝ 1/ln(2eM)`) is cheaper once
/// `mem_bits` is large enough. We use the prototype's balancing heuristic
/// (`research/thm2/thm2.py::run`), `M = max(4, 2·round(2n/mem_bits))`: `C`'s cost scales like
/// `N²/mem_bits`, `B`'s like `M·N`, and `N ∝ n/M` to first order, so equalising them gives
/// `M ∝ n/mem_bits` up to the slowly-varying log factors doc §5.4 spells out exactly. This
/// coarser heuristic is what the prototype measured against MPFR up to `n = 1.28e6`; a tuned
/// constant could shave a further constant factor but isn't needed for correctness.
#[derive(Debug, Clone, Copy)]
pub struct Params2 {
    /// Degree parameter `M`. Always even, `>= 4`.
    pub big_m: u64,
    /// Degree parameter `N`. Always even, `<= n + 2`.
    pub big_n: u64,
    /// The memory budget this run was chosen for (bits of one ART chunk modulus).
    pub mem_bits: u64,
}

impl Params2 {
    /// Chooses `M` from `mem_bits`, then `N` exactly as [`crate::nthdigit::Params::new`] does
    /// from `n`, `n0`, `M` (doc §2's `N = ⌈(n+n0+1)ln10/ln(2eM)⌉`, rounded up to even). Panics
    /// under the same conditions as Theorem 1's `Params::new` (`N` would exceed `n+2`).
    pub fn new(n: u64, n0: u32, mem_bits: u64) -> Self {
        assert!(n >= 1, "n must be >= 1");
        assert!(mem_bits >= 64, "mem_bits must be >= 64, got {mem_bits}");
        let raw = ((2.0 * n as f64 / mem_bits as f64).round() as u64).max(1);
        let mut big_m = 2 * raw;
        if big_m < 4 {
            big_m = 4;
        }
        let log_2em = (2.0 * std::f64::consts::E * big_m as f64).ln();
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
            "N={big_n} exceeds n+2={} for n={n}, n0={n0}, mem_bits={mem_bits}: n is too small \
             for the requested precision; use the MPFR fallback instead",
            n + 2
        );
        Params2 {
            big_m,
            big_n,
            mem_bits,
        }
    }

    /// `log10` of Gourdon's truncation bound `π/(2eM)^N` (doc §2), identical formula to
    /// Theorem 1's [`crate::nthdigit::Params::error_bound_log10`].
    pub fn error_bound_log10(&self) -> f64 {
        let log_2em = (2.0 * std::f64::consts::E * self.big_m as f64).ln();
        (std::f64::consts::PI.ln() - self.big_n as f64 * log_2em) / 10f64.ln()
    }
}

// ---------------------------------------------------------------------------------
// Sieve of Eratosthenes + segmented factoring of the m_k interval (doc §4.4, §5.3)
// ---------------------------------------------------------------------------------

/// All primes `<= x`, plain sieve of Eratosthenes. `x` is `O(sqrt(mem_bits-scaled N))` in
/// every call site here, so this is cheap and its `O(x)`-bit memory is not the bottleneck.
fn primes_upto(x: u64) -> Vec<u64> {
    if x < 2 {
        return Vec::new();
    }
    let x = x as usize;
    let mut sieve = vec![true; x + 1];
    sieve[0] = false;
    sieve[1] = false;
    let mut i = 2usize;
    while i * i <= x {
        if sieve[i] {
            let mut j = i * i;
            while j <= x {
                sieve[j] = false;
                j += i;
            }
        }
        i += 1;
    }
    (0..=x as u64).filter(|&i| sieve[i as usize]).collect()
}

/// Complete factorisations of `m_k = 2MN + 2k + 1` for every `k < N`, via one segmented sieve
/// pass over primes `<= sqrt(max m_k)` (doc §5.3: `O(N loglog N)` word ops total). Holds
/// `O(N)` words: this is the **reference** (non-streaming) construction, kept only for tests
/// that check the streaming item generation ([`stream_needs_and_chunks`] et al.) produces the
/// exact same multiset of items as this original approach (see the module docs' history: this
/// used to be the production path).
#[cfg(test)]
fn factor_interval(big_m: u64, big_n: u64) -> Vec<Vec<(u64, u32)>> {
    let base = 2 * big_m * big_n + 1;
    let top = base + 2 * (big_n - 1);
    let mut rem: Vec<u64> = (0..big_n).map(|k| base + 2 * k).collect();
    let mut fac: Vec<Vec<(u64, u32)>> = vec![Vec::new(); big_n as usize];
    let limit = (top as f64).sqrt() as u64 + 2;
    for p in primes_upto(limit) {
        if p == 2 {
            continue; // every m_k is odd (base odd + even 2k)
        }
        let inv2 = mod_inverse(2 % p, p);
        let base_mod_p = base % p;
        let neg_base_mod_p = (p - base_mod_p) % p;
        let mut k = mulmod(neg_base_mod_p, inv2, p);
        while k < big_n {
            let mut e = 0u32;
            while rem[k as usize].is_multiple_of(p) {
                rem[k as usize] /= p;
                e += 1;
            }
            if e > 0 {
                fac[k as usize].push((p, e));
            }
            k += p;
        }
    }
    for k in 0..big_n as usize {
        if rem[k] > 1 {
            fac[k].push((rem[k], 1));
        }
    }
    fac
}

// ---------------------------------------------------------------------------------
// Streaming factorisation (doc §4.4 last paragraph, doc §7's fix): factor a *window* of `k`
// values at a time (`O(len)` memory, freed when the window is dropped) instead of the whole
// `[0, N)` range up front. Same segmented-sieve technique as [`crate::nthdigit::factor_segment`]
// (residue-class jump per small prime), but keeps exponents (needed to tell Lucas `e==1` from
// p-adic `e>=2`) and also returns each `k`'s leftover cofactor.
// ---------------------------------------------------------------------------------

/// For `k` in `[k0, k0+len)`, `m_k = base + 2k`'s prime-power factors `p^e` with `p` in
/// `small_primes`, plus each `k`'s leftover cofactor (`1`, or a single prime `> max(small_primes)`
/// — see [`factor_interval`]'s docs for why there's at most one). `O(len)` memory for the
/// duration of the call; `small_primes` is shared read-only across every window (built once per
/// [`c_part`] run, `O(pi(sqrt(max m_k)))` words — explicitly within budget, doc §4.4/§7).
fn factor_window(
    base: u64,
    k0: u64,
    len: u64,
    small_primes: &[u64],
) -> (Vec<u64>, Vec<Vec<(u64, u32)>>) {
    let len_usize = len as usize;
    let mut rem: Vec<u64> = (0..len).map(|i| base + 2 * (k0 + i)).collect();
    let mut fac: Vec<Vec<(u64, u32)>> = vec![Vec::new(); len_usize];
    for &p in small_primes {
        if p == 2 {
            continue; // m_k is always odd
        }
        let inv2 = p.div_ceil(2); // inverse of 2 mod odd p, same trick as factor_segment
        let target_k_mod_p = ((p - base % p) % p) * inv2 % p;
        let mut i = ((target_k_mod_p + p - k0 % p) % p) as usize;
        while i < len_usize {
            let mut e = 0u32;
            while rem[i].is_multiple_of(p) {
                rem[i] /= p;
                e += 1;
            }
            if e > 0 {
                fac[i].push((p, e));
            }
            i += p as usize;
        }
    }
    (rem, fac)
}

/// The `Main` item's modulus for target `t`, given `k`'s small-prime factorisation and leftover
/// cofactor (doc §4.3: the "large primes" part, `p > t`). `None` if every prime factor is `<= t`
/// (nothing for the ART to do at this `k` — it's entirely Lucas/p-adic).
fn main_modulus(t: u64, rem: u64, fac: &[(u64, u32)]) -> Option<u64> {
    let mut good = 1u64;
    for &(p, e) in fac {
        if p > t {
            good *= p.pow(e);
        }
    }
    if rem > 1 && rem > t {
        good *= rem;
    }
    (good > 1).then_some(good)
}

/// Compact stand-in for "every `k` in `[0, N)` with a given small-prime property, split by
/// which side of the mirror it's on" (doc §4.4/§7's fix for the memory this used to cost).
///
/// For a *fixed* prime `p` and a fixed side, the set of `k` sharing a small-prime property
/// (`p | m_k`, or `p^2 | m_k`) is a single arithmetic progression in `k` (`m_k = base + 2k` is
/// linear, so `m_k ≡ 0 (mod q)` picks out one residue class of `k` mod `q`), and — because
/// [`stream_needs_and_chunks`] scans `t` in increasing order and, within one side, `t` is a
/// strictly monotonic linear function of `k` (`t = k` on the low side, `t = N-1-k` on the high
/// side) — the matching `t`s are visited in strictly increasing order too. So instead of
/// collecting every matching `(k, t)` (what the old `Vec<(k, t)>` per key did: `O(N loglog
/// sqrt(max m_k))` words total, measured as the dominant term at `n = 1e6..3e6`, see the module
/// docs), each side only needs its first and last `t` seen: the step (`p` for Lucas, `p^2` for
/// p-adic) is recoverable from the key, so `(t_min, t_max)` plus that step describes the *exact*
/// same membership with two `u64`s instead of up to `O(N/p)` pairs. This collapses both
/// `lucas_small` and `padic_small` to `O(pi(sqrt(max m_k)))` words, matching what the module
/// docs originally (incorrectly) claimed the `HashMap`s already cost.
#[derive(Default, Clone, Copy)]
struct SideRanges {
    /// Low side (`k = t`): inclusive `(t_min, t_max)` of every eligible `t` seen.
    lo: Option<(u64, u64)>,
    /// High side (`k = N-1-t`): same, in `t`-space.
    hi: Option<(u64, u64)>,
}

impl SideRanges {
    /// Records one more eligible `t` on the given side. Relies on the caller visiting `t` in
    /// non-decreasing order per side (true for both [`stream_needs_and_chunks`]'s scan and,
    /// since it's built from that same scan's output, the cofactor regrouping in [`c_part`]).
    fn extend(&mut self, is_high: bool, t: u64) {
        let slot = if is_high { &mut self.hi } else { &mut self.lo };
        match slot {
            Some((_, t_max)) => *t_max = t,
            None => *slot = Some((t, t)),
        }
    }
}

/// Every non-`Main` contribution for one `k` (doc §4.3): small-prime Lucas needs (`p <= t`,
/// `e == 1`), small-prime p-adic needs (`p <= t`, `e >= 2`), and the one-off "cofactor Lucas"
/// case (leftover cofactor `<= t`) that doc §7 calls out as the part that doesn't fit neatly
/// into a per-prime arithmetic progression (the cofactor is essentially unique to this `k`, not
/// shared by a residue class of other `k`s the way a small sieve prime is). Pushed into the
/// caller's collectors rather than returned, so a hot per-`k` loop doesn't allocate.
fn classify_extras(
    t: u64,
    k: u64,
    rem: u64,
    fac: &[(u64, u32)],
    lucas_small: &mut HashMap<u64, SideRanges>,
    padic_small: &mut HashMap<u64, SideRanges>,
    cofactor: &mut Vec<(u64, u64, u64)>, // (p, k, t)
) {
    let is_high = k != t;
    for &(p, e) in fac {
        if p > t {
            continue;
        }
        if e == 1 {
            lucas_small.entry(p).or_default().extend(is_high, t);
        } else {
            padic_small.entry(p).or_default().extend(is_high, t);
        }
    }
    if rem > 1 && rem <= t {
        cofactor.push((rem, k, t));
    }
}

/// Everything [`c_part`] needs from a single streaming pass over `k in [0, N)` (via `t in [0,
/// N/2)`, both `k = t` and `k = N-1-t` per doc §4.4's "come out naturally in target order"):
/// where to cut the Main-item chunks (target windows, `O(N log n / mem_bits)` of them, doc
/// §5.1), and the small-prime/cofactor Lucas and p-adic needs. Everything here is `O(pi(sqrt(max
/// m_k)))` or `O(chunks)` sized *except* `cofactor` (doc §7: the one remaining piece whose size
/// isn't bounded independent of `N` — see the module docs).
struct StreamNeeds {
    /// Main-item chunk boundaries, as target windows `[t0, t1)`; [`art_chunk_by_range`]
    /// regenerates each chunk's items from scratch (cheap re-factoring) rather than this pass
    /// keeping them around.
    chunk_bounds: Vec<(u64, u64)>,
    /// `p -> ranges` for small primes (`p <= sqrt(max m_k)`, `e == 1`): see [`SideRanges`].
    lucas_small: HashMap<u64, SideRanges>,
    /// `p -> ranges` for small primes with `e >= 2` (step `p^2`, exact `e` recomputed per `k` at
    /// consumption time — see [`padic_consumers`] — rather than stored, since within one `p` it
    /// varies member-to-member and storing it was exactly the `Vec<(k,t,e)>` this replaces).
    padic_small: HashMap<u64, SideRanges>,
    /// `(p, k, t)` triples for the cofactor-Lucas case (doc §7 caveat: `O(fraction * N)`, not
    /// bounded independent of `N`).
    cofactor: Vec<(u64, u64, u64)>,
}

/// Batch size for [`stream_needs_and_chunks`]'s factoring windows: large enough to amortise
/// per-call overhead, small enough that its `O(batch)` temporary memory never approaches `N`
/// (it's sized off `sqrt(half)`, not `half` itself).
fn stream_batch_size(half: u64) -> u64 {
    ((half as f64).sqrt() as u64).clamp(1024, 1 << 16)
}

/// The single sequential pass over `t in [0, N/2)` (doc §4.4): factors both `k = t` and `k =
/// N-1-t` in `O(batch)`-sized windows via [`factor_window`], decides Main-item chunk boundaries
/// by running bit total (doc's hint: "a cheap first pass counting bits per window without
/// storing items" — Main items themselves are discarded here, only their bit-length counts),
/// and collects the small-prime/cofactor Lucas and p-adic needs. Sequential because chunk-cutting
/// is inherently a running accumulation; factoring itself is a small fraction of total cost
/// (doc §5.3, and the profiling note in `nthdigit.rs`'s module docs), so this doesn't cost
/// parallelism where it matters.
///
/// Requires `big_n` even (always true for a real run — [`Params2::new`] guarantees it — so the
/// `t in [0, N/2)` pairing `(k = t, k = N-1-t)` covers every `k` exactly once with no leftover
/// middle element).
fn stream_needs_and_chunks(
    base: u64,
    big_n: u64,
    small_primes: &[u64],
    mem_bits: u64,
) -> StreamNeeds {
    debug_assert!(big_n.is_multiple_of(2), "big_n must be even, got {big_n}");
    let half = big_n / 2;
    let batch = stream_batch_size(half);
    let mut chunk_bounds = Vec::new();
    let mut lucas_small: HashMap<u64, SideRanges> = HashMap::new();
    let mut padic_small: HashMap<u64, SideRanges> = HashMap::new();
    let mut cofactor = Vec::new();
    let mut bits_acc = 0u64;
    let mut chunk_start = 0u64;
    let mut t0 = 0u64;
    while t0 < half {
        let t1 = (t0 + batch).min(half);
        let len = t1 - t0;
        let (rem_lo, fac_lo) = factor_window(base, t0, len, small_primes);
        let (rem_hi, fac_hi) = factor_window(base, big_n - t1, len, small_primes);
        for t in t0..t1 {
            let i_lo = (t - t0) as usize;
            let i_hi = (t1 - 1 - t) as usize;
            let k_lo = t;
            let k_hi = big_n - 1 - t;

            classify_extras(
                t,
                k_lo,
                rem_lo[i_lo],
                &fac_lo[i_lo],
                &mut lucas_small,
                &mut padic_small,
                &mut cofactor,
            );
            if let Some(q) = main_modulus(t, rem_lo[i_lo], &fac_lo[i_lo]) {
                bits_acc += 64 - q.leading_zeros() as u64;
            }
            classify_extras(
                t,
                k_hi,
                rem_hi[i_hi],
                &fac_hi[i_hi],
                &mut lucas_small,
                &mut padic_small,
                &mut cofactor,
            );
            if let Some(q) = main_modulus(t, rem_hi[i_hi], &fac_hi[i_hi]) {
                bits_acc += 64 - q.leading_zeros() as u64;
            }

            if bits_acc >= mem_bits {
                chunk_bounds.push((chunk_start, t + 1));
                chunk_start = t + 1;
                bits_acc = 0;
            }
        }
        t0 = t1;
    }
    if chunk_start < half {
        chunk_bounds.push((chunk_start, half));
    }
    StreamNeeds {
        chunk_bounds,
        lucas_small,
        padic_small,
        cofactor,
    }
}

// ---------------------------------------------------------------------------------
// The recurrence and binary splitting (doc §4.2, §4.4)
// ---------------------------------------------------------------------------------

/// `(alpha, delta, tau)` for the recurrence block `(j1, j2]`: `P' = alpha*P`,
/// `T' = delta*T + tau*P`, `D' = delta*D` (doc §4.2). Below `LEAF`, computed directly by the
/// step recurrence `P <- (N-j+1)P; T <- jT + (N-j+1)P_old; D <- jD`; above it, split in half
/// and composed via `(a1*a2, d1*d2, d2*t1 + t2*a1)` (standard binary-splitting composition).
const LEAF: u64 = 24;

fn bs(big_n: u64, j1: u64, j2: u64) -> (Big, Big, Big) {
    if j2 - j1 <= LEAF {
        let (mut a, mut d, mut t) = (Big::one(), Big::one(), Big::zero());
        for j in (j1 + 1)..=j2 {
            let u = big_n - j + 1;
            t = t.mul_u64(j).add(&a.mul_u64(u));
            a = a.mul_u64(u);
            d = d.mul_u64(j);
        }
        return (a, d, t);
    }
    let mid = (j1 + j2) >> 1;
    let (a1, d1, t1) = bs(big_n, j1, mid);
    let (a2, d2, t2) = bs(big_n, mid, j2);
    let tau = d2.mul(&t1).add(&t2.mul(&a1));
    (a1.mul(&a2), d1.mul(&d2), tau)
}

/// Applies steps `(j1, j2]` to `state` modulo `Q`, in groups sized so each group's binary
/// splitting stays close to `bits(Q)` bits (doc §4.4 step 2/3: "groups of `g = m/log2 N`").
fn advance(
    big_n: u64,
    state: (Big, Big, Big),
    j1: u64,
    j2: u64,
    q: &Big,
    lg_n: u32,
) -> (Big, Big, Big) {
    if j2 <= j1 {
        return state;
    }
    let qb = q.bits().max(64);
    let g = (qb / lg_n as u64).max(8);
    let (mut p, mut t, mut d) = state;
    let mut j = j1;
    while j < j2 {
        let e = (j + g).min(j2);
        let (a, dd, tt) = bs(big_n, j, e);
        let new_p = a.mul(&p).rem(q);
        let new_t = dd.mul(&t).add(&tt.mul(&p)).rem(q);
        let new_d = dd.mul(&d).rem(q);
        p = new_p;
        t = new_t;
        d = new_d;
        j = e;
    }
    (p, t, d)
}

// ---------------------------------------------------------------------------------
// ART items (doc §4.3, §4.4): partial-fraction parts of each m_k
// ---------------------------------------------------------------------------------

/// What an ART item's leaf result feeds into.
#[derive(Clone, Copy, Debug)]
enum Tag {
    /// `m_k`'s "large primes" part (all `p > t_k`): result folds directly into the
    /// accumulator via [`add_contribution`].
    Main(u64),
    /// A Lucas item shared by every `k` with `p || m_k` (exponent 1) and `t_k mod p == r`:
    /// result is `(s_r, C(N,r)) mod p`, consumed later by [`lucas_consumers`].
    Lucas(u64, u64),
}

struct Item {
    /// Target `t` this item resolves the recurrence state at.
    t: u64,
    /// Modulus (fits `u64`, see module docs).
    q: u64,
    tag: Tag,
}

/// `t! mod m` used by the Lucas/p-adic identities is never computed directly; instead every
/// part's contribution is folded in here: `(-1)^k * frac(X_k * s_mod_q * ((m_k/q)^-1 mod q) /
/// q)`, the partial-fraction reassembly from doc §3 "obstacles 2/4". `s_mod_q` must already be
/// `s_{t_k} mod q` (i.e. the *un-mirrored* value at the item's own target `t`); the `k != t`
/// symmetry flip (doc §4.2 "Symmetry") happens inside.
fn add_contribution(n: u64, big_n: u64, base: u64, k: u64, t: u64, q: u64, s_mod_q: u64) -> u128 {
    let mut s = s_mod_q % q;
    if k != t {
        let two_n = powmod(2, big_n, q);
        s = (two_n + q - s) % q;
    }
    let m_k = base + 2 * k;
    let cofactor = (m_k / q) % q;
    let a = mod_inverse(cofactor, q);
    let e1 = big_n - 2;
    let e2 = (n + 2) - big_n;
    let pow5 = powmod(5, e1, q);
    let pow10 = powmod(10, e2, q);
    let xk = mulmod(pow5, pow10, q);
    let y = mulmod(mulmod(xk, s, q), a, q);
    signed(frac_fixed_point(y, q), k)
}

// ---------------------------------------------------------------------------------
// Lucas' theorem (doc §4.3 "Lucas formula")
// ---------------------------------------------------------------------------------

/// `s_t mod p` from `s_r, C(N,r) mod p` (`r = t mod p`, from an ART item) via Lucas' theorem
/// (doc §4.3). `cache` memoises, per prime `p`, the digit-`i` binomial row/prefix-sum table
/// (`(p, i) -> (C(N_i, ·), prefix sums)`) shared by every `k` that needs prime `p`.
fn lucas_s(big_n: u64, t: u64, p: u64, sr: u64, cr: u64, cache: &mut LucasRowCache) -> u64 {
    let mut n_digits = Vec::new();
    let mut t_digits = Vec::new();
    let (mut x, mut y) = (big_n, t);
    while x > 0 {
        n_digits.push(x % p);
        t_digits.push(y % p);
        x /= p;
        y /= p;
    }
    let depth = n_digits.len();
    if depth == 0 {
        // N == 0: s_t is trivially 0 or 1, handled by the t >= N shortcut upstream. Defensive.
        return 0;
    }
    let mut pow2low = vec![1u64; depth + 1];
    for i in 0..depth {
        pow2low[i + 1] = mulmod(pow2low[i], powmod(2, n_digits[i], p), p);
    }
    let mut hi_prod = 1u64;
    let mut acc = 0u64;
    for i in (0..depth).rev() {
        let (g, cb) = if i == 0 {
            ((sr + p - cr % p) % p, cr % p)
        } else {
            let ni = n_digits[i];
            let (cs, pref) = cache.entry((p, i as u32)).or_insert_with(|| {
                let mut c = 1u64;
                let mut cs = Vec::with_capacity(ni as usize + 1);
                let mut pref = vec![0u64];
                for yv in 0..=ni {
                    if yv > 0 {
                        let inv_y = mod_inverse(yv, p);
                        c = mulmod(mulmod(c, (ni - yv + 1) % p, p), inv_y, p);
                    }
                    cs.push(c);
                    let last = *pref.last().unwrap();
                    pref.push((last + c) % p);
                }
                (cs, pref)
            });
            let ti = t_digits[i] as usize;
            let g = pref[ti.min(cs.len())];
            let cb = if ti < cs.len() { cs[ti] } else { 0 };
            (g, cb)
        };
        acc = (acc + mulmod(mulmod(hi_prod, g, p), pow2low[i], p)) % p;
        hi_prod = mulmod(hi_prod, cb, p);
        if hi_prod == 0 {
            break;
        }
    }
    (acc + hi_prod) % p
}

// ---------------------------------------------------------------------------------
// p-adic Lucas recursion for p^e, e >= 2 (doc §4.3 "p-adic recursion")
// ---------------------------------------------------------------------------------

/// `s_t(N) mod p^e` and `C(N,t) mod p^e` for a fixed prime `p`, via the recursion
/// `(1+x)^N = (1+x)^N0 · ((1+x^p) + p·g(x))^K` (`N = pK+N0`, doc §4.3). Builds `O(p)`-word
/// tables (`gp`: powers of `g` up to `g^(emax-1)`) once per prime; callers build one of these
/// per prime and drop it before moving to the next, so peak memory is `O(max p)` over the
/// primes actually used (doc §5.2), not `O(Σp)`.
struct PadicBinom {
    p: u64,
    pow_p: Vec<u64>,                                 // pow_p[i] = p^i, i <= emax
    gp: Vec<Vec<u64>>, // gp[i] = coefficients of g^i mod p^emax, i < emax
    rows: HashMap<(u64, u64), (Vec<u64>, Vec<u64>)>, // (N0, mod) -> (coeffs, prefix sums)
    /// `(N, t, e) -> s_t(N) mod p^e` / `C(N,t) mod p^e`, memoising *within one query's*
    /// recursion tree (bounded, `O(e^2 log_p N)` per doc §5.2's per-query cost). **Must be
    /// cleared between top-level queries** ([`PadicBinom::clear_query_memo`]) — see that
    /// method's docs for why: left to grow across an entire prime's `klist`, this was measured
    /// to reach `O(klist.len())` entries for small `p` (most of one profiled 447 MiB peak at
    /// `n=1e7`), because a top-level call's `t` is essentially unique per `k`, so nothing above
    /// the first recursion level ever gets reused across different `k`s anyway.
    memo_s: HashMap<(u64, i64, u32), u64>,
    memo_c: HashMap<(u64, i64, u32), u64>,
}

impl PadicBinom {
    fn new(p: u64, emax: u32) -> Self {
        let mut pow_p = vec![1u64; emax as usize + 1];
        for i in 1..=emax as usize {
            pow_p[i] = pow_p[i - 1] * p;
        }
        let pe = pow_p[emax as usize];
        // g[i] = C(p-1, i-1) / i mod pe, for i in 1..p, via the running-product trick.
        let mut g = vec![0u64; p as usize];
        let mut c = 1u64 % pe.max(1);
        for i in 1..p {
            let inv_i = mod_inverse(i % pe, pe);
            g[i as usize] = mulmod(c, inv_i, pe);
            c = mulmod(mulmod(c, (p - i) % pe, pe), inv_i, pe);
        }
        let mut gp: Vec<Vec<u64>> = vec![vec![1u64]];
        for _ in 1..emax {
            let prev = gp.last().unwrap();
            let mut nxt = vec![0u64; prev.len() + p as usize - 1];
            for (a_, &x) in prev.iter().enumerate() {
                if x == 0 {
                    continue;
                }
                for (b_, &gb) in g.iter().enumerate().skip(1) {
                    if gb == 0 {
                        continue;
                    }
                    let idx = a_ + b_;
                    nxt[idx] = (nxt[idx] + mulmod(x, gb, pe)) % pe;
                }
            }
            gp.push(nxt);
        }
        PadicBinom {
            p,
            pow_p,
            gp,
            rows: HashMap::new(),
            memo_s: HashMap::new(),
            memo_c: HashMap::new(),
        }
    }

    /// Coefficients and prefix sums of `(1+x)^N0` mod `modulus`, `N0 < p`. Memoised per
    /// `(N0, modulus)` pair (several `modulus = p^(e-i)` values are used across recursion
    /// levels).
    fn row(&mut self, n0: u64, modulus: u64) -> &(Vec<u64>, Vec<u64>) {
        self.rows.entry((n0, modulus)).or_insert_with(|| {
            let mut a = Vec::with_capacity(n0 as usize + 1);
            let mut c = 1u64 % modulus.max(1);
            for y in 0..=n0 {
                if y > 0 {
                    let inv_y = mod_inverse(y % modulus, modulus);
                    c = mulmod(mulmod(c, (n0 - y + 1) % modulus, modulus), inv_y, modulus);
                }
                a.push(c);
            }
            let mut pref = Vec::with_capacity(a.len());
            let mut s = 0u64;
            for &v in &a {
                s = (s + v) % modulus;
                pref.push(s);
            }
            (a, pref)
        })
    }

    /// Coefficient `u` (or, if `prefix`, the prefix sum up to `u`) of `h_i = (1+x)^N0 * g^i`,
    /// mod `modulus`.
    fn h(&mut self, n0: u64, i: usize, u: u64, modulus: u64, prefix: bool) -> u64 {
        let gi_len;
        {
            let gi = &self.gp[i];
            gi_len = gi.len();
        }
        let (a, pref) = self.row(n0, modulus).clone(); // small (<= p entries); clone avoids
        // holding an immutable borrow of self.rows across the mutable gp access below.
        let gi = &self.gp[i];
        let lo = if prefix { 0 } else { u.saturating_sub(n0) };
        let hi_z = u.min(gi_len as u64 - 1);
        let mut tot = 0u128;
        let mut z = lo;
        while z <= hi_z {
            let gz = gi[z as usize];
            if gz != 0 {
                let w = u - z;
                let term = if prefix {
                    gz as u128 * a_prefix_at(&pref, w, n0) as u128
                } else {
                    gz as u128 * a[w as usize] as u128
                };
                tot += term;
            }
            z += 1;
        }
        (tot % modulus as u128) as u64
    }

    fn s(&mut self, big_n: u64, t: i64, e: u32) -> u64 {
        if e == 0 || t < 0 {
            return 0;
        }
        let modulus = self.pow_p[e as usize];
        if t as u64 >= big_n {
            return powmod(2, big_n, modulus);
        }
        let t = t as u64;
        if let Some(&v) = self.memo_s.get(&(big_n, t as i64, e)) {
            return v;
        }
        let p = self.p;
        let r = if big_n < p {
            self.row(big_n, modulus).1[t as usize]
        } else {
            let (k, n0) = (big_n / p, big_n % p);
            let (tt, rr) = (t / p, t % p);
            let mut r = 0u128;
            for i in 0..e {
                if k < i as u64 {
                    break;
                }
                let coef_mod = Big::binomial(k, i)
                    .mul_u64(self.pow_p[i as usize])
                    .rem_u64(modulus);
                if coef_mod == 0 {
                    continue;
                }
                let m2 = self.pow_p[(e - i) as usize];
                let full = mulmod(
                    powmod(2, n0, m2),
                    self.gp[i as usize]
                        .iter()
                        .fold(0u64, |acc, &v| (acc + v) % m2),
                    m2,
                );
                let mut part = mulmod(
                    full,
                    self.s(k - i as u64, tt as i64 - i as i64 - 1, e - i),
                    m2,
                ) as u128;
                for d in 0..=i as u64 {
                    let hv = self.h(n0, i as usize, d * p + rr, m2, true);
                    let cv = self.c(k - i as u64, tt as i64 - d as i64, e - i);
                    part += hv as u128 * cv as u128;
                }
                r += coef_mod as u128 * (part % m2 as u128);
            }
            (r % modulus as u128) as u64
        };
        self.memo_s.insert((big_n, t as i64, e), r);
        r
    }

    fn c(&mut self, big_n: u64, t: i64, e: u32) -> u64 {
        if e == 0 || t < 0 || t as u64 > big_n {
            return 0;
        }
        let modulus = self.pow_p[e as usize];
        let t = t as u64;
        if let Some(&v) = self.memo_c.get(&(big_n, t as i64, e)) {
            return v;
        }
        let p = self.p;
        let r = if big_n < p {
            self.row(big_n, modulus).0[t as usize]
        } else {
            let (k, n0) = (big_n / p, big_n % p);
            let (tt, rr) = (t / p, t % p);
            let mut r = 0u128;
            for i in 0..e {
                if k < i as u64 {
                    break;
                }
                let coef_mod = Big::binomial(k, i)
                    .mul_u64(self.pow_p[i as usize])
                    .rem_u64(modulus);
                if coef_mod == 0 {
                    continue;
                }
                let m2 = self.pow_p[(e - i) as usize];
                let mut part = 0u128;
                for d in 0..=i as u64 {
                    let hv = self.h(n0, i as usize, d * p + rr, m2, false);
                    let cv = self.c(k - i as u64, tt as i64 - d as i64, e - i);
                    part += hv as u128 * cv as u128;
                }
                r += coef_mod as u128 * (part % m2 as u128);
            }
            (r % modulus as u128) as u64
        };
        self.memo_c.insert((big_n, t as i64, e), r);
        r
    }

    /// Clears `memo_s`/`memo_c` between top-level queries (see their field docs). Cheap
    /// (`HashMap::clear` keeps the allocation, so the capacity settles at the largest single
    /// query's recursion tree, not the number of queries) and correctness-neutral: `s`/`c` are
    /// pure functions of `(self.p, big_n, t, e)`, so discarding memo entries only forces
    /// recomputation, never a different answer.
    fn clear_query_memo(&mut self) {
        self.memo_s.clear();
        self.memo_c.clear();
    }
}

/// Prefix-sum lookup clamped to `N0` (the degree of the row), matching Python's
/// `A[min(w, N0)]`: once past the polynomial's degree the prefix sum is constant.
fn a_prefix_at(pref: &[u64], w: u64, n0: u64) -> u64 {
    pref[(w.min(n0)) as usize]
}

// ---------------------------------------------------------------------------------
// The C part: assembling items, running the ART, Lucas/p-adic consumers
// ---------------------------------------------------------------------------------

/// Sorted-by-target ART items are cut into chunks whose modulus product has `~mem_bits` bits
/// (doc §4.4), each processed independently (they only share the *shape* of the recurrence,
/// not any live state) so chunks can run in parallel.
fn chunk_bounds(items: &[Item], mem_bits: u64) -> Vec<(usize, usize)> {
    let mut bounds = Vec::new();
    let mut i = 0usize;
    while i < items.len() {
        let mut bits = 0u64;
        let mut j = i;
        while j < items.len() && (bits < mem_bits || j == i) {
            bits += 64 - items[j].q.leading_zeros() as u64;
            j += 1;
        }
        bounds.push((i, j));
        i = j;
    }
    bounds
}

/// Runs the ART over one chunk of items (already sorted by target `t`): builds the product
/// tree of their moduli, advances the recurrence from step 0 to the chunk's first target mod
/// `Q`, then descends the tree, reducing mod each subtree's modulus and re-advancing between
/// siblings (doc §4.4 steps 1-4). `Main` leaves fold straight into the accumulator; `Lucas`
/// leaves are returned for the caller to feed into [`lucas_consumers`].
fn art_chunk(
    n: u64,
    big_n: u64,
    base: u64,
    items: &[Item],
    lg_n: u32,
) -> (u128, u64, Vec<LucasLeaf>) {
    let l = items.len();
    let mut tree: HashMap<(usize, usize), Big> = HashMap::new();

    fn build(lo: usize, hi: usize, items: &[Item], tree: &mut HashMap<(usize, usize), Big>) -> Big {
        if hi - lo == 1 {
            let v = Big::from_u64(items[lo].q);
            tree.insert((lo, hi), v.clone());
            return v;
        }
        let mid = (lo + hi) / 2;
        let l = build(lo, mid, items, tree);
        let r = build(mid, hi, items, tree);
        let prod = l.mul(&r);
        tree.insert((lo, hi), prod.clone());
        prod
    }

    let q_total = build(0, l, items, &mut tree);
    let one = (Big::one(), Big::one(), Big::one());
    let x0 = advance(big_n, one, 0, items[0].t, &q_total, lg_n);

    #[allow(clippy::too_many_arguments)]
    fn rec(
        n: u64,
        big_n: u64,
        base: u64,
        lo: usize,
        hi: usize,
        x: (Big, Big, Big),
        items: &[Item],
        tree: &HashMap<(usize, usize), Big>,
        lg_n: u32,
        acc: &mut u128,
        terms: &mut u64,
        lucas_out: &mut Vec<LucasLeaf>,
    ) {
        if hi - lo == 1 {
            let it = &items[lo];
            let q_big = Big::from_u64(it.q);
            let pm = x.0.rem(&q_big).to_u64().unwrap();
            let tm = x.1.rem(&q_big).to_u64().unwrap();
            let dm = x.2.rem(&q_big).to_u64().unwrap();
            match it.tag {
                Tag::Main(k) => {
                    let dinv = mod_inverse(dm, it.q);
                    let s = mulmod(tm, dinv, it.q);
                    *acc = acc.wrapping_add(add_contribution(n, big_n, base, k, it.t, it.q, s));
                    *terms += 1;
                }
                Tag::Lucas(p, r) => {
                    let dinv = mod_inverse(dm, it.q);
                    let sr = mulmod(tm, dinv, it.q);
                    let cr = mulmod(pm, dinv, it.q);
                    lucas_out.push(((p, r), (sr, cr)));
                }
            }
            return;
        }
        let mid = (lo + hi) / 2;
        let ql = tree.get(&(lo, mid)).unwrap();
        let qr = tree.get(&(mid, hi)).unwrap();
        let xl = (x.0.rem(ql), x.1.rem(ql), x.2.rem(ql));
        rec(
            n, big_n, base, lo, mid, xl, items, tree, lg_n, acc, terms, lucas_out,
        );
        let xr0 = (x.0.rem(qr), x.1.rem(qr), x.2.rem(qr));
        let xr = advance(big_n, xr0, items[lo].t, items[mid].t, qr, lg_n);
        rec(
            n, big_n, base, mid, hi, xr, items, tree, lg_n, acc, terms, lucas_out,
        );
    }

    let mut acc = 0u128;
    let mut terms = 0u64;
    let mut lucas_out = Vec::new();
    rec(
        n,
        big_n,
        base,
        0,
        l,
        x0,
        items,
        &tree,
        lg_n,
        &mut acc,
        &mut terms,
        &mut lucas_out,
    );
    (acc, terms, lucas_out)
}

/// One Main-item chunk's worth of the ART, for the streaming path: regenerates the chunk's items
/// from scratch (target window `[t0, t1)`, both `k = t` and `k = N-1-t`, doc §4.4's "generate
/// each chunk's items on the fly by segment-factoring just those k-windows") instead of slicing
/// a pre-built global item list, then defers to the unchanged [`art_chunk`]. Doubles the
/// small-prime factoring work for this window (it was already done once by
/// [`stream_needs_and_chunks`] to size the chunk) in exchange for never holding more than one
/// chunk's items (`O(mem_bits / log n)`, doc §5.1) at a time. Every item here is `Tag::Main`
/// (small-prime and cofactor Lucas/p-adic needs were already pulled out by
/// [`classify_extras`]), so there's nothing for the caller to do with `art_chunk`'s `lucas_out`.
fn art_chunk_by_range(
    n: u64,
    big_n: u64,
    base: u64,
    small_primes: &[u64],
    t0: u64,
    t1: u64,
    lg_n: u32,
) -> (u128, u64) {
    let len = t1 - t0;
    let (rem_lo, fac_lo) = factor_window(base, t0, len, small_primes);
    let (rem_hi, fac_hi) = factor_window(base, big_n - t1, len, small_primes);
    let mut items: Vec<Item> = Vec::new();
    for t in t0..t1 {
        let i_lo = (t - t0) as usize;
        let i_hi = (t1 - 1 - t) as usize;
        if let Some(q) = main_modulus(t, rem_lo[i_lo], &fac_lo[i_lo]) {
            items.push(Item {
                t,
                q,
                tag: Tag::Main(t),
            });
        }
        if let Some(q) = main_modulus(t, rem_hi[i_hi], &fac_hi[i_hi]) {
            items.push(Item {
                t,
                q,
                tag: Tag::Main(big_n - 1 - t),
            });
        }
    }
    if items.is_empty() {
        // Can happen for a target window where every k's prime factors are all <= t (nothing
        // for Main to do) — art_chunk's product-tree `build` needs at least one item.
        return (0, 0);
    }
    let (acc, terms, lucas_out) = art_chunk(n, big_n, base, &items, lg_n);
    debug_assert!(
        lucas_out.is_empty(),
        "art_chunk_by_range's items are all Tag::Main"
    );
    (acc, terms)
}

/// Runs a batch of `Tag::Lucas` items (already built, sorted by target) through the same
/// chunked ART as Main items (`O(mem_bits)` per chunk, doc §4.4), returning the resolved
/// `(s_r, C(N,r))` values keyed by `(p, r)`. Shared by the small-prime pass and the cofactor
/// pass in [`c_part`] — the ART's cost model doesn't care *why* an item's target is what it is,
/// only that items arrive sorted (doc's "an ART item's target is a function of k" hint applies
/// equally to a Lucas item's `t = r`).
fn resolve_lucas_items(
    n: u64,
    big_n: u64,
    base: u64,
    items: &mut [Item],
    mem_bits: u64,
    lg_n: u32,
) -> HashMap<(u64, u64), LucasVal> {
    items.sort_by_key(|it| it.t);
    let bounds = chunk_bounds(items, mem_bits);
    let mapped =
        maybe_par_iter!(&bounds).map(|&(lo, hi)| art_chunk(n, big_n, base, &items[lo..hi], lg_n).2);
    #[cfg(feature = "parallel")]
    {
        mapped
            .fold(HashMap::new, |mut map, pairs| {
                map.extend(pairs);
                map
            })
            .reduce(HashMap::new, |mut m1, m2| {
                m1.extend(m2);
                m1
            })
    }
    #[cfg(not(feature = "parallel"))]
    {
        mapped.fold(HashMap::new(), |mut map, pairs| {
            map.extend(pairs);
            map
        })
    }
}

/// One [`Tag::Lucas`] ART item per populated side of every prime in `need` (doc §4.4: "all `k`
/// with `p | m_k` share one or two such items"). Target `r = t_min % p` — well-defined because
/// every `t` on one side of one prime's [`SideRanges`] shares the same residue mod `p` (that's
/// exactly why the range compresses to two numbers in the first place).
fn lucas_art_items(need: &HashMap<u64, SideRanges>) -> Vec<Item> {
    let mut items = Vec::new();
    for (&p, ranges) in need {
        for (t_min, _) in [ranges.lo, ranges.hi].into_iter().flatten() {
            let r = t_min % p;
            items.push(Item {
                t: r,
                q: p,
                tag: Tag::Lucas(p, r),
            });
        }
    }
    items
}

/// Resolves every `k` that needed a Lucas item (prime `p`, exponent 1) from the `(s_r, C(N,r))`
/// values the ART produced, and folds each into the accumulator. Walks each side's `(t_min,
/// t_max)` range directly (step `p`, doc §4.4/§7: the whole point of [`SideRanges`] is that this
/// is the only place the full `k`/`t` list needs to exist, generated on the fly instead of held).
fn lucas_consumers(
    n: u64,
    big_n: u64,
    base: u64,
    lucas_need: &HashMap<u64, SideRanges>,
    lucas_val: &HashMap<(u64, u64), LucasVal>,
) -> (u128, u64) {
    let iter = maybe_par_iter!(lucas_need).map(|(&p, ranges)| {
        let mut cache = HashMap::new();
        let mut acc = 0u128;
        let mut terms = 0u64;
        for (is_high, range) in [(false, ranges.lo), (true, ranges.hi)] {
            let Some((t_min, t_max)) = range else {
                continue;
            };
            let r = t_min % p;
            let (sr, cr) = lucas_val[&(p, r)];
            let mut t = t_min;
            loop {
                let k = if is_high { big_n - 1 - t } else { t };
                // The "exactly e == 1" set has holes at the p^2 sub-progression (those k
                // belong to padic_small instead, doc §4.3's e>=2 case): a SideRanges range
                // only pins down min/max t, not that every step in between is a member, so
                // skip anything p^2 also divides here rather than double- or mis-counting.
                if !(base + 2 * k).is_multiple_of(p * p) {
                    let s = lucas_s(big_n, t, p, sr, cr, &mut cache);
                    acc = acc.wrapping_add(add_contribution(n, big_n, base, k, t, p, s));
                    terms += 1;
                }
                if t == t_max {
                    break;
                }
                t += p;
            }
        }
        (acc, terms)
    });
    maybe_reduce!(iter, || (0u128, 0u64), |(a1, c1), (a2, c2)| (
        a1.wrapping_add(a2),
        c1 + c2
    ))
}

/// `v_p(base + 2k)`: the exact p-adic valuation of `m_k`, by trial division. Cheap (`e` is
/// small in practice — this is only ever called for `k`s [`stream_needs_and_chunks`] already
/// determined have `p^2 | m_k`, so `e >= 2`, and `p^e <= m_max` bounds it to `O(log_p m_max)`
/// divisions). Recomputed here rather than stored because, unlike `p` itself, `e` is *not*
/// constant along a [`SideRanges`] range (different members can have different exact powers of
/// `p`) — storing it per-`k` is exactly the `Vec<(k,t,e)>` this module no longer keeps.
fn padic_valuation(base: u64, k: u64, p: u64) -> u32 {
    let mut m = base + 2 * k;
    let mut e = 0u32;
    while m.is_multiple_of(p) {
        m /= p;
        e += 1;
    }
    e
}

/// Resolves every `k` that needed the p-adic recursion (`p^e | m_k`, `e >= 2`), one
/// [`PadicBinom`] table per prime (built and dropped independently, so peak memory across
/// primes run in parallel is `threads * O(max p)`, not `O(Σp)`). Walks each side's `(t_min,
/// t_max)` range directly (step `p^2`: `p^2 | m_k` is what put this prime's `t`s in
/// `padic_small` at all, doc §5.2's "arithmetic progression of difference p²").
fn padic_consumers(
    n: u64,
    big_n: u64,
    base: u64,
    padic_need: &HashMap<u64, SideRanges>,
    m_max: u64,
) -> (u128, u64) {
    let iter = maybe_par_iter!(padic_need).map(|(&p, ranges)| {
        // Safe emax bound: the largest e with p^e <= m_max (m_k never exceeds m_max), since
        // we no longer track each k's actual e up front (see padic_valuation's docs). p >=
        // 3 always (m_k is odd), so this loop is O(log_3 m_max) at worst, negligible.
        let mut emax = 2u32;
        while (p as u128).pow(emax + 1) <= m_max as u128 {
            emax += 1;
        }
        let mut pb = PadicBinom::new(p, emax);
        let mut acc = 0u128;
        let mut terms = 0u64;
        let step = p * p;
        for (is_high, range) in [(false, ranges.lo), (true, ranges.hi)] {
            let Some((t_min, t_max)) = range else {
                continue;
            };
            let mut t = t_min;
            loop {
                let k = if is_high { big_n - 1 - t } else { t };
                let e = padic_valuation(base, k, p);
                // Bounds memo_s/memo_c to one query's recursion tree instead of letting
                // them grow with the range's length — see PadicBinom's field docs.
                pb.clear_query_memo();
                let s = pb.s(big_n, t as i64, e);
                let q = pb.pow_p[e as usize];
                acc = acc.wrapping_add(add_contribution(n, big_n, base, k, t, q, s));
                terms += 1;
                if t == t_max {
                    break;
                }
                t += step;
            }
        }
        (acc, terms)
    });
    maybe_reduce!(iter, || (0u128, 0u64), |(a1, c1), (a2, c2)| (
        a1.wrapping_add(a2),
        c1 + c2
    ))
}

/// The full C part (streaming, doc §4.4/§7): one sequential pass ([`stream_needs_and_chunks`])
/// decides Main-item chunk boundaries and collects small-prime/cofactor Lucas and p-adic needs
/// without ever holding a global factor table or item list; Main-item chunks are then
/// regenerated and run through the ART in parallel ([`art_chunk_by_range`]); small-prime and
/// cofactor Lucas items each get their own bounded ART sub-pass ([`resolve_lucas_items`]); then
/// the usual Lucas/p-adic consumer passes. Returns `(C, term_count)`, `term_count` being the
/// exact number of rounded fixed-point terms folded in (for [`crate::nthdigit::error_units`]).
fn c_part(n: u64, p: Params2) -> (u128, u64) {
    let (big_m, big_n) = (p.big_m, p.big_n);
    let base = 2 * big_m * big_n + 1;
    let lg_n = 64 - big_n.leading_zeros();
    let m_max = base + 2 * (big_n - 1);
    let small_primes = primes_upto((m_max as f64).sqrt() as u64 + 2);
    crate::mem_profile::checkpoint("c_part: start");

    let needs = stream_needs_and_chunks(base, big_n, &small_primes, p.mem_bits);
    #[cfg(feature = "mem-profile")]
    log_needs_sizes(&needs, &small_primes);
    crate::mem_profile::checkpoint("c_part: after stream_needs_and_chunks");

    // Main items: regenerate + run each chunk independently, in parallel.
    let main_iter = maybe_par_iter!(&needs.chunk_bounds)
        .map(|&(t0, t1)| art_chunk_by_range(n, big_n, base, &small_primes, t0, t1, lg_n));
    let (main_acc, main_terms) = maybe_reduce!(
        main_iter,
        || (0u128, 0u64),
        |(a1, t1), (a2, t2)| (a1.wrapping_add(a2), t1 + t2)
    );
    crate::mem_profile::checkpoint("c_part: after Main chunks");

    // Small-prime Lucas items: their own bounded ART sub-pass (targets <= sqrt(max m_k), doc
    // §4.4 "process small primes in their own pass with their own small ART").
    let mut small_items = lucas_art_items(&needs.lucas_small);
    let small_lucas_val = resolve_lucas_items(n, big_n, base, &mut small_items, p.mem_bits, lg_n);
    let (lucas_acc, lucas_terms) =
        lucas_consumers(n, big_n, base, &needs.lucas_small, &small_lucas_val);
    crate::mem_profile::checkpoint("c_part: after small-prime Lucas");
    let (padic_acc, padic_terms) = padic_consumers(n, big_n, base, &needs.padic_small, m_max);
    crate::mem_profile::checkpoint("c_part: after p-adic");

    // Cofactor Lucas items (doc §7: the one piece whose count isn't bounded independent of N —
    // regroup the raw (p, k, t) triples into the same SideRanges shape as lucas_small (doc §7's
    // fix: k with q | m_k is an AP of step q here too, since q is prime — see SideRanges' docs),
    // then reuse the exact same machinery. This doesn't shrink `cofactor` itself (still O(N) —
    // most of these primes are each unique to one or two k, so the range rarely has more than
    // one member), but it does mean there's only one Lucas-item/consumer implementation to get
    // right, and it costs no more than the flat Vec did.
    let mut cofactor_need: HashMap<u64, SideRanges> = HashMap::new();
    for (p_, k, t) in needs.cofactor {
        cofactor_need.entry(p_).or_default().extend(k != t, t);
    }
    let mut cofactor_items = lucas_art_items(&cofactor_need);
    let cofactor_lucas_val =
        resolve_lucas_items(n, big_n, base, &mut cofactor_items, p.mem_bits, lg_n);
    let (cofactor_acc, cofactor_terms) =
        lucas_consumers(n, big_n, base, &cofactor_need, &cofactor_lucas_val);
    crate::mem_profile::checkpoint("c_part: after cofactor Lucas (end)");

    (
        main_acc
            .wrapping_add(lucas_acc)
            .wrapping_add(padic_acc)
            .wrapping_add(cofactor_acc),
        main_terms + lucas_terms + padic_terms + cofactor_terms,
    )
}

/// Prints an explicit heap-footprint breakdown of [`StreamNeeds`]' structures (task: "profile
/// memory composition" — a global allocator peak alone can't say *which* structure is big).
/// `.capacity()`, not `.len()`, since that's what's actually resident.
#[cfg(feature = "mem-profile")]
fn log_needs_sizes(needs: &StreamNeeds, small_primes: &[u64]) {
    use std::mem::{size_of, size_of_val};
    let lucas_bytes = needs.lucas_small.capacity() * (size_of::<u64>() + size_of::<SideRanges>());
    let padic_bytes = needs.padic_small.capacity() * (size_of::<u64>() + size_of::<SideRanges>());
    let cofactor_bytes = needs.cofactor.capacity() * size_of::<(u64, u64, u64)>();
    let chunk_bounds_bytes = needs.chunk_bounds.capacity() * size_of::<(u64, u64)>();
    let small_primes_bytes = size_of_val(small_primes);
    eprintln!(
        "[mem-profile] StreamNeeds breakdown: lucas_small={:.2} MiB ({} keys/primes) \
         padic_small={:.2} MiB ({} keys/primes) cofactor={:.2} MiB ({} entries) \
         chunk_bounds={:.2} MiB ({} chunks) small_primes={:.2} MiB ({} primes)",
        lucas_bytes as f64 / (1024.0 * 1024.0),
        needs.lucas_small.len(),
        padic_bytes as f64 / (1024.0 * 1024.0),
        needs.padic_small.len(),
        cofactor_bytes as f64 / (1024.0 * 1024.0),
        needs.cofactor.len(),
        chunk_bounds_bytes as f64 / (1024.0 * 1024.0),
        needs.chunk_bounds.len(),
        small_primes_bytes as f64 / (1024.0 * 1024.0),
        small_primes.len(),
    );
}

// ---------------------------------------------------------------------------------
// frac(10^n pi) and digit extraction
// ---------------------------------------------------------------------------------

/// `frac(10^n π)` to within `10^-n0`, using Theorem 2 (chunked ART with `~mem_bits`-bit
/// chunks) for the `C` part and Theorem 1's `B` part unchanged. Returns `(value, term_count)`
/// so callers can certify with [`crate::nthdigit::error_units`] using the true rounding count.
pub fn frac_10n_pi_m_with_terms(n: u64, n0: u32, mem_bits: u64) -> (u128, u64) {
    let p = Params2::new(n, n0, mem_bits);
    let b_terms = (p.big_m + 1) * p.big_n;
    let (b, (c, c_terms)) = maybe_join!(|| nthdigit::b_sum(n, b_terms), || c_part(n, p));
    (b.wrapping_sub(c), b_terms + c_terms)
}

/// `frac(10^n π)` to within `10^-n0`, as a 128-bit fixed-point fraction (value `= x / 2^128`).
/// Same fixed-point convention as [`crate::nthdigit::frac_10n_pi`]. See [`digits`] for
/// certified digit extraction with an MPFR fallback for small `n`.
pub fn frac_10n_pi_m(n: u64, n0: u32, mem_bits: u64) -> u128 {
    frac_10n_pi_m_with_terms(n, n0, mem_bits).0
}

/// A sensible default memory budget: `mem_bits ≈ 4·sqrt(n)·log2(10)` bits, i.e. about `4·sqrt(n)`
/// decimal digits — the "headline case" the reconstruction doc measures (`m ∝ √n` gives
/// Gourdon's `n^1.5` special case, doc §6.1).
pub fn default_mem_bits(n: u64) -> u64 {
    let bits = 4.0 * (n as f64).sqrt() * 10f64.log2();
    (bits.ceil() as u64).max(256)
}

/// `count` decimal digits of π at positions `n+1..=n+count`, computed with Theorem 2 at memory
/// budget `mem_bits`. Mirrors [`crate::nthdigit::digits`] exactly (guard-doubling
/// certification loop, [`nthdigit::digits_fallback`] for small/ill-conditioned `n`, chunking
/// beyond `MAX_N0`-certifiable length) — see that module's docs for the certification argument, which
/// applies unchanged since both algorithms use the same fixed-point encoding and the same
/// "one ulp per rounded term" error accounting (just with a different, exactly-counted, term
/// count here).
pub fn digits(n: u64, count: usize, mem_bits: u64) -> String {
    if count == 0 {
        return String::new();
    }
    if count > nthdigit::CHUNK {
        return (0..count)
            .step_by(nthdigit::CHUNK)
            .map(|i| digits(n + i as u64, nthdigit::CHUNK.min(count - i), mem_bits))
            .collect();
    }
    let mut guard: u64 = 4;
    loop {
        let n0 = count as u64 + guard;
        if n < nthdigit::SMALL_N_THRESHOLD || n < 4 * n0 || n0 > MAX_N0 as u64 {
            return nthdigit::digits_fallback(n, count);
        }
        let n0 = n0 as u32;
        let (x, terms) = frac_10n_pi_m_with_terms(n, n0, mem_bits);
        let err = nthdigit::error_units(n0, terms);
        let base = extract_digits(x, count);
        let plus = extract_digits(x.wrapping_add(err), count);
        let minus = extract_digits(x.wrapping_sub(err), count);
        if base == plus && base == minus {
            return base;
        }
        guard *= 2;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rug::Integer;
    use rug::ops::Pow;

    #[test]
    fn params_are_sane() {
        for &(n, n0, mem_bits) in &[
            (2000u64, 20u32, 512u64),
            (10_000, 30, 2048),
            (100_000, 40, 8192),
        ] {
            let p = Params2::new(n, n0, mem_bits);
            assert!(p.big_m.is_multiple_of(2));
            assert!(p.big_m >= 4);
            assert!(p.big_n.is_multiple_of(2));
            assert!(p.big_n <= n + 2);
            assert!(p.error_bound_log10() < -((n + n0 as u64) as f64));
        }
    }

    #[test]
    fn default_mem_bits_grows_like_sqrt_n() {
        let small = default_mem_bits(10_000);
        let big = default_mem_bits(1_000_000);
        // sqrt(1e6)/sqrt(1e4) = 10, so big should be roughly 10x small.
        let ratio = big as f64 / small as f64;
        assert!((8.0..12.0).contains(&ratio), "ratio={ratio}");
    }

    /// A tiny deterministic xorshift64 PRNG (matches `tests/nthdigit.rs`'s), so these tests are
    /// reproducible without a `rand` dependency.
    struct Xorshift(u64);
    impl Xorshift {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
    }

    /// Exact `Σ_{j=0}^{k} C(N,j)`, brute force via `rug::Integer::binomial`.
    fn s_k_exact(big_n: u64, k: u64) -> Integer {
        let mut sum = Integer::from(0);
        for j in 0..=k {
            sum += Integer::from(big_n).binomial(j as u32);
        }
        sum
    }

    // -----------------------------------------------------------------------------
    // Recurrence state (P_t, T_t, D_t) vs exact binomial sums (doc §4.2)
    // -----------------------------------------------------------------------------

    #[test]
    fn recurrence_state_matches_exact_binomial_sums() {
        for &(big_n, t) in &[
            (10u64, 0u64),
            (10, 5),
            (10, 10),
            (37, 20),
            (100, 50),
            (200, 1),
        ] {
            let (alpha, delta, tau) = bs(big_n, 0, t);
            // Folding the block onto the identity state (1,1,1): P_t=alpha, D_t=delta,
            // T_t=delta+tau (see module test docs / the derivation in the PR).
            let p_t = alpha;
            let d_t = delta;
            let t_t = d_t.add(&tau);

            let p_exact = {
                let mut prod = Integer::from(1);
                for j in 1..=t {
                    prod *= big_n - j + 1;
                }
                prod
            };
            let d_exact = {
                let mut prod = Integer::from(1);
                for j in 1..=t {
                    prod *= j;
                }
                prod
            };
            let t_exact = Integer::from(&d_exact * &s_k_exact(big_n, t));

            assert_eq!(
                p_t.to_decimal_string(),
                p_exact.to_string(),
                "N={big_n} t={t}: P_t mismatch"
            );
            assert_eq!(
                d_t.to_decimal_string(),
                d_exact.to_string(),
                "N={big_n} t={t}: D_t mismatch"
            );
            assert_eq!(
                t_t.to_decimal_string(),
                t_exact.to_string(),
                "N={big_n} t={t}: T_t mismatch"
            );
        }
    }

    #[test]
    fn advance_mod_q_matches_direct_reduction() {
        let lg_n = 10u32;
        for &(big_n, j1, j2, q) in &[
            (50u64, 0u64, 30u64, 97u64),
            (200, 10, 150, 9973),
            (500, 0, 500, 10_007),
        ] {
            let qb = Big::from_u64(q);
            let one_mod_q = Big::one().rem(&qb);
            let start = (one_mod_q.clone(), one_mod_q.clone(), one_mod_q);
            let (p, t, d) = advance(big_n, start, j1, j2, &qb, lg_n);

            // direct step-by-step reference, reduced mod q at every step
            let (mut pr, mut tr, mut dr) = (1u64 % q, 1u64 % q, 1u64 % q);
            for j in (j1 + 1)..=j2 {
                let u = big_n - j + 1;
                let new_p = mulmod(u % q, pr, q);
                let new_t = (mulmod(j % q, tr, q) + mulmod(u % q, pr, q)) % q;
                let new_d = mulmod(j % q, dr, q);
                pr = new_p;
                tr = new_t;
                dr = new_d;
            }
            assert_eq!(
                p.to_u64().unwrap(),
                pr,
                "N={big_n} j1={j1} j2={j2} q={q}: P"
            );
            assert_eq!(
                t.to_u64().unwrap(),
                tr,
                "N={big_n} j1={j1} j2={j2} q={q}: T"
            );
            assert_eq!(
                d.to_u64().unwrap(),
                dr,
                "N={big_n} j1={j1} j2={j2} q={q}: D"
            );
        }
    }

    // -----------------------------------------------------------------------------
    // Lucas path vs exact binomial sums
    // -----------------------------------------------------------------------------

    #[test]
    fn lucas_matches_exact_binomial_sums() {
        let mut rng = Xorshift(0x9E37_79B9_7F4A_7C15);
        let mut checked = 0;
        for _ in 0..400 {
            // Fresh cache per case: it's keyed on (p, digit-position) only, which is only
            // valid for a fixed N (true in the real algorithm, where N never changes mid-run;
            // here N is randomised per case, so a shared cache across cases would serve stale
            // rows from a *different* N's digits).
            let mut cache = HashMap::new();
            let big_n = 20 + rng.next() % 2000;
            let p = *[3u64, 5, 7, 11, 13, 17, 19, 23, 29, 31]
                .get((rng.next() % 10) as usize)
                .unwrap();
            if p > big_n {
                continue;
            }
            let t = rng.next() % (big_n + 1);
            let r = t % p;
            let sr: u64 = (s_k_exact(big_n, r) % Integer::from(p))
                .to_string()
                .parse()
                .unwrap();
            let cr: u64 = (Integer::from(big_n).binomial(r as u32) % Integer::from(p))
                .to_string()
                .parse()
                .unwrap();
            let got = lucas_s(big_n, t, p, sr, cr, &mut cache);
            let expected: u64 = (s_k_exact(big_n, t) % Integer::from(p))
                .to_string()
                .parse()
                .unwrap();
            assert_eq!(got, expected, "N={big_n} t={t} p={p} r={r}");
            checked += 1;
        }
        assert!(
            checked > 100,
            "too few valid (N,t,p) cases generated: {checked}"
        );
    }

    // -----------------------------------------------------------------------------
    // p-adic Lucas recursion vs exact binomial sums, e >= 1 (e=1 degenerates to Lucas)
    // -----------------------------------------------------------------------------

    #[test]
    fn padic_matches_exact_binomial_sums() {
        let mut rng = Xorshift(0xD1B5_4A32_D192_ED03);
        let primes = [2u64, 3, 5, 7, 11, 13, 17, 23, 31, 37, 53, 101];
        let mut checked = 0;
        for _ in 0..400 {
            let p = primes[(rng.next() % primes.len() as u64) as usize];
            let e = 1 + (rng.next() % 4) as u32; // 1..=4
            let big_n = 5 + rng.next() % 3000;
            let t = rng.next() % (big_n + 1);
            let modulus = p.pow(e);

            let mut pb = PadicBinom::new(p, e);
            let got_s = pb.s(big_n, t as i64, e);
            let got_c = pb.c(big_n, t as i64, e);

            let exp_s: u64 = (s_k_exact(big_n, t) % Integer::from(modulus))
                .to_string()
                .parse()
                .unwrap();
            let exp_c: u64 = (Integer::from(big_n).binomial(t as u32) % Integer::from(modulus))
                .to_string()
                .parse()
                .unwrap();
            assert_eq!(got_s, exp_s, "S: N={big_n} t={t} p={p} e={e}");
            assert_eq!(got_c, exp_c, "C: N={big_n} t={t} p={p} e={e}");
            checked += 1;
        }
        assert_eq!(checked, 400);
    }

    #[test]
    fn padic_matches_exact_binomial_sums_boundary_cases() {
        // t >= N and t == 0 edge cases, plus every e up to 5 for a couple of primes.
        for &p in &[3u64, 5, 7] {
            for e in 1..=5u32 {
                let modulus = p.pow(e);
                for &(big_n, t) in &[
                    (0u64, 0u64),
                    (1, 0),
                    (1, 1),
                    (50, 0),
                    (50, 50),
                    (50, 51),
                    (200, 200),
                ] {
                    let mut pb = PadicBinom::new(p, e);
                    let got = pb.s(big_n, t as i64, e);
                    let expected: u64 = if t > big_n {
                        (Integer::from(2).pow(big_n as u32) % Integer::from(modulus))
                            .to_string()
                            .parse()
                            .unwrap()
                    } else {
                        (s_k_exact(big_n, t) % Integer::from(modulus))
                            .to_string()
                            .parse()
                            .unwrap()
                    };
                    assert_eq!(got, expected, "p={p} e={e} N={big_n} t={t}");
                }
            }
        }
    }

    // -----------------------------------------------------------------------------
    // Partial-fraction split reassembles: contribution via split parts == contribution via
    // the unsplit modulus (doc §3 "obstacles 2/4").
    // -----------------------------------------------------------------------------

    #[test]
    fn partial_fraction_split_reassembles() {
        // m = p1^e1 * p2^e2 * p3, all coprime; compute frac(X*s/m) directly, then via the
        // partial-fraction split (each part computed mod its own prime power and combined the
        // way add_contribution does), and check they match to within the fixed-point ulp.
        let big_n = 40u64;
        let n = 1000u64;
        let base = 987_654_321u64; // arbitrary base so m_k = base + 2k is realistic-looking
        for &t in &[3u64, 10, 20] {
            let k = t; // unmirrored (k == t)
            let m = base + 2 * k;
            // Factor m by trial division (m is small enough here).
            let mut rem = m;
            let mut factors = Vec::new();
            let mut d = 2u64;
            while d * d <= rem {
                if rem.is_multiple_of(d) {
                    let mut e = 0;
                    while rem.is_multiple_of(d) {
                        rem /= d;
                        e += 1;
                    }
                    factors.push((d, e));
                }
                d += 1;
            }
            if rem > 1 {
                factors.push((rem, 1));
            }

            let s_exact: u64 = (s_k_exact(big_n, t) % Integer::from(m))
                .to_string()
                .parse()
                .unwrap();

            // Direct (unsplit) contribution: same formula as add_contribution but mod the
            // whole m (cofactor = 1, so the modular inverse step is a no-op).
            let direct = add_contribution(n, big_n, base, k, t, m, s_exact);

            // Split: for each prime power part q=p^e, s mod q via Lucas (e==1) or PadicBinom
            // (e>=2), reassembled via add_contribution's partial-fraction logic.
            let mut split_acc = 0u128;
            for &(p, e) in &factors {
                let q = p.pow(e);
                let s_part: u64 = (s_k_exact(big_n, t) % Integer::from(q))
                    .to_string()
                    .parse()
                    .unwrap();
                split_acc =
                    split_acc.wrapping_add(add_contribution(n, big_n, base, k, t, q, s_part));
            }
            // Each part rounds independently (one floor per add_contribution call), so the
            // split sum can differ from the single unsplit computation by a few ulps of the
            // u128 fixed-point encoding (this is exactly why error_units is charged one ulp
            // per *actual* rounded term, not one per k — see the module docs). Bound the
            // difference by the number of parts, not zero.
            let diff = direct.wrapping_sub(split_acc);
            let diff = diff.min(diff.wrapping_neg());
            assert!(
                diff <= factors.len() as u128 + 1,
                "t={t} m={m} factors={factors:?}: direct={direct} split={split_acc} diff={diff}"
            );
        }
    }

    // -----------------------------------------------------------------------------
    // Remainder-tree (art_chunk) results vs direct mod
    // -----------------------------------------------------------------------------

    #[test]
    fn art_chunk_matches_direct_mod() {
        let big_n = 60u64;
        let n = 2000u64;
        let base = 123_456_789u64;
        let lg_n = 64 - big_n.leading_zeros();
        // A handful of synthetic Main items at various targets and moduli.
        let items = vec![
            Item {
                t: 5,
                q: 97,
                tag: Tag::Main(5),
            },
            Item {
                t: 12,
                q: 9973,
                tag: Tag::Main(12),
            },
            Item {
                t: 12,
                q: 10_007,
                tag: Tag::Main(47),
            }, // k=N-1-t=47 mirrored
            Item {
                t: 30,
                q: 104_729,
                tag: Tag::Main(30),
            },
        ];
        let (acc, terms, lucas_out) = art_chunk(n, big_n, base, &items, lg_n);
        assert_eq!(terms, items.len() as u64);
        assert!(lucas_out.is_empty());

        let mut expected = 0u128;
        for it in &items {
            let s: u64 = (s_k_exact(big_n, it.t) % Integer::from(it.q))
                .to_string()
                .parse()
                .unwrap();
            let Tag::Main(k) = it.tag else { unreachable!() };
            expected = expected.wrapping_add(add_contribution(n, big_n, base, k, it.t, it.q, s));
        }
        assert_eq!(acc, expected);
    }

    // -----------------------------------------------------------------------------
    // Segmented-sieve factorisation vs trial division
    // -----------------------------------------------------------------------------

    #[test]
    fn factor_interval_matches_trial_division() {
        let (big_m, big_n) = (12u64, 200u64);
        let fac = factor_interval(big_m, big_n);
        let base = 2 * big_m * big_n + 1;
        for k in 0..big_n {
            let mut rem = base + 2 * k;
            let mut expected = Vec::new();
            let mut d = 3u64; // m_k is always odd
            while d * d <= rem {
                if rem.is_multiple_of(d) {
                    let mut e = 0u32;
                    while rem.is_multiple_of(d) {
                        rem /= d;
                        e += 1;
                    }
                    expected.push((d, e));
                }
                d += 2;
            }
            if rem > 1 {
                expected.push((rem, 1));
            }
            let mut got = fac[k as usize].clone();
            got.sort();
            expected.sort();
            assert_eq!(got, expected, "k={k}");
        }
    }

    // -----------------------------------------------------------------------------
    // Full frac_10n_pi_m vs Theorem 1 and vs MPFR, across mem_bits values
    // -----------------------------------------------------------------------------

    #[test]
    fn frac_10n_pi_m_matches_theorem1_and_mpfr() {
        use crate::digits_to_bits;
        use rug::{Float, float::Constant};

        for &n in &[2500u64, 5000, 20_000] {
            for &mem_bits in &[256u64, 1024, 8192] {
                let n0 = 20u32;
                let (thm2, terms) = frac_10n_pi_m_with_terms(n, n0, mem_bits);
                let thm1 = nthdigit::frac_10n_pi(n, n0);
                let err2 = nthdigit::error_units(n0, terms);
                let diff = thm2.wrapping_sub(thm1);
                let diff = diff.min(diff.wrapping_neg());
                assert!(
                    diff <= 2 * err2,
                    "n={n} mem_bits={mem_bits}: thm1/thm2 disagree by {diff} > 2*{err2}"
                );

                // Cross-check against MPFR directly too.
                let bits = digits_to_bits((n + 40) as u32) + 8;
                let pi = Float::with_val(bits, Constant::Pi);
                let scale = Float::with_val(bits, Integer::from(10).pow((n + 20) as u32));
                let scaled = Float::with_val(bits, &pi * &scale);
                let int_part = scaled.to_integer().unwrap();
                let s = int_part.to_string();
                // s = "3" followed by (n+20) decimal digits, so s[i] is the digit at position
                // i (1-based); we want positions n+1..=n+10.
                let start = (n + 1) as usize;
                let digits10: String = s[start..start + 10].to_string();
                let got10 = extract_digits(thm2, 10);
                assert_eq!(got10, digits10, "n={n} mem_bits={mem_bits}: vs MPFR");
            }
        }
    }

    #[test]
    fn digits_matches_theorem1_across_mem_bits() {
        for &n in &[3000u64, 15_000, 50_000] {
            for &mem_bits in &[512u64, default_mem_bits(n), 16384] {
                let count = 10;
                let got1 = nthdigit::digits(n, count);
                let got2 = digits(n, count, mem_bits);
                assert_eq!(got1, got2, "n={n} mem_bits={mem_bits}");
            }
        }
    }

    // -----------------------------------------------------------------------------
    // Streamed item generation vs the old global (reference) construction: must produce
    // exactly the same multiset of (target, modulus, tag) items, per the task's correctness
    // gate. `factor_interval` above is kept #[cfg(test)] specifically to make this comparison
    // possible.
    // -----------------------------------------------------------------------------

    /// `(Main items as (t,q,k), lucas_need, padic_need)`, the shape both
    /// [`reference_construction`] and [`streamed_construction`] return for comparison.
    type ConstructionResult = (
        Vec<(u64, u64, u64)>,
        HashMap<(u64, u64), Vec<(u64, u64)>>,
        HashMap<u64, Vec<(u64, u64, u32)>>,
    );

    /// The old (global, non-streaming) item/need construction, i.e. exactly what `c_part` did
    /// before streaming: every `k`'s complete factorisation via `factor_interval`, classified
    /// into Main items / Lucas needs / p-adic needs in one pass with no windowing.
    fn reference_construction(big_m: u64, big_n: u64) -> ConstructionResult {
        let fac = factor_interval(big_m, big_n);
        let mut main_items = Vec::new();
        let mut lucas_need: HashMap<(u64, u64), Vec<(u64, u64)>> = HashMap::new();
        let mut padic_need: HashMap<u64, Vec<(u64, u64, u32)>> = HashMap::new();
        for k in 0..big_n {
            let t = k.min(big_n - 1 - k);
            let mut good = 1u64;
            for &(pp, e) in &fac[k as usize] {
                if pp > t {
                    good *= pp.pow(e);
                } else if e == 1 {
                    lucas_need.entry((pp, t % pp)).or_default().push((k, t));
                } else {
                    padic_need.entry(pp).or_default().push((k, t, e));
                }
            }
            if good > 1 {
                main_items.push((t, good, k));
            }
        }
        (main_items, lucas_need, padic_need)
    }

    /// Expands a [`SideRanges`] map back into the explicit `(p,r) -> [(k,t)]` shape the
    /// pre-streaming reference construction used (test-only: production code never
    /// materialises this, that's the whole point of `SideRanges` — see its docs).
    fn materialize_lucas(
        need: &HashMap<u64, SideRanges>,
        big_n: u64,
        base: u64,
    ) -> HashMap<(u64, u64), Vec<(u64, u64)>> {
        let mut out: HashMap<(u64, u64), Vec<(u64, u64)>> = HashMap::new();
        for (&p, ranges) in need {
            for (is_high, range) in [(false, ranges.lo), (true, ranges.hi)] {
                let Some((t_min, t_max)) = range else {
                    continue;
                };
                let r = t_min % p;
                let mut t = t_min;
                loop {
                    let k = if is_high { big_n - 1 - t } else { t };
                    if !(base + 2 * k).is_multiple_of(p * p) {
                        out.entry((p, r)).or_default().push((k, t));
                    }
                    if t == t_max {
                        break;
                    }
                    t += p;
                }
            }
        }
        out
    }

    /// Same idea for `padic_small`, recomputing each member's exact `e` the same way
    /// [`padic_consumers`] does (step `p^2`, not `p`).
    fn materialize_padic(
        need: &HashMap<u64, SideRanges>,
        big_n: u64,
        base: u64,
    ) -> HashMap<u64, Vec<(u64, u64, u32)>> {
        let mut out: HashMap<u64, Vec<(u64, u64, u32)>> = HashMap::new();
        for (&p, ranges) in need {
            let step = p * p;
            for (is_high, range) in [(false, ranges.lo), (true, ranges.hi)] {
                let Some((t_min, t_max)) = range else {
                    continue;
                };
                let mut t = t_min;
                loop {
                    let k = if is_high { big_n - 1 - t } else { t };
                    let e = padic_valuation(base, k, p);
                    out.entry(p).or_default().push((k, t, e));
                    if t == t_max {
                        break;
                    }
                    t += step;
                }
            }
        }
        out
    }

    /// The new streaming construction, gathering every chunk's Main items (regenerated exactly
    /// as [`art_chunk_by_range`] would) plus the small-prime and (regrouped) cofactor Lucas
    /// needs and the p-adic needs, into the same shape [`reference_construction`] returns (via
    /// [`materialize_lucas`]/[`materialize_padic`]), so the two can be compared directly.
    fn streamed_construction(big_m: u64, big_n: u64, mem_bits: u64) -> ConstructionResult {
        let base = 2 * big_m * big_n + 1;
        let m_max = base + 2 * (big_n - 1);
        let small_primes = primes_upto((m_max as f64).sqrt() as u64 + 2);
        let needs = stream_needs_and_chunks(base, big_n, &small_primes, mem_bits);

        let mut main_items = Vec::new();
        for &(t0, t1) in &needs.chunk_bounds {
            let len = t1 - t0;
            let (rem_lo, fac_lo) = factor_window(base, t0, len, &small_primes);
            let (rem_hi, fac_hi) = factor_window(base, big_n - t1, len, &small_primes);
            for t in t0..t1 {
                let i_lo = (t - t0) as usize;
                let i_hi = (t1 - 1 - t) as usize;
                if let Some(q) = main_modulus(t, rem_lo[i_lo], &fac_lo[i_lo]) {
                    main_items.push((t, q, t));
                }
                if let Some(q) = main_modulus(t, rem_hi[i_hi], &fac_hi[i_hi]) {
                    main_items.push((t, q, big_n - 1 - t));
                }
            }
        }

        let mut combined_lucas = materialize_lucas(&needs.lucas_small, big_n, base);
        for (p, k, t) in needs.cofactor {
            combined_lucas.entry((p, t % p)).or_default().push((k, t));
        }
        let padic = materialize_padic(&needs.padic_small, big_n, base);
        (main_items, combined_lucas, padic)
    }

    #[test]
    fn streamed_items_match_reference_construction() {
        for &(big_m, big_n, mem_bits) in &[
            (4u64, 60u64, 128u64),
            (6, 200, 256),
            (12, 734, 512),
            (30, 4000, 1024),
            (2, 52, 64), // small case
        ] {
            let (ref_main, ref_lucas, ref_padic) = reference_construction(big_m, big_n);
            let (str_main, str_lucas, str_padic) = streamed_construction(big_m, big_n, mem_bits);

            let mut a = ref_main.clone();
            a.sort();
            let mut b = str_main.clone();
            b.sort();
            assert_eq!(a, b, "Main (t,q,k) multiset differs: M={big_m} N={big_n}");

            let mut la: Vec<_> = ref_lucas.keys().copied().collect();
            la.sort();
            let mut lb: Vec<_> = str_lucas.keys().copied().collect();
            lb.sort();
            assert_eq!(la, lb, "Lucas (p,r) key set differs: M={big_m} N={big_n}");
            for key in &la {
                let mut rv = ref_lucas[key].clone();
                rv.sort();
                let mut sv = str_lucas[key].clone();
                sv.sort();
                assert_eq!(
                    rv, sv,
                    "Lucas need-list differs for key={key:?}: M={big_m} N={big_n}"
                );
            }

            let mut pa: Vec<_> = ref_padic.keys().copied().collect();
            pa.sort();
            let mut pb: Vec<_> = str_padic.keys().copied().collect();
            pb.sort();
            assert_eq!(pa, pb, "p-adic prime set differs: M={big_m} N={big_n}");
            for p in &pa {
                let mut rv = ref_padic[p].clone();
                rv.sort();
                let mut sv = str_padic[p].clone();
                sv.sort();
                assert_eq!(
                    rv, sv,
                    "p-adic need-list differs for p={p}: M={big_m} N={big_n}"
                );
            }

            // Term count (every (k,t) pair, across Main + Lucas + p-adic, contributes exactly
            // one add_contribution call) must match too.
            let ref_terms = ref_main.len()
                + ref_lucas.values().map(Vec::len).sum::<usize>()
                + ref_padic.values().map(Vec::len).sum::<usize>();
            let str_terms = str_main.len()
                + str_lucas.values().map(Vec::len).sum::<usize>()
                + str_padic.values().map(Vec::len).sum::<usize>();
            assert_eq!(
                ref_terms, str_terms,
                "term count differs: M={big_m} N={big_n}"
            );
        }
    }

    #[test]
    fn chunk_bounds_cover_target_range_exactly_once() {
        for &(big_m, big_n, mem_bits) in &[
            (4u64, 60u64, 128u64),
            (6, 200, 256),
            (12, 734, 512),
            (30, 4000, 1024),
            (96, 36_806, 4202), // n=1e5-scale M,N,mem_bits
        ] {
            let base = 2 * big_m * big_n + 1;
            let m_max = base + 2 * (big_n - 1);
            let small_primes = primes_upto((m_max as f64).sqrt() as u64 + 2);
            let needs = stream_needs_and_chunks(base, big_n, &small_primes, mem_bits);
            let half = big_n / 2;
            let mut expected_start = 0u64;
            for &(t0, t1) in &needs.chunk_bounds {
                assert_eq!(
                    t0, expected_start,
                    "gap or overlap before chunk starting at {t0}: M={big_m} N={big_n}"
                );
                assert!(
                    t1 > t0,
                    "empty chunk window ({t0},{t1}): M={big_m} N={big_n}"
                );
                expected_start = t1;
            }
            assert_eq!(
                expected_start, half,
                "chunk windows don't cover [0, N/2) exactly: M={big_m} N={big_n}"
            );
        }
    }
}
