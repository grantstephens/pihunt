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
//! memory is that times the (small, fixed) thread count — **with the default
//! [`USE_SIEVED_FACTORING`] path**, `c_sum` additionally holds one process-wide list of
//! primes `<= sqrt(2MN)` (`O(pi(sqrt(2MN)))` words, a few KiB through n=10^6, an estimated
//! low single-digit MiB by n=10^9) to replace `O(k)` trial division per `k` with a segmented
//! sieve; see `docs/nthdigit.md` for why (short version: trial division was never actually
//! the dominant cost, so this alone was a small win — Montgomery multiplication in
//! [`s_k_mod_with_primes`] is what the real speedup came from) and for the strict-memory
//! fallback that keeps the bound above exactly.
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

/// Montgomery multiplication context for a fixed *odd* modulus `m`. Every modulus in this
/// module is odd (`m_k = 2MN + 2k + 1` and `2k + 1` in `b_sum` are both of the form `2x+1`),
/// so this applies everywhere `mulmod`'s repeated `u128` divisions show up as the hot loop's
/// bottleneck — see [`s_k_mod_with_primes`], by far the dominant cost (profiling: the
/// binomial loop is ~25x the cost of factoring itself, sieved or not — see
/// `docs/nthdigit.md`). Montgomery reduction ([`Self::redc`]) replaces each `u128` division
/// with a handful of multiplications and a shift, at the cost of keeping values in
/// "Montgomery form" (`x*2^64 mod m`) for the duration, converting in/out at the boundary.
#[derive(Clone, Copy)]
struct Montgomery {
    m: u64,
    /// `-m^-1 mod 2^64`.
    m_inv_neg: u64,
    /// `2^64 mod m`: the Montgomery form of `1`.
    r_mod_m: u64,
    /// `(2^64)^2 mod m`: multiplying a plain residue by this (then reducing) converts it to
    /// Montgomery form.
    r2_mod_m: u64,
}

impl Montgomery {
    /// `m` must be odd (so it's invertible mod `2^64`) — checked with a `debug_assert`, not
    /// enforced at runtime, since every caller in this module already guarantees it structurally.
    fn new(m: u64) -> Self {
        debug_assert!(
            !m.is_multiple_of(2),
            "Montgomery reduction needs an odd modulus"
        );
        // `redc`'s intermediate `t + mprime*m` can overflow `u128` once `m` gets within a
        // small constant factor of `2^64` (worst case around `m ≈ 1.24 * 2^63`); every real
        // `m` in this module (`2*M*N + 2k + 1`, `2k + 1`) stays many orders of magnitude
        // below this for any `n` the crate could plausibly run, so this is a documented
        // assumption, not a limitation that bites in practice.
        debug_assert!(
            m < (1u64 << 62),
            "Montgomery reduction here assumes m < 2^62"
        );
        // Newton's method for `m^-1 mod 2^64`: an odd `m` is its own inverse mod 8 (3 correct
        // bits), and each iteration below doubles the number of correct bits (3, 6, 12, 24,
        // 48, 96), so 5 iterations comfortably clears 64.
        let mut inv = m;
        for _ in 0..5 {
            inv = inv.wrapping_mul(2u64.wrapping_sub(m.wrapping_mul(inv)));
        }
        debug_assert_eq!(
            m.wrapping_mul(inv),
            1,
            "inv must be the exact inverse of m mod 2^64"
        );
        let m_inv_neg = inv.wrapping_neg();
        let r_mod_m = ((1u128 << 64) % m as u128) as u64;
        let r2_mod_m = mulmod(r_mod_m, r_mod_m, m);
        Montgomery {
            m,
            m_inv_neg,
            r_mod_m,
            r2_mod_m,
        }
    }

    /// `t * 2^-64 mod m`, for `t < m * 2^64` (always true here: every `t` passed in is a
    /// product of two values `< m`, or `< m` outright, so `t < m^2 <= m * 2^64` trivially
    /// for the `m` this module ever sees). The low 64 bits of `t + (t*m_inv_neg mod 2^64)*m`
    /// are exactly zero by construction of `m_inv_neg`, so the shift below is an exact
    /// division, and the result needs at most one conditional subtraction to land in `[0,
    /// m)` (it's `< 2m` beforehand).
    #[inline]
    fn redc(&self, t: u128) -> u64 {
        let mprime = (t as u64).wrapping_mul(self.m_inv_neg);
        let sum = t + (mprime as u128) * (self.m as u128);
        let mut result = (sum >> 64) as u64;
        if result >= self.m {
            result -= self.m;
        }
        result
    }

    /// The Montgomery form of `1` (i.e. `2^64 mod m`) — the correct initial accumulator
    /// value for a running Montgomery-form product (in place of plain `1 % m`).
    #[inline]
    fn one(&self) -> u64 {
        self.r_mod_m
    }

    /// Converts a plain residue `a mod m` (`a < m`) into Montgomery form (`a * 2^64 mod m`).
    #[inline]
    fn encode(&self, a: u64) -> u64 {
        self.redc(a as u128 * self.r2_mod_m as u128)
    }

    /// Converts a Montgomery-form value back to a plain residue.
    #[inline]
    fn decode(&self, a: u64) -> u64 {
        self.redc(a as u128)
    }

    /// Multiplies two Montgomery-form values, giving a Montgomery-form result (i.e. this is
    /// `mulmod` on the underlying plain values, without ever leaving Montgomery form).
    #[inline]
    fn mul(&self, a: u64, b: u64) -> u64 {
        self.redc(a as u128 * b as u128)
    }
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
/// early once `m`'s cofactor is fully spent). `O(k)` time, `O(log m)` memory. This is the
/// strict-memory fallback (see [`USE_SIEVED_FACTORING`]): correct and `O(log^2 n)` memory
/// on the nose, but the dominant cost of the whole C-sum once `k` gets into the hundreds of
/// thousands (measured: ~100 s at n=10^6, ~2 h at n=10^7 on 6 cores). [`factor_segment`]
/// below replaces this for the hot path.
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
/// Finds the primes `<= k` dividing `m` by trial division (`O(k)`, the strict-memory path;
/// see [`USE_SIEVED_FACTORING`]), then defers to [`s_k_mod_with_primes`] for the actual
/// Algorithm 2 loop.
fn s_k_mod_direct(big_n: u64, k: u64, m: u64) -> u64 {
    if m == 1 {
        return 0;
    }
    #[cfg(feature = "nthdigit-profile")]
    let t0 = std::time::Instant::now();
    let primes = primes_le_k_dividing(m, k);
    #[cfg(feature = "nthdigit-profile")]
    {
        FACTOR_NS.fetch_add(
            t0.elapsed().as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
    }
    #[cfg(feature = "nthdigit-profile")]
    let t1 = std::time::Instant::now();
    let result = s_k_mod_with_primes(big_n, k, m, &primes);
    #[cfg(feature = "nthdigit-profile")]
    {
        LOOP_NS.fetch_add(
            t1.elapsed().as_nanos() as u64,
            std::sync::atomic::Ordering::Relaxed,
        );
    }
    result
}

/// The actual Algorithm 2 loop: `Σ_{j=0}^{k} C(N,j) mod m`, given the (distinct) primes
/// `<= k` dividing `m` (however they were found — trial division or the segmented sieve).
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
fn s_k_mod_with_primes(big_n: u64, k: u64, m: u64, primes: &[u64]) -> u64 {
    // `m` is always odd here (every caller's `m` is `2*something + 1`), so Montgomery
    // reduction applies; see [`Montgomery`] for why this loop is where it matters (measured
    // as ~96% of the whole computation's cost, dwarfing prime-factoring by ~25x whichever
    // way the factors are found).
    let mont = Montgomery::new(m);
    let mut r_vals = vec![1u64; primes.len()];
    // Montgomery form of each `r_vals[idx]`, cached and only refreshed on the (relatively
    // rare — most `j` share no factor with most `p_i`) iterations where `r_vals[idx]`
    // actually changes: re-`encode`-ing every `r_val` on every `j` regardless, measured,
    // cost more than the division-based `mulmod` it was replacing (every extra `encode` is
    // itself a `redc`), which is exactly the kind of "looks free, isn't" trap this whole
    // exercise is about avoiding — see `docs/nthdigit.md`.
    let mut r_vals_mont = vec![mont.one(); primes.len()];
    let (mut a_acc, mut b_acc, mut c_acc) = (mont.one(), mont.one(), mont.one());

    for j in 1..=k {
        let a = big_n - j + 1;
        let mut a_star = a;
        let mut b_star = j;
        for (idx, &p) in primes.iter().enumerate() {
            let alpha = strip_power(&mut a_star, p);
            let beta = strip_power(&mut b_star, p);
            if alpha > 0 || beta > 0 {
                for _ in 0..alpha {
                    r_vals[idx] *= p;
                }
                for _ in 0..beta {
                    r_vals[idx] /= p;
                }
                r_vals_mont[idx] = mont.encode(r_vals[idx] % m);
            }
        }
        debug_assert!(
            primes.iter().enumerate().all(|(i, _)| r_vals[i] <= big_n),
            "R_i must stay <= N (Gourdon lemma 1)"
        );

        let a_star_m = mont.encode(a_star % m);
        let b_star_m = mont.encode(b_star % m);
        a_acc = mont.mul(a_acc, a_star_m);
        b_acc = mont.mul(b_acc, b_star_m);

        let mut r_prod = mont.one();
        for &r in &r_vals_mont {
            r_prod = mont.mul(r_prod, r);
        }
        let term1 = mont.mul(c_acc, b_star_m);
        let term2 = mont.mul(a_acc, r_prod);
        c_acc = (term1 + term2) % m;
    }

    let b_inv = mod_inverse(mont.decode(b_acc), m);
    mulmod(mont.decode(c_acc), b_inv, m)
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

/// Same as [`sum_binomials_mod`], but given `m`'s *full* distinct prime factorisation
/// (every prime factor of `m`, not just those `<= k`) instead of computing it internally.
/// Used by the segmented-sieve path: [`factor_segment`] finds all prime factors of `m` (up
/// to and including a possible large cofactor) once per segment, and this filters that list
/// down to `<= k` (or `<= comp_k` under the complement trick — the same factorisation of
/// `m` serves either threshold, since it doesn't depend on `k` at all).
fn sum_binomials_mod_sieved(big_n: u64, k: u64, m: u64, all_factors: &[u64]) -> u64 {
    if m == 1 {
        return 0;
    }
    if 2 * k > big_n {
        let comp_k = big_n - k - 1;
        let primes: Vec<u64> = all_factors
            .iter()
            .copied()
            .filter(|&p| p <= comp_k)
            .collect();
        let s_comp = s_k_mod_with_primes(big_n, comp_k, m, &primes);
        let two_pow_n = powmod(2, big_n, m);
        (two_pow_n + m - s_comp % m) % m
    } else {
        let primes: Vec<u64> = all_factors.iter().copied().filter(|&p| p <= k).collect();
        s_k_mod_with_primes(big_n, k, m, &primes)
    }
}

// ---------------------------------------------------------------------------------
// Segmented sieve over the arithmetic progression m_k = c0 + 2k (Gourdon §3.1-3.2 leaves
// the choice of factoring method open; this is the practical replacement for
// `primes_le_k_dividing`'s O(k) trial division in the hot C-sum loop).
// ---------------------------------------------------------------------------------

/// Selects the factoring strategy [`c_sum`] uses to find, for each `k`, the primes `<= k`
/// dividing `m_k = 2MN + 2k + 1`:
///
/// - `true` (default): the segmented sieve below. Needs an extra `O(pi(sqrt(2MN)))` words
///   for the small-prime list (a few thousand to a few hundred thousand `u64`s across the
///   whole n=10^4..10^9 range we care about — see "Memory bound" in `docs/nthdigit.md`),
///   plus small per-segment scratch. In exchange it turns the dominant cost of the whole
///   computation (O(k) trial division per k, i.e. O(N^2) over the C-sum) into an O(N log
///   log N)-ish sieve pass, which is why n=10^6 goes from ~100 s to low single-digit
///   seconds (see the benchmark table).
/// - `false`: falls back to `primes_le_k_dividing`'s O(k) trial division per k (the original
///   v1 path), giving back the strict `O(log^2 n)` per-thread memory bound Gourdon's
///   theorem promises, at the cost of the 60x-100x slowdown this change fixes. Flip this if
///   the sieve's extra memory is ever actually a problem (it isn't at any n we've run: peak
///   RSS stays in the single-digit MiB through 10^7).
const USE_SIEVED_FACTORING: bool = true;

/// Target number of segments *per thread* for [`c_sum_sieved`]'s rayon split. Algorithm 2's
/// per-`k` cost is `O(k)`, i.e. wildly skewed across a segment range (the last segment costs
/// ~`N` times more than the first), so we deliberately oversubscribe well past the thread
/// count and let rayon's work-stealing even out the imbalance; 32x was enough in practice
/// (see `docs/nthdigit.md` — an earlier fixed `SEGMENT_SIZE = 2^15` regressed n=10^4..10^6
/// by 4-5x because it collapsed to only 1-2 segments total at those `N`, i.e. no
/// parallelism at all, before the per-`k` sieve savings had a chance to matter).
const SEGMENTS_PER_THREAD: u64 = 32;

/// Floor on segment length, so the fixed per-segment cost of scanning `small_primes` (an
/// `O(1)` modular check per prime, to find that prime's first hit in the segment) doesn't
/// dominate at very small `N` (where `big_n / (threads * SEGMENTS_PER_THREAD)` would
/// otherwise round down to single digits or zero).
const MIN_SEGMENT_SIZE: u64 = 64;

/// Segment length for [`c_sum_sieved`]: `big_n` split into roughly
/// `available_threads * SEGMENTS_PER_THREAD` pieces (never shorter than
/// [`MIN_SEGMENT_SIZE`]). Scales with `N`, not fixed, precisely so it keeps giving rayon
/// enough independent tasks to load-balance regardless of how big the run is.
fn segment_size_for(big_n: u64) -> u64 {
    let threads = rayon::current_num_threads().max(1) as u64;
    (big_n / (threads * SEGMENTS_PER_THREAD)).max(MIN_SEGMENT_SIZE)
}

/// Integer square root of `x` (largest `r` with `r*r <= x`), via `f64::sqrt` plus an exact
/// correction (needed since `f64` can't represent every `u64` exactly).
fn isqrt(x: u64) -> u64 {
    if x == 0 {
        return 0;
    }
    let mut r = (x as f64).sqrt() as u64;
    while r > 0 && r.checked_mul(r).is_none_or(|sq| sq > x) {
        r -= 1;
    }
    while (r + 1).checked_mul(r + 1).is_some_and(|sq| sq <= x) {
        r += 1;
    }
    r
}

/// All primes `<= bound`, via a plain sieve of Eratosthenes. `O(bound)` time and space,
/// called once per [`c_sum_sieved`] invocation (not per segment, not per k) — `bound` here
/// is `sqrt(2MN)`, not `N` itself, so this stays cheap even though it's the one place in
/// this module whose memory does scale with `n` (see [`USE_SIEVED_FACTORING`]).
fn primes_up_to(bound: u64) -> Vec<u64> {
    if bound < 2 {
        return Vec::new();
    }
    let bound = bound as usize;
    let mut is_composite = vec![false; bound + 1];
    let mut p = 2usize;
    while p * p <= bound {
        if !is_composite[p] {
            let mut j = p * p;
            while j <= bound {
                is_composite[j] = true;
                j += p;
            }
        }
        p += 1;
    }
    (2..=bound)
        .filter(|&i| !is_composite[i])
        .map(|i| i as u64)
        .collect()
}

/// For each `k` in `[k0, k0+len)`, the complete set of distinct prime factors of
/// `m_k = c0 + 2*k` (which is always odd), found by:
///
/// 1. Dividing every `m_k` in the segment by every prime in `small_primes` that hits it
///    (located via modular arithmetic: `m_k ≡ 0 (mod p)` iff `k ≡ -c0 * inv2(p) (mod p)`,
///    and `inv2(p) = (p+1)/2` for odd `p`, so no extended-gcd is needed).
/// 2. Whatever's left of each `m_k` after that is `1` or a single prime `> sqrt(max m_k)`
///    (that's what `small_primes` covering everything up to `sqrt(max m_k)` guarantees) —
///    included as-is if it's `> 1`.
///
/// `small_primes` must contain every odd prime `<= sqrt(c0 + 2*(k0+len-1))`; it's the
/// caller's job to size it for the whole k-range once (see [`c_sum_sieved`]), not just this
/// segment, so it can be shared read-only across every segment/thread.
///
/// Callers filter the returned lists down to primes `<= k` (or `<= comp_k`) themselves —
/// see [`sum_binomials_mod_sieved`] — since Algorithm 2 only cares about primes that small,
/// but the same factorisation serves both the direct and complement thresholds.
fn factor_segment(c0: u64, k0: u64, len: u64, small_primes: &[u64]) -> Vec<Vec<u64>> {
    let len_usize = len as usize;
    let mut cofactor: Vec<u64> = (0..len).map(|i| c0 + 2 * (k0 + i)).collect();
    let mut factors: Vec<Vec<u64>> = vec![Vec::new(); len_usize];

    for &p in small_primes {
        if p == 2 {
            // m_k = c0 + 2k is always odd (c0 = 2MN+1 is odd); 2 never divides it.
            continue;
        }
        let inv2 = p.div_ceil(2);
        // k (absolute) with (c0 + 2k) ≡ 0 (mod p)  <=>  k ≡ -c0 * inv2 (mod p).
        let target_k_mod_p = ((p - c0 % p) % p) * inv2 % p;
        // Segment-relative index i (k0+i ≡ target_k_mod_p mod p) of the first hit.
        let mut i = ((target_k_mod_p + p - k0 % p) % p) as usize;
        while i < len_usize {
            debug_assert!(cofactor[i].is_multiple_of(p));
            while cofactor[i].is_multiple_of(p) {
                cofactor[i] /= p;
            }
            factors[i].push(p);
            i += p as usize;
        }
    }

    for i in 0..len_usize {
        if cofactor[i] > 1 {
            factors[i].push(cofactor[i]);
        }
    }
    factors
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
    if USE_SIEVED_FACTORING {
        c_sum_sieved(n, p)
    } else {
        c_sum_trial_division(n, p)
    }
}

/// The strict-`O(log^2 n)`-memory C-sum: one rayon task per `k`, `O(k)` trial division to
/// find `m_k`'s prime factors `<= k` (see [`USE_SIEVED_FACTORING`]).
fn c_sum_trial_division(n: u64, p: Params) -> u128 {
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

/// The segmented-sieve C-sum: `k in 0..N` is chunked into [`segment_size_for`]-sized ranges,
/// one rayon task per chunk. Each task sieves its own segment's `m_k` values against a
/// once-per-call list of primes `<= sqrt(max m_k)` (via [`factor_segment`]), then runs
/// Algorithm 2's `O(k)` binomial loop for every `k` in the segment using the factors just
/// found (via [`sum_binomials_mod_sieved`]) instead of re-deriving them by trial division.
fn c_sum_sieved(n: u64, p: Params) -> u128 {
    use rayon::prelude::*;
    let (big_n, big_m) = (p.big_n, p.m);
    let c0 = 2 * big_m * big_n + 1; // m_k = c0 + 2k, always odd
    let m_max = c0 + 2 * (big_n - 1);
    let small_primes = primes_up_to(isqrt(m_max));
    let segment_size = segment_size_for(big_n);

    let segment_starts: Vec<u64> = (0..big_n).step_by(segment_size as usize).collect();
    segment_starts
        .into_par_iter()
        .map(|k0| {
            let len = segment_size.min(big_n - k0);
            #[cfg(feature = "nthdigit-profile")]
            let t0 = std::time::Instant::now();
            let factors = factor_segment(c0, k0, len, &small_primes);
            #[cfg(feature = "nthdigit-profile")]
            {
                FACTOR_NS.fetch_add(
                    t0.elapsed().as_nanos() as u64,
                    std::sync::atomic::Ordering::Relaxed,
                );
            }
            #[cfg(feature = "nthdigit-profile")]
            let t1 = std::time::Instant::now();
            let e1 = big_n - 2;
            let e2 = n - big_n + 2;
            let mut partial = 0u128;
            for (i, prime_factors) in factors.iter().enumerate() {
                let k = k0 + i as u64;
                let m = c0 + 2 * k;
                let s = sum_binomials_mod_sieved(big_n, k, m, prime_factors);
                let pow5 = powmod(5, e1, m);
                let pow10 = powmod(10, e2, m);
                let y = mulmod(mulmod(pow5, pow10, m), s, m);
                partial = partial.wrapping_add(signed(frac_fixed_point(y, m), k));
            }
            #[cfg(feature = "nthdigit-profile")]
            {
                LOOP_NS.fetch_add(
                    t1.elapsed().as_nanos() as u64,
                    std::sync::atomic::Ordering::Relaxed,
                );
            }
            partial
        })
        .reduce(|| 0u128, u128::wrapping_add)
}

#[cfg(feature = "nthdigit-profile")]
static FACTOR_NS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);
#[cfg(feature = "nthdigit-profile")]
static LOOP_NS: std::sync::atomic::AtomicU64 = std::sync::atomic::AtomicU64::new(0);

/// Throwaway profiling hook (see `nthdigit-profile` feature): total time spent inside
/// [`factor_segment`] and total time spent in the per-segment binomial loop, summed across
/// every thread, since the process started (nanoseconds).
#[cfg(feature = "nthdigit-profile")]
pub fn profile_totals_ns() -> (u64, u64) {
    (
        FACTOR_NS.load(std::sync::atomic::Ordering::Relaxed),
        LOOP_NS.load(std::sync::atomic::Ordering::Relaxed),
    )
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

    /// A tiny xorshift64 PRNG (deterministic, no external `rand` dependency), mirroring the
    /// one in `tests/nthdigit.rs`.
    struct Xorshift(u64);
    impl Xorshift {
        fn next(&mut self) -> u64 {
            self.0 ^= self.0 << 13;
            self.0 ^= self.0 >> 7;
            self.0 ^= self.0 << 17;
            self.0
        }
    }

    /// [`Montgomery`]'s multiplication must agree with plain [`mulmod`] for random odd `m`
    /// and random `a, b < m` — the correctness oracle for the hot-loop rewrite in
    /// [`s_k_mod_with_primes`]. Also checks `encode`/`decode` round-trip and that `one()`
    /// really is the Montgomery form of `1`.
    #[test]
    fn montgomery_mul_matches_mulmod() {
        let mut rng = Xorshift(0x9E37_79B9_7F4A_7C15);
        // A handful of small, fixed edge cases first (m=1 is deliberately excluded: every
        // caller guards `m == 1` before ever constructing a `Montgomery`).
        for &m in &[3u64, 5, 7, 9, 2 * 3 * 5 * 7 + 1, (1u64 << 61) - 1] {
            check_montgomery(m, 0, 0);
            check_montgomery(m, 1 % m, 1 % m);
            check_montgomery(m, m - 1, m - 1);
        }
        for _ in 0..5_000 {
            // Kept `< 2^61` (see `Montgomery::new`'s overflow note) and odd.
            let mut m = (rng.next() % (1u64 << 61)) | 1;
            if m == 1 {
                m = 3;
            }
            let a = rng.next() % m;
            let b = rng.next() % m;
            check_montgomery(m, a, b);
        }
    }

    fn check_montgomery(m: u64, a: u64, b: u64) {
        let mont = Montgomery::new(m);
        assert_eq!(
            mont.decode(mont.one()),
            1 % m,
            "m={m}: one() must be Montgomery-form 1"
        );
        assert_eq!(
            mont.decode(mont.encode(a)),
            a,
            "m={m} a={a}: encode/decode round-trip"
        );
        let expected = mulmod(a, b, m);
        let got = mont.decode(mont.mul(mont.encode(a), mont.encode(b)));
        assert_eq!(
            got, expected,
            "m={m} a={a} b={b}: Montgomery mul disagrees with mulmod"
        );
    }

    /// The full distinct-prime factorisation of `m`, by trial division against every
    /// integer up to `m` itself (i.e. [`primes_le_k_dividing`] with no `k` cutoff) — an
    /// independent (if slower) way to factor `m`, used below as ground truth for
    /// [`factor_segment`].
    fn full_factorisation(m: u64) -> Vec<u64> {
        let mut v = primes_le_k_dividing(m, m);
        v.sort_unstable();
        v
    }

    /// [`factor_segment`]'s output (a full factorisation per `k`, before any `<= k`
    /// filtering) must agree with trial division for every `k` in `0..big_n`, across several
    /// `(M, N)` pairs and several segment sizes (including ones that don't evenly divide
    /// `big_n`, to exercise the last-segment-is-short case).
    #[test]
    fn sieved_factorisation_matches_trial_division() {
        let cases: &[(u64, u64)] = &[(4, 20), (10, 50), (6, 200), (2, 733)];
        for &(big_m, big_n) in cases {
            let c0 = 2 * big_m * big_n + 1;
            let m_max = c0 + 2 * (big_n - 1);
            let small_primes = primes_up_to(isqrt(m_max));
            for &seg_len in &[big_n, 7, 32] {
                let mut k0 = 0u64;
                while k0 < big_n {
                    let len = seg_len.min(big_n - k0);
                    let factors = factor_segment(c0, k0, len, &small_primes);
                    assert_eq!(factors.len(), len as usize);
                    for (i, sieved) in factors.iter().enumerate() {
                        let k = k0 + i as u64;
                        let m = c0 + 2 * k;
                        let mut got = sieved.clone();
                        got.sort_unstable();
                        assert_eq!(
                            got,
                            full_factorisation(m),
                            "M={big_m} N={big_n} k={k} m={m} seg_len={seg_len}"
                        );
                    }
                    k0 += len;
                }
            }
        }
    }

    /// A hand-picked `m` with a repeated small prime factor (`3^2`) and a large prime
    /// cofactor (`9973`, well above `sqrt(m) ≈ 669.9`) that the sieve must recover as a
    /// single prime even though it never appears in `small_primes`.
    #[test]
    fn sieved_factorisation_handles_repeated_factor_and_large_cofactor() {
        let m: u64 = 3 * 3 * 5 * 9973; // = 448785, odd, sqrt ≈ 669.9
        assert!(!m.is_multiple_of(2), "m must be odd");
        let small_primes = primes_up_to(isqrt(m));
        assert!(small_primes.contains(&5));
        assert!(
            !small_primes.contains(&9973),
            "9973 should be the leftover cofactor, not sieved"
        );
        let factors = factor_segment(m, 0, 1, &small_primes);
        let mut got = factors[0].clone();
        got.sort_unstable();
        assert_eq!(got, vec![3, 5, 9973]);
    }

    /// `sum_binomials_mod` (trial division) and the sieved path must agree end-to-end for
    /// every `k < N`, given the same `(M, N)` — i.e. the sieve doesn't just factor correctly
    /// in isolation, it also feeds Algorithm 2 the right `<= k` subset (both directly and
    /// under the `k > N/2` complement trick).
    #[test]
    fn sieved_algorithm2_matches_trial_division_algorithm2() {
        let cases: &[(u64, u64)] = &[(4, 20), (10, 50), (6, 200)];
        for &(big_m, big_n) in cases {
            let c0 = 2 * big_m * big_n + 1;
            let m_max = c0 + 2 * (big_n - 1);
            let small_primes = primes_up_to(isqrt(m_max));
            let factors = factor_segment(c0, 0, big_n, &small_primes);
            for k in 0..big_n {
                let m = c0 + 2 * k;
                let expected = sum_binomials_mod(big_n, k, m);
                let got = sum_binomials_mod_sieved(big_n, k, m, &factors[k as usize]);
                assert_eq!(got, expected, "M={big_m} N={big_n} k={k} m={m}");
            }
        }
    }
}
