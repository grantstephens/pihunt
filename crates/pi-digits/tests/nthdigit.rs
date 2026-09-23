//! Tests for `pi_digits::nthdigit` (Gourdon's low-memory n-th decimal digit algorithm).
//! See `docs/nthdigit.md` and the module docs in `src/nthdigit.rs`.

use pi_digits::digits_to_bits;
use pi_digits::nthdigit::{self, Params};
use rug::{Float, Integer, float::Constant, ops::Pow};

/// A tiny deterministic xorshift64 PRNG, so tests are reproducible without a `rand` dep.
struct Xorshift(u64);
impl Xorshift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

/// Exact `Σ_{j=0}^{k} C(N,j)`, via `rug::Integer::binomial` (brute force, only used for
/// small N/k in these tests).
fn s_k_exact(big_n: u64, k: u64) -> Integer {
    let mut sum = Integer::from(0);
    for j in 0..=k {
        sum += Integer::from(big_n).binomial(j as u32);
    }
    sum
}

// ---------------------------------------------------------------------------------
// Parameters
// ---------------------------------------------------------------------------------

#[test]
fn params_are_sane() {
    for &(n, n0) in &[
        (2000u64, 20u32),
        (10_000, 30),
        (100_000, 40),
        (1_000_000, 50),
    ] {
        let p = Params::new(n, n0);
        assert!(p.m.is_multiple_of(2), "M must be even, got {}", p.m);
        assert!(p.m >= 4, "M must be >= 4, got {}", p.m);
        assert!(p.big_n.is_multiple_of(2), "N must be even, got {}", p.big_n);
        assert!(p.big_n <= n + 2, "N={} must be <= n+2={}", p.big_n, n + 2);
        assert!(
            p.error_bound_log10() < -((n + n0 as u64) as f64),
            "n={n} n0={n0}: error bound 1e{} not below 1e-{}",
            p.error_bound_log10(),
            n + n0 as u64
        );
    }
}

#[test]
#[should_panic(expected = "exceeds n+2")]
fn params_reject_n_too_small_for_n0() {
    // n0 way too large relative to n: N would have to exceed n+2.
    Params::new(100, 1000);
}

/// The identity `5^(N-2) * 10^(n-N+2) * 2^N == 4 * 10^n` that makes the C-sum numerators
/// integers (Gourdon's key observation, just above Proposition 1).
#[test]
fn power_of_ten_identity_holds() {
    for &(n, n0) in &[(50u64, 5u32), (200u64, 20u32), (2000u64, 30u32)] {
        let p = Params::new(n, n0);
        let lhs = Integer::from(5).pow((p.big_n - 2) as u32)
            * Integer::from(10).pow((n - p.big_n + 2) as u32)
            * Integer::from(2).pow(p.big_n as u32);
        let rhs = Integer::from(4) * Integer::from(10).pow(n as u32);
        assert_eq!(lhs, rhs, "n={n} n0={n0}: identity fails, N={}", p.big_n);
    }
}

/// Validates the Cohen-Villegas-Zagier acceleration itself (formula 7), independent of the
/// modular tricks used to evaluate it cheaply: compute `S` exactly (well, to high MPFR
/// precision) for small `n` and check it lands within Gourdon's own bound of `π`.
#[test]
fn accelerated_series_converges_to_pi() {
    for &(n, n0) in &[(50u64, 10u32), (200u64, 20u32)] {
        let p = Params::new(n, n0);
        let (big_m, big_n) = (p.m, p.big_n);
        let bits = digits_to_bits((n + n0 as u64 + 50) as u32);

        let mut s = Float::with_val(bits, 0);
        let terms1 = (big_m + 1) * big_n;
        for k in 0..terms1 {
            let term = Float::with_val(bits, 4) / Float::with_val(bits, 2 * k + 1);
            if k % 2 == 0 {
                s += &term;
            } else {
                s -= &term;
            }
        }

        let two_pow_n = Float::with_val(bits, Integer::from(2).pow(big_n as u32));
        for k in 0..big_n {
            let sk = s_k_exact(big_n, k);
            let m = 2 * big_m * big_n + 2 * k + 1;
            let num = Float::with_val(bits, 4) * Float::with_val(bits, &sk);
            let den = Float::with_val(bits, &two_pow_n) * Float::with_val(bits, m);
            let term = num / den;
            if k % 2 == 0 {
                s -= &term;
            } else {
                s += &term;
            }
        }

        let pi = Float::with_val(bits, Constant::Pi);
        let diff = Float::with_val(bits, &s - &pi).abs();
        let bound = Float::with_val(bits, 1)
            / Float::with_val(bits, Integer::from(10).pow((n + n0 as u64) as u32));
        assert!(
            diff < bound,
            "n={n}: |S-pi|={diff} exceeds bound 1e-{}",
            n + n0 as u64
        );
    }
}

// ---------------------------------------------------------------------------------
// Algorithm 2 (sum of binomials mod m)
// ---------------------------------------------------------------------------------

#[test]
fn algorithm2_matches_brute_force() {
    let cases: &[(u64, u64, u64)] = &[
        (10, 3, 7),
        (10, 7, 13),
        // composite m with a repeated small prime factor, k below and above N/2
        (20, 5, 3u64.pow(4) * 5 * 7),
        (20, 15, 3u64.pow(4) * 5 * 7),
        (50, 10, 3u64.pow(3) * 5u64.pow(2) * 11),
        (50, 40, 3u64.pow(3) * 5u64.pow(2) * 11),
        (64, 30, 3u64.pow(2) * 7u64.pow(2) * 13),
        (64, 50, 3u64.pow(2) * 7u64.pow(2) * 13),
        // odd composite with several distinct small primes
        (100, 20, 3 * 5 * 7 * 11 * 13),
        (100, 80, 3 * 5 * 7 * 11 * 13),
        // k exactly at N/2
        (40, 20, 3u64.pow(3) * 5 * 7 * 11),
        // m odd prime larger than k
        (30, 10, 9973),
    ];
    for &(big_n, k, m) in cases {
        let expected: u64 = (s_k_exact(big_n, k) % Integer::from(m))
            .to_string()
            .parse()
            .unwrap();
        let got = nthdigit::sum_binomials_mod(big_n, k, m);
        assert_eq!(got, expected, "N={big_n} k={k} m={m}");
    }
}

// ---------------------------------------------------------------------------------
// 128-by-64 fixed-point division
// ---------------------------------------------------------------------------------

#[test]
fn fixed_point_division_matches_rug_for_random_inputs() {
    let mut rng = Xorshift(0x243F_6A88_85A3_08D3);
    // A few fixed edge cases first.
    for &(x, m) in &[(0u64, 1u64), (0, 2), (1, 2), (u64::MAX - 1, u64::MAX)] {
        check_fixed_point(x, m);
    }
    for _ in 0..5000 {
        let m = (rng.next() % (u64::MAX - 1)) + 2; // m >= 2
        let x = rng.next() % m; // 0 <= x < m
        check_fixed_point(x, m);
    }
}

fn check_fixed_point(x: u64, m: u64) {
    let num = Integer::from(x) << 128u32;
    let den = Integer::from(m);
    let expected = Integer::from(&num / &den); // truncating == floor, both non-negative
    let expected: u128 = expected.to_string().parse().unwrap();
    assert_eq!(nthdigit::frac_fixed_point(x, m), expected, "x={x} m={m}");
}

// ---------------------------------------------------------------------------------
// Digits vs MPFR
// ---------------------------------------------------------------------------------

/// First `count` decimal digits of pi (`reference[i]` = digit at position `i+1`), computed
/// once via MPFR directly (independent of `nthdigit`'s own MPFR fallback path in the sense
/// that it doesn't share any code with `digits_via_mpfr`, only the underlying `Constant::Pi`
/// value).
fn reference_pi_digits(count: u64) -> String {
    let guard = 20u64;
    let total = count + guard;
    let bits = digits_to_bits(total as u32) + 8;
    let pi = Float::with_val(bits, Constant::Pi);
    let scale = Float::with_val(bits, Integer::from(10).pow(total as u32));
    let scaled = Float::with_val(bits, &pi * &scale);
    let int_part = scaled.to_integer().unwrap();
    let s = int_part.to_string(); // "3" + `total` digits
    s[1..1 + count as usize].to_string()
}

fn assert_digits_match(reference: &str, n: u64, count: usize) {
    let got = nthdigit::digits(n, count);
    let start = n as usize;
    let expected = &reference[start..start + count];
    assert_eq!(got, expected, "n={n} count={count}");
}

#[test]
fn digits_match_mpfr_reference() {
    let count_ref: u64 = 20_000;
    let reference = reference_pi_digits(count_ref);
    let count = 10usize;

    for n in 0..200u64 {
        assert_digits_match(&reference, n, count);
    }

    let mut rng = Xorshift(0xD1B5_4A32_D192_ED03);
    let bound = count_ref - count as u64;
    for _ in 0..200 {
        let n = rng.next() % bound;
        assert_digits_match(&reference, n, count);
    }

    // The Feynman point: six consecutive 9s starting at decimal position 762.
    assert_digits_match(&reference, 761, 10);
    assert_eq!(&reference[761..767], "999999");
}

#[test]
#[ignore] // slow: MPFR at ~1e5 decimal digits of precision, plus the full Gourdon run
fn digit_at_1e5_matches_mpfr() {
    check_large_digit(100_000);
}

// Capped at 2*10^5: 10^6 takes ~75 s here. The 10^6 and 10^7 digits are MPFR-verified and
// recorded in docs/nthdigit.md.
#[test]
#[ignore] // slow: MPFR at ~2e5 decimal digits of precision, plus the full Gourdon run
fn digit_at_2e5_matches_mpfr() {
    check_large_digit(200_000);
}

/// The large-position checks each use every core; run them one at a time even under
/// `--include-ignored`.
static HEAVY: std::sync::Mutex<()> = std::sync::Mutex::new(());

fn check_large_digit(n: u64) {
    let _one_at_a_time = HEAVY.lock().unwrap_or_else(|e| e.into_inner());
    let count = 10usize;
    let got = nthdigit::digits(n, count);

    let guard = 20u64;
    let total = n + count as u64 + guard;
    let bits = digits_to_bits(total as u32) + 8;
    let pi = Float::with_val(bits, Constant::Pi);
    let scale = Float::with_val(bits, Integer::from(10).pow(total as u32));
    let scaled = Float::with_val(bits, &pi * &scale);
    let int_part = scaled.to_integer().unwrap();
    let s = int_part.to_string();
    let expected = &s[(n as usize + 1)..(n as usize + 1 + count)];
    assert_eq!(got, expected, "n={n}");
}

/// Requests longer than one certified chunk are stitched from independent chunks and must
/// still match MPFR exactly across the chunk seams.
#[test]
fn long_requests_are_chunked_correctly() {
    let reference = reference_pi_digits(6_000);
    assert_digits_match(&reference, 5_000, 50);
}
