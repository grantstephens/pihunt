//! Integration tests for `pi_digits::nthdigit2` (the Theorem-2 port). Building-block tests (the
//! recurrence, Lucas, p-adic recursion, the remainder tree, partial-fraction reassembly) live
//! as unit tests inside `src/nthdigit2.rs` itself, since they exercise private helpers; this
//! file covers what's reachable through the public API: full digit equality against Theorem 1
//! and against MPFR, across several `mem_bits` values, plus the CLI. See
//! `docs/nthdigit-theorem2.md` and `tests/nthdigit.rs` (the Theorem-1 equivalent this mirrors).

use pi_digits::digits_to_bits;
use pi_digits::nthdigit;
use pi_digits::nthdigit2::{self, default_mem_bits};
use rug::{Float, Integer, float::Constant, ops::Pow};

/// A tiny deterministic xorshift64 PRNG (matches `tests/nthdigit.rs`), reproducible without a
/// `rand` dependency.
struct Xorshift(u64);
impl Xorshift {
    fn next(&mut self) -> u64 {
        self.0 ^= self.0 << 13;
        self.0 ^= self.0 >> 7;
        self.0 ^= self.0 << 17;
        self.0
    }
}

/// First `count` decimal digits of pi (`reference[i]` = digit at position `i+1`), via MPFR
/// directly. Same approach as `tests/nthdigit.rs::reference_pi_digits`.
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

fn assert_digits_match(reference: &str, n: u64, count: usize, mem_bits: u64) {
    let got = nthdigit2::digits(n, count, mem_bits);
    let start = n as usize;
    let expected = &reference[start..start + count];
    assert_eq!(got, expected, "n={n} count={count} mem_bits={mem_bits}");
}

// ---------------------------------------------------------------------------------
// digits() vs MPFR, across several mem_bits values (tiny, ~sqrt(n)-ish, large)
// ---------------------------------------------------------------------------------

#[test]
fn digits_match_mpfr_reference_across_mem_bits() {
    let count_ref: u64 = 20_000;
    let reference = reference_pi_digits(count_ref);
    let count = 10usize;

    for &mem_bits in &[192u64, 1024, 65536] {
        for n in (0..2000u64).step_by(37) {
            assert_digits_match(&reference, n, count, mem_bits);
        }
        // The Feynman point: six consecutive 9s starting at decimal position 762.
        assert_digits_match(&reference, 761, 10, mem_bits);
    }

    let mut rng = Xorshift(0xD1B5_4A32_D192_ED03);
    let bound = count_ref - count as u64;
    for _ in 0..100 {
        let n = rng.next() % bound;
        assert_digits_match(&reference, n, count, default_mem_bits(n.max(1)));
    }
}

/// Requests longer than one certified chunk are stitched from independent chunks (same
/// chunking logic as Theorem 1's `digits`) and must still match MPFR exactly across the seams.
#[test]
fn long_requests_are_chunked_correctly() {
    let reference = reference_pi_digits(6_000);
    assert_eq!(
        nthdigit2::digits(5_000, 50, default_mem_bits(5_000)),
        reference[5_000..5_050]
    );
}

// ---------------------------------------------------------------------------------
// digits() vs Theorem 1, across several mem_bits values and many positions
// ---------------------------------------------------------------------------------

#[test]
fn digits_match_theorem1_across_positions_and_mem_bits() {
    let mut rng = Xorshift(0x243F_6A88_85A3_08D3);
    for &mem_bits_kind in &["tiny", "sqrt", "large"] {
        for _ in 0..15 {
            let n = 2_000 + rng.next() % 40_000;
            let mem_bits = match mem_bits_kind {
                "tiny" => 192,
                "sqrt" => default_mem_bits(n),
                _ => 131_072,
            };
            let count = 10;
            let thm1 = nthdigit::digits(n, count);
            let thm2 = nthdigit2::digits(n, count, mem_bits);
            assert_eq!(thm1, thm2, "n={n} mem_bits={mem_bits} ({mem_bits_kind})");
        }
    }
}

// ---------------------------------------------------------------------------------
// Large-n sanity checks (slow: real work, not just MPFR reference generation)
// ---------------------------------------------------------------------------------

#[test]
#[ignore] // slow: full Theorem-2 run + MPFR at ~1e5 digits of precision
fn digit_at_1e5_matches_mpfr() {
    check_large_digit(100_000);
}

#[test]
#[ignore] // slow: full Theorem-2 run + MPFR at ~1e6 digits of precision
fn digit_at_1e6_matches_mpfr() {
    check_large_digit(1_000_000);
}

// No 10^7 test: ~2 minutes even on its own (peak RSS is ~70 MiB since the 2026-09-23 memory
// fix; it was ~450 MiB, which alongside the other ignored checks helped OOM an 8 GiB machine).
// Its digits are MPFR-verified in docs/nthdigit.md ("Position convention"):
// `nthdigit2::digits(10_000_000, 10, ..)` = `pihunt digit 10000001` = `2591513361`, and
// `pihunt digit 10000000` = `7259151336`.

/// The large-position checks each use every core and up to ~100 MiB; run them one at a time
/// even when the test harness runs `--include-ignored` in parallel.
static HEAVY: std::sync::Mutex<()> = std::sync::Mutex::new(());

fn check_large_digit(n: u64) {
    let _one_at_a_time = HEAVY.lock().unwrap_or_else(|e| e.into_inner());
    let count = 10usize;
    let got = nthdigit2::digits(n, count, default_mem_bits(n));

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
