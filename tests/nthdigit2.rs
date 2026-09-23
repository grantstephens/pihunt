//! Integration tests for `pihunt::nthdigit2` (the Theorem-2 port). Building-block tests (the
//! recurrence, Lucas, p-adic recursion, the remainder tree, partial-fraction reassembly) live
//! as unit tests inside `src/nthdigit2.rs` itself, since they exercise private helpers; this
//! file covers what's reachable through the public API: full digit equality against Theorem 1
//! and against MPFR, across several `mem_bits` values, plus the CLI. See
//! `docs/nthdigit-theorem2.md` and `tests/nthdigit.rs` (the Theorem-1 equivalent this mirrors).

use pihunt::nthdigit;
use pihunt::nthdigit2::{self, default_mem_bits};
use pihunt::pslq::digits_to_bits;
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

fn check_large_digit(n: u64) {
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

// ---------------------------------------------------------------------------------
// CLI
// ---------------------------------------------------------------------------------

#[test]
fn cli_digit_thm2_matches_thm1() {
    let bin = env!("CARGO_BIN_EXE_pihunt");
    let out1 = std::process::Command::new(bin)
        .args(["digit", "1", "--count", "20", "--method", "thm1"])
        .output()
        .unwrap();
    let out2 = std::process::Command::new(bin)
        .args([
            "digit", "1", "--count", "20", "--method", "thm2", "--mem", "512",
        ])
        .output()
        .unwrap();
    assert!(out1.status.success());
    assert!(
        out2.status.success(),
        "{}",
        String::from_utf8_lossy(&out2.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&out1.stdout).trim(),
        String::from_utf8_lossy(&out2.stdout).trim()
    );
    assert_eq!(
        String::from_utf8_lossy(&out2.stdout).trim(),
        "14159265358979323846"
    );
}

#[test]
fn cli_stream_thm2_prints_independent_blocks() {
    let bin = env!("CARGO_BIN_EXE_pihunt");
    let out = std::process::Command::new(bin)
        .args([
            "stream", "--from", "1", "--block", "5", "--blocks", "2", "--method", "thm2", "--mem",
            "512",
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines, vec!["14159", "26535"]);
}
