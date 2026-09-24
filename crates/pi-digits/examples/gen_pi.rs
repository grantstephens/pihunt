//! Generates `site/data/pi-200k.txt`: the first 200,000 decimal digits of π after the point
//! (positions 1..=200000), as a single line of ASCII digits with no leading "3." and no
//! trailing newline.
//!
//! Uses [`pi_digits::pi_ref::digits`] (the exact Machin-formula fallback, certified via its own
//! guard-band-and-retry error bound) in one call, then cross-checks the result against an
//! independent MPFR computation (`rug::Float` at generous precision) before writing anything —
//! belt and suspenders for the file every other task in this plan treats as ground truth.
//!
//! Deliberately not `pihunt digit 1 --count 200000`: that CLI path chunks any request over 16
//! digits, and chunks past position ~2000 fall through to Gourdon's Theorem 1, which is
//! correct but would take on the order of hours for a request this size. Run via
//! `site/scripts/gen-pi.sh`, which builds this example under the `gmp` (default) feature.

use rug::{Float, Integer, float::Constant, ops::Pow};
use std::io::Write;
use std::path::Path;

const COUNT: usize = 200_000;

/// Independent MPFR reference for positions `1..=n+count` (1-based), used only to cross-check
/// [`pi_digits::pi_ref::digits`]'s output before trusting it. Mirrors the helper in
/// `crates/pi-digits/tests/pi_ref.rs`.
fn mpfr_digits(n: u64, count: usize) -> String {
    let total = n + count as u64 + 25;
    let bits = pi_digits::digits_to_bits(total as u32) + 16;
    let pi = Float::with_val(bits, Constant::Pi);
    let scaled = Float::with_val(
        bits,
        &pi * &Float::with_val(bits, Integer::from(10).pow(total as u32)),
    );
    let s = scaled.to_integer().unwrap().to_string(); // "3" + digits
    s[(n as usize + 1)..(n as usize + 1 + count)].to_string()
}

fn main() {
    eprintln!("computing {COUNT} digits via pi_ref::digits(0, {COUNT})...");
    let got = pi_digits::pi_ref::digits(0, COUNT);
    assert_eq!(
        got.len(),
        COUNT,
        "pi_ref::digits returned {} digits, expected {COUNT}",
        got.len()
    );
    assert!(
        got.bytes().all(|b| b.is_ascii_digit()),
        "pi_ref::digits returned a non-digit byte"
    );

    eprintln!("cross-checking against independent MPFR computation...");
    let want = mpfr_digits(0, COUNT);
    assert_eq!(got, want, "pi_ref::digits disagrees with MPFR");

    // Sanity checks a human can eyeball in the diff/log.
    assert_eq!(&got[..5], "14159", "first 5 digits should be 14159");
    assert_eq!(
        &got[761..767],
        "999999",
        "Feynman point (positions 762..=767) should be 999999"
    );

    let out_path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../site/data/pi-200k.txt");
    let mut f = std::fs::File::create(&out_path)
        .unwrap_or_else(|e| panic!("creating {}: {e}", out_path.display()));
    f.write_all(got.as_bytes())
        .unwrap_or_else(|e| panic!("writing {}: {e}", out_path.display()));

    eprintln!(
        "wrote {} bytes to {} (MPFR cross-check passed)",
        got.len(),
        out_path.display()
    );
}
