//! Low-memory decimal digit extraction for π: Gourdon's Theorem 1 (`nthdigit`) and the
//! reconstructed Theorem 2 (`nthdigit2`). Shared by the `pihunt` CLI and the site's WASM demo.

pub mod mem_profile;
pub mod nthdigit;
pub mod nthdigit2;

/// Decimal digits → bits, rounded up.
pub fn digits_to_bits(digits: u32) -> u32 {
    (digits as f64 * std::f64::consts::LOG2_10).ceil() as u32
}
