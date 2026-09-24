//! `wasm-bindgen` wrapper over `pi-digits` (built with `--no-default-features --features pure`,
//! since WASM has no system GMP to link). Exposes both nthdigit algorithms plus the current
//! WASM linear memory size, for the site's live demo.

use wasm_bindgen::prelude::*;

const MAX_POS: f64 = 200_000.0;

/// Validates a 1-based `pos` and converts it to `n = pos - 1` (the 0-based offset the underlying
/// `pi_digits` functions take). Rejects non-finite, non-integer, or out-of-range values with a
/// JS `Error` rather than panicking or silently truncating.
fn check_pos(pos: f64) -> Result<u64, JsError> {
    if !(pos.is_finite() && pos.fract() == 0.0 && pos >= 1.0 && pos <= MAX_POS) {
        return Err(JsError::new("position must be an integer in 1..=200000"));
    }
    Ok(pos as u64)
}

/// Digits at 1-based positions `pos..pos+count` via Gourdon's Theorem 1.
#[wasm_bindgen]
pub fn digits_thm1(pos: f64, count: u32) -> Result<String, JsError> {
    let p = check_pos(pos)?;
    Ok(pi_digits::nthdigit::digits(p - 1, count as usize))
}

/// Same, via the reconstructed Theorem 2 at its default memory budget.
#[wasm_bindgen]
pub fn digits_thm2(pos: f64, count: u32) -> Result<String, JsError> {
    let p = check_pos(pos)?;
    let n = p - 1;
    Ok(pi_digits::nthdigit2::digits(
        n,
        count as usize,
        pi_digits::nthdigit2::default_mem_bits(n.max(1)),
    ))
}

/// Current size of this module's linear memory, in bytes (it only ever grows). Returns `0.0`
/// outside a `wasm32` target (e.g. under native `cargo test`/`clippy` of this wrapper), where
/// there's no WASM linear memory to report.
#[wasm_bindgen]
pub fn wasm_memory_bytes() -> f64 {
    #[cfg(target_arch = "wasm32")]
    {
        (core::arch::wasm32::memory_size(0) * 65536) as f64
    }
    #[cfg(not(target_arch = "wasm32"))]
    {
        0.0
    }
}
