//! Checks a candidate relation against columns computed at twice the search precision.

use rug::{Float, Integer};

/// log10 |a · x|, or -inf if the residual is exactly zero.
pub fn residual_log10(relation: &[Integer], cols: &[Float]) -> f64 {
    let bits = cols[0].prec();
    let mut sum = Float::with_val(bits, 0);
    for (a, x) in relation.iter().zip(cols) {
        sum += Float::with_val(bits, x * a);
    }
    if sum.is_zero() {
        return f64::NEG_INFINITY;
    }
    // log10 via MPFR so tiny residuals don't underflow f64.
    Float::with_val(64, sum.abs_ref()).log10().to_f64()
}

/// A true relation's residual shrinks with the doubled precision; a spurious one stays put.
/// `hi` must be accurate to 2 * `digits` decimal digits.
pub fn passes(relation: &[Integer], hi: &[Float], digits: u32) -> (bool, f64) {
    let r = residual_log10(relation, hi);
    (r < -(2.0 * digits as f64 - 60.0), r)
}
