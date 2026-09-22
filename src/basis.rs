//! Builds the vector of constants PSLQ searches over.

use rug::{Float, Integer, float::Constant, ops::Pow};
use serde::{Deserialize, Serialize};

/// Extra constants that may appear alongside π and the series columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Extra {
    Pi2,
    Log2,
    Log3,
    Log5,
    Catalan,
    Zeta3,
}

impl Extra {
    pub fn name(self) -> &'static str {
        match self {
            Extra::Pi2 => "pi2",
            Extra::Log2 => "log2",
            Extra::Log3 => "log3",
            Extra::Log5 => "log5",
            Extra::Catalan => "catalan",
            Extra::Zeta3 => "zeta3",
        }
    }

    pub fn value(self, bits: u32) -> Float {
        match self {
            Extra::Pi2 => Float::with_val(bits, Constant::Pi).square(),
            Extra::Log2 => Float::with_val(bits, Constant::Log2),
            Extra::Log3 => Float::with_val(bits, 3).ln(),
            Extra::Log5 => Float::with_val(bits, 5).ln(),
            Extra::Catalan => Float::with_val(bits, Constant::Catalan),
            Extra::Zeta3 => Float::with_val(bits, Float::zeta_u(3)),
        }
    }
}

/// The formula shape for one job: base b, period m, degrees s_lo..=s_hi, extras.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Shape {
    pub base: u32,
    pub period: u32,
    pub s_lo: u32,
    pub s_hi: u32,
    /// Kept sorted and deduplicated.
    pub extras: Vec<Extra>,
}

impl Shape {
    pub fn new(base: u32, period: u32, s_lo: u32, s_hi: u32, mut extras: Vec<Extra>) -> Self {
        extras.sort();
        extras.dedup();
        Shape {
            base,
            period,
            s_lo,
            s_hi,
            extras,
        }
    }

    /// Total column count including π.
    pub fn columns(&self) -> usize {
        1 + (self.period * (self.s_hi - self.s_lo + 1)) as usize + self.extras.len()
    }

    /// Column labels, in the same order as `build`.
    pub fn column_names(&self) -> Vec<String> {
        let mut names = vec!["pi".to_string()];
        for s in self.s_lo..=self.s_hi {
            for j in 1..=self.period {
                names.push(format!("S(j={j},s={s})"));
            }
        }
        names.extend(self.extras.iter().map(|e| e.name().to_string()));
        names
    }

    /// All column values at `bits` of precision: π, series (degree-major), extras.
    pub fn build(&self, bits: u32) -> Vec<Float> {
        let mut cols = vec![Float::with_val(bits, Constant::Pi)];
        for s in self.s_lo..=self.s_hi {
            for j in 1..=self.period {
                cols.push(series(self.base, self.period, j, s, bits));
            }
        }
        cols.extend(self.extras.iter().map(|e| e.value(bits)));
        cols
    }
}

/// S(j,s) = sum_{k>=0} 1 / (b^k (mk+j)^s), accurate to `bits`.
pub fn series(base: u32, period: u32, j: u32, s: u32, bits: u32) -> Float {
    let log2b = (base as f64).log2();
    let rough_terms = (bits as f64 / log2b).ceil() + 2.0;
    let guard = rough_terms.log2().ceil() as u32 + 32;
    let wp = bits + guard;
    let terms = ((wp as f64) / log2b).ceil() as u64 + 2;

    let mut sum = Float::with_val(wp, 0);
    let mut pk = Float::with_val(wp, 1); // b^-k
    for k in 0..terms {
        let d = Integer::from(period as u64 * k + j as u64).pow(s);
        sum += Float::with_val(wp, &pk / &d);
        pk /= base;
    }
    Float::with_val(bits, &sum)
}

/// Decimal digits needed so PSLQ can find relations with max |coeff| <= c among n columns.
/// Factor 1.5 measured in prototyping: 1.25 let spurious 10^7-size relations through at n ~ 45.
pub fn auto_digits(n: usize, coeff_bound: u64) -> u32 {
    (n as f64 * (coeff_bound as f64).log10() * 1.5).ceil() as u32 + 50
}

/// Column values at working precision (`lo`) and at twice it (`hi`, for verification).
pub struct Columns {
    pub names: Vec<String>,
    pub digits: u32,
    pub lo: Vec<Float>,
    pub hi: Vec<Float>,
}

impl Columns {
    pub fn build(shape: &Shape, digits: u32) -> Self {
        let hi = shape.build(crate::pslq::digits_to_bits(2 * digits));
        let bits = crate::pslq::digits_to_bits(digits);
        let lo = hi.iter().map(|v| Float::with_val(bits, v)).collect();
        Columns {
            names: shape.column_names(),
            digits,
            lo,
            hi,
        }
    }
}

impl std::str::FromStr for Extra {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        [
            Extra::Pi2,
            Extra::Log2,
            Extra::Log3,
            Extra::Log5,
            Extra::Catalan,
            Extra::Zeta3,
        ]
        .into_iter()
        .find(|e| e.name() == s)
        .ok_or_else(|| format!("unknown extra {s:?}"))
    }
}
