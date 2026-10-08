//! Integer relation finding.

pub mod classic;
pub mod multilevel;
mod state;

use rug::{Float, Integer};

/// Tuning and stopping parameters for one relation search.
#[derive(Debug, Clone)]
pub struct PslqParams {
    /// Must be > sqrt(4/3).
    pub gamma: f64,
    /// Largest |coefficient| we care about; drives the exclusion bound.
    pub coeff_bound: u64,
    pub max_iterations: u64,
    /// Working precision in decimal digits. Inputs must be accurate to this.
    pub digits: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    /// `coeffs · x ≈ 0`, one coefficient per input.
    Relation {
        coeffs: Vec<Integer>,
        iterations: u64,
        bound: f64,
        closeness: f64,
    },
    /// No relation with max |coefficient| ≤ coeff_bound exists.
    Excluded {
        bound: f64,
        iterations: u64,
        closeness: f64,
    },
    PrecisionExhausted {
        bound: f64,
        iterations: u64,
        closeness: f64,
    },
    IterationCap {
        bound: f64,
        iterations: u64,
        closeness: f64,
    },
}

impl Outcome {
    pub fn iterations(&self) -> u64 {
        match self {
            Outcome::Relation { iterations, .. }
            | Outcome::Excluded { iterations, .. }
            | Outcome::PrecisionExhausted { iterations, .. }
            | Outcome::IterationCap { iterations, .. } => *iterations,
        }
    }

    pub fn bound(&self) -> f64 {
        match self {
            Outcome::Relation { bound, .. }
            | Outcome::Excluded { bound, .. }
            | Outcome::PrecisionExhausted { bound, .. }
            | Outcome::IterationCap { bound, .. } => *bound,
        }
    }

    /// How close the reduction's best candidate got to a relation: the smallest `|y_i|` at
    /// termination, divided by the detection threshold `eps`. `< 1` means a relation was found
    /// (`Outcome::Relation`); `>= 1` for every other outcome, and the closer to 1, the nearer
    /// this shape came to one without quite getting there. Purely observational — never
    /// affects exclusion/hit logic, just richer per-job data for later analysis (e.g. picking
    /// where to sample next). See `docs/findings-bbp-hunt.md`.
    pub fn closeness(&self) -> f64 {
        match self {
            Outcome::Relation { closeness, .. }
            | Outcome::Excluded { closeness, .. }
            | Outcome::PrecisionExhausted { closeness, .. }
            | Outcome::IterationCap { closeness, .. } => *closeness,
        }
    }
}

pub trait RelationFinder: Sync {
    fn name(&self) -> &'static str;
    fn find(&self, x: &[Float], params: &PslqParams) -> Outcome;
}

/// Decimal digits → MPFR bits.
pub fn digits_to_bits(digits: u32) -> u32 {
    (digits as f64 * std::f64::consts::LOG2_10).ceil() as u32
}

/// Divide by the gcd and flip sign so the first nonzero coefficient is positive.
pub fn primitive(coeffs: &[Integer]) -> Vec<Integer> {
    let g = coeffs.iter().fold(Integer::new(), |g, c| g.gcd(c));
    let sign = match coeffs.iter().find(|c| **c != 0) {
        Some(c) if *c < 0 => -1,
        _ => 1,
    };
    if g == 0 {
        return coeffs.to_vec();
    }
    coeffs
        .iter()
        .map(|c| Integer::from(c / &g) * sign)
        .collect()
}

/// True if any |coefficient| is larger than `bound`.
pub fn exceeds(coeffs: &[Integer], bound: u64) -> bool {
    coeffs.iter().any(|c| *c.as_abs() > bound)
}
