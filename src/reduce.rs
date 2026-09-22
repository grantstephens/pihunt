//! Strips rational linear dependencies out of the non-π columns.

use crate::basis::Columns;
use crate::pslq::{Outcome, PslqParams, RelationFinder, exceeds, primitive};
use crate::verify;
use rug::{Float, Integer};

/// A column removed because it is a rational combination of the others.
#[derive(Debug, Clone, PartialEq)]
pub struct Dropped {
    pub column: usize,
    /// Primitive relation, one coefficient per column of the full basis.
    pub relation: Vec<Integer>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Reduced {
    /// `keep` holds the surviving non-π column indices, ascending.
    Done {
        keep: Vec<usize>,
        dropped: Vec<Dropped>,
    },
    /// A step could neither find nor exclude a relation, or found one that failed checks.
    Inconclusive {
        outcome: Outcome,
        dropped: Vec<Dropped>,
    },
}

/// Repeatedly search columns 1.. (everything except π) and drop dependent ones.
/// A relation is only trusted if max |coeff| <= coeff_bound and it verifies at 2x precision.
pub fn reduce(cols: &Columns, finder: &dyn RelationFinder, params: &PslqParams) -> Reduced {
    let n = cols.lo.len();
    let mut keep: Vec<usize> = (1..n).collect();
    let mut dropped = Vec::new();
    while keep.len() >= 2 {
        let x: Vec<Float> = keep.iter().map(|&i| cols.lo[i].clone()).collect();
        let outcome = finder.find(&x, params);
        let Outcome::Relation { coeffs, .. } = &outcome else {
            if matches!(outcome, Outcome::Excluded { .. }) {
                break;
            }
            return Reduced::Inconclusive { outcome, dropped };
        };
        let mut full = vec![Integer::new(); n];
        for (&i, c) in keep.iter().zip(coeffs) {
            full[i] = c.clone();
        }
        let relation = primitive(&full);
        let too_big = exceeds(&relation, params.coeff_bound);
        if too_big || !verify::passes(&relation, &cols.hi, cols.digits).0 {
            return Reduced::Inconclusive { outcome, dropped };
        }
        let column = *keep
            .iter()
            .rev()
            .find(|&&i| relation[i] != 0)
            .expect("nonzero relation");
        keep.retain(|&i| i != column);
        dropped.push(Dropped { column, relation });
    }
    Reduced::Done { keep, dropped }
}
