//! Turns a main-search PSLQ outcome into a log verdict.

use crate::basis::{Columns, Shape};
use crate::known;
use crate::log::{Kind, Verify};
use crate::pslq::{Outcome, exceeds};
use crate::verify;
use rug::Integer;

#[derive(Debug, Clone, PartialEq)]
pub struct Verdict {
    pub kind: Kind,
    /// Aligned with the full column list; normalised so the π coefficient is positive.
    pub relation: Option<Vec<Integer>>,
    pub verify: Option<Verify>,
    pub tag: Option<String>,
    pub note: Option<String>,
}

/// `relation` from a main-search `Outcome::Relation` must already be expanded to full columns.
pub fn classify(
    shape: &Shape,
    cols: &Columns,
    outcome: &Outcome,
    relation: Option<Vec<Integer>>,
    coeff_bound: u64,
) -> Verdict {
    let plain = |kind, note: Option<&str>| Verdict {
        kind,
        relation: None,
        verify: None,
        tag: None,
        note: note.map(String::from),
    };
    match outcome {
        Outcome::Excluded { .. } => plain(Kind::Excluded, None),
        Outcome::PrecisionExhausted { .. } => {
            plain(Kind::Inconclusive, Some("precision exhausted"))
        }
        Outcome::IterationCap { .. } => plain(Kind::Inconclusive, Some("iteration cap")),
        Outcome::Relation { .. } => {
            let rel = normalise(relation.expect("relation outcome carries coefficients"));
            let mut v = Verdict {
                kind: Kind::Junk,
                relation: None,
                verify: None,
                tag: None,
                note: None,
            };
            if rel[0] == 0 {
                v.note = Some("relation without pi after reduction — bug signal".into());
            } else if exceeds(&rel, coeff_bound) {
                v.kind = Kind::Suspicious;
                v.note = Some("coefficient exceeds coeff_bound".into());
            } else {
                let (passed, residual_log10) = verify::passes(&rel, &cols.hi, cols.digits);
                v.verify = Some(Verify {
                    passed,
                    residual_log10,
                });
                v.kind = if passed { Kind::Hit } else { Kind::Spurious };
                if passed {
                    let known = known::is_known(shape, &cols.names, &rel);
                    v.tag = Some(if known { "known" } else { "NEW" }.into());
                }
            }
            v.relation = Some(rel);
            v
        }
    }
}

/// Divide by gcd; make the π coefficient positive (or the first nonzero one if π is absent).
fn normalise(rel: Vec<Integer>) -> Vec<Integer> {
    let mut p = crate::pslq::primitive(&rel);
    if p[0] < 0 {
        p.iter_mut().for_each(|c| *c = Integer::from(-&*c));
    }
    p
}
