//! Relations already in the literature, so rediscoveries don't get flagged as NEW.

use crate::basis::Shape;
use crate::pslq::primitive;
use rug::Integer;

struct Known {
    base: u32,
    period: u32,
    /// (column name, coefficient) for every nonzero coefficient, primitive form.
    terms: &'static [(&'static str, i64)],
}

const KNOWN: &[Known] = &[
    // Bailey–Borwein–Plouffe (1995).
    Known {
        base: 16,
        period: 8,
        terms: &[
            ("pi", 1),
            ("S(j=1,s=1)", -4),
            ("S(j=4,s=1)", 2),
            ("S(j=5,s=1)", 1),
            ("S(j=6,s=1)", 1),
        ],
    },
    // The base-16 "zero relation" that makes BBP non-unique.
    Known {
        base: 16,
        period: 8,
        terms: &[
            ("S(j=1,s=1)", 8),
            ("S(j=2,s=1)", -8),
            ("S(j=3,s=1)", -4),
            ("S(j=4,s=1)", -8),
            ("S(j=5,s=1)", -2),
            ("S(j=6,s=1)", -2),
            ("S(j=7,s=1)", 1),
        ],
    },
    // Bailey's base-64 formula for pi^2.
    Known {
        base: 64,
        period: 6,
        terms: &[
            ("S(j=1,s=2)", 144),
            ("S(j=2,s=2)", -216),
            ("S(j=3,s=2)", -72),
            ("S(j=4,s=2)", -54),
            ("S(j=5,s=2)", 9),
            ("pi2", -8),
        ],
    },
    // log 2 = sum 1 / (2^(k+1) (k+1)).
    Known {
        base: 2,
        period: 1,
        terms: &[("S(j=1,s=1)", 1), ("log2", -2)],
    },
];

/// True if `relation` (aligned with `names`) is a known formula for this shape.
pub fn is_known(shape: &Shape, names: &[String], relation: &[Integer]) -> bool {
    let rel = primitive(relation);
    let mut found: Vec<(&str, Integer)> = names
        .iter()
        .zip(rel)
        .filter(|(_, c)| *c != 0)
        .map(|(n, c)| (n.as_str(), c))
        .collect();
    found.sort_by(|a, b| a.0.cmp(b.0));
    KNOWN
        .iter()
        .filter(|k| k.base == shape.base && k.period == shape.period)
        .any(|k| {
            let mut want: Vec<(&str, i64)> = k.terms.to_vec();
            want.sort_by(|a, b| a.0.cmp(b.0));
            want.len() == found.len()
                && want
                    .iter()
                    .zip(&found)
                    .all(|(w, f)| w.0 == f.0 && f.1 == w.1)
        })
}
