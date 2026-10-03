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
    // Bellard (1997): pi = 2^-6 sum_k (-1)^k 2^(-10k) (-2^5/(4k+1) - 1/(4k+3) + 2^8/(10k+1)
    // - 2^6/(10k+3) - 2^2/(10k+5) - 2^2/(10k+7) + 1/(10k+9)). The (-1)^k/2^(10k) weight is
    // folded into a positive base by pairing k = 2K, 2K+1: base 2^20, period 40, with
    // 1/(4k+a) = 5/(40K+20r+5a) and 1/(10k+c) = 2/(40K+20r+2c). The first entry is that
    // form written out directly; the second is what PSLQ returns after the reduction step
    // drops dependent columns. They differ by zero relations (checked in tests/known.rs).
    Known {
        base: 1 << 20,
        period: 40,
        terms: &[
            ("pi", 65536),
            ("S(j=2,s=1)", -524288),
            ("S(j=5,s=1)", 163840),
            ("S(j=6,s=1)", 131072),
            ("S(j=10,s=1)", 8192),
            ("S(j=14,s=1)", 8192),
            ("S(j=15,s=1)", 5120),
            ("S(j=18,s=1)", -2048),
            ("S(j=22,s=1)", 512),
            ("S(j=25,s=1)", -160),
            ("S(j=26,s=1)", -128),
            ("S(j=30,s=1)", -8),
            ("S(j=34,s=1)", -8),
            ("S(j=35,s=1)", -5),
            ("S(j=38,s=1)", 2),
        ],
    },
    Known {
        base: 1 << 20,
        period: 40,
        terms: &[
            ("pi", 65536),
            ("S(j=2,s=1)", -1048576),
            ("S(j=4,s=1)", 786432),
            ("S(j=5,s=1)", 819200),
            ("S(j=8,s=1)", -327680),
            ("S(j=10,s=1)", -24576),
            ("S(j=12,s=1)", 49152),
            ("S(j=15,s=1)", -15360),
            ("S(j=16,s=1)", -20480),
            ("S(j=18,s=1)", -4096),
            ("S(j=20,s=1)", -2048),
            ("S(j=24,s=1)", -1280),
            ("S(j=25,s=1)", -800),
            ("S(j=26,s=1)", -256),
            ("S(j=28,s=1)", 192),
            ("S(j=30,s=1)", -40),
            ("S(j=32,s=1)", -80),
            ("S(j=34,s=1)", -16),
            ("S(j=35,s=1)", 15),
            ("S(j=36,s=1)", 12),
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
