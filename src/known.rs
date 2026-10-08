//! Relations already in the literature, so rediscoveries don't get flagged as NEW.

use crate::basis::Shape;
use crate::pslq::primitive;
use rug::{Integer, ops::Pow};
use std::collections::HashMap;

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

/// True if `relation` is a mechanical corollary of some entry in [`KNOWN`]: a "base power +
/// period multiple" descendant, not an independent new formula.
///
/// Two trivial moves generate corollaries from any base-`b0`, period-`m0` relation among
/// `S(b0,m0,j,s)` columns (plus base-independent terms like π, π², log 2, ...):
///
/// - **Base power** (any `t ≥ 1`): splitting the index `k = t·k' + r` gives
///   `S(b0,m0,j,s) = Σ_{r=0}^{t-1} b0^(-r)·S(b0^t, m0·t, m0·r+j, s)`.
/// - **Period multiple** (any `d ≥ 1`, same base): `S(B,M,J,s) = d^s·S(B, M·d, J·d, s)`.
///
/// Composing both and clearing the `b0^(-r)` fractions by scaling the whole equation by
/// `b0^(t-1)` gives an exact, checkable relation at base `B = b0^t`, period `M = m0·t·d`. If the
/// candidate's base/period fit that pattern for some known entry and some `(t, d)`, and the
/// predicted relation matches the candidate's up to the usual gcd/sign normalisation, it's that
/// descendant, not a discovery. See `docs/findings-bbp-hunt.md` ("a NEW hit that wasn't") for the
/// worked example (base 256, period 32 — a `t=2, d=2` descendant of the 1996 BBP formula) that
/// motivated this check.
pub fn is_derived(shape: &Shape, names: &[String], relation: &[Integer]) -> bool {
    let rel = primitive(relation);
    let idx: HashMap<&str, usize> = names
        .iter()
        .enumerate()
        .map(|(i, n)| (n.as_str(), i))
        .collect();
    KNOWN.iter().any(|k| {
        let Some(t) = exact_power(shape.base, k.base) else {
            return false;
        };
        let denom = k.period * t;
        if denom == 0 || !shape.period.is_multiple_of(denom) {
            return false;
        }
        let d = shape.period / denom;
        match predicted_relation(k, t, d, &idx, names.len()) {
            Some(predicted) => primitive(&predicted) == rel,
            None => false,
        }
    })
}

/// Smallest `t` with `root^t == base`, or `None` if `base` isn't an exact power of `root`.
fn exact_power(base: u32, root: u32) -> Option<u32> {
    if root < 2 {
        return None;
    }
    let mut val = root as u64;
    let mut t = 1u32;
    while val < base as u64 {
        val *= root as u64;
        t += 1;
    }
    (val == base as u64).then_some(t)
}

/// Parses `"S(j=J,s=S)"` into `(J, S)`; `None` for non-series columns (π, extras).
fn parse_series_name(name: &str) -> Option<(u32, u32)> {
    let rest = name.strip_prefix("S(j=")?;
    let (j, rest) = rest.split_once(",s=")?;
    let s = rest.strip_suffix(')')?;
    Some((j.parse().ok()?, s.parse().ok()?))
}

/// The base-power (`t`) + period-multiple (`d`) descendant of `k`, aligned to `names`'
/// column order. `None` if the descendant shape doesn't have a matching column for some term
/// (e.g. the candidate's degree range doesn't cover a term's `s`).
fn predicted_relation(
    k: &Known,
    t: u32,
    d: u32,
    idx: &HashMap<&str, usize>,
    len: usize,
) -> Option<Vec<Integer>> {
    let mut predicted = vec![Integer::new(); len];
    for &(name, coeff) in k.terms {
        if let Some((j0, s0)) = parse_series_name(name) {
            for r in 0..t {
                let j = d * (k.period * r + j0);
                let col = format!("S(j={j},s={s0})");
                let i = *idx.get(col.as_str())?;
                let weight = Integer::from(k.base).pow(t - 1 - r) * Integer::from(d).pow(s0);
                predicted[i] += Integer::from(coeff) * weight;
            }
        } else {
            let i = *idx.get(name)?;
            let weight = Integer::from(k.base).pow(t - 1);
            predicted[i] += Integer::from(coeff) * weight;
        }
    }
    Some(predicted)
}
