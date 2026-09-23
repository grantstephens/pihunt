//! Measures how far `auto_digits`'s factor `f` (in `ceil(n · log10(C) · f) + 50`) needs to
//! grow with `n` to avoid `Inconclusive` outcomes on real shapes. Feeds `docs/precision-rule.md`.
//!
//! cargo test --release --test precision_rule -- --ignored --nocapture

use pihunt::basis::{Extra, Shape};
use pihunt::job::run_job;
use pihunt::log::Kind;
use pihunt::plan::Job;
use pihunt::pslq::multilevel::MultilevelPslq;
use rayon::prelude::*;
use std::time::Instant;

/// Candidates searched in increasing order; the first that genuinely resolves wins.
const FACTORS: [f64; 7] = [1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5];

/// A finished sweep at the old fixed factor 1.5 showed that `Inconclusive` isn't the only
/// precision artifact: some jobs' *main* search finds a "relation" with coefficients far past
/// `coeff_bound` (a precision-floor near-miss the checks correctly reject as `Suspicious`).
/// That's just as much a sign of too few digits as `Inconclusive`, so only `Excluded` or `Hit`
/// count as a genuine resolution here — anything else keeps searching higher `f`.
fn resolved(kind: Kind) -> bool {
    matches!(kind, Kind::Excluded | Kind::Hit)
}

/// Cap iterations so a badly-underprecisioned attempt fails fast (`IterationCap`) instead of
/// grinding for minutes. The timing baseline shows genuine excludes at n ~ 67 settle in well
/// under this; it is an experiment-only cap, not the production default (1_000_000).
const MAX_ITERATIONS: u64 = 300_000;

fn digits_for(n: usize, coeff_bound: u64, f: f64) -> u32 {
    (n as f64 * (coeff_bound as f64).log10() * f).ceil() as u32 + 50
}

#[test]
#[ignore]
fn measure_precision_rule() {
    let all = [
        Extra::Pi2,
        Extra::Log2,
        Extra::Log3,
        Extra::Log5,
        Extra::Catalan,
        Extra::Zeta3,
    ];
    // Spread of real shapes: bases 10/100/1000, n from ~7 to 100, mirroring batches/*.toml.
    let shapes: Vec<Shape> = vec![
        Shape::new(10, 6, 1, 1, vec![]),              // n=7
        Shape::new(10, 10, 1, 1, all[..2].to_vec()),  // n=13
        Shape::new(100, 10, 1, 1, all[..3].to_vec()), // n=14
        Shape::new(10, 8, 1, 2, vec![]),              // n=17
        Shape::new(10, 12, 1, 1, all.to_vec()),       // n=19
        Shape::new(1000, 14, 1, 1, vec![]),           // n=15
        Shape::new(10, 12, 1, 2, all[..2].to_vec()),  // n=27
        Shape::new(100, 16, 1, 2, all[..3].to_vec()), // n=36
        Shape::new(10, 16, 1, 2, all.to_vec()),       // n=39
        Shape::new(10, 20, 1, 2, all.to_vec()),       // n=47
        Shape::new(1000, 20, 1, 2, all.to_vec()),     // n=47
        Shape::new(100, 24, 1, 2, all.to_vec()),      // n=55
        Shape::new(10, 16, 1, 3, all[..2].to_vec()),  // n=51
        Shape::new(10, 20, 1, 3, all.to_vec()),       // n=67
        Shape::new(10, 21, 1, 3, all.to_vec()),       // n=70 - flagged Suspicious upstream
        Shape::new(10, 30, 1, 2, all.to_vec()),       // n=67 - flagged Suspicious upstream
        Shape::new(10, 31, 1, 2, all.to_vec()),       // n=69 - flagged Suspicious upstream
        Shape::new(10, 32, 1, 2, all.to_vec()),       // n=71 - flagged Suspicious upstream
        Shape::new(10, 24, 1, 3, all.to_vec()),       // n=79
        Shape::new(10, 28, 1, 3, all.to_vec()),       // n=91
        Shape::new(1000, 28, 1, 3, all.to_vec()),     // n=91
        Shape::new(10, 31, 1, 3, all.to_vec()),       // n=100, the reported example
    ];
    let coeff_bounds = [1_000u64, 10_000, 100_000];

    let mut cases: Vec<(Shape, u64)> = Vec::new();
    for s in &shapes {
        for &c in &coeff_bounds {
            cases.push((s.clone(), c));
        }
    }

    // (shape, coeff_bound, n, Some((minimal f, seconds, outcome)) or None if unresolved).
    type CaseResult = (Shape, u64, usize, Option<(f64, f64, Kind)>);
    let mut results: Vec<CaseResult> = cases
        .par_iter()
        .map(|(shape, coeff_bound)| {
            let n = shape.columns();
            let mut found = None;
            for &f in &FACTORS {
                let digits = digits_for(n, *coeff_bound, f);
                let job = Job {
                    shape: shape.clone(),
                    coeff_bound: *coeff_bound,
                    digits,
                    gamma: 1.16,
                    max_iterations: MAX_ITERATIONS,
                };
                let t = Instant::now();
                let rec = run_job(&job, "precision_rule", &MultilevelPslq, 200);
                let secs = t.elapsed().as_secs_f64();
                if resolved(rec.outcome) {
                    found = Some((f, secs, rec.outcome));
                    break;
                }
            }
            (shape.clone(), *coeff_bound, n, found)
        })
        .collect();

    results.sort_by_key(|(_, c, n, _)| (*n, *c));

    println!("| base | period | degrees | #extras | n | C | minimal f | outcome | seconds |");
    println!("|---|---|---|---|---|---|---|---|---|");
    for (shape, c, n, found) in &results {
        let degrees = format!("{}..{}", shape.s_lo, shape.s_hi);
        match found {
            Some((f, secs, kind)) => println!(
                "| {} | {} | {} | {} | {} | {:.0e} | {} | {:?} | {:.2} |",
                shape.base,
                shape.period,
                degrees,
                shape.extras.len(),
                n,
                c,
                f,
                kind,
                secs
            ),
            None => println!(
                "| {} | {} | {} | {} | {} | {:.0e} | >{} | - | - |",
                shape.base,
                shape.period,
                degrees,
                shape.extras.len(),
                n,
                c,
                FACTORS.last().unwrap()
            ),
        }
    }

    let unresolved: Vec<_> = results.iter().filter(|(_, _, _, f)| f.is_none()).collect();
    assert!(
        unresolved.is_empty(),
        "{} case(s) unresolved even at f = {}: {unresolved:?}",
        unresolved.len(),
        FACTORS.last().unwrap()
    );
}

/// The exact case from the problem report: base 10, period 31, degrees 1..3, all six extras
/// (n = 100), C = 1e5. Old flat-1.5 rule gives 800 digits and comes back `Inconclusive`; the
/// new rule gives 1600 and should settle `Excluded` on the first attempt, no escalation.
#[test]
#[ignore]
fn sanity_check_n100_example() {
    use pihunt::basis::auto_digits;

    let shape = Shape::new(
        10,
        31,
        1,
        3,
        vec![
            Extra::Pi2,
            Extra::Log2,
            Extra::Log3,
            Extra::Log5,
            Extra::Catalan,
            Extra::Zeta3,
        ],
    );
    let n = shape.columns();
    assert_eq!(n, 100);
    let coeff_bound = 100_000u64;

    for (label, digits) in [
        ("old rule (flat 1.5)", 800u32),
        ("new rule", auto_digits(n, coeff_bound)),
    ] {
        let job = Job {
            shape: shape.clone(),
            coeff_bound,
            digits,
            gamma: 1.16,
            max_iterations: 300_000,
        };
        let t = Instant::now();
        let rec = run_job(&job, "sanity_n100", &MultilevelPslq, 200);
        println!(
            "{label}: digits={digits} outcome={:?} iterations={} seconds={:.2} note={:?}",
            rec.outcome,
            rec.iterations,
            t.elapsed().as_secs_f64(),
            rec.note
        );
    }
}
