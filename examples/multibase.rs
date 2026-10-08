//! Does π need columns from several different base-2^k towers in the same relation?
//!
//! Every search so far (`batches/*.toml`, `src/basis.rs::Shape`) fixes one base per relation.
//! BBP and Bellard are each just "π = a Q-linear combination of several weight-1 periods" at a
//! *single* base — nothing about that requires the periods to share a base. Mixing columns from
//! several different clean bases (see `batches/base2-degree1.toml` — the 12 bases that aren't a
//! derived corollary of anything in `known::KNOWN`) is still weight-1 + weight-1 + ..., so the
//! weight-conservation conjecture that ruled out mixing degree/extras doesn't apply here. It's a
//! real, untried generalization, not a stretch — see docs/findings-bbp-hunt.md ("multibase").
//!
//! `ARITY` is how many bases combine in one relation. Bigger isn't free: columns = 1 + ARITY*m,
//! and `basis::auto_digits`'s precision factor grows with column count, so going past ARITY=3
//! starts pushing past n ~ 100 (the range `docs/precision-rule.md` actually measured) for no
//! principled reason — the weight argument doesn't favour more bases over fewer, so higher
//! arity is just more dice rolled in a regime the precision formula isn't validated for, not a
//! more promising search. See docs/findings-bbp-hunt.md for the arity-7-vs-3 reasoning.
//!
//! Standalone rather than wired into `Shape`/the batch-TOML machinery, which is one-base-per-shape
//! throughout (`Shape::new(base: u32, ...)`) — retrofitting that would touch dozens of call
//! sites for what's currently a one-off experiment. This reuses `basis::series` (the same
//! per-column series the real pipeline builds) and the multilevel PSLQ engine directly.
//!
//! Deliberately excludes any base already known to carry a formula (16, 64, 2^20, and their
//! `known::is_derived` descendants): grouping a known base with anything would trivially
//! "succeed" by rediscovering it with zero coefficients on the other bases, telling us nothing.
//!
//! Run: `cargo run --release --example multibase`

use pihunt::basis::{auto_digits, series};
use pihunt::pslq::multilevel::MultilevelPslq;
use pihunt::pslq::{Outcome, PslqParams, RelationFinder, digits_to_bits};
use pihunt::verify;
use rayon::prelude::*;
use rug::float::Constant;
use rug::{Float, Integer};
use std::fs::OpenOptions;
use std::io::Write;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::sync::mpsc::channel;
use std::time::Instant;

const CLEAN_BASES: [u32; 12] = [
    32, 128, 512, 2048, 8192, 16384, 32768, 131072, 524288, 2097152, 4194304, 8388608,
];
const ARITY: usize = 3;
const PERIODS: std::ops::RangeInclusive<u32> = 2..=24;
const COEFF_BOUND: u64 = 10_000_000;
const GAMMA: f64 = 1.16;
const MAX_ITERATIONS: u64 = 1_000_000;
const OUTPUT: &str = "results/multibase3.jsonl";

/// All `k`-element subsets of `items`, as sorted `Vec<u32>`s.
fn combinations(items: &[u32], k: usize) -> Vec<Vec<u32>> {
    if k == 0 {
        return vec![vec![]];
    }
    let Some((&first, rest)) = items.split_first() else {
        return vec![];
    };
    let mut out = combinations(rest, k - 1);
    for c in &mut out {
        c.insert(0, first);
    }
    out.extend(combinations(rest, k));
    out
}

fn main() {
    let groups = combinations(&CLEAN_BASES, ARITY);
    let jobs: Vec<(&Vec<u32>, u32)> = groups
        .iter()
        .flat_map(|g| PERIODS.map(move |m| (g, m)))
        .collect();
    let total = jobs.len();
    eprintln!(
        "multibase: {} base {ARITY}-tuples x periods {:?} = {total} jobs → {OUTPUT}",
        groups.len(),
        PERIODS
    );

    if let Some(dir) = std::path::Path::new(OUTPUT).parent() {
        std::fs::create_dir_all(dir).expect("create results dir");
    }
    let (tx, rx) = channel::<serde_json::Value>();
    let writer = std::thread::spawn(move || {
        let mut file = OpenOptions::new()
            .create(true)
            .append(true)
            .open(OUTPUT)
            .expect("open output");
        for v in rx {
            writeln!(file, "{v}").expect("write result");
        }
    });

    let finished = AtomicUsize::new(0);
    let excluded = AtomicUsize::new(0);
    let hits = AtomicUsize::new(0);
    let basis_only = AtomicUsize::new(0);
    let unresolved = AtomicUsize::new(0);

    jobs.par_iter().for_each_with(tx, |tx, &(bases, m)| {
        let n = 1 + ARITY * m as usize;
        let digits = auto_digits(n, COEFF_BOUND);
        let hi: Vec<Float> = columns(bases, m, digits_to_bits(2 * digits));
        let lo: Vec<Float> = hi
            .iter()
            .map(|v| Float::with_val(digits_to_bits(digits), v))
            .collect();
        let params = PslqParams {
            gamma: GAMMA,
            coeff_bound: COEFF_BOUND,
            max_iterations: MAX_ITERATIONS,
            digits,
        };

        let t = Instant::now();
        let outcome = MultilevelPslq.find(&lo, &params);
        let elapsed_ms = t.elapsed().as_millis() as u64;

        let (kind, relation, verified) = match &outcome {
            Outcome::Relation { coeffs, .. } if coeffs[0] == 0 => {
                basis_only.fetch_add(1, Ordering::Relaxed);
                ("basis_relation_no_pi", Some(coeffs.clone()), None)
            }
            Outcome::Relation { coeffs, .. } => {
                hits.fetch_add(1, Ordering::Relaxed);
                let (passed, residual_log10) = verify::passes(coeffs, &hi, digits);
                ("hit", Some(coeffs.clone()), Some((passed, residual_log10)))
            }
            Outcome::Excluded { .. } => {
                excluded.fetch_add(1, Ordering::Relaxed);
                ("excluded", None, None)
            }
            _ => {
                unresolved.fetch_add(1, Ordering::Relaxed);
                ("inconclusive", None, None)
            }
        };

        let names = column_names(bases, m);
        let rec = serde_json::json!({
            "bases": bases, "period": m, "columns": names,
            "digits": digits, "coeff_bound": COEFF_BOUND,
            "outcome": kind,
            "bound": format!("{:.3e}", outcome.bound()),
            "closeness": outcome.closeness(),
            "iterations": outcome.iterations(),
            "elapsed_ms": elapsed_ms,
            "relation": relation.map(|r| r.iter().map(Integer::to_string).collect::<Vec<_>>()),
            "verify": verified.map(|(passed, residual_log10)| serde_json::json!({
                "passed": passed, "residual_log10": residual_log10,
            })),
        });
        tx.send(rec).expect("writer thread alive");

        let i = finished.fetch_add(1, Ordering::Relaxed) + 1;
        let bases_str = bases
            .iter()
            .map(u32::to_string)
            .collect::<Vec<_>>()
            .join(",");
        eprintln!(
            "[{i}/{total}] bases=[{bases_str}] m={m} → {kind} ({elapsed_ms} ms, {digits} digits)"
        );
    });

    // `tx` was moved into `for_each_with` above; rayon clones it per split and drops each
    // clone when that split finishes, so once the parallel iteration completes, the last
    // clone drops, the channel closes, and the writer thread's `for v in rx` loop ends.
    writer.join().expect("writer thread panicked");

    eprintln!("\n== summary ==");
    eprintln!("Excluded: {}", excluded.load(Ordering::Relaxed));
    eprintln!("Hit: {}", hits.load(Ordering::Relaxed));
    eprintln!(
        "Basis relation (no π): {}",
        basis_only.load(Ordering::Relaxed)
    );
    eprintln!("Inconclusive: {}", unresolved.load(Ordering::Relaxed));
}

fn columns(bases: &[u32], m: u32, bits: u32) -> Vec<Float> {
    let mut cols = vec![Float::with_val(bits, Constant::Pi)];
    for &b in bases {
        cols.extend((1..=m).map(|j| series(b, m, j, 1, bits)));
    }
    cols
}

fn column_names(bases: &[u32], m: u32) -> Vec<String> {
    let mut names = vec!["pi".to_string()];
    for &b in bases {
        names.extend((1..=m).map(|j| format!("S(b={b},j={j},s=1)")));
    }
    names
}
