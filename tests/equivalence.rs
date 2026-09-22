//! Multilevel PSLQ must reach the same verdict as classic on every job in the corpus.
//! The full scout corpus is slow in debug builds, so it's ignored by default:
//! cargo test --release --test equivalence -- --ignored

use pihunt::config;
use pihunt::job::run_job;
use pihunt::log::Record;
use pihunt::plan::plan;
use pihunt::pslq::{classic::ClassicPslq, multilevel::MultilevelPslq};
use rayon::prelude::*;
use std::path::Path;

/// The parts of a record that must match; timings, iterations and bounds may differ.
fn verdict(r: &Record) -> impl PartialEq + std::fmt::Debug + use<> {
    let dropped: Vec<(String, Vec<String>)> = r
        .dropped
        .iter()
        .map(|d| (d.column.clone(), d.relation.clone()))
        .collect();
    (r.outcome, r.relation.clone(), r.tag.clone(), dropped)
}

fn check_batch(path: &str) {
    let batch = config::load(Path::new(path)).unwrap();
    let jobs = plan(&batch);
    let mismatches: Vec<String> = jobs
        .par_iter()
        .filter_map(|job| {
            let max = batch.defaults.max_columns;
            let c = run_job(job, "eq", &ClassicPslq, max);
            let m = run_job(job, "eq", &MultilevelPslq, max);
            (verdict(&c) != verdict(&m)).then(|| {
                format!(
                    "{:?}\n  classic:    {:?}\n  multilevel: {:?}",
                    job.shape,
                    verdict(&c),
                    verdict(&m)
                )
            })
        })
        .collect();
    assert!(
        mismatches.is_empty(),
        "{} of {} jobs differ:\n{}",
        mismatches.len(),
        jobs.len(),
        mismatches.join("\n")
    );
}

#[test]
fn known_batch_matches() {
    check_batch("batches/known.toml");
}

#[test]
#[ignore]
fn scout_batch_matches() {
    check_batch("batches/base10-scout.toml");
}
