//! Runs one job end to end: basis → reduce → main PSLQ → classify → Record.

use crate::basis::Columns;
use crate::classify::classify;
use crate::known;
use crate::log::{DroppedRecord, Kind, Params, Record};
use crate::plan::{ALGO_VERSION, Job};
use crate::pslq::{Outcome, RelationFinder};
use crate::reduce::{Dropped, Reduced, reduce};
use rug::{Float, Integer};
use std::time::Instant;

pub fn run_job(job: &Job, batch: &str, finder: &dyn RelationFinder, max_columns: usize) -> Record {
    let started = jiff::Timestamp::now().to_string();
    let t = Instant::now();
    let names = job.shape.column_names();
    let mut rec = Record {
        job_id: job.id(finder.name()),
        batch: batch.to_string(),
        pihunt_version: env!("CARGO_PKG_VERSION").to_string(),
        started,
        elapsed_ms: 0,
        params: params(job, finder),
        columns: names.clone(),
        dropped: vec![],
        outcome: Kind::Skipped,
        bound: None,
        iterations: 0,
        relation: None,
        verify: None,
        tag: None,
        note: None,
        escalated_from: None,
    };
    if names.len() > max_columns {
        rec.note = Some(format!(
            "{} columns > max_columns {max_columns}",
            names.len()
        ));
        return rec;
    }

    let cols = Columns::build(&job.shape, job.digits);
    let p = job.params();
    let (keep, dropped) = match reduce(&cols, finder, &p) {
        Reduced::Done { keep, dropped } => (keep, dropped),
        Reduced::Inconclusive { outcome, dropped } => {
            rec.dropped = dropped_records(job, &names, &dropped);
            rec.outcome = Kind::Inconclusive;
            rec.bound = Some(fmt_bound(outcome.bound()));
            rec.iterations = outcome.iterations();
            rec.note = Some(match &outcome {
                Outcome::Relation { coeffs, .. } => {
                    let max = coeffs
                        .iter()
                        .map(|c| c.as_abs().to_string())
                        .max_by_key(|s| (s.len(), s.clone()));
                    format!(
                        "reduction found a relation that failed checks (max |coeff| {})",
                        max.unwrap_or_default()
                    )
                }
                Outcome::PrecisionExhausted { .. } => "reduction: precision exhausted".into(),
                _ => "reduction: iteration cap".into(),
            });
            rec.elapsed_ms = t.elapsed().as_millis() as u64;
            return rec;
        }
    };
    rec.dropped = dropped_records(job, &names, &dropped);

    let mut idx = vec![0];
    idx.extend(keep);
    let x: Vec<Float> = idx.iter().map(|&i| cols.lo[i].clone()).collect();
    let outcome = finder.find(&x, &p);
    let relation = match &outcome {
        Outcome::Relation { coeffs, .. } => {
            let mut full = vec![Integer::new(); names.len()];
            for (&i, c) in idx.iter().zip(coeffs) {
                full[i] = c.clone();
            }
            Some(full)
        }
        _ => None,
    };
    let v = classify(&job.shape, &cols, &outcome, relation, job.coeff_bound);
    rec.outcome = v.kind;
    rec.bound = Some(fmt_bound(outcome.bound()));
    rec.iterations = outcome.iterations();
    rec.relation = v
        .relation
        .map(|r| r.iter().map(|c| c.to_string()).collect());
    rec.verify = v.verify;
    rec.tag = v.tag;
    rec.note = v.note;
    rec.elapsed_ms = t.elapsed().as_millis() as u64;
    rec
}

fn params(job: &Job, finder: &dyn RelationFinder) -> Params {
    let s = &job.shape;
    Params {
        base: s.base,
        period: s.period,
        degrees: [s.s_lo, s.s_hi],
        extras: s.extras.iter().map(|e| e.name().to_string()).collect(),
        coeff_bound: job.coeff_bound,
        precision_digits: job.digits,
        gamma: job.gamma,
        max_iterations: job.max_iterations,
        finder: finder.name().to_string(),
        algo_version: ALGO_VERSION,
    }
}

fn dropped_records(job: &Job, names: &[String], dropped: &[Dropped]) -> Vec<DroppedRecord> {
    dropped
        .iter()
        .map(|d| DroppedRecord {
            column: names[d.column].clone(),
            relation: d.relation.iter().map(|c| c.to_string()).collect(),
            tag: known::is_known(&job.shape, names, &d.relation).then(|| "known".to_string()),
        })
        .collect()
}

pub fn fmt_bound(b: f64) -> String {
    format!("{b:.3e}")
}
