//! Runs a job as a chain of attempts: resume from the log, escalate precision when inconclusive.

use crate::job::run_job;
use crate::log::{self, Kind, Record};
use crate::plan::Job;
use crate::pslq::RelationFinder;
use std::collections::HashMap;
use std::path::Path;

/// Settings shared by every job in a batch.
#[derive(Clone, Copy)]
pub struct Ctx<'a> {
    pub batch: &'a str,
    pub finder: &'a dyn RelationFinder,
    pub max_columns: usize,
    /// Retries at 2x, 4x, ... digits for inconclusive attempts.
    pub escalate: u32,
}

pub struct Chain {
    /// The attempt that settled the job (or the last one tried).
    pub last: Record,
    /// Attempts taken from the existing log instead of being re-run.
    pub resumed: usize,
}

/// Attempt `job` at its own digits, then escalate while inconclusive. Attempts already in
/// `done` (keyed by job ID) are reused, not re-run. Every freshly run attempt goes to
/// `on_record` as soon as it finishes, so a crash loses at most the attempt in flight.
pub fn run_chain(
    job: &Job,
    ctx: &Ctx,
    done: &HashMap<String, Record>,
    on_record: &mut dyn FnMut(&Record),
) -> Chain {
    let mut resumed = 0;
    let mut prev: Option<String> = None;
    let mut last = None;
    for k in 0..=ctx.escalate {
        let attempt = job.escalated(k);
        let rec = match done.get(&attempt.id(ctx.finder.name())) {
            Some(rec) => {
                resumed += 1;
                rec.clone()
            }
            None => {
                let mut rec = run_job(&attempt, ctx.batch, ctx.finder, ctx.max_columns);
                rec.escalated_from = prev.clone();
                on_record(&rec);
                rec
            }
        };
        let settled = rec.outcome != Kind::Inconclusive;
        prev = Some(rec.job_id.clone());
        last = Some(rec);
        if settled {
            break;
        }
    }
    Chain {
        last: last.expect("at least one attempt"),
        resumed,
    }
}

/// Every record already in `path`, keyed by job ID. A missing file means a fresh run.
pub fn load_done(path: &Path) -> Result<HashMap<String, Record>, String> {
    if !path.exists() {
        return Ok(HashMap::new());
    }
    Ok(log::read(path)?
        .into_iter()
        .map(|r| (r.job_id.clone(), r))
        .collect())
}
