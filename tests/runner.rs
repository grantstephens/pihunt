use pihunt::basis::{Extra, Shape};
use pihunt::log::{Kind, Record};
use pihunt::plan::Job;
use pihunt::pslq::{RelationFinder, multilevel::MultilevelPslq};
use pihunt::runner::{Ctx, load_done, run_chain};
use std::collections::HashMap;

fn ctx(escalate: u32) -> Ctx<'static> {
    Ctx {
        batch: "t",
        finder: &MultilevelPslq,
        max_columns: 80,
        escalate,
    }
}

/// Inconclusive at 60 digits (reduction relation fails checks), excluded at 120.
fn needs_escalation() -> Job {
    let all = vec![
        Extra::Pi2,
        Extra::Log2,
        Extra::Log3,
        Extra::Log5,
        Extra::Catalan,
        Extra::Zeta3,
    ];
    Job {
        shape: Shape::new(10, 6, 1, 1, all),
        coeff_bound: 1000,
        digits: 60,
        gamma: 1.16,
        max_iterations: 1_000_000,
    }
}

fn bbp() -> Job {
    Job {
        shape: Shape::new(16, 8, 1, 1, vec![]),
        coeff_bound: 1000,
        digits: 100,
        gamma: 1.16,
        max_iterations: 1_000_000,
    }
}

fn run(job: &Job, c: &Ctx, done: &HashMap<String, Record>) -> (Record, Vec<Record>, usize) {
    let mut logged = vec![];
    let chain = run_chain(job, c, done, &mut |r| logged.push(r.clone()));
    (chain.last, logged, chain.resumed)
}

#[test]
fn escalates_inconclusive_until_resolved() {
    let job = needs_escalation();
    let (last, logged, resumed) = run(&job, &ctx(2), &HashMap::new());
    assert_eq!(resumed, 0);
    assert_eq!(logged.len(), 2);
    assert_eq!(logged[0].outcome, Kind::Inconclusive);
    assert_eq!(logged[0].escalated_from, None);
    assert_eq!(logged[1].params.precision_digits, 120);
    assert_eq!(logged[1].escalated_from.as_ref(), Some(&logged[0].job_id));
    assert_eq!(last, logged[1]);
    assert_eq!(last.outcome, Kind::Excluded);
}

#[test]
fn escalation_can_be_disabled() {
    let (last, logged, _) = run(&needs_escalation(), &ctx(0), &HashMap::new());
    assert_eq!(logged.len(), 1);
    assert_eq!(last.outcome, Kind::Inconclusive);
}

#[test]
fn resolved_jobs_are_not_escalated() {
    let (last, logged, _) = run(&bbp(), &ctx(2), &HashMap::new());
    assert_eq!(logged.len(), 1);
    assert_eq!(last.outcome, Kind::Hit);
}

#[test]
fn skipped_jobs_are_not_escalated() {
    let c = Ctx {
        max_columns: 5,
        ..ctx(2)
    };
    let (last, logged, _) = run(&bbp(), &c, &HashMap::new());
    assert_eq!(logged.len(), 1);
    assert_eq!(last.outcome, Kind::Skipped);
}

#[test]
fn resume_reuses_logged_attempts() {
    let job = needs_escalation();
    let (_, first, _) = run(&job, &ctx(2), &HashMap::new());
    let done: HashMap<String, Record> = first
        .iter()
        .map(|r| (r.job_id.clone(), r.clone()))
        .collect();
    let (last, logged, resumed) = run(&job, &ctx(2), &done);
    assert!(logged.is_empty(), "nothing should be re-run");
    assert_eq!(resumed, 2);
    assert_eq!(last, first[1]);

    // Only the first attempt logged (e.g. killed mid-chain): resume runs just the second.
    let partial: HashMap<String, Record> = [(first[0].job_id.clone(), first[0].clone())].into();
    let (last, logged, resumed) = run(&job, &ctx(2), &partial);
    assert_eq!((logged.len(), resumed), (1, 1));
    assert_eq!(last.outcome, Kind::Excluded);
    assert_eq!(logged[0].escalated_from.as_ref(), Some(&first[0].job_id));
}

#[test]
fn job_ids_depend_on_finder() {
    let job = bbp();
    assert_ne!(job.id(MultilevelPslq.name()), job.id("classic"));
}

#[test]
fn load_done_handles_missing_and_existing_logs() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("r.jsonl");
    assert!(load_done(&path).unwrap().is_empty());
    let (_, logged, _) = run(&bbp(), &ctx(0), &HashMap::new());
    let (tx, h) = pihunt::log::spawn_writer(&path).unwrap();
    tx.send(logged[0].clone()).unwrap();
    drop(tx);
    h.join().unwrap().unwrap();
    let done = load_done(&path).unwrap();
    assert_eq!(done.get(&logged[0].job_id), Some(&logged[0]));
}

#[test]
fn precision_limited_outcomes_escalate() {
    use pihunt::runner::needs_more_precision;
    for k in [Kind::Inconclusive, Kind::Suspicious, Kind::Spurious] {
        assert!(needs_more_precision(k), "{k:?} should escalate");
    }
    for k in [Kind::Hit, Kind::Excluded, Kind::Skipped, Kind::Junk] {
        assert!(!needs_more_precision(k), "{k:?} is settled");
    }
}
