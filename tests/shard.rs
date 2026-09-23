use pihunt::basis::Shape;
use pihunt::config;
use pihunt::plan::{Job, plan};
use pihunt::shard::Shard;
use std::path::{Path, PathBuf};

fn job() -> Job {
    Job {
        shape: Shape::new(10, 6, 1, 1, vec![]),
        coeff_bound: 1000,
        digits: 100,
        gamma: 1.16,
        max_iterations: 1_000_000,
    }
}

#[test]
fn parses_valid_specs() {
    assert_eq!("1/4".parse::<Shard>().unwrap(), Shard { k: 1, n: 4 });
    assert_eq!("4/4".parse::<Shard>().unwrap(), Shard { k: 4, n: 4 });
    assert_eq!("1/1".parse::<Shard>().unwrap(), Shard { k: 1, n: 1 });
}

#[test]
fn rejects_malformed_specs() {
    for bad in [
        "0/4", "5/4", "2/0", "abc", "1/", "/4", "1/4/5", "-1/4", "1.5/4", "",
    ] {
        assert!(bad.parse::<Shard>().is_err(), "expected error for {bad:?}");
    }
}

#[test]
fn every_job_lands_in_exactly_one_shard() {
    let batch = config::load(Path::new("batches/known.toml")).unwrap();
    let jobs = plan(&batch);
    let finder = batch.defaults.finder.finder().name();
    assert!(!jobs.is_empty());
    for n in 1..=5u32 {
        for job in &jobs {
            let owners: Vec<u32> = (1..=n)
                .filter(|&k| Shard { k, n }.owns(job, finder))
                .collect();
            assert_eq!(
                owners.len(),
                1,
                "job {:?} owned by shards {owners:?} for n={n}",
                job.shape
            );
        }
    }
}

#[test]
fn assignment_is_stable_across_calls() {
    let jobs = plan(&config::load(Path::new("batches/known.toml")).unwrap());
    let shard = Shard { k: 2, n: 3 };
    for job in &jobs {
        let first = shard.owns(job, "classic");
        let second = shard.owns(job, "classic");
        assert_eq!(first, second);
    }
}

#[test]
fn escalation_does_not_change_shard() {
    let job = job();
    let finder = "classic";
    // Escalation attempt 0 is the base job itself: hashing it is exactly "the base job's ID".
    assert_eq!(job.escalated(0), job);
    for n in 1..=5u32 {
        for k in 1..=n {
            let shard = Shard { k, n };
            assert_eq!(
                shard.owns(&job, finder),
                shard.owns(&job.escalated(0), finder)
            );
        }
    }
    // Escalated attempts have different digits, and therefore a different raw ID; sharding on
    // the base job (rather than whichever attempt happens to be running) is what keeps a whole
    // escalation chain together.
    assert_ne!(job.id(finder), job.escalated(1).id(finder));
    assert_ne!(job.id(finder), job.escalated(2).id(finder));
}

#[test]
fn output_path_inserts_shard_suffix() {
    let shard = Shard { k: 2, n: 4 };
    assert_eq!(
        shard.output_path(Path::new("results/base10-wide.jsonl")),
        PathBuf::from("results/base10-wide.shard-2-of-4.jsonl")
    );
    // No extension.
    assert_eq!(
        shard.output_path(Path::new("results/base10-wide")),
        PathBuf::from("results/base10-wide.shard-2-of-4")
    );
    // No directory.
    assert_eq!(
        shard.output_path(Path::new("out.jsonl")),
        PathBuf::from("out.shard-2-of-4.jsonl")
    );
}
