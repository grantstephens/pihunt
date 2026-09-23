use std::collections::HashSet;
use std::process::Command;

const BATCH: &str = r#"
name   = "cli"
output = "OUT"

[defaults]
coeff_bound = 1000

[search]
mode    = "grid"
bases   = [16]
periods = { from = 7, to = 8 }
degrees = [[1, 1]]
extras  = [[]]
"#;

#[test]
fn run_then_verify() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("r.jsonl");
    let batch = dir.path().join("b.toml");
    std::fs::write(&batch, BATCH.replace("OUT", out.to_str().unwrap())).unwrap();
    let bin = env!("CARGO_BIN_EXE_pihunt");

    let plan = Command::new(bin)
        .args(["plan", batch.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(plan.status.success());
    assert!(String::from_utf8_lossy(&plan.stdout).contains("2 jobs"));

    let run = Command::new(bin)
        .args(["run", batch.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let stdout = String::from_utf8_lossy(&run.stdout);
    assert!(stdout.contains("hit [known] b=16 m=8"), "{stdout}");
    assert_eq!(std::fs::read_to_string(&out).unwrap().lines().count(), 2);

    let verify = Command::new(bin)
        .args(["verify", out.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(verify.status.success());
    assert!(String::from_utf8_lossy(&verify.stdout).contains("1 hits checked, 0 failed"));
}

#[test]
fn bad_batch_fails_cleanly() {
    let dir = tempfile::tempdir().unwrap();
    let batch = dir.path().join("b.toml");
    std::fs::write(&batch, "name = 1").unwrap();
    let run = Command::new(env!("CARGO_BIN_EXE_pihunt"))
        .args(["run", batch.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!run.status.success());
    assert!(String::from_utf8_lossy(&run.stderr).starts_with("error:"));
}

#[test]
fn rerun_resumes_instead_of_repeating() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("r.jsonl");
    let batch = dir.path().join("b.toml");
    let text = BATCH.replace("OUT", out.to_str().unwrap()).replace(
        "coeff_bound = 1000",
        "coeff_bound = 1000\nfinder = \"multilevel\"",
    );
    std::fs::write(&batch, text).unwrap();
    let bin = env!("CARGO_BIN_EXE_pihunt");
    let run = || {
        Command::new(bin)
            .args(["run", batch.to_str().unwrap()])
            .output()
            .unwrap()
    };

    assert!(run().status.success());
    let first = std::fs::read_to_string(&out).unwrap();
    let again = run();
    assert!(again.status.success());
    assert_eq!(
        std::fs::read_to_string(&out).unwrap(),
        first,
        "resume must not append"
    );
    let stdout = String::from_utf8_lossy(&again.stdout);
    assert!(
        stdout.contains("resumed 2 attempts from the log"),
        "{stdout}"
    );
    assert!(stdout.contains("hit [known] b=16 m=8"), "{stdout}");
}

#[test]
fn report_reads_logs() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("r.jsonl");
    let batch = dir.path().join("b.toml");
    std::fs::write(&batch, BATCH.replace("OUT", out.to_str().unwrap())).unwrap();
    let bin = env!("CARGO_BIN_EXE_pihunt");
    assert!(
        Command::new(bin)
            .args(["run", batch.to_str().unwrap()])
            .output()
            .unwrap()
            .status
            .success()
    );

    let report = Command::new(bin)
        .args(["report", out.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(report.status.success());
    let md = String::from_utf8_lossy(&report.stdout);
    assert!(md.starts_with("# pihunt exclusion report"), "{md}");
    assert!(md.contains("| 16 | 7 | 1..1 | - | 1000 |"), "{md}");

    let missing = Command::new(bin)
        .args(["report", "nope.jsonl"])
        .output()
        .unwrap();
    assert!(!missing.status.success());
}

fn job_ids(path: &std::path::Path) -> HashSet<String> {
    pihunt::log::read(path)
        .unwrap()
        .into_iter()
        .map(|r| r.job_id)
        .collect()
}

#[test]
fn sharding_splits_the_batch_without_dropping_or_duplicating_jobs() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("r.jsonl");
    let batch = dir.path().join("b.toml");
    std::fs::write(&batch, BATCH.replace("OUT", out.to_str().unwrap())).unwrap();
    let bin = env!("CARGO_BIN_EXE_pihunt");

    let unsharded = Command::new(bin)
        .args(["run", batch.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(unsharded.status.success());
    let unsharded_ids = job_ids(&out);

    let shard1 = out.with_file_name("r.shard-1-of-2.jsonl");
    let shard2 = out.with_file_name("r.shard-2-of-2.jsonl");
    for (k, path) in [(1, &shard1), (2, &shard2)] {
        assert!(!path.exists());
        let run = Command::new(bin)
            .args(["run", batch.to_str().unwrap(), "--shard", &format!("{k}/2")])
            .output()
            .unwrap();
        assert!(
            run.status.success(),
            "{}",
            String::from_utf8_lossy(&run.stderr)
        );
        assert_ne!(path, &out);
    }

    let combined: HashSet<String> = job_ids(&shard1).union(&job_ids(&shard2)).cloned().collect();
    assert_eq!(combined, unsharded_ids);
    // The two shards must not overlap.
    assert_eq!(job_ids(&shard1).intersection(&job_ids(&shard2)).count(), 0);
}

#[test]
fn plan_reports_shard_job_count() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("r.jsonl");
    let batch = dir.path().join("b.toml");
    std::fs::write(&batch, BATCH.replace("OUT", out.to_str().unwrap())).unwrap();
    let plan = Command::new(env!("CARGO_BIN_EXE_pihunt"))
        .args(["plan", batch.to_str().unwrap(), "--shard", "1/2"])
        .output()
        .unwrap();
    assert!(plan.status.success());
    assert!(
        String::from_utf8_lossy(&plan.stdout).contains("shard 1/2:"),
        "{}",
        String::from_utf8_lossy(&plan.stdout)
    );
}

#[test]
fn bad_shard_spec_fails_cleanly() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("r.jsonl");
    let batch = dir.path().join("b.toml");
    std::fs::write(&batch, BATCH.replace("OUT", out.to_str().unwrap())).unwrap();
    for spec in ["0/4", "5/4", "2/0", "abc"] {
        let run = Command::new(env!("CARGO_BIN_EXE_pihunt"))
            .args(["run", batch.to_str().unwrap(), "--shard", spec])
            .output()
            .unwrap();
        assert!(!run.status.success(), "{spec}");
        assert!(
            String::from_utf8_lossy(&run.stderr).starts_with("error:"),
            "{spec}: {}",
            String::from_utf8_lossy(&run.stderr)
        );
    }
}
