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
