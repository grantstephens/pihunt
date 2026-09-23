//! CLI tests for the `digit`/`stream` subcommands (Theorem 1 and Theorem 2), covering both
//! `nthdigit` and `nthdigit2` end-to-end through the `pihunt` binary. Moved out of
//! `crates/pi-digits/tests/{nthdigit,nthdigit2}.rs`: those live in the `pi-digits` crate, which
//! has no `pihunt` binary of its own, so `env!("CARGO_BIN_EXE_pihunt")` only resolves here.

#[test]
fn cli_digit_prints_expected_digits() {
    let bin = env!("CARGO_BIN_EXE_pihunt");
    let out = std::process::Command::new(bin)
        .args(["digit", "1", "--count", "5"])
        .output()
        .unwrap();
    assert!(out.status.success());
    assert_eq!(String::from_utf8_lossy(&out.stdout).trim(), "14159");
}

#[test]
fn cli_stream_prints_independent_blocks() {
    let bin = env!("CARGO_BIN_EXE_pihunt");
    let out = std::process::Command::new(bin)
        .args(["stream", "--from", "1", "--block", "5", "--blocks", "2"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines, vec!["14159", "26535"]);
}

#[test]
fn cli_digit_thm2_matches_thm1() {
    let bin = env!("CARGO_BIN_EXE_pihunt");
    let out1 = std::process::Command::new(bin)
        .args(["digit", "1", "--count", "20", "--method", "thm1"])
        .output()
        .unwrap();
    let out2 = std::process::Command::new(bin)
        .args([
            "digit", "1", "--count", "20", "--method", "thm2", "--mem", "512",
        ])
        .output()
        .unwrap();
    assert!(out1.status.success());
    assert!(
        out2.status.success(),
        "{}",
        String::from_utf8_lossy(&out2.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&out1.stdout).trim(),
        String::from_utf8_lossy(&out2.stdout).trim()
    );
    assert_eq!(
        String::from_utf8_lossy(&out2.stdout).trim(),
        "14159265358979323846"
    );
}

#[test]
fn cli_stream_thm2_prints_independent_blocks() {
    let bin = env!("CARGO_BIN_EXE_pihunt");
    let out = std::process::Command::new(bin)
        .args([
            "stream", "--from", "1", "--block", "5", "--blocks", "2", "--method", "thm2", "--mem",
            "512",
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stdout = String::from_utf8_lossy(&out.stdout);
    let lines: Vec<&str> = stdout.lines().collect();
    assert_eq!(lines, vec!["14159", "26535"]);
}
