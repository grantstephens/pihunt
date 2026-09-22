use pihunt::log::{Kind, Params, Record, read, spawn_writer};

fn record(id: &str) -> Record {
    Record {
        job_id: id.into(),
        batch: "t".into(),
        pihunt_version: "0.1.0".into(),
        started: "2026-09-22T00:00:00Z".into(),
        elapsed_ms: 5,
        params: Params {
            base: 10,
            period: 2,
            degrees: [1, 1],
            extras: vec!["log2".into()],
            coeff_bound: 1000,
            precision_digits: 80,
            gamma: 1.16,
            max_iterations: 10,
            finder: "classic".into(),
            algo_version: 1,
        },
        columns: vec!["pi".into()],
        dropped: vec![],
        outcome: Kind::Excluded,
        bound: Some("1.000e4".into()),
        iterations: 3,
        relation: None,
        verify: None,
        tag: None,
        note: None,
    }
}

#[test]
fn writes_appends_and_reads_back() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("sub/out.jsonl");
    for ids in [["a", "b"], ["c", "d"]] {
        let (tx, h) = spawn_writer(&path).unwrap();
        for id in ids {
            tx.send(record(id)).unwrap();
        }
        drop(tx);
        h.join().unwrap().unwrap();
    }
    let back = read(&path).unwrap();
    assert_eq!(
        back.iter().map(|r| r.job_id.as_str()).collect::<Vec<_>>(),
        ["a", "b", "c", "d"]
    );
    assert_eq!(back[0], record("a"));
    let line = std::fs::read_to_string(&path).unwrap();
    assert!(line.contains("\"outcome\":\"excluded\""));
}
