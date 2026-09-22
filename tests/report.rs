use pihunt::log::{DroppedRecord, Kind, Params, Record};
use pihunt::report::render;

fn rec(base: u32, period: u32, extras: &[&str], outcome: Kind, c: u64, digits: u32) -> Record {
    let mut columns = vec!["pi".to_string()];
    columns.extend((1..=period).map(|j| format!("S(j={j},s=1)")));
    columns.extend(extras.iter().map(|e| e.to_string()));
    Record {
        job_id: format!("{base}-{period}-{}-{c}-{digits}", extras.join(",")),
        batch: "t".into(),
        pihunt_version: "0.2.0".into(),
        started: "2026-09-22T00:00:00Z".into(),
        elapsed_ms: 1,
        params: Params {
            base,
            period,
            degrees: [1, 1],
            extras: extras.iter().map(|e| e.to_string()).collect(),
            coeff_bound: c,
            precision_digits: digits,
            gamma: 1.16,
            max_iterations: 10,
            finder: "multilevel".into(),
            algo_version: 1,
        },
        columns,
        dropped: vec![],
        outcome,
        bound: Some("1.0e4".into()),
        iterations: 1,
        relation: None,
        verify: None,
        tag: None,
        note: None,
        escalated_from: None,
    }
}

fn hit(tag: &str) -> Record {
    let mut r = rec(16, 8, &[], Kind::Hit, 1000, 100);
    r.relation = Some(
        ["1", "-4", "0", "0", "2", "1", "1", "0", "0"]
            .map(String::from)
            .to_vec(),
    );
    r.tag = Some(tag.into());
    r
}

#[test]
fn empty_input() {
    assert!(render(&[]).contains("No records."));
}

#[test]
fn strongest_exclusion_per_shape() {
    let out = render(&[
        rec(10, 4, &["log2"], Kind::Excluded, 1000, 80),
        rec(10, 4, &["log2"], Kind::Excluded, 100000, 120),
        rec(10, 5, &[], Kind::Excluded, 1000, 80),
    ]);
    assert!(
        out.contains("| 10 | 4 | 1..1 | log2 | 100000 | 120 |"),
        "{out}"
    );
    assert!(!out.contains("| 10 | 4 | 1..1 | log2 | 1000 |"), "{out}");
    assert!(out.contains("| 10 | 5 | 1..1 | - | 1000 | 80 |"), "{out}");
    assert!(out.contains("Shapes excluded: 2"), "{out}");
}

#[test]
fn hits_are_listed_and_new_ones_shout() {
    let out = render(&[hit("known")]);
    assert!(out.contains("1·pi + -4·S(j=1,s=1) + 2·S(j=4,s=1) + 1·S(j=5,s=1) + 1·S(j=6,s=1) = 0"));
    assert!(!out.contains("NEW HITS"));
    let out = render(&[hit("NEW")]);
    assert!(out.contains("## 🚨 NEW HITS"), "{out}");
}

#[test]
fn unresolved_only_when_never_settled() {
    let mut stuck = rec(10, 7, &[], Kind::Inconclusive, 1000, 240);
    stuck.note = Some("reduction: precision exhausted".into());
    let early = rec(10, 6, &[], Kind::Inconclusive, 1000, 60);
    let settled = rec(10, 6, &[], Kind::Excluded, 1000, 120);
    let out = render(&[stuck, early, settled]);
    let unresolved = out.split("## Unresolved").nth(1).expect("section");
    assert!(
        unresolved.contains(
            "| 10 | 7 | 1..1 | - | inconclusive | 240 | reduction: precision exhausted |"
        ),
        "{out}"
    );
    assert!(!unresolved.contains("| 10 | 6 |"), "{out}");
}

#[test]
fn basis_relations_are_deduplicated() {
    let mut a = rec(10, 6, &["log2", "log3", "log5"], Kind::Excluded, 1000, 80);
    a.dropped = vec![DroppedRecord {
        column: "log5".into(),
        relation: ["0", "0", "0", "0", "0", "0", "3", "-5", "10", "-5"]
            .map(String::from)
            .to_vec(),
        tag: None,
    }];
    let mut b = a.clone();
    b.params.coeff_bound = 10000;
    b.job_id = "other".into();
    let out = render(&[a, b]);
    let rel = "3·S(j=6,s=1) + -5·log2 + 10·log3 + -5·log5 = 0";
    assert_eq!(out.matches(rel).count(), 1, "{out}");
    assert!(
        out.contains("| 10 | 6 | 1..1 | log2,log3,log5 | 10000 | 80 | log5 |"),
        "{out}"
    );
}
