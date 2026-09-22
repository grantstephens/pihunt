use pihunt::basis::{Extra, Shape, auto_digits};
use pihunt::job::run_job;
use pihunt::log::Kind;
use pihunt::plan::Job;
use pihunt::pslq::classic::ClassicPslq;

fn job(shape: Shape) -> Job {
    let digits = auto_digits(shape.columns(), 1000);
    Job {
        shape,
        coeff_bound: 1000,
        digits,
        gamma: 1.16,
        max_iterations: 1_000_000,
    }
}

#[test]
fn rediscovers_bbp() {
    let rec = run_job(&job(Shape::new(16, 8, 1, 1, vec![])), "t", &ClassicPslq, 80);
    assert_eq!(rec.outcome, Kind::Hit);
    assert_eq!(rec.tag.as_deref(), Some("known"));
    let rel: Vec<&str> = rec
        .relation
        .as_ref()
        .unwrap()
        .iter()
        .map(|s| s.as_str())
        .collect();
    assert_eq!(rel, ["1", "-4", "0", "0", "2", "1", "1", "0", "0"]);
    assert_eq!(rec.dropped[0].column, "S(j=7,s=1)");
    assert_eq!(rec.dropped[0].tag.as_deref(), Some("known"));
}

#[test]
fn rediscovers_bailey_pi_squared_as_basis_relation() {
    let rec = run_job(
        &job(Shape::new(64, 6, 2, 2, vec![Extra::Pi2])),
        "t",
        &ClassicPslq,
        80,
    );
    let d = rec
        .dropped
        .iter()
        .find(|d| d.column == "pi2")
        .expect("pi2 dropped");
    assert_eq!(
        d.relation,
        ["0", "144", "-216", "-72", "-54", "9", "0", "-8"]
    );
    assert_eq!(d.tag.as_deref(), Some("known"));
    assert_eq!(rec.outcome, Kind::Excluded);
}

#[test]
fn skips_oversized_jobs() {
    let rec = run_job(&job(Shape::new(10, 8, 1, 1, vec![])), "t", &ClassicPslq, 5);
    assert_eq!(rec.outcome, Kind::Skipped);
    assert!(rec.note.unwrap().contains("max_columns"));
}
