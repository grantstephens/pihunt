use pihunt::basis::{Columns, Shape};
use pihunt::classify::classify;
use pihunt::log::Kind;
use pihunt::pslq::Outcome;
use rug::Integer;

fn ints(v: &[i64]) -> Vec<Integer> {
    v.iter().map(|&c| Integer::from(c)).collect()
}

fn setup() -> (Shape, Columns) {
    let shape = Shape::new(16, 8, 1, 1, vec![]);
    let cols = Columns::build(&shape, 100);
    (shape, cols)
}

fn rel_outcome(v: &[i64]) -> Outcome {
    Outcome::Relation {
        coeffs: ints(v),
        iterations: 1,
        bound: 1.0,
    }
}

#[test]
fn non_relations() {
    let (shape, cols) = setup();
    let ex = Outcome::Excluded {
        bound: 5e4,
        iterations: 9,
    };
    assert_eq!(
        classify(&shape, &cols, &ex, None, 1000).kind,
        Kind::Excluded
    );
    let pe = Outcome::PrecisionExhausted {
        bound: 5.0,
        iterations: 9,
    };
    assert_eq!(
        classify(&shape, &cols, &pe, None, 1000).kind,
        Kind::Inconclusive
    );
    let ic = Outcome::IterationCap {
        bound: 5.0,
        iterations: 9,
    };
    assert_eq!(
        classify(&shape, &cols, &ic, None, 1000).kind,
        Kind::Inconclusive
    );
}

#[test]
fn bbp_is_known_hit_with_positive_pi() {
    let (shape, cols) = setup();
    let neg = [-1, 4, 0, 0, -2, -1, -1, 0, 0];
    let v = classify(&shape, &cols, &rel_outcome(&neg), Some(ints(&neg)), 1000);
    assert_eq!(v.kind, Kind::Hit);
    assert_eq!(v.tag.as_deref(), Some("known"));
    assert_eq!(v.relation, Some(ints(&[1, -4, 0, 0, 2, 1, 1, 0, 0])));
    assert!(v.verify.unwrap().passed);
}

#[test]
fn relation_without_pi_is_junk() {
    let (shape, cols) = setup();
    let r = [0, 8, -8, -4, -8, -2, -2, 1, 0];
    assert_eq!(
        classify(&shape, &cols, &rel_outcome(&r), Some(ints(&r)), 1000).kind,
        Kind::Junk
    );
}

#[test]
fn oversized_is_suspicious() {
    let (shape, cols) = setup();
    let r = [1, -4, 0, 0, 2, 1, 1, 0, 0];
    assert_eq!(
        classify(&shape, &cols, &rel_outcome(&r), Some(ints(&r)), 3).kind,
        Kind::Suspicious
    );
}

#[test]
fn false_relation_is_spurious() {
    let (shape, cols) = setup();
    let r = [1, -4, 0, 0, 2, 1, 2, 0, 0];
    let v = classify(&shape, &cols, &rel_outcome(&r), Some(ints(&r)), 1000);
    assert_eq!(v.kind, Kind::Spurious);
    assert!(!v.verify.unwrap().passed);
    assert_eq!(v.tag, None);
}
