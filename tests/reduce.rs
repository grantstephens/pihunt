use pihunt::basis::{Columns, Extra, Shape};
use pihunt::pslq::{PslqParams, classic::ClassicPslq};
use pihunt::reduce::{Reduced, reduce};
use rug::Integer;

fn params(digits: u32) -> PslqParams {
    PslqParams {
        gamma: 1.16,
        coeff_bound: 1000,
        max_iterations: 100_000,
        digits,
    }
}

fn ints(v: &[i64]) -> Vec<Integer> {
    v.iter().map(|&c| Integer::from(c)).collect()
}

#[test]
fn drops_planted_duplicate() {
    // Base 2, period 1: S(j=1,s=1) = 2 log 2, so log2 is redundant.
    let cols = Columns::build(&Shape::new(2, 1, 1, 1, vec![Extra::Log2]), 80);
    match reduce(&cols, &ClassicPslq, &params(80)) {
        Reduced::Done { keep, dropped } => {
            assert_eq!(keep, vec![1]);
            assert_eq!(dropped.len(), 1);
            assert_eq!(dropped[0].column, 2);
            assert_eq!(dropped[0].relation, ints(&[0, 1, -2]));
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn drops_bbp_zero_relation() {
    let cols = Columns::build(&Shape::new(16, 8, 1, 1, vec![]), 100);
    match reduce(&cols, &ClassicPslq, &params(100)) {
        Reduced::Done { keep, dropped } => {
            assert_eq!(keep, vec![1, 2, 3, 4, 5, 6, 8]);
            assert_eq!(dropped[0].relation, ints(&[0, 8, -8, -4, -8, -2, -2, 1, 0]));
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn independent_columns_untouched() {
    let cols = Columns::build(
        &Shape::new(10, 1, 1, 1, vec![Extra::Catalan, Extra::Zeta3]),
        80,
    );
    match reduce(&cols, &ClassicPslq, &params(80)) {
        Reduced::Done { keep, dropped } => {
            assert_eq!(keep, vec![1, 2, 3]);
            assert!(dropped.is_empty());
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn oversized_relation_is_inconclusive_not_dropped() {
    // Same duplicate as above, but a bound of 1 makes the true relation (1, -2) too big to trust.
    let cols = Columns::build(&Shape::new(2, 1, 1, 1, vec![Extra::Log2]), 80);
    let p = PslqParams {
        coeff_bound: 1,
        ..params(80)
    };
    assert!(matches!(
        reduce(&cols, &ClassicPslq, &p),
        Reduced::Inconclusive { .. }
    ));
}
