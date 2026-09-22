use pihunt::basis::{Columns, Shape};
use pihunt::verify::passes;
use rug::Integer;

fn ints(v: &[i64]) -> Vec<Integer> {
    v.iter().map(|&c| Integer::from(c)).collect()
}

#[test]
fn true_relation_passes() {
    let cols = Columns::build(&Shape::new(16, 8, 1, 1, vec![]), 100);
    let (ok, r) = passes(&ints(&[1, -4, 0, 0, 2, 1, 1, 0, 0]), &cols.hi, 100);
    assert!(ok, "residual 1e{r}");
}

#[test]
fn wrong_relation_fails() {
    let cols = Columns::build(&Shape::new(16, 8, 1, 1, vec![]), 100);
    let (ok, r) = passes(&ints(&[1, -4, 0, 0, 2, 1, 2, 0, 0]), &cols.hi, 100);
    assert!(!ok);
    assert!(r > -5.0);
}
