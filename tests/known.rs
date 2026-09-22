use pihunt::basis::Shape;
use pihunt::known::is_known;
use rug::Integer;

fn ints(v: &[i64]) -> Vec<Integer> {
    v.iter().map(|&c| Integer::from(c)).collect()
}

#[test]
fn recognises_bbp_in_any_scaling_and_sign() {
    let shape = Shape::new(16, 8, 1, 1, vec![]);
    let names = shape.column_names();
    assert!(is_known(
        &shape,
        &names,
        &ints(&[1, -4, 0, 0, 2, 1, 1, 0, 0])
    ));
    assert!(is_known(
        &shape,
        &names,
        &ints(&[-3, 12, 0, 0, -6, -3, -3, 0, 0])
    ));
}

#[test]
fn rejects_near_misses() {
    let shape = Shape::new(16, 8, 1, 1, vec![]);
    let names = shape.column_names();
    assert!(!is_known(
        &shape,
        &names,
        &ints(&[1, -4, 0, 0, 2, 1, 2, 0, 0])
    ));
    let other_base = Shape::new(10, 8, 1, 1, vec![]);
    assert!(!is_known(
        &other_base,
        &names,
        &ints(&[1, -4, 0, 0, 2, 1, 1, 0, 0])
    ));
}
