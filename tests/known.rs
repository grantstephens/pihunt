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

/// Bellard (1997), written out directly in the base-2^20, period-40 columns.
const BELLARD_DIRECT: &[(usize, i64)] = &[
    (2, -524288),
    (5, 163840),
    (6, 131072),
    (10, 8192),
    (14, 8192),
    (15, 5120),
    (18, -2048),
    (22, 512),
    (25, -160),
    (26, -128),
    (30, -8),
    (34, -8),
    (35, -5),
    (38, 2),
];

/// The same relation as PSLQ returns it, with dependent columns dropped first.
const BELLARD_REDUCED: &[(usize, i64)] = &[
    (2, -1048576),
    (4, 786432),
    (5, 819200),
    (8, -327680),
    (10, -24576),
    (12, 49152),
    (15, -15360),
    (16, -20480),
    (18, -4096),
    (20, -2048),
    (24, -1280),
    (25, -800),
    (26, -256),
    (28, 192),
    (30, -40),
    (32, -80),
    (34, -16),
    (35, 15),
    (36, 12),
];

fn bellard_relation(shape: &Shape, terms: &[(usize, i64)]) -> Vec<Integer> {
    let mut rel = vec![Integer::new(); shape.columns()];
    rel[0] = Integer::from(65536);
    for &(j, c) in terms {
        rel[j] = Integer::from(c); // column j is S(j, s=1)
    }
    rel
}

#[test]
fn bellard_forms_are_genuine_relations_and_known() {
    let shape = Shape::new(1 << 20, 40, 1, 1, vec![]);
    let names = shape.column_names();
    let cols = shape.build(1800);
    for terms in [BELLARD_DIRECT, BELLARD_REDUCED] {
        let rel = bellard_relation(&shape, terms);
        // Independent of the PSLQ code: the dot product with the columns vanishes to ~500 digits.
        let mut sum = rug::Float::with_val(1800, 0);
        for (c, x) in rel.iter().zip(&cols) {
            sum += rug::Float::with_val(1800, x * c);
        }
        assert!(sum.is_zero() || sum.get_exp().unwrap() < -1700, "residual {sum:e}");
        assert!(is_known(&shape, &names, &rel));
    }
}

#[test]
fn bellard_rejects_near_miss() {
    let shape = Shape::new(1 << 20, 40, 1, 1, vec![]);
    let names = shape.column_names();
    let mut rel = bellard_relation(&shape, BELLARD_DIRECT);
    rel[38] = Integer::from(3);
    assert!(!is_known(&shape, &names, &rel));
}
