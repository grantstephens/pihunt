use pihunt::basis::{Extra, Shape};
use pihunt::known::{is_derived, is_known};
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
        assert!(
            sum.is_zero() || sum.get_exp().unwrap() < -1700,
            "residual {sum:e}"
        );
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

/// `batches/base2-focused.toml`'s real hit at base=256, period=32: base 256 = 16^2 (t=2 radix
/// split of BBP's base 16, period 8), then period-multiplied by d=2 (16 -> 32). Not a new
/// formula -- see docs/findings-bbp-hunt.md.
#[test]
fn recognises_bbp_base_power_and_period_multiple_descendant() {
    let shape = Shape::new(256, 32, 1, 1, vec![]);
    let names = shape.column_names();
    let mut rel = vec![0i64; shape.columns()];
    rel[0] = 8;
    for &(j, c) in &[
        (2, -64),
        (8, 32),
        (10, 16),
        (12, 16),
        (18, -4),
        (24, 2),
        (26, 1),
        (28, 1),
    ] {
        rel[j] = c;
    }
    let rel = ints(&rel);
    assert!(
        !is_known(&shape, &names, &rel),
        "not a verbatim known entry"
    );
    assert!(is_derived(&shape, &names, &rel));
}

/// Same base/period as above, but an unrelated relation -- must not match just because the
/// shape happens to be a BBP descendant.
#[test]
fn rejects_unrelated_relation_at_a_descendant_shape() {
    let shape = Shape::new(256, 32, 1, 1, vec![]);
    let names = shape.column_names();
    let mut rel = vec![0i64; shape.columns()];
    rel[0] = 1;
    rel[2] = -1;
    let rel = ints(&rel);
    assert!(!is_derived(&shape, &names, &rel));
}

/// Bailey's base-64 pi^2 formula (period 6, degree 2, extra pi2) radix-split with t=2 to base
/// 4096 = 64^2, period 12, d=1 (no period multiply this time).
#[test]
fn recognises_bailey_pi2_base_power_descendant() {
    let shape = Shape::new(4096, 12, 2, 2, vec![Extra::Pi2]);
    let names = shape.column_names();
    let mut rel = vec![0i64; shape.columns()];
    for &(j, c) in &[
        (1, 9216),
        (2, -13824),
        (3, -4608),
        (4, -3456),
        (5, 576),
        (7, 144),
        (8, -216),
        (9, -72),
        (10, -54),
        (11, 9),
    ] {
        rel[j] = c;
    }
    *rel.last_mut().unwrap() = -512; // pi2, the last column
    let rel = ints(&rel);
    assert!(is_derived(&shape, &names, &rel));
}

/// 512 = 2^9 isn't a power of 16, 64, or 2^20 at a period any known formula's radix split +
/// period-multiply reaches -- a hit here can't be waved away by this mechanism.
#[test]
fn clean_base_is_not_derived_from_anything() {
    let shape = Shape::new(512, 8, 1, 1, vec![]);
    let names = shape.column_names();
    let mut rel = vec![0i64; shape.columns()];
    rel[0] = 1;
    rel[1] = -4;
    rel[4] = 2;
    rel[5] = 1;
    rel[6] = 1;
    let rel = ints(&rel);
    assert!(!is_derived(&shape, &names, &rel));
}
