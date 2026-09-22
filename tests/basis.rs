use pihunt::basis::{Columns, Extra, Shape, auto_digits, series};
use pihunt::pslq::digits_to_bits;
use rug::{Float, float::Constant};

fn close(a: &Float, b: &Float, digits: u32) -> bool {
    let d = Float::with_val(a.prec(), a - b).abs();
    d < Float::with_val(a.prec(), Float::u_pow_u(10, digits)).recip()
}

#[test]
fn bbp_identity_holds() {
    let bits = digits_to_bits(100);
    let s = |j| series(16, 8, j, 1, bits);
    let rhs = Float::with_val(bits, 4 * s(1) - 2 * s(4)) - s(5) - s(6);
    assert!(close(&rhs, &Float::with_val(bits, Constant::Pi), 98));
}

#[test]
fn base2_period1_is_two_log2() {
    let bits = digits_to_bits(100);
    let two_log2 = Float::with_val(bits, Constant::Log2) * 2;
    assert!(close(&series(2, 1, 1, 1, bits), &two_log2, 98));
}

#[test]
fn extras_match_mpfr() {
    let bits = digits_to_bits(80);
    assert!(close(
        &Extra::Log5.value(bits),
        &Float::with_val(bits, 5).ln(),
        78
    ));
    assert!(close(
        &Extra::Pi2.value(bits),
        &Float::with_val(bits, Constant::Pi).square(),
        78
    ));
    assert!(close(
        &Extra::Zeta3.value(bits),
        &Float::with_val(bits, Float::zeta_u(3)),
        78
    ));
}

#[test]
fn column_layout() {
    let shape = Shape::new(10, 2, 1, 2, vec![Extra::Log5, Extra::Log2, Extra::Log5]);
    assert_eq!(shape.extras, vec![Extra::Log2, Extra::Log5]);
    assert_eq!(shape.columns(), 7);
    assert_eq!(
        shape.column_names(),
        [
            "pi",
            "S(j=1,s=1)",
            "S(j=2,s=1)",
            "S(j=1,s=2)",
            "S(j=2,s=2)",
            "log2",
            "log5"
        ]
    );
    assert_eq!(shape.build(128).len(), 7);
}

#[test]
fn auto_digits_monotone() {
    assert!(auto_digits(10, 1000) < auto_digits(11, 1000));
    assert!(auto_digits(10, 1000) < auto_digits(10, 10_000));
    assert_eq!(auto_digits(9, 1000), 91); // ceil(9 * 3 * 1.5) + 50
}

#[test]
fn columns_have_two_precisions() {
    let cols = Columns::build(&Shape::new(16, 8, 1, 1, vec![]), 100);
    assert_eq!(cols.lo[0].prec(), digits_to_bits(100));
    assert_eq!(cols.hi[0].prec(), digits_to_bits(200));
    assert!(close(
        &Float::with_val(cols.hi[0].prec(), &cols.hi[0]),
        &Float::with_val(cols.hi[0].prec(), Constant::Pi),
        198
    ));
}

#[test]
fn extra_names_round_trip() {
    for e in [
        Extra::Pi2,
        Extra::Log2,
        Extra::Log3,
        Extra::Log5,
        Extra::Catalan,
        Extra::Zeta3,
    ] {
        assert_eq!(e.name().parse::<Extra>(), Ok(e));
    }
    assert!("log7".parse::<Extra>().is_err());
}
