use rug::{Float, Integer, float::Constant, ops::Pow};

fn mpfr_digits(n: u64, count: usize) -> String {
    let total = n + count as u64 + 25;
    let bits = pi_digits::digits_to_bits(total as u32) + 16;
    let pi = Float::with_val(bits, Constant::Pi);
    let scaled = Float::with_val(bits, &pi * &Float::with_val(bits, Integer::from(10).pow(total as u32)));
    let s = scaled.to_integer().unwrap().to_string(); // "3" + digits
    s[(n as usize + 1)..(n as usize + 1 + count)].to_string()
}

#[test]
fn pi_ref_matches_mpfr_everywhere_it_is_used() {
    for n in (0..2100).step_by(7) {
        assert_eq!(pi_digits::pi_ref::digits(n, 10), mpfr_digits(n, 10), "n={n}");
    }
    // Feynman point: positions 762..=767 are 999999.
    assert_eq!(&pi_digits::pi_ref::digits(761, 6), "999999");
    assert_eq!(pi_digits::pi_ref::digits(0, 5), "14159");
}
