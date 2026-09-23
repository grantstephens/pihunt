use pi_digits::bignum::Big;

fn lcg(seed: &mut u64) -> u64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    *seed
}

/// A pseudo-random Big with roughly `words` 64-bit words, built through the public API only.
fn rand_big(seed: &mut u64, words: usize) -> Big {
    let mut x = Big::zero();
    for _ in 0..words {
        x = x.shl(64).add(&Big::from_u64(lcg(seed)));
    }
    x
}

#[test]
fn ring_identities_hold_on_large_values() {
    let mut s = 1u64;
    for words in [1usize, 3, 17, 64, 200] {
        let (a, b, c) = (
            rand_big(&mut s, words),
            rand_big(&mut s, words),
            rand_big(&mut s, words / 2 + 1),
        );
        assert_eq!(
            a.mul(&b).add(&a.mul(&c)),
            a.mul(&b.add(&c)),
            "distributive, words={words}"
        );
        let r = a.rem(&c);
        let q_times_c = a.sub(&r); // a = q*c + r, so a - r is divisible by c
        assert_eq!(q_times_c.rem(&c), Big::zero());
        assert!(r.bits() <= c.bits());
        let m = lcg(&mut s) | 1;
        assert_eq!(a.rem_u64(m), a.rem(&Big::from_u64(m)).to_u64().unwrap());
    }
}

#[test]
fn exact_division_and_small_helpers() {
    let x = Big::from_u64(12345).mul(&Big::pow_u64(10, 30));
    assert_eq!(x.div_u64_exact(12345), Big::pow_u64(10, 30));
    assert_eq!(Big::from_u64(17).div_u64(5), Big::from_u64(3)); // truncates, doesn't round
    assert_eq!(Big::pow_u64(10, 3).to_u64(), Some(1000));
    assert_eq!(Big::binomial(10, 3).to_u64(), Some(120));
    assert_eq!(
        Big::binomial(100_000_000, 4).to_decimal_string(),
        "4166666416666671249999975000000"
    );
    assert_eq!(
        Big::pow_u64(2, 70).to_decimal_string(),
        "1180591620717411303424"
    );
    assert!(Big::zero().is_zero());
    assert_eq!(Big::one().mul_u64(7), Big::from_u64(7));
}

#[test]
fn matches_rug_on_random_inputs() {
    use rug::Integer;
    let to_rug = |b: &Big| b.to_decimal_string().parse::<Integer>().unwrap();
    let mut s = 7u64;
    for words in [1usize, 2, 9, 40, 150] {
        let (a, b) = (rand_big(&mut s, words), rand_big(&mut s, words / 3 + 1));
        let (ra, rb) = (to_rug(&a), to_rug(&b));
        assert_eq!(to_rug(&a.mul(&b)), Integer::from(&ra * &rb));
        assert_eq!(to_rug(&a.add(&b)), Integer::from(&ra + &rb));
        assert_eq!(to_rug(&a.rem(&b)), Integer::from(&ra % &rb));
        let m = lcg(&mut s) | 1;
        assert_eq!(a.rem_u64(m), Integer::from(&ra % m).to_u64().unwrap());
        assert_eq!(to_rug(&a.div_u64(m)), Integer::from(&ra / m));
    }
}
