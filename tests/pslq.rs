use pihunt::pslq::{Outcome, PslqParams, RelationFinder, classic::ClassicPslq, digits_to_bits};
use rug::{Float, Integer, rand::RandState};

fn params(digits: u32, coeff_bound: u64) -> PslqParams {
    PslqParams {
        gamma: 1.16,
        coeff_bound,
        max_iterations: 100_000,
        digits,
    }
}

fn relation(o: Outcome) -> Vec<i64> {
    match o {
        Outcome::Relation { coeffs, .. } => {
            let mut v: Vec<i64> = coeffs.iter().map(|c| c.to_i64().unwrap()).collect();
            if v.iter().find(|c| **c != 0).unwrap() < &0 {
                v.iter_mut().for_each(|c| *c = -*c);
            }
            v
        }
        other => panic!("expected relation, got {other:?}"),
    }
}

#[test]
fn sqrt2_sqrt8() {
    let p = params(60, 100);
    let b = digits_to_bits(60);
    let x = [Float::with_val(b, 2).sqrt(), Float::with_val(b, 8).sqrt()];
    assert_eq!(relation(ClassicPslq.find(&x, &p)), vec![2, -1]);
}

#[test]
fn logs() {
    let p = params(60, 100);
    let b = digits_to_bits(60);
    let x = [
        Float::with_val(b, 2).ln(),
        Float::with_val(b, 3).ln(),
        Float::with_val(b, 6).ln(),
    ];
    assert_eq!(relation(ClassicPslq.find(&x, &p)), vec![1, 1, -1]);
}

fn random_reals(n: usize, bits: u32, seed: u64) -> Vec<Float> {
    let mut rng = RandState::new();
    rng.seed(&Integer::from(seed));
    (0..n)
        .map(|_| Float::with_val(bits, Float::random_bits(&mut rng)))
        .collect()
}

#[test]
fn planted_relations_recovered() {
    let digits = 80;
    let bits = digits_to_bits(digits);
    for seed in 0..50u64 {
        let mut x = random_reals(5, bits, seed);
        let a: Vec<i64> = (0..5)
            .map(|i| ((seed * 7 + i * 13) % 19) as i64 - 9)
            .collect();
        let mut last = Float::with_val(bits, 0);
        for (ai, xi) in a.iter().zip(&x) {
            last -= Float::with_val(bits, xi * *ai);
        }
        x.push(last);
        let mut want: Vec<i64> = a.clone();
        want.push(1);
        if want.iter().find(|c| **c != 0).unwrap() < &0 {
            want.iter_mut().for_each(|c| *c = -*c);
        }
        let got = relation(ClassicPslq.find(&x, &params(digits, 1000)));
        assert_eq!(got, want, "seed {seed}");
    }
}

#[test]
fn no_false_relations_on_random_input() {
    let digits = 60;
    let bits = digits_to_bits(digits);
    for seed in 0..1000u64 {
        let x = random_reals(6, bits, seed);
        if let Outcome::Relation { coeffs, .. } = ClassicPslq.find(&x, &params(digits, 1000)) {
            panic!("seed {seed}: false relation {coeffs:?}");
        }
    }
}

#[test]
fn bound_is_monotone() {
    let digits = 60;
    let bits = digits_to_bits(digits);
    let x = random_reals(8, bits, 7);
    let mut prev = 0.0;
    for cap in 0..60 {
        let mut p = params(digits, 1_000_000);
        p.max_iterations = cap;
        let b = ClassicPslq.find(&x, &p).bound();
        assert!(
            b >= prev,
            "bound fell from {prev} to {b} at iteration {cap}"
        );
        prev = b;
    }
}
