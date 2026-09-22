use pihunt::pslq::{
    Outcome, PslqParams, RelationFinder, classic::ClassicPslq, digits_to_bits,
    multilevel::MultilevelPslq,
};
use rug::{Float, Integer, rand::RandState};

fn params(digits: u32, coeff_bound: u64) -> PslqParams {
    PslqParams {
        gamma: 1.16,
        coeff_bound,
        max_iterations: 100_000,
        digits,
    }
}

/// Every finder must pass every test in this file.
fn finders() -> [&'static dyn RelationFinder; 2] {
    [&ClassicPslq, &MultilevelPslq]
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
    for f in finders() {
        assert_eq!(relation(f.find(&x, &p)), vec![2, -1], "{}", f.name());
    }
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
    for f in finders() {
        assert_eq!(relation(f.find(&x, &p)), vec![1, 1, -1], "{}", f.name());
    }
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
        for f in finders() {
            let got = relation(f.find(&x, &params(digits, 1000)));
            assert_eq!(got, want, "{} seed {seed}", f.name());
        }
    }
}

#[test]
fn no_false_relations_on_random_input() {
    let digits = 60;
    let bits = digits_to_bits(digits);
    for seed in 0..1000u64 {
        let x = random_reals(6, bits, seed);
        for f in finders() {
            if let Outcome::Relation { coeffs, .. } = f.find(&x, &params(digits, 1000)) {
                panic!("{} seed {seed}: false relation {coeffs:?}", f.name());
            }
        }
    }
}

#[test]
fn bound_is_monotone() {
    let digits = 60;
    let bits = digits_to_bits(digits);
    let x = random_reals(8, bits, 7);
    for f in finders() {
        let mut prev = 0.0;
        for cap in 0..60 {
            let mut p = params(digits, 1_000_000);
            p.max_iterations = cap;
            let b = f.find(&x, &p).bound();
            assert!(
                b >= prev,
                "{}: bound fell from {prev} to {b} at iteration {cap}",
                f.name()
            );
            prev = b;
        }
    }
}

/// Big enough that multilevel needs many f64 inner runs and full-precision syncs.
#[test]
fn planted_relation_needing_many_syncs() {
    let digits = 250;
    let bits = digits_to_bits(digits);
    let n = 24;
    let mut x = random_reals(n - 1, bits, 99);
    let a: Vec<i64> = (0..n as i64 - 1)
        .map(|i| (i * 7919 % 2001) - 1000)
        .collect();
    let mut last = Float::with_val(bits, 0);
    for (ai, xi) in a.iter().zip(&x) {
        last -= Float::with_val(bits, xi * *ai);
    }
    x.push(last);
    let mut want = a.clone();
    want.push(1);
    if want.iter().find(|c| **c != 0).unwrap() < &0 {
        want.iter_mut().for_each(|c| *c = -*c);
    }
    for f in finders() {
        let got = relation(f.find(&x, &params(digits, 10_000)));
        assert_eq!(got, want, "{}", f.name());
    }
}

#[test]
fn multilevel_excludes_random_input_like_classic() {
    let digits = 200;
    let bits = digits_to_bits(digits);
    for seed in 0..5u64 {
        let x = random_reals(20, bits, seed);
        let c = ClassicPslq.find(&x, &params(digits, 1000));
        let m = MultilevelPslq.find(&x, &params(digits, 1000));
        assert!(
            matches!(c, Outcome::Excluded { .. }),
            "classic seed {seed}: {c:?}"
        );
        assert!(
            matches!(m, Outcome::Excluded { .. }),
            "multilevel seed {seed}: {m:?}"
        );
    }
}
