//! Textbook Ferguson–Bailey PSLQ, everything at full MPFR precision.

use super::{Outcome, PslqParams, RelationFinder, digits_to_bits};
use rug::{Float, Integer};

pub struct ClassicPslq;

impl RelationFinder for ClassicPslq {
    fn name(&self) -> &'static str {
        "classic"
    }

    fn find(&self, x: &[Float], p: &PslqParams) -> Outcome {
        pslq(x, p)
    }
}

fn pslq(x: &[Float], p: &PslqParams) -> Outcome {
    let n = x.len();
    assert!(n >= 2, "PSLQ needs at least 2 inputs");
    assert!(
        p.gamma > (4.0f64 / 3.0).sqrt(),
        "gamma must exceed sqrt(4/3)"
    );
    let prec = digits_to_bits(p.digits);
    let f = |v: f64| Float::with_val(prec, v);

    // Detection threshold and precision-exhaustion limit.
    let margin = p.digits.saturating_sub(30);
    let eps = Float::with_val(prec, Float::u_pow_u(10, margin)).recip();
    let a_limit = Integer::from(Integer::u_pow_u(10, margin));
    let exclude_at = p.coeff_bound as f64 * (n as f64).sqrt();

    // Partial norms s_k = sqrt(sum_{j>=k} x_j^2), then normalise.
    let mut s: Vec<Float> = vec![f(0.0); n];
    let mut acc = f(0.0);
    for k in (0..n).rev() {
        acc += Float::with_val(prec, x[k].square_ref());
        s[k] = Float::with_val(prec, acc.sqrt_ref());
    }
    let t0 = s[0].clone();
    let mut y: Vec<Float> = x.iter().map(|v| Float::with_val(prec, v) / &t0).collect();
    for sk in s.iter_mut() {
        *sk /= &t0;
    }

    // H is n x (n-1), lower trapezoidal.
    let mut h: Vec<Vec<Float>> = vec![vec![f(0.0); n - 1]; n];
    for i in 0..n {
        for j in 0..(n - 1).min(i + 1) {
            h[i][j] = if i == j {
                Float::with_val(prec, &s[j + 1] / &s[j])
            } else {
                let num = Float::with_val(prec, &y[i] * &y[j]);
                let den = Float::with_val(prec, &s[j] * &s[j + 1]);
                -(num / den)
            };
        }
    }

    let mut a: Vec<Vec<Integer>> = identity(n);
    let mut b: Vec<Vec<Integer>> = identity(n);

    hermite_reduce(&mut h, &mut y, &mut a, &mut b, 1, n - 1, prec);

    let gpow: Vec<Float> = (0..n - 1).map(|r| f(p.gamma.powi(r as i32 + 1))).collect();

    let mut iterations = 0u64;
    loop {
        let bound = norm_bound(&h);
        if let Some(i) = tiny_y(&y, &eps) {
            let coeffs = (0..n).map(|k| b[k][i].clone()).collect();
            return Outcome::Relation {
                coeffs,
                iterations,
                bound,
            };
        }
        if bound > exclude_at {
            return Outcome::Excluded { bound, iterations };
        }
        if max_abs(&a) > a_limit {
            return Outcome::PrecisionExhausted { bound, iterations };
        }
        if iterations >= p.max_iterations {
            return Outcome::IterationCap { bound, iterations };
        }
        iterations += 1;

        // Pick r maximising gamma^(r+1) |H_rr|.
        let mut r = 0;
        let mut best = f(-1.0);
        for j in 0..n - 1 {
            let v = Float::with_val(prec, h[j][j].abs_ref()) * &gpow[j];
            if v > best {
                best = v;
                r = j;
            }
        }

        // Swap r and r+1 in y, A, H (rows) and B (columns).
        y.swap(r, r + 1);
        a.swap(r, r + 1);
        h.swap(r, r + 1);
        for row in b.iter_mut() {
            row.swap(r, r + 1);
        }

        // Corner fix-up restores lower-trapezoidal form.
        if r < n - 2 {
            let t1 = h[r][r].clone();
            let t2 = h[r][r + 1].clone();
            let t3 = Float::with_val(prec, t1.hypot_ref(&t2));
            for row in h.iter_mut().skip(r) {
                let u = row[r].clone();
                let v = row[r + 1].clone();
                row[r] = Float::with_val(prec, &t1 * &u + &t2 * &v) / &t3;
                row[r + 1] = Float::with_val(prec, &t1 * &v - &t2 * &u) / &t3;
            }
        }

        hermite_reduce(&mut h, &mut y, &mut a, &mut b, r + 1, r + 1, prec);
    }
}

/// Reduce rows `from..n` against columns `min(i-1, jmax) ..= 0`.
#[allow(clippy::needless_range_loop)] // rows i and j of H are read and written together
fn hermite_reduce(
    h: &mut [Vec<Float>],
    y: &mut [Float],
    a: &mut [Vec<Integer>],
    b: &mut [Vec<Integer>],
    from: usize,
    jmax: usize,
    prec: u32,
) {
    let n = y.len();
    for i in from..n {
        for j in (0..=(i - 1).min(jmax)).rev() {
            if h[j][j].is_zero() {
                continue;
            }
            let q = Float::with_val(prec, &h[i][j] / &h[j][j]).round();
            let t = q.to_integer().expect("finite H entries");
            if t == 0 {
                continue;
            }
            let yi = Float::with_val(prec, &y[i] * &t);
            y[j] += yi;
            for k in 0..=j {
                let d = Float::with_val(prec, &h[j][k] * &t);
                h[i][k] -= d;
            }
            for k in 0..n {
                let d = Integer::from(&t * &a[j][k]);
                a[i][k] -= d;
                let e = Integer::from(&t * &b[k][i]);
                b[k][j] += e;
            }
        }
    }
}

fn identity(n: usize) -> Vec<Vec<Integer>> {
    (0..n)
        .map(|i| (0..n).map(|j| Integer::from((i == j) as u8)).collect())
        .collect()
}

/// Any integer relation has Euclidean norm >= 1 / max_j |H_jj|.
fn norm_bound(h: &[Vec<Float>]) -> f64 {
    let m = (0..h[0].len())
        .map(|j| h[j][j].to_f64().abs())
        .fold(0.0, f64::max);
    if m == 0.0 { f64::INFINITY } else { 1.0 / m }
}

fn tiny_y(y: &[Float], eps: &Float) -> Option<usize> {
    let (i, v) = y
        .iter()
        .enumerate()
        .map(|(i, v)| (i, Float::with_val(v.prec(), v.abs_ref())))
        .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap())?;
    (v < *eps).then_some(i)
}

fn max_abs(a: &[Vec<Integer>]) -> Integer {
    a.iter()
        .flatten()
        .map(|v| Integer::from(v.abs_ref()))
        .max()
        .unwrap_or_default()
}
