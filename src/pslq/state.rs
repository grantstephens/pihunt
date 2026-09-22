//! Full-precision PSLQ state shared by every finder.
//!
//! Invariants: `a · b = I`, `y = x_normalised · b`, and `h` is (up to a right orthogonal
//! factor) `a` applied to the initial H. That last one is what makes `norm_bound` a valid
//! lower bound on any relation's norm, so only this state is allowed to make claims.

use super::{Outcome, PslqParams, digits_to_bits};
use rug::{Float, Integer};

pub(super) struct State {
    pub n: usize,
    pub prec: u32,
    pub y: Vec<Float>,
    /// n x (n-1), lower trapezoidal between steps.
    pub h: Vec<Vec<Float>>,
    pub a: Vec<Vec<Integer>>,
    pub b: Vec<Vec<Integer>>,
    gamma: f64,
    gpow: Vec<Float>,
    eps: Float,
    a_limit: Integer,
    exclude_at: f64,
    max_iterations: u64,
}

impl State {
    pub fn new(x: &[Float], p: &PslqParams) -> Self {
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
        let y: Vec<Float> = x.iter().map(|v| Float::with_val(prec, v) / &t0).collect();
        for sk in s.iter_mut() {
            *sk /= &t0;
        }

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

        let gpow = (0..n - 1).map(|r| f(p.gamma.powi(r as i32 + 1))).collect();
        let mut st = State {
            n,
            prec,
            y,
            h,
            a: identity(n),
            b: identity(n),
            gamma: p.gamma,
            gpow,
            eps,
            a_limit,
            exclude_at,
            max_iterations: p.max_iterations,
        };
        st.hermite_reduce(1, n - 1);
        st
    }

    pub fn gamma(&self) -> f64 {
        self.gamma
    }

    /// Bound past which `check` reports `Excluded`.
    pub fn exclude_at(&self) -> f64 {
        self.exclude_at
    }

    /// The termination checks, in order. `None` means keep going.
    pub fn check(&self, iterations: u64) -> Option<Outcome> {
        let bound = norm_bound(&self.h);
        if let Some(i) = tiny_y(&self.y, &self.eps) {
            let coeffs = (0..self.n).map(|k| self.b[k][i].clone()).collect();
            return Some(Outcome::Relation {
                coeffs,
                iterations,
                bound,
            });
        }
        if bound > self.exclude_at {
            return Some(Outcome::Excluded { bound, iterations });
        }
        if max_abs(&self.a) > self.a_limit {
            return Some(Outcome::PrecisionExhausted { bound, iterations });
        }
        if iterations >= self.max_iterations {
            return Some(Outcome::IterationCap { bound, iterations });
        }
        None
    }

    /// One textbook PSLQ iteration at full precision.
    #[allow(clippy::needless_range_loop)] // H's diagonal and gpow are indexed together
    pub fn step(&mut self) {
        let (n, prec) = (self.n, self.prec);
        let h = &mut self.h;

        // Pick r maximising gamma^(r+1) |H_rr|.
        let mut r = 0;
        let mut best = Float::with_val(prec, -1);
        for j in 0..n - 1 {
            let v = Float::with_val(prec, h[j][j].abs_ref()) * &self.gpow[j];
            if v > best {
                best = v;
                r = j;
            }
        }

        // Swap r and r+1 in y, A, H (rows) and B (columns).
        self.y.swap(r, r + 1);
        self.a.swap(r, r + 1);
        h.swap(r, r + 1);
        for row in self.b.iter_mut() {
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

        self.hermite_reduce(r + 1, r + 1);
    }

    /// Reduce rows `from..n` against columns `min(i-1, jmax) ..= 0`.
    #[allow(clippy::needless_range_loop)] // rows i and j of H are read and written together
    pub fn hermite_reduce(&mut self, from: usize, jmax: usize) {
        let (n, prec) = (self.n, self.prec);
        let (h, y, a, b) = (&mut self.h, &mut self.y, &mut self.a, &mut self.b);
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
