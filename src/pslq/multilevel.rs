//! Two-level PSLQ (Bailey's pslqm2): cheap f64 iterations, periodic full-precision syncs.
//!
//! The f64 inner loop only steers. It accumulates an exact integer transform (entries kept
//! below 2^52, so f64 holds them exactly), which is then applied to the full-precision
//! `State`. All termination checks, and therefore all claims, come from that `State`.

use super::state::State;
use super::{Outcome, PslqParams, RelationFinder};
use rayon::prelude::*;
use rug::ops::NegAssign;
use rug::{Assign, Float, Integer};

pub struct MultilevelPslq;

/// f64 represents every integer below 2^53 exactly; stay a factor of 2 clear.
const EXACT: f64 = 4_503_599_627_370_496.0; // 2^52
/// Inner loop stops once the scaled y spans this many orders of magnitude: f64 can't resolve it.
const Y_FLOOR: f64 = 1e-12;

impl RelationFinder for MultilevelPslq {
    fn name(&self) -> &'static str {
        "multilevel"
    }

    fn find(&self, x: &[Float], p: &PslqParams) -> Outcome {
        let mut st = State::new(x, p);
        let mut iterations = 0u64;
        loop {
            if let Some(outcome) = st.check(iterations) {
                return outcome;
            }
            let mut inner = Inner::from_state(&st);
            let used = inner.run(st.gamma(), st.exclude_at(), p.max_iterations - iterations);
            if used == 0 {
                // f64 can't make progress here (e.g. y spans too many magnitudes).
                iterations += 1;
                st.step();
                continue;
            }
            iterations += used;
            sync(&mut st, &inner);
        }
    }
}

/// How far one f64 iteration got.
enum Step {
    /// Swap, corner and every reduction applied.
    Full,
    /// Swap and corner applied, reduction stopped before leaving the exact range.
    Partial,
    /// Nothing applied.
    Refused,
}

/// f64 copy of the state plus the integer transform accumulated since the last sync.
struct Inner {
    n: usize,
    y: Vec<f64>,
    h: Vec<Vec<f64>>,
    /// Left transform applied to rows (A-side).
    a: Vec<Vec<f64>>,
    /// Transpose of the right transform (B-side, the inverse of `a`). Stored transposed so
    /// B's column operations are contiguous row operations here.
    bt: Vec<Vec<f64>>,
}

impl Inner {
    fn from_state(st: &State) -> Self {
        let n = st.n;
        let ymax =
            st.y.iter()
                .map(|v| Float::with_val(st.prec, v.abs_ref()))
                .max_by(|a, b| a.partial_cmp(b).unwrap())
                .expect("n >= 2");
        let y =
            st.y.iter()
                .map(|v| Float::with_val(st.prec, v / &ymax).to_f64())
                .collect();
        let h =
            st.h.iter()
                .map(|row| row.iter().map(Float::to_f64).collect())
                .collect();
        let eye = |i: usize| (0..n).map(|j| (i == j) as u8 as f64).collect();
        Inner {
            n,
            y,
            h,
            a: (0..n).map(eye).collect(),
            bt: (0..n).map(eye).collect(),
        }
    }

    /// Run up to `budget` f64 PSLQ iterations. Returns how many were applied.
    /// Stops early once the f64 bound estimate reaches `exclude_at`, so the full-precision
    /// state gets to make the exclusion call instead of the search running past it.
    fn run(&mut self, gamma: f64, exclude_at: f64, budget: u64) -> u64 {
        let gpow: Vec<f64> = (0..self.n - 1).map(|r| gamma.powi(r as i32 + 1)).collect();
        let mut used = 0;
        while used < budget && !self.y_collapsed() && self.bound_estimate() <= exclude_at {
            match self.step(&gpow) {
                Step::Full => used += 1,
                Step::Partial => return used + 1,
                Step::Refused => return used,
            }
        }
        used
    }

    fn bound_estimate(&self) -> f64 {
        let m = (0..self.n - 1).fold(0.0f64, |m, j| m.max(self.h[j][j].abs()));
        1.0 / m
    }

    fn y_collapsed(&self) -> bool {
        let max = self.y.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let min = self.y.iter().fold(f64::INFINITY, |m, v| m.min(v.abs()));
        min <= Y_FLOOR * max
    }

    /// Same iteration as `State::step`, in f64. Every reduction is checked *before* it is
    /// applied, so the transform stays exact without snapshots: if one would leave the exact
    /// integer range, the step stops there. What was applied is still a valid (partially
    /// reduced) PSLQ state, and the full reduction at the next sync finishes the job.
    #[allow(clippy::needless_range_loop)] // rows i and j of H are read and written together
    fn step(&mut self, gpow: &[f64]) -> Step {
        let n = self.n;
        let h = &mut self.h;
        let mut r = 0;
        let mut best = -1.0;
        for j in 0..n - 1 {
            let v = h[j][j].abs() * gpow[j];
            if v > best {
                best = v;
                r = j;
            }
        }

        // The corner rotation divides by this; refuse the step if f64 can't do it.
        let corner = r < n - 2;
        let (t1, t2) = if corner {
            (h[r + 1][r], h[r + 1][r + 1])
        } else {
            (0.0, 0.0)
        };
        let t3 = t1.hypot(t2);
        if corner && !(t3 > 0.0 && t3.is_finite()) {
            return Step::Refused;
        }

        self.y.swap(r, r + 1);
        self.a.swap(r, r + 1);
        h.swap(r, r + 1);
        self.bt.swap(r, r + 1);

        if corner {
            for row in h.iter_mut().skip(r) {
                let (u, v) = (row[r], row[r + 1]);
                row[r] = (t1 * u + t2 * v) / t3;
                row[r + 1] = (t1 * v - t2 * u) / t3;
            }
        }

        let (y, a, bt) = (&mut self.y, &mut self.a, &mut self.bt);
        for i in r + 1..n {
            for j in (0..=(i - 1).min(r + 1)).rev() {
                if h[j][j] == 0.0 {
                    continue;
                }
                let t = (h[i][j] / h[j][j]).round();
                if t == 0.0 {
                    continue;
                }
                // j < i, so split to hold row j (read) and row i (write) at once.
                let (a_lo, a_hi) = a.split_at_mut(i);
                let (a_j, a_i) = (&a_lo[j], &mut a_hi[0]);
                let (b_lo, b_hi) = bt.split_at_mut(i);
                let (b_j, b_i) = (&mut b_lo[j], &b_hi[0]);
                let fits = t.abs() < EXACT
                    && a_i.iter().zip(a_j).all(|(x, z)| (x - t * z).abs() < EXACT)
                    && b_j.iter().zip(b_i).all(|(x, z)| (x + t * z).abs() < EXACT);
                if !fits {
                    return Step::Partial;
                }
                y[j] += t * y[i];
                let (h_lo, h_hi) = h.split_at_mut(i);
                for (x, z) in h_hi[0][..=j].iter_mut().zip(&h_lo[j][..=j]) {
                    *x -= t * z;
                }
                for (x, z) in a_i.iter_mut().zip(a_j) {
                    *x -= t * z;
                }
                for (x, z) in b_j.iter_mut().zip(b_i) {
                    *x += t * z;
                }
            }
        }
        Step::Full
    }
}

/// Apply the inner loop's exact integer transform to the full-precision state,
/// restore H's lower-trapezoidal shape, then do a full Hermite reduction.
#[allow(clippy::needless_range_loop)] // matrix products index both operands by k
fn sync(st: &mut State, inner: &Inner) {
    let (n, prec) = (st.n, st.prec);
    let int = |m: &[Vec<f64>]| -> Vec<Vec<Integer>> {
        m.iter()
            .map(|row| {
                row.iter()
                    .map(|v| Integer::from_f64(*v).expect("exact integer"))
                    .collect()
            })
            .collect()
    };
    // The same entries as 64-bit Floats (exact: all below 2^52), for fused multiply-adds.
    let float = |m: &[Vec<f64>]| -> Vec<Vec<Float>> {
        m.iter()
            .map(|row| row.iter().map(|v| Float::with_val(64, *v)).collect())
            .collect()
    };
    let (ta, ta_f) = (int(&inner.a), float(&inner.a));
    let (tbt, tbt_f) = (int(&inner.bt), float(&inner.bt));

    // Rows are independent in every product below, so they run in parallel. Inside a busy
    // batch the extra tasks mostly stay on the current thread; a lone big job gets every core.
    // y <- y · TB, i.e. y_j = y · (TB^T)_j
    let y = &st.y;
    st.y = tbt_f
        .par_iter()
        .map(|col| {
            let mut s = Float::with_val(prec, 0);
            for (yk, c) in y.iter().zip(col) {
                if !c.is_zero() {
                    s += yk * c;
                }
            }
            s
        })
        .collect();
    // B <- B · TB: B_ij = B_i · (TB^T)_j
    let b = &st.b;
    st.b = b
        .par_iter()
        .map(|row| {
            tbt.iter()
                .map(|col| {
                    let mut s = Integer::new();
                    for (x, c) in row.iter().zip(col) {
                        if *c != 0 {
                            s += x * c;
                        }
                    }
                    s
                })
                .collect()
        })
        .collect();
    // A <- TA · A and H <- TA · H: row i is a combination of rows k weighted by TA_ik.
    let a = &st.a;
    st.a = ta
        .par_iter()
        .map(|t_row| {
            let mut out = vec![Integer::new(); n];
            for (c, a_row) in t_row.iter().zip(a) {
                if *c != 0 {
                    for (o, x) in out.iter_mut().zip(a_row) {
                        *o += c * x;
                    }
                }
            }
            out
        })
        .collect();
    let h = &st.h;
    st.h = ta_f
        .par_iter()
        .map(|t_row| {
            let mut out = vec![Float::with_val(prec, 0); n - 1];
            for (c, h_row) in t_row.iter().zip(h) {
                if !c.is_zero() {
                    for (o, x) in out.iter_mut().zip(h_row) {
                        *o += c * x;
                    }
                }
            }
            out
        })
        .collect();

    lq(&mut st.h, prec);
    st.hermite_reduce(1, n - 1);
}

/// Householder reflections from the right until H is lower trapezoidal again (H <- H·Q).
/// This is the hot spot of a sync: roughly half the multiplies of Givens rotations, and
/// every Float is updated in place (no allocation inside the loops).
fn lq(h: &mut [Vec<Float>], prec: u32) {
    let m = h[0].len();
    let mut u = vec![Float::new(prec); m];
    let mut norm2 = Float::new(prec);
    let mut alpha = Float::new(prec);
    let mut beta = Float::new(prec);
    for i in 0..m {
        if h[i][i + 1..].iter().all(|v| v.is_zero()) {
            continue;
        }
        // Reflect row i's tail v = h[i][i..] onto (alpha, 0, ..., 0) with |alpha| = |v|.
        norm2.assign(0);
        for v in &h[i][i..] {
            norm2 += v * v;
        }
        alpha.assign(norm2.sqrt_ref());
        if h[i][i].is_sign_positive() {
            alpha.neg_assign();
        }
        // u = v - alpha e1, and beta = 2 / (u · u) = 1 / (norm2 - alpha v0).
        u[i].assign(&h[i][i] - &alpha);
        for j in i + 1..m {
            u[j].assign(&h[i][j]);
        }
        beta.assign(&alpha * &h[i][i]);
        beta = Float::with_val(prec, &norm2 - &beta).recip();
        let (u, beta) = (&u, &beta);
        h[i..].par_iter_mut().for_each(|row| {
            let mut dot = Float::with_val(prec, 0);
            for (x, uj) in row[i..].iter().zip(&u[i..]) {
                dot += x * uj;
            }
            dot *= beta;
            for (x, uj) in row[i..].iter_mut().zip(&u[i..]) {
                *x -= &dot * uj;
            }
        });
        // Exact zeros where the reflection annihilated row i.
        h[i][i].assign(&alpha);
        for v in &mut h[i][i + 1..] {
            v.assign(0);
        }
    }
}
