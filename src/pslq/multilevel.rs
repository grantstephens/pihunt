//! Two-level PSLQ (Bailey's pslqm2): cheap f64 iterations, periodic full-precision syncs.
//!
//! The f64 inner loop only steers. It accumulates an exact integer transform (entries kept
//! below 2^52, so f64 holds them exactly), which is then applied to the full-precision
//! `State`. All termination checks, and therefore all claims, come from that `State`.

use super::state::State;
use super::{Outcome, PslqParams, RelationFinder};
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

/// f64 copy of the state plus the integer transform accumulated since the last sync.
#[derive(Clone)]
struct Inner {
    n: usize,
    y: Vec<f64>,
    h: Vec<Vec<f64>>,
    /// Left transform applied to rows (A-side).
    a: Vec<Vec<f64>>,
    /// Right transform applied to columns (B-side), the inverse of `a`.
    b: Vec<Vec<f64>>,
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
            b: (0..n).map(eye).collect(),
        }
    }

    /// Run up to `budget` f64 PSLQ iterations. Returns how many were kept.
    /// Stops early once the f64 bound estimate reaches `exclude_at`, so the full-precision
    /// state gets to make the exclusion call instead of the search running past it.
    fn run(&mut self, gamma: f64, exclude_at: f64, budget: u64) -> u64 {
        let gpow: Vec<f64> = (0..self.n - 1).map(|r| gamma.powi(r as i32 + 1)).collect();
        let mut used = 0;
        let mut snapshot = self.clone();
        while used < budget && !self.y_collapsed() && self.bound_estimate() <= exclude_at {
            snapshot.copy_from(self);
            self.step(&gpow);
            if !self.exact() {
                self.copy_from(&snapshot);
                break;
            }
            used += 1;
        }
        used
    }

    /// `*self = other.clone()` without reallocating (this runs every inner iteration).
    fn copy_from(&mut self, other: &Inner) {
        self.y.copy_from_slice(&other.y);
        for (dst, src) in [
            (&mut self.h, &other.h),
            (&mut self.a, &other.a),
            (&mut self.b, &other.b),
        ] {
            for (d, s) in dst.iter_mut().zip(src) {
                d.copy_from_slice(s);
            }
        }
    }

    fn bound_estimate(&self) -> f64 {
        let m = (0..self.n - 1).fold(0.0f64, |m, j| m.max(self.h[j][j].abs()));
        1.0 / m
    }

    fn y_collapsed(&self) -> bool {
        let max = self.y.iter().fold(0.0f64, |m, v| m.max(v.abs()));
        let min = self.y.iter().fold(f64::INFINITY, |m, v| m.min(v.abs()));
        min <= Y_FLOOR * max // y is finite here: `exact` rejects anything else
    }

    /// Transform entries are exact integers and H/y are finite.
    fn exact(&self) -> bool {
        let ok = |v: &f64| v.is_finite() && v.abs() < EXACT;
        self.a.iter().flatten().all(ok)
            && self.b.iter().flatten().all(ok)
            && self.y.iter().all(|v| v.is_finite())
            && self.h.iter().flatten().all(|v| v.is_finite())
    }

    /// Same iteration as `State::step`, in f64.
    #[allow(clippy::needless_range_loop)] // rows i and j of H are read and written together
    fn step(&mut self, gpow: &[f64]) {
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

        self.y.swap(r, r + 1);
        self.a.swap(r, r + 1);
        h.swap(r, r + 1);
        for row in self.b.iter_mut() {
            row.swap(r, r + 1);
        }

        if r < n - 2 {
            let (t1, t2) = (h[r][r], h[r][r + 1]);
            let t3 = t1.hypot(t2);
            for row in h.iter_mut().skip(r) {
                let (u, v) = (row[r], row[r + 1]);
                row[r] = (t1 * u + t2 * v) / t3;
                row[r + 1] = (t1 * v - t2 * u) / t3;
            }
        }

        let (y, a, b) = (&mut self.y, &mut self.a, &mut self.b);
        for i in r + 1..n {
            for j in (0..=(i - 1).min(r + 1)).rev() {
                if h[j][j] == 0.0 {
                    continue;
                }
                let t = (h[i][j] / h[j][j]).round();
                if t == 0.0 {
                    continue;
                }
                y[j] += t * y[i];
                for k in 0..=j {
                    h[i][k] -= t * h[j][k];
                }
                for k in 0..n {
                    a[i][k] -= t * a[j][k];
                    b[k][j] += t * b[k][i];
                }
            }
        }
    }
}

/// Apply the inner loop's exact integer transform to the full-precision state,
/// restore H's lower-trapezoidal shape, then do a full Hermite reduction.
#[allow(clippy::needless_range_loop)] // matrix products index both operands by k
fn sync(st: &mut State, inner: &Inner) {
    let (n, prec) = (st.n, st.prec);
    let to_int = |m: &[Vec<f64>]| -> Vec<Vec<Integer>> {
        m.iter()
            .map(|row| {
                row.iter()
                    .map(|v| Integer::from_f64(*v).expect("exact integer"))
                    .collect()
            })
            .collect()
    };
    let ta = to_int(&inner.a);
    let tb = to_int(&inner.b);
    // The same entries as 64-bit Floats (exact: all below 2^52), for fused multiply-adds.
    let to_float = |m: &[Vec<f64>]| -> Vec<Vec<Float>> {
        m.iter()
            .map(|row| row.iter().map(|v| Float::with_val(64, *v)).collect())
            .collect()
    };
    let ta_f = to_float(&inner.a);
    let tb_f = to_float(&inner.b);

    // y <- y · TB
    st.y = (0..n)
        .map(|j| {
            let mut s = Float::with_val(prec, 0);
            for k in 0..n {
                if tb[k][j] != 0 {
                    s += &st.y[k] * &tb_f[k][j];
                }
            }
            s
        })
        .collect();
    // B <- B · TB
    st.b = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let mut s = Integer::new();
                    for k in 0..n {
                        if tb[k][j] != 0 {
                            s += &st.b[i][k] * &tb[k][j];
                        }
                    }
                    s
                })
                .collect()
        })
        .collect();
    // A <- TA · A
    st.a = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| {
                    let mut s = Integer::new();
                    for k in 0..n {
                        if ta[i][k] != 0 {
                            s += &ta[i][k] * &st.a[k][j];
                        }
                    }
                    s
                })
                .collect()
        })
        .collect();
    // H <- TA · H
    st.h = (0..n)
        .map(|i| {
            (0..n - 1)
                .map(|j| {
                    let mut s = Float::with_val(prec, 0);
                    for k in 0..n {
                        if ta[i][k] != 0 {
                            s += &st.h[k][j] * &ta_f[i][k];
                        }
                    }
                    s
                })
                .collect()
        })
        .collect();

    lq(&mut st.h, prec);
    st.hermite_reduce(1, n - 1);
}

/// Givens rotations from the right until H is lower trapezoidal again (H <- H·Q).
/// This is the hot spot of a sync, so the inner loop reuses two scratch Floats instead of
/// allocating four per row.
fn lq(h: &mut [Vec<Float>], prec: u32) {
    let m = h[0].len();
    let mut r = Float::new(prec);
    let mut c = Float::new(prec);
    let mut s = Float::new(prec);
    let mut t1 = Float::new(prec);
    let mut t2 = Float::new(prec);
    for i in 0..m {
        for j in i + 1..m {
            if h[i][j].is_zero() {
                continue;
            }
            r.assign(h[i][i].hypot_ref(&h[i][j]));
            c.assign(&h[i][i] / &r);
            s.assign(&h[i][j] / &r);
            for row in h.iter_mut().skip(i) {
                let (left, right) = row.split_at_mut(j);
                let (u, v) = (&mut left[i], &mut right[0]);
                t1.assign(&c * &*u);
                t1 += &s * &*v; // c·u + s·v
                t2.assign(&c * &*v);
                t2 -= &s * &*u; // c·v − s·u
                std::mem::swap(u, &mut t1);
                std::mem::swap(v, &mut t2);
            }
        }
    }
}
