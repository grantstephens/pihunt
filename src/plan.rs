//! Expands a batch into concrete jobs.

use crate::basis::{Shape, auto_digits};
use crate::config::{Batch, Method, Mode, Precision};
use crate::pslq::PslqParams;

/// Bump whenever a change could alter any job's outcome.
pub const ALGO_VERSION: u32 = 1;

#[derive(Debug, Clone, PartialEq)]
pub struct Job {
    pub shape: Shape,
    pub coeff_bound: u64,
    /// Resolved working precision, decimal digits.
    pub digits: u32,
    pub gamma: f64,
    pub max_iterations: u64,
}

impl Job {
    pub fn params(&self) -> PslqParams {
        PslqParams {
            gamma: self.gamma,
            coeff_bound: self.coeff_bound,
            max_iterations: self.max_iterations,
            digits: self.digits,
        }
    }

    /// Stable 32-hex-char ID over everything that determines the outcome.
    pub fn id(&self, finder: &str) -> String {
        let s = &self.shape;
        let extras: Vec<&str> = s.extras.iter().map(|e| e.name()).collect();
        let canon = format!(
            "v{ALGO_VERSION}|{finder}|b{}|m{}|s{}-{}|x{}|C{}|D{}|g{:?}|i{}",
            s.base,
            s.period,
            s.s_lo,
            s.s_hi,
            extras.join(","),
            self.coeff_bound,
            self.digits,
            self.gamma,
            self.max_iterations
        );
        blake3::hash(canon.as_bytes()).to_hex()[..32].to_string()
    }
}

/// All jobs for a batch, deduplicated, in a deterministic order.
pub fn plan(batch: &Batch) -> Vec<Job> {
    let s = &batch.search;
    let periods: Vec<u32> = (s.periods.from..=s.periods.to).collect();
    let axes = [
        s.bases.len(),
        periods.len(),
        s.degrees.len(),
        s.extras.len(),
    ];
    let picks: Vec<[usize; 4]> = match s.mode {
        Mode::Grid => grid(axes),
        Mode::Sample => {
            let smp = batch.sample.as_ref().expect("validated");
            match smp.method {
                Method::Sobol => sobol(axes, smp.count, smp.seed),
                Method::Lhs => lhs(axes, smp.count, smp.seed),
            }
        }
    };
    let mut jobs: Vec<Job> = Vec::new();
    for [bi, pi, di, ei] in picks {
        let [lo, hi] = s.degrees[di];
        let shape = Shape::new(s.bases[bi], periods[pi], lo, hi, s.extras[ei].clone());
        let job = make_job(batch, shape);
        if !jobs.contains(&job) {
            jobs.push(job);
        }
    }
    jobs
}

fn make_job(batch: &Batch, shape: Shape) -> Job {
    let d = &batch.defaults;
    let digits = match d.precision_digits {
        Precision::Auto => auto_digits(shape.columns(), d.coeff_bound),
        Precision::Fixed(p) => p,
    };
    Job {
        shape,
        coeff_bound: d.coeff_bound,
        digits,
        gamma: d.gamma,
        max_iterations: d.max_iterations,
    }
}

fn grid(axes: [usize; 4]) -> Vec<[usize; 4]> {
    let mut out = Vec::new();
    for a in 0..axes[0] {
        for b in 0..axes[1] {
            for c in 0..axes[2] {
                for d in 0..axes[3] {
                    out.push([a, b, c, d]);
                }
            }
        }
    }
    out
}

fn to_index(u: f64, len: usize) -> usize {
    ((u * len as f64) as usize).min(len - 1)
}

fn sobol(axes: [usize; 4], count: u32, seed: u32) -> Vec<[usize; 4]> {
    (0..count)
        .map(|i| {
            std::array::from_fn(|d| {
                to_index(sobol_burley::sample(i, d as u32, seed) as f64, axes[d])
            })
        })
        .collect()
}

fn lhs(axes: [usize; 4], count: u32, seed: u32) -> Vec<[usize; 4]> {
    let n = count as usize;
    let mut rng = SplitMix64(seed as u64);
    let strata: Vec<Vec<usize>> = (0..4)
        .map(|_| {
            let mut p: Vec<usize> = (0..n).collect();
            for i in (1..n).rev() {
                p.swap(i, (rng.next() % (i as u64 + 1)) as usize);
            }
            p
        })
        .collect();
    (0..n)
        .map(|i| {
            std::array::from_fn(|d| {
                let u = (strata[d][i] as f64 + rng.unit()) / n as f64;
                to_index(u, axes[d])
            })
        })
        .collect()
}

struct SplitMix64(u64);

impl SplitMix64 {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / (1u64 << 53) as f64
    }
}
