# pihunt Stage 1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `pihunt` stage 1: a Rust CLI that runs batches of PSLQ searches for BBP-type π formulas, rediscovers the known ones, and logs every result (including exclusion bounds) to JSONL.

**Architecture:** One binary crate with a library core. `basis` computes constants at P and 2P digits. `reduce` strips dependent columns using PSLQ plus verification. The main PSLQ search runs over π plus the surviving columns, and `classify` turns the outcome into a verdict. `job` wires one job end to end, `plan` expands TOML batches into jobs, `log` appends JSONL, and `main` drives it all with rayon.

**Tech Stack:** Rust 2024 (stable 1.98 via mise), `rug` 1.30 on system GMP/MPFR, `rayon`, `serde`/`serde_json`/`toml`, `clap`, `blake3`, `sobol_burley`, `jiff`, and `tempfile` (tests).

**Spec:** `docs/superpowers/specs/2026-09-22-pihunt-design.md`

**Provenance:** every code block in this plan was compiled and tested in a scratch prototype before the plan was written: 37 tests pass, and clippy is clean. Type it in exactly as given. If something doesn't compile, the transcription is wrong, not the code.

## Global Constraints

- Toolchain via mise: `mise.toml` pins `rust = "stable"` and sets `CARGO_TARGET_DIR = "{{env.HOME}}/.cache/pihunt/target"`. **Never** build into the repo, because `~/sync` is synced.
- Run cargo through mise so the env applies: `mise exec -- cargo …` (or have mise shell activation enabled).
- `gmp-mpfr-sys` must use the `use-system-libs` feature (system GMP/MPFR are installed).
- Auto precision: `ceil(n · log10(C) · 1.5) + 50` decimal digits. Fixed precision must be ≥ 60.
- PSLQ `gamma` must be > sqrt(4/3). Default 1.16.
- Every relation, whether a basis relation or a main-search hit, is trusted only if max |coeff| ≤ `coeff_bound` **and** it verifies against columns built at 2× precision (`residual < 10^-(2·digits − 60)`).
- Integers in the JSONL log are serialised as strings.
- Bump `ALGO_VERSION` in `src/plan.rs` whenever outcome-affecting math changes.
- Commit messages end with `Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>`.

## File Structure

```
Cargo.toml, Cargo.lock
src/lib.rs            module list
src/main.rs           CLI: run / plan / verify, progress + summary
src/pslq/mod.rs       PslqParams, Outcome, RelationFinder, digits_to_bits, primitive, exceeds
src/pslq/classic.rs   ClassicPslq (Ferguson–Bailey PSLQ at full MPFR precision)
src/basis.rs          Extra, Shape, series, auto_digits, Columns (P and 2P values)
src/verify.rs         residual_log10, passes
src/known.rs          table of literature formulas, is_known
src/reduce.rs         reduce → Reduced::{Done, Inconclusive}, Dropped
src/config.rs         Batch TOML parsing + validation
src/plan.rs           Job, ALGO_VERSION, plan (grid / sobol / lhs)
src/log.rs            Record, Kind, spawn_writer, read
src/classify.rs       classify → Verdict
src/job.rs            run_job → Record
tests/*.rs            one integration test file per module + cli.rs + timing.rs
batches/*.toml        known.toml, base10-scout.toml
results/              committed JSONL logs
docs/timing-baseline.md
```

---

### Task 1: Crate scaffold and classic PSLQ

**Files:**
- Create: `Cargo.toml`, `src/lib.rs`, `src/main.rs` (stub), `src/pslq/mod.rs`, `src/pslq/classic.rs`
- Test: `tests/pslq.rs`

**Interfaces:**
- Produces: `pihunt::pslq::{PslqParams { gamma: f64, coeff_bound: u64, max_iterations: u64, digits: u32 }, Outcome::{Relation { coeffs: Vec<Integer>, iterations: u64, bound: f64 }, Excluded { bound, iterations }, PrecisionExhausted { bound, iterations }, IterationCap { bound, iterations }}, Outcome::{iterations(), bound()}, trait RelationFinder { fn name(&self) -> &'static str; fn find(&self, x: &[Float], params: &PslqParams) -> Outcome }, digits_to_bits(u32) -> u32, primitive(&[Integer]) -> Vec<Integer>, exceeds(&[Integer], u64) -> bool}` and `pihunt::pslq::classic::ClassicPslq`.

- [ ] **Step 1: Create `Cargo.toml`**

`Cargo.toml`:

```toml
[package]
name = "pihunt"
version = "0.1.0"
edition = "2024"

[dependencies]
blake3 = "1.8.7"
clap = { version = "4.6.7", features = ["derive"] }
# Link the system GMP/MPFR instead of building them from source.
gmp-mpfr-sys = { version = "1.7.1", features = ["use-system-libs"] }
jiff = "0.2.37"
rayon = "1.12.0"
rug = "1.30.0"
serde = { version = "1.0.229", features = ["derive"] }
serde_json = "1.0.151"
sobol_burley = "0.5.0"
toml = "1.1.6"

[dev-dependencies]
tempfile = "3.27.0"
```

- [ ] **Step 2: Create `src/lib.rs` with only the PSLQ module, and a stub `src/main.rs`**

```rust
pub mod pslq;
```

`src/main.rs` (replaced in Task 10):

```rust
fn main() {}
```

- [ ] **Step 3: Write the failing tests `tests/pslq.rs`**

These cover the toy relations, planted relations across 50 seeds, **zero false relations over 1000 random seeds**, and a monotone exclusion bound.

`tests/pslq.rs`:

```rust
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
```

- [ ] **Step 4: Run and confirm it fails**

Run: `mise exec -- cargo test --test pslq`
Expected: compile error, `file not found for module pslq` (the module files don't exist yet).

- [ ] **Step 5: Write `src/pslq/mod.rs`**

`src/pslq/mod.rs`:

```rust
//! Integer relation finding.

pub mod classic;

use rug::{Float, Integer};

/// Tuning and stopping parameters for one relation search.
#[derive(Debug, Clone)]
pub struct PslqParams {
    /// Must be > sqrt(4/3).
    pub gamma: f64,
    /// Largest |coefficient| we care about; drives the exclusion bound.
    pub coeff_bound: u64,
    pub max_iterations: u64,
    /// Working precision in decimal digits. Inputs must be accurate to this.
    pub digits: u32,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    /// `coeffs · x ≈ 0`, one coefficient per input.
    Relation {
        coeffs: Vec<Integer>,
        iterations: u64,
        bound: f64,
    },
    /// No relation with max |coefficient| ≤ coeff_bound exists.
    Excluded {
        bound: f64,
        iterations: u64,
    },
    PrecisionExhausted {
        bound: f64,
        iterations: u64,
    },
    IterationCap {
        bound: f64,
        iterations: u64,
    },
}

impl Outcome {
    pub fn iterations(&self) -> u64 {
        match self {
            Outcome::Relation { iterations, .. }
            | Outcome::Excluded { iterations, .. }
            | Outcome::PrecisionExhausted { iterations, .. }
            | Outcome::IterationCap { iterations, .. } => *iterations,
        }
    }

    pub fn bound(&self) -> f64 {
        match self {
            Outcome::Relation { bound, .. }
            | Outcome::Excluded { bound, .. }
            | Outcome::PrecisionExhausted { bound, .. }
            | Outcome::IterationCap { bound, .. } => *bound,
        }
    }
}

pub trait RelationFinder: Sync {
    fn name(&self) -> &'static str;
    fn find(&self, x: &[Float], params: &PslqParams) -> Outcome;
}

/// Decimal digits → MPFR bits.
pub fn digits_to_bits(digits: u32) -> u32 {
    (digits as f64 * std::f64::consts::LOG2_10).ceil() as u32
}

/// Divide by the gcd and flip sign so the first nonzero coefficient is positive.
pub fn primitive(coeffs: &[Integer]) -> Vec<Integer> {
    let g = coeffs.iter().fold(Integer::new(), |g, c| g.gcd(c));
    let sign = match coeffs.iter().find(|c| **c != 0) {
        Some(c) if *c < 0 => -1,
        _ => 1,
    };
    if g == 0 {
        return coeffs.to_vec();
    }
    coeffs
        .iter()
        .map(|c| Integer::from(c / &g) * sign)
        .collect()
}

/// True if any |coefficient| is larger than `bound`.
pub fn exceeds(coeffs: &[Integer], bound: u64) -> bool {
    coeffs.iter().any(|c| *c.as_abs() > bound)
}
```

- [ ] **Step 6: Write `src/pslq/classic.rs`**

Notes for the implementer: H is n×(n−1). The initial Hermite reduction is full. Per-iteration reduction only touches rows `r+1..n` against columns `min(i−1, r+1)..=0` (Bailey's reduced variant). `1/max|H_jj|` is a lower bound on any relation's norm, and it never decreases.

`src/pslq/classic.rs`:

```rust
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
```

- [ ] **Step 7: Run tests and confirm they pass**

Run: `mise exec -- cargo test --release --test pslq`
Expected: `5 passed`. (A debug build also passes, in about 2 s.)

- [ ] **Step 8: Commit**

```bash
git add Cargo.toml Cargo.lock src tests/pslq.rs
git commit -m "feat: crate scaffold and classic PSLQ

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 2: Basis columns

**Files:**
- Create: `src/basis.rs`
- Modify: `src/lib.rs` (add `pub mod basis;`)
- Test: `tests/basis.rs`

**Interfaces:**
- Consumes: `pslq::digits_to_bits`.
- Produces: `basis::{Extra::{Pi2, Log2, Log3, Log5, Catalan, Zeta3}` (serde lowercase, `name()`, `value(bits)`, `FromStr`), `Shape { base: u32, period: u32, s_lo: u32, s_hi: u32, extras: Vec<Extra> }` with `Shape::new(..)` (sorts and dedups extras), `columns() -> usize`, `column_names() -> Vec<String>` (`"pi"`, `"S(j=1,s=1)"`…, extras names), `build(bits) -> Vec<Float>`, `series(base, period, j, s, bits) -> Float`, `auto_digits(n, coeff_bound) -> u32`, `Columns { names, digits, lo, hi }` with `Columns::build(&Shape, digits)`}.

- [ ] **Step 1: Write the failing tests `tests/basis.rs`**

`tests/basis.rs`:

```rust
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
```

- [ ] **Step 2: Run and confirm it fails**

Run: `mise exec -- cargo test --test basis`
Expected: compile error, `unresolved import pihunt::basis`.

- [ ] **Step 3: Write `src/basis.rs` and add `pub mod basis;` to `src/lib.rs`**

`src/basis.rs`:

```rust
//! Builds the vector of constants PSLQ searches over.

use rug::{Float, Integer, float::Constant, ops::Pow};
use serde::{Deserialize, Serialize};

/// Extra constants that may appear alongside π and the series columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Extra {
    Pi2,
    Log2,
    Log3,
    Log5,
    Catalan,
    Zeta3,
}

impl Extra {
    pub fn name(self) -> &'static str {
        match self {
            Extra::Pi2 => "pi2",
            Extra::Log2 => "log2",
            Extra::Log3 => "log3",
            Extra::Log5 => "log5",
            Extra::Catalan => "catalan",
            Extra::Zeta3 => "zeta3",
        }
    }

    pub fn value(self, bits: u32) -> Float {
        match self {
            Extra::Pi2 => Float::with_val(bits, Constant::Pi).square(),
            Extra::Log2 => Float::with_val(bits, Constant::Log2),
            Extra::Log3 => Float::with_val(bits, 3).ln(),
            Extra::Log5 => Float::with_val(bits, 5).ln(),
            Extra::Catalan => Float::with_val(bits, Constant::Catalan),
            Extra::Zeta3 => Float::with_val(bits, Float::zeta_u(3)),
        }
    }
}

/// The formula shape for one job: base b, period m, degrees s_lo..=s_hi, extras.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Shape {
    pub base: u32,
    pub period: u32,
    pub s_lo: u32,
    pub s_hi: u32,
    /// Kept sorted and deduplicated.
    pub extras: Vec<Extra>,
}

impl Shape {
    pub fn new(base: u32, period: u32, s_lo: u32, s_hi: u32, mut extras: Vec<Extra>) -> Self {
        extras.sort();
        extras.dedup();
        Shape {
            base,
            period,
            s_lo,
            s_hi,
            extras,
        }
    }

    /// Total column count including π.
    pub fn columns(&self) -> usize {
        1 + (self.period * (self.s_hi - self.s_lo + 1)) as usize + self.extras.len()
    }

    /// Column labels, in the same order as `build`.
    pub fn column_names(&self) -> Vec<String> {
        let mut names = vec!["pi".to_string()];
        for s in self.s_lo..=self.s_hi {
            for j in 1..=self.period {
                names.push(format!("S(j={j},s={s})"));
            }
        }
        names.extend(self.extras.iter().map(|e| e.name().to_string()));
        names
    }

    /// All column values at `bits` of precision: π, series (degree-major), extras.
    pub fn build(&self, bits: u32) -> Vec<Float> {
        let mut cols = vec![Float::with_val(bits, Constant::Pi)];
        for s in self.s_lo..=self.s_hi {
            for j in 1..=self.period {
                cols.push(series(self.base, self.period, j, s, bits));
            }
        }
        cols.extend(self.extras.iter().map(|e| e.value(bits)));
        cols
    }
}

/// S(j,s) = sum_{k>=0} 1 / (b^k (mk+j)^s), accurate to `bits`.
pub fn series(base: u32, period: u32, j: u32, s: u32, bits: u32) -> Float {
    let log2b = (base as f64).log2();
    let rough_terms = (bits as f64 / log2b).ceil() + 2.0;
    let guard = rough_terms.log2().ceil() as u32 + 32;
    let wp = bits + guard;
    let terms = ((wp as f64) / log2b).ceil() as u64 + 2;

    let mut sum = Float::with_val(wp, 0);
    let mut pk = Float::with_val(wp, 1); // b^-k
    for k in 0..terms {
        let d = Integer::from(period as u64 * k + j as u64).pow(s);
        sum += Float::with_val(wp, &pk / &d);
        pk /= base;
    }
    Float::with_val(bits, &sum)
}

/// Decimal digits needed so PSLQ can find relations with max |coeff| <= c among n columns.
/// Factor 1.5 measured in prototyping: 1.25 let spurious 10^7-size relations through at n ~ 45.
pub fn auto_digits(n: usize, coeff_bound: u64) -> u32 {
    (n as f64 * (coeff_bound as f64).log10() * 1.5).ceil() as u32 + 50
}

/// Column values at working precision (`lo`) and at twice it (`hi`, for verification).
pub struct Columns {
    pub names: Vec<String>,
    pub digits: u32,
    pub lo: Vec<Float>,
    pub hi: Vec<Float>,
}

impl Columns {
    pub fn build(shape: &Shape, digits: u32) -> Self {
        let hi = shape.build(crate::pslq::digits_to_bits(2 * digits));
        let bits = crate::pslq::digits_to_bits(digits);
        let lo = hi.iter().map(|v| Float::with_val(bits, v)).collect();
        Columns {
            names: shape.column_names(),
            digits,
            lo,
            hi,
        }
    }
}

impl std::str::FromStr for Extra {
    type Err = String;
    fn from_str(s: &str) -> Result<Self, String> {
        [
            Extra::Pi2,
            Extra::Log2,
            Extra::Log3,
            Extra::Log5,
            Extra::Catalan,
            Extra::Zeta3,
        ]
        .into_iter()
        .find(|e| e.name() == s)
        .ok_or_else(|| format!("unknown extra {s:?}"))
    }
}
```

- [ ] **Step 4: Run tests and confirm they pass**

Run: `mise exec -- cargo test --release --test basis`
Expected: `7 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/basis.rs src/lib.rs tests/basis.rs
git commit -m "feat: basis columns with P and 2P precision

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 3: Verification and known-formula table

**Files:**
- Create: `src/verify.rs`, `src/known.rs`
- Modify: `src/lib.rs` (add `pub mod known;` and `pub mod verify;`)
- Test: `tests/verify.rs`, `tests/known.rs`

**Interfaces:**
- Consumes: `basis::{Shape, Columns}`, `pslq::primitive`.
- Produces: `verify::{residual_log10(&[Integer], &[Float]) -> f64, passes(relation: &[Integer], hi: &[Float], digits: u32) -> (bool, f64)}` and `known::is_known(&Shape, names: &[String], relation: &[Integer]) -> bool` (scale- and sign-insensitive; matches base + period + exact nonzero terms by column name).

- [ ] **Step 1: Write the failing tests `tests/verify.rs` and `tests/known.rs`**

`tests/verify.rs`:

```rust
use pihunt::basis::{Columns, Shape};
use pihunt::verify::passes;
use rug::Integer;

fn ints(v: &[i64]) -> Vec<Integer> {
    v.iter().map(|&c| Integer::from(c)).collect()
}

#[test]
fn true_relation_passes() {
    let cols = Columns::build(&Shape::new(16, 8, 1, 1, vec![]), 100);
    let (ok, r) = passes(&ints(&[1, -4, 0, 0, 2, 1, 1, 0, 0]), &cols.hi, 100);
    assert!(ok, "residual 1e{r}");
}

#[test]
fn wrong_relation_fails() {
    let cols = Columns::build(&Shape::new(16, 8, 1, 1, vec![]), 100);
    let (ok, r) = passes(&ints(&[1, -4, 0, 0, 2, 1, 2, 0, 0]), &cols.hi, 100);
    assert!(!ok);
    assert!(r > -5.0);
}
```

`tests/known.rs`:

```rust
use pihunt::basis::Shape;
use pihunt::known::is_known;
use rug::Integer;

fn ints(v: &[i64]) -> Vec<Integer> {
    v.iter().map(|&c| Integer::from(c)).collect()
}

#[test]
fn recognises_bbp_in_any_scaling_and_sign() {
    let shape = Shape::new(16, 8, 1, 1, vec![]);
    let names = shape.column_names();
    assert!(is_known(
        &shape,
        &names,
        &ints(&[1, -4, 0, 0, 2, 1, 1, 0, 0])
    ));
    assert!(is_known(
        &shape,
        &names,
        &ints(&[-3, 12, 0, 0, -6, -3, -3, 0, 0])
    ));
}

#[test]
fn rejects_near_misses() {
    let shape = Shape::new(16, 8, 1, 1, vec![]);
    let names = shape.column_names();
    assert!(!is_known(
        &shape,
        &names,
        &ints(&[1, -4, 0, 0, 2, 1, 2, 0, 0])
    ));
    let other_base = Shape::new(10, 8, 1, 1, vec![]);
    assert!(!is_known(
        &other_base,
        &names,
        &ints(&[1, -4, 0, 0, 2, 1, 1, 0, 0])
    ));
}
```

- [ ] **Step 2: Run and confirm they fail**

Run: `mise exec -- cargo test --test verify --test known`
Expected: compile errors, `unresolved import pihunt::verify` / `pihunt::known`.

- [ ] **Step 3: Write `src/verify.rs`**

`src/verify.rs`:

```rust
//! Checks a candidate relation against columns computed at twice the search precision.

use rug::{Float, Integer};

/// log10 |a · x|, or -inf if the residual is exactly zero.
pub fn residual_log10(relation: &[Integer], cols: &[Float]) -> f64 {
    let bits = cols[0].prec();
    let mut sum = Float::with_val(bits, 0);
    for (a, x) in relation.iter().zip(cols) {
        sum += Float::with_val(bits, x * a);
    }
    if sum.is_zero() {
        return f64::NEG_INFINITY;
    }
    // log10 via MPFR so tiny residuals don't underflow f64.
    Float::with_val(64, sum.abs_ref()).log10().to_f64()
}

/// A true relation's residual shrinks with the doubled precision; a spurious one stays put.
/// `hi` must be accurate to 2 * `digits` decimal digits.
pub fn passes(relation: &[Integer], hi: &[Float], digits: u32) -> (bool, f64) {
    let r = residual_log10(relation, hi);
    (r < -(2.0 * digits as f64 - 60.0), r)
}
```

- [ ] **Step 4: Write `src/known.rs` and add both modules to `src/lib.rs`**

`src/known.rs`:

```rust
//! Relations already in the literature, so rediscoveries don't get flagged as NEW.

use crate::basis::Shape;
use crate::pslq::primitive;
use rug::Integer;

struct Known {
    base: u32,
    period: u32,
    /// (column name, coefficient) for every nonzero coefficient, primitive form.
    terms: &'static [(&'static str, i64)],
}

const KNOWN: &[Known] = &[
    // Bailey–Borwein–Plouffe (1995).
    Known {
        base: 16,
        period: 8,
        terms: &[
            ("pi", 1),
            ("S(j=1,s=1)", -4),
            ("S(j=4,s=1)", 2),
            ("S(j=5,s=1)", 1),
            ("S(j=6,s=1)", 1),
        ],
    },
    // The base-16 "zero relation" that makes BBP non-unique.
    Known {
        base: 16,
        period: 8,
        terms: &[
            ("S(j=1,s=1)", 8),
            ("S(j=2,s=1)", -8),
            ("S(j=3,s=1)", -4),
            ("S(j=4,s=1)", -8),
            ("S(j=5,s=1)", -2),
            ("S(j=6,s=1)", -2),
            ("S(j=7,s=1)", 1),
        ],
    },
    // Bailey's base-64 formula for pi^2.
    Known {
        base: 64,
        period: 6,
        terms: &[
            ("S(j=1,s=2)", 144),
            ("S(j=2,s=2)", -216),
            ("S(j=3,s=2)", -72),
            ("S(j=4,s=2)", -54),
            ("S(j=5,s=2)", 9),
            ("pi2", -8),
        ],
    },
    // log 2 = sum 1 / (2^(k+1) (k+1)).
    Known {
        base: 2,
        period: 1,
        terms: &[("S(j=1,s=1)", 1), ("log2", -2)],
    },
];

/// True if `relation` (aligned with `names`) is a known formula for this shape.
pub fn is_known(shape: &Shape, names: &[String], relation: &[Integer]) -> bool {
    let rel = primitive(relation);
    let mut found: Vec<(&str, Integer)> = names
        .iter()
        .zip(rel)
        .filter(|(_, c)| *c != 0)
        .map(|(n, c)| (n.as_str(), c))
        .collect();
    found.sort_by(|a, b| a.0.cmp(b.0));
    KNOWN
        .iter()
        .filter(|k| k.base == shape.base && k.period == shape.period)
        .any(|k| {
            let mut want: Vec<(&str, i64)> = k.terms.to_vec();
            want.sort_by(|a, b| a.0.cmp(b.0));
            want.len() == found.len()
                && want
                    .iter()
                    .zip(&found)
                    .all(|(w, f)| w.0 == f.0 && f.1 == w.1)
        })
}
```

- [ ] **Step 5: Run tests and confirm they pass**

Run: `mise exec -- cargo test --release --test verify --test known`
Expected: `2 passed` for each file.

- [ ] **Step 6: Commit**

```bash
git add src/verify.rs src/known.rs src/lib.rs tests/verify.rs tests/known.rs
git commit -m "feat: 2x-precision verification and known-formula table

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 4: Basis reduction

**Files:**
- Create: `src/reduce.rs`
- Modify: `src/lib.rs` (add `pub mod reduce;`)
- Test: `tests/reduce.rs`

**Interfaces:**
- Consumes: `basis::Columns`, `pslq::{exceeds, primitive, Outcome, PslqParams, RelationFinder}`, `verify::passes`.
- Produces: `reduce::{Dropped { column: usize, relation: Vec<Integer> }, Reduced::{Done { keep: Vec<usize>, dropped: Vec<Dropped> }, Inconclusive { outcome: Outcome, dropped: Vec<Dropped> }}, reduce(&Columns, &dyn RelationFinder, &PslqParams) -> Reduced}`. `keep` holds non-π column indices, ascending, and never includes 0.

Why the checks matter: in prototyping, unchecked reduction at n ≈ 60 accepted spurious relations with 10⁶–10⁸ coefficients and dropped genuine columns (zeta3, catalan, log3). A relation that fails the coefficient bound or 2P verification makes the job **inconclusive**. It never drops a column.

- [ ] **Step 1: Write the failing tests `tests/reduce.rs`**

`tests/reduce.rs`:

```rust
use pihunt::basis::{Columns, Extra, Shape};
use pihunt::pslq::{PslqParams, classic::ClassicPslq};
use pihunt::reduce::{Reduced, reduce};
use rug::Integer;

fn params(digits: u32) -> PslqParams {
    PslqParams {
        gamma: 1.16,
        coeff_bound: 1000,
        max_iterations: 100_000,
        digits,
    }
}

fn ints(v: &[i64]) -> Vec<Integer> {
    v.iter().map(|&c| Integer::from(c)).collect()
}

#[test]
fn drops_planted_duplicate() {
    // Base 2, period 1: S(j=1,s=1) = 2 log 2, so log2 is redundant.
    let cols = Columns::build(&Shape::new(2, 1, 1, 1, vec![Extra::Log2]), 80);
    match reduce(&cols, &ClassicPslq, &params(80)) {
        Reduced::Done { keep, dropped } => {
            assert_eq!(keep, vec![1]);
            assert_eq!(dropped.len(), 1);
            assert_eq!(dropped[0].column, 2);
            assert_eq!(dropped[0].relation, ints(&[0, 1, -2]));
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn drops_bbp_zero_relation() {
    let cols = Columns::build(&Shape::new(16, 8, 1, 1, vec![]), 100);
    match reduce(&cols, &ClassicPslq, &params(100)) {
        Reduced::Done { keep, dropped } => {
            assert_eq!(keep, vec![1, 2, 3, 4, 5, 6, 8]);
            assert_eq!(dropped[0].relation, ints(&[0, 8, -8, -4, -8, -2, -2, 1, 0]));
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn independent_columns_untouched() {
    let cols = Columns::build(
        &Shape::new(10, 1, 1, 1, vec![Extra::Catalan, Extra::Zeta3]),
        80,
    );
    match reduce(&cols, &ClassicPslq, &params(80)) {
        Reduced::Done { keep, dropped } => {
            assert_eq!(keep, vec![1, 2, 3]);
            assert!(dropped.is_empty());
        }
        other => panic!("{other:?}"),
    }
}

#[test]
fn oversized_relation_is_inconclusive_not_dropped() {
    // Same duplicate as above, but a bound of 1 makes the true relation (1, -2) too big to trust.
    let cols = Columns::build(&Shape::new(2, 1, 1, 1, vec![Extra::Log2]), 80);
    let p = PslqParams {
        coeff_bound: 1,
        ..params(80)
    };
    assert!(matches!(
        reduce(&cols, &ClassicPslq, &p),
        Reduced::Inconclusive { .. }
    ));
}
```

- [ ] **Step 2: Run and confirm it fails**

Run: `mise exec -- cargo test --test reduce`
Expected: compile error, `unresolved import pihunt::reduce`.

- [ ] **Step 3: Write `src/reduce.rs` and add `pub mod reduce;` to `src/lib.rs`**

`src/reduce.rs`:

```rust
//! Strips rational linear dependencies out of the non-π columns.

use crate::basis::Columns;
use crate::pslq::{Outcome, PslqParams, RelationFinder, exceeds, primitive};
use crate::verify;
use rug::{Float, Integer};

/// A column removed because it is a rational combination of the others.
#[derive(Debug, Clone, PartialEq)]
pub struct Dropped {
    pub column: usize,
    /// Primitive relation, one coefficient per column of the full basis.
    pub relation: Vec<Integer>,
}

#[derive(Debug, Clone, PartialEq)]
pub enum Reduced {
    /// `keep` holds the surviving non-π column indices, ascending.
    Done {
        keep: Vec<usize>,
        dropped: Vec<Dropped>,
    },
    /// A step could neither find nor exclude a relation, or found one that failed checks.
    Inconclusive {
        outcome: Outcome,
        dropped: Vec<Dropped>,
    },
}

/// Repeatedly search columns 1.. (everything except π) and drop dependent ones.
/// A relation is only trusted if max |coeff| <= coeff_bound and it verifies at 2x precision.
pub fn reduce(cols: &Columns, finder: &dyn RelationFinder, params: &PslqParams) -> Reduced {
    let n = cols.lo.len();
    let mut keep: Vec<usize> = (1..n).collect();
    let mut dropped = Vec::new();
    while keep.len() >= 2 {
        let x: Vec<Float> = keep.iter().map(|&i| cols.lo[i].clone()).collect();
        let outcome = finder.find(&x, params);
        let Outcome::Relation { coeffs, .. } = &outcome else {
            if matches!(outcome, Outcome::Excluded { .. }) {
                break;
            }
            return Reduced::Inconclusive { outcome, dropped };
        };
        let mut full = vec![Integer::new(); n];
        for (&i, c) in keep.iter().zip(coeffs) {
            full[i] = c.clone();
        }
        let relation = primitive(&full);
        let too_big = exceeds(&relation, params.coeff_bound);
        if too_big || !verify::passes(&relation, &cols.hi, cols.digits).0 {
            return Reduced::Inconclusive { outcome, dropped };
        }
        let column = *keep
            .iter()
            .rev()
            .find(|&&i| relation[i] != 0)
            .expect("nonzero relation");
        keep.retain(|&i| i != column);
        dropped.push(Dropped { column, relation });
    }
    Reduced::Done { keep, dropped }
}
```

- [ ] **Step 4: Run tests and confirm they pass**

Run: `mise exec -- cargo test --release --test reduce`
Expected: `4 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/reduce.rs src/lib.rs tests/reduce.rs
git commit -m "feat: verified basis reduction

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 5: Batch config

**Files:**
- Create: `src/config.rs`
- Modify: `src/lib.rs` (add `pub mod config;`)
- Test: `tests/config.rs`

**Interfaces:**
- Consumes: `basis::Extra` (deserialised from lowercase names).
- Produces: `config::{Batch { name: String, output: PathBuf, threads: Option<usize>, defaults: Defaults, search: Search, sample: Option<Sample> }, Defaults { coeff_bound: u64, precision_digits: Precision, gamma: f64, max_iterations: u64, max_columns: usize }, Precision::{Auto, Fixed(u32)}, Mode::{Grid, Sample}, Range { from: u32, to: u32 }, Search { mode, bases: Vec<u32>, periods: Range, degrees: Vec<[u32; 2]>, extras: Vec<Vec<Extra>> }, Method::{Sobol, Lhs}, Sample { method, count: u32, seed: u32 }, load(&Path) -> Result<Batch, String>, parse(&str) -> Result<Batch, String>}`.

- [ ] **Step 1: Write the failing tests `tests/config.rs`**

`tests/config.rs`:

```rust
use pihunt::basis::Extra;
use pihunt::config::{Method, Mode, Precision, parse};

const GOOD: &str = r#"
name   = "t"
output = "results/t.jsonl"

[defaults]
coeff_bound = 1000

[search]
mode    = "grid"
bases   = [10, 16]
periods = { from = 2, to = 4 }
degrees = [[1, 1], [1, 2]]
extras  = [[], ["log2", "log5"]]
"#;

#[test]
fn parses_with_defaults() {
    let b = parse(GOOD).unwrap();
    assert_eq!(b.defaults.gamma, 1.16);
    assert_eq!(b.defaults.max_iterations, 1_000_000);
    assert_eq!(b.defaults.max_columns, 80);
    assert_eq!(b.defaults.precision_digits, Precision::Auto);
    assert_eq!(b.search.mode, Mode::Grid);
    assert_eq!(b.search.extras[1], vec![Extra::Log2, Extra::Log5]);
    assert_eq!(b.threads, None);
}

#[test]
fn parses_fixed_precision_and_sample() {
    let text = GOOD
        .replace(
            "coeff_bound = 1000",
            "coeff_bound = 1000\nprecision_digits = 300",
        )
        .replace("mode    = \"grid\"", "mode    = \"sample\"")
        + "\n[sample]\nmethod = \"lhs\"\ncount = 10\nseed = 7\n";
    let b = parse(&text).unwrap();
    assert_eq!(b.defaults.precision_digits, Precision::Fixed(300));
    assert_eq!(b.sample.unwrap().method, Method::Lhs);
}

#[test]
fn rejects_bad_configs() {
    let cases = [
        GOOD.replace("coeff_bound = 1000", "coeff_bound = 1000\ngamma = 1.1"),
        GOOD.replace("\"log5\"", "\"log7\""),
        GOOD.replace("from = 2, to = 4", "from = 5, to = 4"),
        GOOD.replace("mode    = \"grid\"", "mode    = \"sample\""),
        GOOD.replace(
            "coeff_bound = 1000",
            "coeff_bound = 1000\nprecision_digits = 30",
        ),
        GOOD.replace(
            "coeff_bound = 1000",
            "coeff_bound = 1000\nprecision_digits = \"lots\"",
        ),
        GOOD.replace("[[1, 1], [1, 2]]", "[[2, 1]]"),
        GOOD.replace("bases   = [10, 16]", "bases   = []"),
        GOOD.replace("name   = \"t\"", "name   = \"t\"\ntypo = 1"),
    ];
    for (i, c) in cases.iter().enumerate() {
        assert!(parse(c).is_err(), "case {i} should fail:\n{c}");
    }
}
```

- [ ] **Step 2: Run and confirm it fails**

Run: `mise exec -- cargo test --test config`
Expected: compile error, `unresolved import pihunt::config`.

- [ ] **Step 3: Write `src/config.rs` and add `pub mod config;` to `src/lib.rs`**

`src/config.rs`:

```rust
//! Batch file parsing and validation.

use crate::basis::Extra;
use serde::Deserialize;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Batch {
    pub name: String,
    pub output: PathBuf,
    /// Defaults to all cores.
    pub threads: Option<usize>,
    pub defaults: Defaults,
    pub search: Search,
    pub sample: Option<Sample>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Defaults {
    pub coeff_bound: u64,
    #[serde(default)]
    pub precision_digits: Precision,
    #[serde(default = "default_gamma")]
    pub gamma: f64,
    #[serde(default = "default_max_iterations")]
    pub max_iterations: u64,
    #[serde(default = "default_max_columns")]
    pub max_columns: usize,
}

fn default_gamma() -> f64 {
    1.16
}
fn default_max_iterations() -> u64 {
    1_000_000
}
fn default_max_columns() -> usize {
    80
}

/// `"auto"` or a fixed number of decimal digits.
#[derive(Debug, Clone, Copy, PartialEq, Default, Deserialize)]
#[serde(try_from = "toml::Value")]
pub enum Precision {
    #[default]
    Auto,
    Fixed(u32),
}

impl TryFrom<toml::Value> for Precision {
    type Error = String;
    fn try_from(v: toml::Value) -> Result<Self, String> {
        match v {
            toml::Value::String(s) if s == "auto" => Ok(Precision::Auto),
            toml::Value::Integer(i) if i > 0 && i <= u32::MAX as i64 => {
                Ok(Precision::Fixed(i as u32))
            }
            other => Err(format!(
                "precision_digits must be \"auto\" or a positive integer, got {other}"
            )),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Mode {
    Grid,
    Sample,
}

#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Range {
    pub from: u32,
    pub to: u32,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Search {
    pub mode: Mode,
    pub bases: Vec<u32>,
    pub periods: Range,
    /// Each entry is an inclusive degree range [s_lo, s_hi].
    pub degrees: Vec<[u32; 2]>,
    /// Each entry is one set of extras; use `[[]]` for none.
    pub extras: Vec<Vec<Extra>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Method {
    Sobol,
    Lhs,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sample {
    pub method: Method,
    pub count: u32,
    pub seed: u32,
}

pub fn load(path: &Path) -> Result<Batch, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    parse(&text)
}

pub fn parse(text: &str) -> Result<Batch, String> {
    let batch: Batch = toml::from_str(text).map_err(|e| e.to_string())?;
    validate(&batch)?;
    Ok(batch)
}

fn validate(b: &Batch) -> Result<(), String> {
    let d = &b.defaults;
    let s = &b.search;
    if d.gamma <= (4.0f64 / 3.0).sqrt() {
        return Err(format!(
            "gamma must be > sqrt(4/3) ≈ 1.1547, got {}",
            d.gamma
        ));
    }
    if d.coeff_bound < 2 {
        return Err("coeff_bound must be >= 2".into());
    }
    if let Precision::Fixed(p) = d.precision_digits
        && p < 60
    {
        return Err(format!("precision_digits must be >= 60, got {p}"));
    }
    if s.bases.is_empty() || s.bases.iter().any(|&b| b < 2) {
        return Err("bases must be non-empty and every base >= 2".into());
    }
    if s.periods.from < 1 || s.periods.from > s.periods.to {
        return Err(format!(
            "periods must satisfy 1 <= from <= to, got {:?}",
            s.periods
        ));
    }
    if s.degrees.is_empty() || s.degrees.iter().any(|[lo, hi]| *lo < 1 || lo > hi) {
        return Err(
            "degrees must be non-empty and every [lo, hi] must satisfy 1 <= lo <= hi".into(),
        );
    }
    if s.extras.is_empty() {
        return Err("extras must be non-empty; use [[]] for no extras".into());
    }
    if b.threads == Some(0) {
        return Err("threads must be >= 1".into());
    }
    match (s.mode, &b.sample) {
        (Mode::Sample, None) => return Err("mode = \"sample\" needs a [sample] section".into()),
        (Mode::Sample, Some(smp)) if smp.count == 0 || smp.count > 65_536 => {
            return Err("sample.count must be in 1..=65536".into());
        }
        _ => {}
    }
    Ok(())
}
```

- [ ] **Step 4: Run tests and confirm they pass**

Run: `mise exec -- cargo test --release --test config`
Expected: `3 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/config.rs src/lib.rs tests/config.rs
git commit -m "feat: batch TOML parsing and validation

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 6: Job planning and IDs

**Files:**
- Create: `src/plan.rs`
- Modify: `src/lib.rs` (add `pub mod plan;`)
- Test: `tests/plan.rs`

**Interfaces:**
- Consumes: `basis::{auto_digits, Shape}`, `config::{Batch, Method, Mode, Precision}`, `pslq::PslqParams`.
- Produces: `plan::{ALGO_VERSION: u32, Job { shape: Shape, coeff_bound: u64, digits: u32, gamma: f64, max_iterations: u64 }, Job::params() -> PslqParams, Job::id(finder: &str) -> String (32 hex chars), plan(&Batch) -> Vec<Job>}`.

- [ ] **Step 1: Write the failing tests `tests/plan.rs`**

`tests/plan.rs`:

```rust
use pihunt::basis::{Extra, Shape};
use pihunt::config::parse;
use pihunt::plan::{Job, plan};

const GRID: &str = r#"
name   = "t"
output = "results/t.jsonl"

[defaults]
coeff_bound = 1000

[search]
mode    = "grid"
bases   = [10, 16]
periods = { from = 2, to = 4 }
degrees = [[1, 1], [1, 2]]
extras  = [[], ["log2", "log5"], ["log5", "log2"]]
"#;

fn sample(method: &str, seed: u32) -> String {
    GRID.replace("mode    = \"grid\"", "mode    = \"sample\"")
        + &format!("\n[sample]\nmethod = \"{method}\"\ncount = 20\nseed = {seed}\n")
}

#[test]
fn grid_is_cartesian_and_deduplicated() {
    // Two extras sets are the same after sorting, so 2 * 3 * 2 * 2 unique jobs.
    let jobs = plan(&parse(GRID).unwrap());
    assert_eq!(jobs.len(), 24);
    assert_eq!(jobs[0].shape, Shape::new(10, 2, 1, 1, vec![]));
    assert_eq!(jobs[0].digits, 50 + (3.0f64 * 3.0 * 1.5).ceil() as u32);
}

#[test]
fn sample_modes_are_deterministic_and_in_range() {
    for method in ["sobol", "lhs"] {
        let a = plan(&parse(&sample(method, 1)).unwrap());
        let b = plan(&parse(&sample(method, 1)).unwrap());
        assert_eq!(a, b, "{method}");
        assert!(!a.is_empty() && a.len() <= 20);
        for j in &a {
            assert!([10, 16].contains(&j.shape.base));
            assert!((2..=4).contains(&j.shape.period));
        }
    }
}

fn job(extras: Vec<Extra>) -> Job {
    Job {
        shape: Shape::new(10, 4, 1, 2, extras),
        coeff_bound: 1000,
        digits: 120,
        gamma: 1.16,
        max_iterations: 10,
    }
}

#[test]
fn job_ids_are_stable_and_precise() {
    let a = job(vec![Extra::Log2, Extra::Log5]);
    let b = job(vec![Extra::Log5, Extra::Log2]);
    assert_eq!(a.id("classic"), b.id("classic"));
    assert_eq!(a.id("classic").len(), 32);
    assert_ne!(a.id("classic"), a.id("multilevel"));
    let variants = [
        Job {
            digits: 121,
            ..a.clone()
        },
        Job {
            coeff_bound: 999,
            ..a.clone()
        },
        Job {
            gamma: 1.2,
            ..a.clone()
        },
        Job {
            max_iterations: 11,
            ..a.clone()
        },
        Job {
            shape: Shape::new(100, 4, 1, 2, vec![Extra::Log2, Extra::Log5]),
            ..a.clone()
        },
        Job {
            shape: Shape::new(10, 5, 1, 2, vec![Extra::Log2, Extra::Log5]),
            ..a.clone()
        },
        Job {
            shape: Shape::new(10, 4, 1, 1, vec![Extra::Log2, Extra::Log5]),
            ..a.clone()
        },
        Job {
            shape: Shape::new(10, 4, 1, 2, vec![Extra::Log2]),
            ..a.clone()
        },
    ];
    for v in variants {
        assert_ne!(a.id("classic"), v.id("classic"), "{v:?}");
    }
}
```

- [ ] **Step 2: Run and confirm it fails**

Run: `mise exec -- cargo test --test plan`
Expected: compile error, `unresolved import pihunt::plan`.

- [ ] **Step 3: Write `src/plan.rs` and add `pub mod plan;` to `src/lib.rs`**

`src/plan.rs`:

```rust
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
```

- [ ] **Step 4: Run tests and confirm they pass**

Run: `mise exec -- cargo test --release --test plan`
Expected: `3 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/plan.rs src/lib.rs tests/plan.rs
git commit -m "feat: grid/sobol/lhs job planning with stable job IDs

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 7: Results log

**Files:**
- Create: `src/log.rs`
- Modify: `src/lib.rs` (add `pub mod log;`)
- Test: `tests/log.rs`

**Interfaces:**
- Produces: `log::{Kind::{Hit, Junk, Suspicious, Spurious, Excluded, Inconclusive, Skipped}` (serde lowercase), `Params { base, period, degrees: [u32; 2], extras: Vec<String>, coeff_bound, precision_digits, gamma, max_iterations, finder: String, algo_version }`, `DroppedRecord { column: String, relation: Vec<String>, tag: Option<String> }`, `Verify { passed: bool, residual_log10: f64 }`, `Record { job_id, batch, pihunt_version, started, elapsed_ms, params, columns, dropped, outcome: Kind, bound: Option<String>, iterations, relation: Option<Vec<String>>, verify: Option<Verify>, tag: Option<String>, note: Option<String> }`, `spawn_writer(&Path) -> io::Result<(Sender<Record>, JoinHandle<io::Result<()>>)>`, `read(&Path) -> Result<Vec<Record>, String>`}.

- [ ] **Step 1: Write the failing tests `tests/log.rs`**

`tests/log.rs`:

```rust
use pihunt::log::{Kind, Params, Record, read, spawn_writer};

fn record(id: &str) -> Record {
    Record {
        job_id: id.into(),
        batch: "t".into(),
        pihunt_version: "0.1.0".into(),
        started: "2026-09-22T00:00:00Z".into(),
        elapsed_ms: 5,
        params: Params {
            base: 10,
            period: 2,
            degrees: [1, 1],
            extras: vec!["log2".into()],
            coeff_bound: 1000,
            precision_digits: 80,
            gamma: 1.16,
            max_iterations: 10,
            finder: "classic".into(),
            algo_version: 1,
        },
        columns: vec!["pi".into()],
        dropped: vec![],
        outcome: Kind::Excluded,
        bound: Some("1.000e4".into()),
        iterations: 3,
        relation: None,
        verify: None,
        tag: None,
        note: None,
    }
}

#[test]
fn writes_appends_and_reads_back() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("sub/out.jsonl");
    for ids in [["a", "b"], ["c", "d"]] {
        let (tx, h) = spawn_writer(&path).unwrap();
        for id in ids {
            tx.send(record(id)).unwrap();
        }
        drop(tx);
        h.join().unwrap().unwrap();
    }
    let back = read(&path).unwrap();
    assert_eq!(
        back.iter().map(|r| r.job_id.as_str()).collect::<Vec<_>>(),
        ["a", "b", "c", "d"]
    );
    assert_eq!(back[0], record("a"));
    let line = std::fs::read_to_string(&path).unwrap();
    assert!(line.contains("\"outcome\":\"excluded\""));
}
```

- [ ] **Step 2: Run and confirm it fails**

Run: `mise exec -- cargo test --test log`
Expected: compile error, `unresolved import pihunt::log`.

- [ ] **Step 3: Write `src/log.rs` and add `pub mod log;` to `src/lib.rs`**

`src/log.rs`:

```rust
//! Append-only JSONL results log.

use serde::{Deserialize, Serialize};
use std::fs::OpenOptions;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::sync::mpsc::{Sender, channel};
use std::thread::JoinHandle;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Kind {
    Hit,
    Junk,
    Suspicious,
    Spurious,
    Excluded,
    Inconclusive,
    Skipped,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Params {
    pub base: u32,
    pub period: u32,
    pub degrees: [u32; 2],
    pub extras: Vec<String>,
    pub coeff_bound: u64,
    pub precision_digits: u32,
    pub gamma: f64,
    pub max_iterations: u64,
    pub finder: String,
    pub algo_version: u32,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DroppedRecord {
    pub column: String,
    pub relation: Vec<String>,
    pub tag: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Verify {
    pub passed: bool,
    pub residual_log10: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Record {
    pub job_id: String,
    pub batch: String,
    pub pihunt_version: String,
    pub started: String,
    pub elapsed_ms: u64,
    pub params: Params,
    pub columns: Vec<String>,
    pub dropped: Vec<DroppedRecord>,
    pub outcome: Kind,
    /// Exclusion bound reached, formatted like "1.2e5". None for skipped jobs.
    pub bound: Option<String>,
    pub iterations: u64,
    /// Aligned with `columns`; present for hit/junk/suspicious/spurious.
    pub relation: Option<Vec<String>>,
    pub verify: Option<Verify>,
    /// "known" or "NEW" for hits.
    pub tag: Option<String>,
    /// Why a job was inconclusive or skipped.
    pub note: Option<String>,
}

/// Spawn the single writer thread. Drop the sender to finish; join to surface IO errors.
pub fn spawn_writer(
    path: &Path,
) -> std::io::Result<(Sender<Record>, JoinHandle<std::io::Result<()>>)> {
    if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)?;
    }
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    let (tx, rx) = channel::<Record>();
    let handle = std::thread::spawn(move || {
        for rec in rx {
            let line = serde_json::to_string(&rec).expect("records always serialise");
            writeln!(file, "{line}")?;
            file.flush()?;
        }
        Ok(())
    });
    Ok((tx, handle))
}

pub fn read(path: &Path) -> Result<Vec<Record>, String> {
    let file = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    BufReader::new(file)
        .lines()
        .enumerate()
        .filter(|(_, l)| !matches!(l, Ok(s) if s.trim().is_empty()))
        .map(|(i, l)| {
            let l = l.map_err(|e| e.to_string())?;
            serde_json::from_str(&l).map_err(|e| format!("{}:{}: {e}", path.display(), i + 1))
        })
        .collect()
}
```

- [ ] **Step 4: Run tests and confirm they pass**

Run: `mise exec -- cargo test --release --test log`
Expected: `1 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/log.rs src/lib.rs tests/log.rs
git commit -m "feat: append-only JSONL results log

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 8: Classification

**Files:**
- Create: `src/classify.rs`
- Modify: `src/lib.rs` (add `pub mod classify;`)
- Test: `tests/classify.rs`

**Interfaces:**
- Consumes: `basis::{Columns, Shape}`, `known::is_known`, `log::{Kind, Verify}`, `pslq::{exceeds, primitive, Outcome}`, `verify::passes`.
- Produces: `classify::{Verdict { kind: Kind, relation: Option<Vec<Integer>>, verify: Option<Verify>, tag: Option<String>, note: Option<String> }, classify(&Shape, &Columns, &Outcome, relation: Option<Vec<Integer>>, coeff_bound: u64) -> Verdict}`. `relation` must already be expanded to the full column list. The returned relation is primitive with a positive π coefficient.

- [ ] **Step 1: Write the failing tests `tests/classify.rs`**

`tests/classify.rs`:

```rust
use pihunt::basis::{Columns, Shape};
use pihunt::classify::classify;
use pihunt::log::Kind;
use pihunt::pslq::Outcome;
use rug::Integer;

fn ints(v: &[i64]) -> Vec<Integer> {
    v.iter().map(|&c| Integer::from(c)).collect()
}

fn setup() -> (Shape, Columns) {
    let shape = Shape::new(16, 8, 1, 1, vec![]);
    let cols = Columns::build(&shape, 100);
    (shape, cols)
}

fn rel_outcome(v: &[i64]) -> Outcome {
    Outcome::Relation {
        coeffs: ints(v),
        iterations: 1,
        bound: 1.0,
    }
}

#[test]
fn non_relations() {
    let (shape, cols) = setup();
    let ex = Outcome::Excluded {
        bound: 5e4,
        iterations: 9,
    };
    assert_eq!(
        classify(&shape, &cols, &ex, None, 1000).kind,
        Kind::Excluded
    );
    let pe = Outcome::PrecisionExhausted {
        bound: 5.0,
        iterations: 9,
    };
    assert_eq!(
        classify(&shape, &cols, &pe, None, 1000).kind,
        Kind::Inconclusive
    );
    let ic = Outcome::IterationCap {
        bound: 5.0,
        iterations: 9,
    };
    assert_eq!(
        classify(&shape, &cols, &ic, None, 1000).kind,
        Kind::Inconclusive
    );
}

#[test]
fn bbp_is_known_hit_with_positive_pi() {
    let (shape, cols) = setup();
    let neg = [-1, 4, 0, 0, -2, -1, -1, 0, 0];
    let v = classify(&shape, &cols, &rel_outcome(&neg), Some(ints(&neg)), 1000);
    assert_eq!(v.kind, Kind::Hit);
    assert_eq!(v.tag.as_deref(), Some("known"));
    assert_eq!(v.relation, Some(ints(&[1, -4, 0, 0, 2, 1, 1, 0, 0])));
    assert!(v.verify.unwrap().passed);
}

#[test]
fn relation_without_pi_is_junk() {
    let (shape, cols) = setup();
    let r = [0, 8, -8, -4, -8, -2, -2, 1, 0];
    assert_eq!(
        classify(&shape, &cols, &rel_outcome(&r), Some(ints(&r)), 1000).kind,
        Kind::Junk
    );
}

#[test]
fn oversized_is_suspicious() {
    let (shape, cols) = setup();
    let r = [1, -4, 0, 0, 2, 1, 1, 0, 0];
    assert_eq!(
        classify(&shape, &cols, &rel_outcome(&r), Some(ints(&r)), 3).kind,
        Kind::Suspicious
    );
}

#[test]
fn false_relation_is_spurious() {
    let (shape, cols) = setup();
    let r = [1, -4, 0, 0, 2, 1, 2, 0, 0];
    let v = classify(&shape, &cols, &rel_outcome(&r), Some(ints(&r)), 1000);
    assert_eq!(v.kind, Kind::Spurious);
    assert!(!v.verify.unwrap().passed);
    assert_eq!(v.tag, None);
}
```

- [ ] **Step 2: Run and confirm it fails**

Run: `mise exec -- cargo test --test classify`
Expected: compile error, `unresolved import pihunt::classify`.

- [ ] **Step 3: Write `src/classify.rs` and add `pub mod classify;` to `src/lib.rs`**

`src/classify.rs`:

```rust
//! Turns a main-search PSLQ outcome into a log verdict.

use crate::basis::{Columns, Shape};
use crate::known;
use crate::log::{Kind, Verify};
use crate::pslq::{Outcome, exceeds};
use crate::verify;
use rug::Integer;

#[derive(Debug, Clone, PartialEq)]
pub struct Verdict {
    pub kind: Kind,
    /// Aligned with the full column list; normalised so the π coefficient is positive.
    pub relation: Option<Vec<Integer>>,
    pub verify: Option<Verify>,
    pub tag: Option<String>,
    pub note: Option<String>,
}

/// `relation` from a main-search `Outcome::Relation` must already be expanded to full columns.
pub fn classify(
    shape: &Shape,
    cols: &Columns,
    outcome: &Outcome,
    relation: Option<Vec<Integer>>,
    coeff_bound: u64,
) -> Verdict {
    let plain = |kind, note: Option<&str>| Verdict {
        kind,
        relation: None,
        verify: None,
        tag: None,
        note: note.map(String::from),
    };
    match outcome {
        Outcome::Excluded { .. } => plain(Kind::Excluded, None),
        Outcome::PrecisionExhausted { .. } => {
            plain(Kind::Inconclusive, Some("precision exhausted"))
        }
        Outcome::IterationCap { .. } => plain(Kind::Inconclusive, Some("iteration cap")),
        Outcome::Relation { .. } => {
            let rel = normalise(relation.expect("relation outcome carries coefficients"));
            let mut v = Verdict {
                kind: Kind::Junk,
                relation: None,
                verify: None,
                tag: None,
                note: None,
            };
            if rel[0] == 0 {
                v.note = Some("relation without pi after reduction — bug signal".into());
            } else if exceeds(&rel, coeff_bound) {
                v.kind = Kind::Suspicious;
                v.note = Some("coefficient exceeds coeff_bound".into());
            } else {
                let (passed, residual_log10) = verify::passes(&rel, &cols.hi, cols.digits);
                v.verify = Some(Verify {
                    passed,
                    residual_log10,
                });
                v.kind = if passed { Kind::Hit } else { Kind::Spurious };
                if passed {
                    let known = known::is_known(shape, &cols.names, &rel);
                    v.tag = Some(if known { "known" } else { "NEW" }.into());
                }
            }
            v.relation = Some(rel);
            v
        }
    }
}

/// Divide by gcd; make the π coefficient positive (or the first nonzero one if π is absent).
fn normalise(rel: Vec<Integer>) -> Vec<Integer> {
    let mut p = crate::pslq::primitive(&rel);
    if p[0] < 0 {
        p.iter_mut().for_each(|c| *c = Integer::from(-&*c));
    }
    p
}
```

- [ ] **Step 4: Run tests and confirm they pass**

Run: `mise exec -- cargo test --release --test classify`
Expected: `5 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/classify.rs src/lib.rs tests/classify.rs
git commit -m "feat: classify PSLQ outcomes into log verdicts

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 9: End-to-end job runner

**Files:**
- Create: `src/job.rs`
- Modify: `src/lib.rs` (add `pub mod job;`)
- Test: `tests/job.rs`

**Interfaces:**
- Consumes: everything above.
- Produces: `job::{run_job(&Job, batch: &str, &dyn RelationFinder, max_columns: usize) -> Record, fmt_bound(f64) -> String}`.

- [ ] **Step 1: Write the failing tests `tests/job.rs`**

This is where BBP and Bailey's π² formula get rediscovered through the real pipeline.

`tests/job.rs`:

```rust
use pihunt::basis::{Extra, Shape, auto_digits};
use pihunt::job::run_job;
use pihunt::log::Kind;
use pihunt::plan::Job;
use pihunt::pslq::classic::ClassicPslq;

fn job(shape: Shape) -> Job {
    let digits = auto_digits(shape.columns(), 1000);
    Job {
        shape,
        coeff_bound: 1000,
        digits,
        gamma: 1.16,
        max_iterations: 1_000_000,
    }
}

#[test]
fn rediscovers_bbp() {
    let rec = run_job(&job(Shape::new(16, 8, 1, 1, vec![])), "t", &ClassicPslq, 80);
    assert_eq!(rec.outcome, Kind::Hit);
    assert_eq!(rec.tag.as_deref(), Some("known"));
    let rel: Vec<&str> = rec
        .relation
        .as_ref()
        .unwrap()
        .iter()
        .map(|s| s.as_str())
        .collect();
    assert_eq!(rel, ["1", "-4", "0", "0", "2", "1", "1", "0", "0"]);
    assert_eq!(rec.dropped[0].column, "S(j=7,s=1)");
    assert_eq!(rec.dropped[0].tag.as_deref(), Some("known"));
}

#[test]
fn rediscovers_bailey_pi_squared_as_basis_relation() {
    let rec = run_job(
        &job(Shape::new(64, 6, 2, 2, vec![Extra::Pi2])),
        "t",
        &ClassicPslq,
        80,
    );
    let d = rec
        .dropped
        .iter()
        .find(|d| d.column == "pi2")
        .expect("pi2 dropped");
    assert_eq!(
        d.relation,
        ["0", "144", "-216", "-72", "-54", "9", "0", "-8"]
    );
    assert_eq!(d.tag.as_deref(), Some("known"));
    assert_eq!(rec.outcome, Kind::Excluded);
}

#[test]
fn skips_oversized_jobs() {
    let rec = run_job(&job(Shape::new(10, 8, 1, 1, vec![])), "t", &ClassicPslq, 5);
    assert_eq!(rec.outcome, Kind::Skipped);
    assert!(rec.note.unwrap().contains("max_columns"));
}
```

- [ ] **Step 2: Run and confirm it fails**

Run: `mise exec -- cargo test --test job`
Expected: compile error, `unresolved import pihunt::job`.

- [ ] **Step 3: Write `src/job.rs` and add `pub mod job;` to `src/lib.rs`**

`src/lib.rs` should now be exactly:

`src/lib.rs`:

```rust
pub mod basis;
pub mod classify;
pub mod config;
pub mod job;
pub mod known;
pub mod log;
pub mod plan;
pub mod pslq;
pub mod reduce;
pub mod verify;
```

`src/job.rs`:

```rust
//! Runs one job end to end: basis → reduce → main PSLQ → classify → Record.

use crate::basis::Columns;
use crate::classify::classify;
use crate::known;
use crate::log::{DroppedRecord, Kind, Params, Record};
use crate::plan::{ALGO_VERSION, Job};
use crate::pslq::{Outcome, RelationFinder};
use crate::reduce::{Dropped, Reduced, reduce};
use rug::{Float, Integer};
use std::time::Instant;

pub fn run_job(job: &Job, batch: &str, finder: &dyn RelationFinder, max_columns: usize) -> Record {
    let started = jiff::Timestamp::now().to_string();
    let t = Instant::now();
    let names = job.shape.column_names();
    let mut rec = Record {
        job_id: job.id(finder.name()),
        batch: batch.to_string(),
        pihunt_version: env!("CARGO_PKG_VERSION").to_string(),
        started,
        elapsed_ms: 0,
        params: params(job, finder),
        columns: names.clone(),
        dropped: vec![],
        outcome: Kind::Skipped,
        bound: None,
        iterations: 0,
        relation: None,
        verify: None,
        tag: None,
        note: None,
    };
    if names.len() > max_columns {
        rec.note = Some(format!(
            "{} columns > max_columns {max_columns}",
            names.len()
        ));
        return rec;
    }

    let cols = Columns::build(&job.shape, job.digits);
    let p = job.params();
    let (keep, dropped) = match reduce(&cols, finder, &p) {
        Reduced::Done { keep, dropped } => (keep, dropped),
        Reduced::Inconclusive { outcome, dropped } => {
            rec.dropped = dropped_records(job, &names, &dropped);
            rec.outcome = Kind::Inconclusive;
            rec.bound = Some(fmt_bound(outcome.bound()));
            rec.iterations = outcome.iterations();
            rec.note = Some(match &outcome {
                Outcome::Relation { coeffs, .. } => {
                    let max = coeffs
                        .iter()
                        .map(|c| c.as_abs().to_string())
                        .max_by_key(|s| (s.len(), s.clone()));
                    format!(
                        "reduction found a relation that failed checks (max |coeff| {})",
                        max.unwrap_or_default()
                    )
                }
                Outcome::PrecisionExhausted { .. } => "reduction: precision exhausted".into(),
                _ => "reduction: iteration cap".into(),
            });
            rec.elapsed_ms = t.elapsed().as_millis() as u64;
            return rec;
        }
    };
    rec.dropped = dropped_records(job, &names, &dropped);

    let mut idx = vec![0];
    idx.extend(keep);
    let x: Vec<Float> = idx.iter().map(|&i| cols.lo[i].clone()).collect();
    let outcome = finder.find(&x, &p);
    let relation = match &outcome {
        Outcome::Relation { coeffs, .. } => {
            let mut full = vec![Integer::new(); names.len()];
            for (&i, c) in idx.iter().zip(coeffs) {
                full[i] = c.clone();
            }
            Some(full)
        }
        _ => None,
    };
    let v = classify(&job.shape, &cols, &outcome, relation, job.coeff_bound);
    rec.outcome = v.kind;
    rec.bound = Some(fmt_bound(outcome.bound()));
    rec.iterations = outcome.iterations();
    rec.relation = v
        .relation
        .map(|r| r.iter().map(|c| c.to_string()).collect());
    rec.verify = v.verify;
    rec.tag = v.tag;
    rec.note = v.note;
    rec.elapsed_ms = t.elapsed().as_millis() as u64;
    rec
}

fn params(job: &Job, finder: &dyn RelationFinder) -> Params {
    let s = &job.shape;
    Params {
        base: s.base,
        period: s.period,
        degrees: [s.s_lo, s.s_hi],
        extras: s.extras.iter().map(|e| e.name().to_string()).collect(),
        coeff_bound: job.coeff_bound,
        precision_digits: job.digits,
        gamma: job.gamma,
        max_iterations: job.max_iterations,
        finder: finder.name().to_string(),
        algo_version: ALGO_VERSION,
    }
}

fn dropped_records(job: &Job, names: &[String], dropped: &[Dropped]) -> Vec<DroppedRecord> {
    dropped
        .iter()
        .map(|d| DroppedRecord {
            column: names[d.column].clone(),
            relation: d.relation.iter().map(|c| c.to_string()).collect(),
            tag: known::is_known(&job.shape, names, &d.relation).then(|| "known".to_string()),
        })
        .collect()
}

pub fn fmt_bound(b: f64) -> String {
    format!("{b:.3e}")
}
```

- [ ] **Step 4: Run tests and confirm they pass**

Run: `mise exec -- cargo test --release --test job`
Expected: `3 passed`.

- [ ] **Step 5: Commit**

```bash
git add src/job.rs src/lib.rs tests/job.rs
git commit -m "feat: end-to-end job runner

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 10: CLI and batch files

**Files:**
- Modify: `src/main.rs` (replace the stub)
- Create: `batches/known.toml`, `batches/base10-scout.toml`
- Test: `tests/cli.rs`

**Interfaces:**
- Consumes: `config::load`, `plan::plan`, `job::run_job`, `log::{spawn_writer, read, Kind, Record}`, `basis::{Columns, Extra, Shape}`, `verify::passes`, `pslq::classic::ClassicPslq`.
- Produces: the binary `pihunt` with `run <batch.toml>`, `plan <batch.toml>` and `verify <results.jsonl>`. Errors print `error: …` to stderr and exit non-zero. `verify` exits non-zero if any hit fails.

- [ ] **Step 1: Write the failing tests `tests/cli.rs`**

`tests/cli.rs`:

```rust
use std::process::Command;

const BATCH: &str = r#"
name   = "cli"
output = "OUT"

[defaults]
coeff_bound = 1000

[search]
mode    = "grid"
bases   = [16]
periods = { from = 7, to = 8 }
degrees = [[1, 1]]
extras  = [[]]
"#;

#[test]
fn run_then_verify() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("r.jsonl");
    let batch = dir.path().join("b.toml");
    std::fs::write(&batch, BATCH.replace("OUT", out.to_str().unwrap())).unwrap();
    let bin = env!("CARGO_BIN_EXE_pihunt");

    let plan = Command::new(bin)
        .args(["plan", batch.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(plan.status.success());
    assert!(String::from_utf8_lossy(&plan.stdout).contains("2 jobs"));

    let run = Command::new(bin)
        .args(["run", batch.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let stdout = String::from_utf8_lossy(&run.stdout);
    assert!(stdout.contains("hit [known] b=16 m=8"), "{stdout}");
    assert_eq!(std::fs::read_to_string(&out).unwrap().lines().count(), 2);

    let verify = Command::new(bin)
        .args(["verify", out.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(verify.status.success());
    assert!(String::from_utf8_lossy(&verify.stdout).contains("1 hits checked, 0 failed"));
}

#[test]
fn bad_batch_fails_cleanly() {
    let dir = tempfile::tempdir().unwrap();
    let batch = dir.path().join("b.toml");
    std::fs::write(&batch, "name = 1").unwrap();
    let run = Command::new(env!("CARGO_BIN_EXE_pihunt"))
        .args(["run", batch.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!run.status.success());
    assert!(String::from_utf8_lossy(&run.stderr).starts_with("error:"));
}
```

- [ ] **Step 2: Run and confirm it fails**

Run: `mise exec -- cargo test --test cli`
Expected: FAIL. The stub binary prints nothing, so the `plan` assertion fails.

- [ ] **Step 3: Replace `src/main.rs`**

`src/main.rs`:

```rust
use clap::{Parser, Subcommand};
use pihunt::basis::{Columns, Extra, Shape};
use pihunt::log::{Kind, Record};
use pihunt::pslq::classic::ClassicPslq;
use pihunt::{config, job, log, plan, verify};
use rayon::prelude::*;
use rug::Integer;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Parser)]
#[command(version, about = "Hunt for BBP-type formulas for pi with PSLQ")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

#[derive(Subcommand)]
enum Cmd {
    /// Run every job in a batch file, appending results to its output log.
    Run { batch: PathBuf },
    /// Dry run: show how big a batch is without running it.
    Plan { batch: PathBuf },
    /// Re-verify every hit in a results log at 2x precision.
    Verify { results: PathBuf },
}

fn main() -> ExitCode {
    let result = match Cli::parse().cmd {
        Cmd::Run { batch } => run(batch),
        Cmd::Plan { batch } => show_plan(batch),
        Cmd::Verify { results } => reverify(results),
    };
    result.unwrap_or_else(|e| {
        eprintln!("error: {e}");
        ExitCode::FAILURE
    })
}

fn run(path: PathBuf) -> Result<ExitCode, String> {
    let batch = config::load(&path)?;
    let jobs = plan::plan(&batch);
    let threads = batch
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| e.to_string())?;
    let (tx, writer) = log::spawn_writer(&batch.output).map_err(|e| e.to_string())?;
    let done = AtomicUsize::new(0);
    let total = jobs.len();
    eprintln!(
        "{}: {total} jobs on {threads} threads → {}",
        batch.name,
        batch.output.display()
    );

    let records: Vec<Record> = pool.install(|| {
        jobs.par_iter()
            .map_with(tx, |tx, j| {
                let rec = job::run_job(j, &batch.name, &ClassicPslq, batch.defaults.max_columns);
                let n = done.fetch_add(1, Ordering::Relaxed) + 1;
                eprintln!(
                    "[{n}/{total}] {} → {:?} ({} ms)",
                    describe(&rec),
                    rec.outcome,
                    rec.elapsed_ms
                );
                tx.send(rec.clone()).expect("writer thread alive");
                rec
            })
            .collect()
    });
    writer
        .join()
        .expect("writer thread panicked")
        .map_err(|e| e.to_string())?;
    summarise(&records);
    Ok(ExitCode::SUCCESS)
}

fn show_plan(path: PathBuf) -> Result<ExitCode, String> {
    let batch = config::load(&path)?;
    let jobs = plan::plan(&batch);
    let max_cols = batch.defaults.max_columns;
    let (skipped, run): (Vec<_>, Vec<_>) = jobs.iter().partition(|j| j.shape.columns() > max_cols);
    println!(
        "{}: {} jobs ({} to run, {} skipped over max_columns {max_cols})",
        batch.name,
        jobs.len(),
        run.len(),
        skipped.len()
    );
    if let (Some(n), Some(lo), Some(hi)) = (
        run.iter().map(|j| j.shape.columns()).max(),
        run.iter().map(|j| j.digits).min(),
        run.iter().map(|j| j.digits).max(),
    ) {
        println!("largest n = {n}, precision {lo}..={hi} digits");
    }
    Ok(ExitCode::SUCCESS)
}

fn reverify(path: PathBuf) -> Result<ExitCode, String> {
    let records = log::read(&path)?;
    let mut failed = 0;
    let hits: Vec<&Record> = records.iter().filter(|r| r.outcome == Kind::Hit).collect();
    for rec in &hits {
        let p = &rec.params;
        let extras = p
            .extras
            .iter()
            .map(|e| e.parse::<Extra>())
            .collect::<Result<Vec<_>, _>>()?;
        let shape = Shape::new(p.base, p.period, p.degrees[0], p.degrees[1], extras);
        let relation = rec.relation.as_ref().ok_or("hit without relation")?;
        let relation = relation
            .iter()
            .map(|c| c.parse::<Integer>().map_err(|e| e.to_string()))
            .collect::<Result<Vec<_>, _>>()?;
        let cols = Columns::build(&shape, p.precision_digits);
        let (ok, r) = verify::passes(&relation, &cols.hi, p.precision_digits);
        println!(
            "{} {} residual 1e{r:.0}",
            if ok { "OK  " } else { "FAIL" },
            describe(rec)
        );
        failed += !ok as usize;
    }
    println!("{} hits checked, {failed} failed", hits.len());
    Ok(if failed == 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    })
}

fn describe(rec: &Record) -> String {
    let p = &rec.params;
    let x = if p.extras.is_empty() {
        "-".to_string()
    } else {
        p.extras.join(",")
    };
    format!(
        "b={} m={} s={}..{} x={x}",
        p.base, p.period, p.degrees[0], p.degrees[1]
    )
}

fn formula(rec: &Record) -> String {
    let Some(rel) = &rec.relation else {
        return String::new();
    };
    rec.columns
        .iter()
        .zip(rel)
        .filter(|(_, c)| c.as_str() != "0")
        .map(|(n, c)| format!("{c}·{n}"))
        .collect::<Vec<_>>()
        .join(" + ")
        + " = 0"
}

fn summarise(records: &[Record]) {
    use Kind::*;
    println!("\n== summary ==");
    for kind in [
        Hit,
        Excluded,
        Inconclusive,
        Skipped,
        Suspicious,
        Spurious,
        Junk,
    ] {
        let n = records.iter().filter(|r| r.outcome == kind).count();
        if n > 0 {
            println!("{kind:?}: {n}");
        }
    }
    for rec in records.iter().filter(|r| r.outcome == Hit) {
        let tag = rec.tag.as_deref().unwrap_or("?");
        let banner = if tag == "NEW" {
            "!!!!!!!!!! NEW HIT !!!!!!!!!!\n"
        } else {
            ""
        };
        println!("{banner}hit [{tag}] {}: {}", describe(rec), formula(rec));
    }
}
```

- [ ] **Step 4: Create `batches/known.toml`**

`batches/known.toml`:

```toml
name   = "known"
output = "results/known.jsonl"

[defaults]
coeff_bound = 1000

[search]
mode    = "grid"
bases   = [16, 64]
periods = { from = 6, to = 8 }
degrees = [[1, 1], [2, 2]]
extras  = [[], ["pi2"]]
```

- [ ] **Step 5: Create `batches/base10-scout.toml`**

`batches/base10-scout.toml`:

```toml
name   = "base10-scout"
output = "results/base10-scout.jsonl"

[defaults]
coeff_bound = 1000
max_columns = 50

[search]
mode    = "sample"
bases   = [10, 100, 1000]
periods = { from = 2, to = 16 }
degrees = [[1, 1], [1, 2], [2, 2]]
extras  = [[], ["log2", "log5"], ["pi2", "log2", "log3", "log5", "catalan", "zeta3"]]

[sample]
method = "sobol"
count  = 200
seed   = 42
```

- [ ] **Step 6: Run the full suite plus clippy**

Run: `mise exec -- cargo test --release && mise exec -- cargo clippy --all-targets`
Expected: 37 tests pass across all test files, and clippy prints no warnings.

- [ ] **Step 7: Commit**

```bash
git add src/main.rs tests/cli.rs batches
git commit -m "feat: pihunt CLI (run/plan/verify) and batch files

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```

---

### Task 11: Stage 1 acceptance run

**Files:**
- Create: `tests/timing.rs`, `docs/timing-baseline.md`
- Create (generated): `results/known.jsonl`, `results/base10-scout.jsonl`

This task checks the spec's Stage 1 done criteria against real runs. There's no new library code.

- [ ] **Step 1: Build the release binary**

Run: `mise exec -- cargo build --release`
Binary: `~/.cache/pihunt/target/release/pihunt` (referred to below as `$PIHUNT`).

- [ ] **Step 2: Known-formula batch (done criteria 2 and 3)**

Run: `$PIHUNT run batches/known.toml`
Expected summary: `Hit: 2`, `Excluded: 22`, and both hits printed as `hit [known] b=16 m=8 …: 1·pi + -4·S(j=1,s=1) + 2·S(j=4,s=1) + 1·S(j=5,s=1) + 1·S(j=6,s=1) = 0`.
Then: `grep '"column":"pi2"' results/known.jsonl` must show the relation `["0","144","-216","-72","-54","9","0","-8"]` with `"tag":"known"`.

- [ ] **Step 3: Base-10 scout (done criterion 4)**

Run: `$PIHUNT plan batches/base10-scout.toml`, then `$PIHUNT run batches/base10-scout.toml`
Expected: `plan` reports 176 jobs (200 samples, deduplicated). `run` completes unattended in about 15 s on 6 cores. The prototype got 175 `Excluded` and 1 `Inconclusive` (with a `note` saying a reduction relation failed checks). **If any `NEW` hit appears, stop and report it to the human before doing anything else.**

- [ ] **Step 4: Re-verify (done criterion 5)**

Run: `$PIHUNT verify results/known.jsonl && $PIHUNT verify results/base10-scout.jsonl`
Expected: `2 hits checked, 0 failed` and `0 hits checked, 0 failed` (unless the scout found hits), both exiting 0.

- [ ] **Step 5: Add the timing harness `tests/timing.rs` (done criterion 6)**

`tests/timing.rs`:

```rust
//! Stage-1 timing baseline. Run with:
//! cargo test --release --test timing -- --ignored --nocapture

use pihunt::basis::{Extra, Shape, auto_digits};
use pihunt::job::run_job;
use pihunt::plan::Job;
use pihunt::pslq::classic::ClassicPslq;

#[test]
#[ignore]
fn timing_baseline() {
    let all = [
        Extra::Pi2,
        Extra::Log2,
        Extra::Log3,
        Extra::Log5,
        Extra::Catalan,
        Extra::Zeta3,
    ];
    let shapes = [
        Shape::new(10, 8, 1, 1, vec![]),
        Shape::new(10, 12, 1, 1, all.to_vec()),
        Shape::new(10, 12, 1, 2, all[..2].to_vec()),
        Shape::new(10, 16, 1, 2, all.to_vec()),
        Shape::new(10, 20, 1, 2, all.to_vec()),
    ];
    println!("| n | digits | outcome | iterations | seconds |");
    println!("|---|---|---|---|---|");
    for shape in shapes {
        let n = shape.columns();
        let digits = auto_digits(n, 1000);
        let job = Job {
            shape,
            coeff_bound: 1000,
            digits,
            gamma: 1.16,
            max_iterations: 10_000_000,
        };
        let rec = run_job(&job, "timing", &ClassicPslq, 200);
        println!(
            "| {n} | {digits} | {:?} | {} | {:.2} |",
            rec.outcome,
            rec.iterations,
            rec.elapsed_ms as f64 / 1000.0
        );
    }
}
```

- [ ] **Step 6: Record the baseline**

Run: `mise exec -- cargo test --release --test timing -- --ignored --nocapture`
Write the printed markdown table into `docs/timing-baseline.md`, under a heading that states the date, CPU (`lscpu | grep 'Model name'`), and the `pihunt` git commit. For reference, the prototype measured 0.01 s at n=9, 0.18 s at n=19, 0.90 s at n=27, 5.3 s at n=39, and 5.6 s at n=47 (inconclusive).

- [ ] **Step 7: Commit**

```bash
git add tests/timing.rs docs/timing-baseline.md results
git commit -m "test: stage 1 acceptance runs and timing baseline

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
```
