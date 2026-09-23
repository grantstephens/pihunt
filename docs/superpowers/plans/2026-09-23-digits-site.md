# Digits Site Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a static site that presents the reconstructed Gourdon Theorem 2 as a paper, together with a live in-browser demo. The demo races Theorem 1 against Theorem 2, both compiled to WASM from the same verified code the CLI uses, and streams digits.

**Architecture:** Move the two digit algorithms into a new workspace crate, `crates/pi-digits`. Its big-number arithmetic sits behind one small `Big` seam with two backends: `gmp` (rug, used by the CLI) and `pure` (dashu-int, used by WASM). Rayon goes behind a `parallel` feature. A `site/wasm` wasm-bindgen crate wraps it. A hand-written static page in `site/src` runs one Web Worker per algorithm. `site/build.sh` produces `site/dist/`, which the user deploys with `wrangler pages deploy`.

**Tech Stack:** Rust 2024 (stable via rustup/mise), rug 1.30, dashu-int, rayon, wasm-bindgen + wasm-pack, binaryen `wasm-opt`, Node 26 (built-ins only, no npm deps), KaTeX (vendored), Inter + JetBrains Mono (self-hosted woff2).

**Spec:** `docs/superpowers/specs/2026-09-23-digits-site-design.md`

## Global Constraints

- **Memory safety (this machine has 7.7 GiB):**
  - Run heavy commands (digit positions ≥ 10⁵, ignored tests, benchmarks) one at a time under `( ulimit -v 3000000; exec <cmd> )`.
  - Never run `cargo test -- --include-ignored` across the whole workspace.
  - Always pass `-j 4` to cargo build/test.
- **Build env:** do NOT use `mise exec`. Use `export PATH=$HOME/.cargo/bin:$PATH CARGO_TARGET_DIR=$HOME/.cache/pihunt/target-site`. The worktree lives under `~/sync`, which is synced, so nothing may build into it.
- **Behaviour must not change:**
  - With the default features, every existing test passes unchanged: 102 non-ignored tests at the start of this plan.
  - `pihunt digit` output is byte-identical at positions 1, 762, 10⁴, 10⁵, 10⁶. The strings are in `docs/nthdigit.md` "Position convention".
  - The certification code must not move behind the seam: u128 fixed-point accumulation, `error_units`, the digit-boundary retry loop, and the exact counting of rounded terms.
- **Features:** `pi-digits` has features `gmp` (default), `pure`, `parallel` (default). Exactly one of `gmp`/`pure`, otherwise `compile_error!`.
- **Positions:** the CLI and WASM API are 1-based (`pos`). The library is `digits(n, count)` with `n = pos − 1`, returning positions `n+1..=n+count`.
- **Site:**
  - No third-party requests at runtime: KaTeX, fonts, WASM and data are all self-hosted.
  - Fonts: Inter (body and headings, fallback `"Helvetica Neue", Helvetica, Arial, system-ui, sans-serif`) and JetBrains Mono (digits and code).
  - Light and dark themes; no horizontal scroll at 360 px width.
- **Demo limits:** presets 1 000 / 10 000 / 50 000 / 100 000; soft warning above 10⁵; hard cap at 2·10⁵. Browser reference file `pi-200k.txt` holds positions 1..=200 000.
- **Deploy:** the user runs `npx wrangler pages deploy site/dist --project-name <their-name>`. `site/dist/` is gitignored.
- **Commits:** every commit message ends with exactly:
  ```
  Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>
  Claude-Session: https://claude.ai/code/session_01ALoLJ3AsymWwAtULqoePob
  ```

## Review Focus

1. **Bad position input** (0, negative, non-integer, empty, above 2·10⁵, `1e5` scientific notation): the demo shows an inline validation message and never starts a worker. *Pinned in Task 6.*
2. **Positions below 2000 take the fallback path** (`pi_ref`, not the Gourdon algorithms) in both backends and in WASM: digits at positions 1..=2000 must match the reference. *Pinned in Tasks 3 and 5.*
3. **Cancel mid-race, then race again:** a late message from a terminated worker must never paint into the new run's lanes. Use a run id (generation counter). *Pinned in Task 6.*
4. **The lanes disagree** (a real bug in one backend): the UI must say so loudly, not show two quiet results. *Pinned in Task 6* via a test hook that forces a mismatch.
5. **Page without JS or without WASM support:** the paper content is fully readable, and the demo area shows a plain explanatory message. *Pinned in Task 8* via the check script plus a manual check.

---

### Task 1: Extract `crates/pi-digits` (pure refactor, gmp backend only)

**Files:**
- Modify: `Cargo.toml` (root). Add `[workspace] members = [".", "crates/pi-digits"]` and `resolver = "3"`. Add dependency `pi-digits = { path = "crates/pi-digits" }`. Move the `nthdigit-profile` and `mem-profile` features so they forward: `nthdigit-profile = ["pi-digits/nthdigit-profile"]`, `mem-profile = ["pi-digits/mem-profile"]`.
- Create: `crates/pi-digits/Cargo.toml`: deps `rug` (same features as root), `rayon`; features `default = ["gmp", "parallel"]`, `gmp = []`, `pure = []`, `parallel = []`, `nthdigit-profile = []`, `mem-profile = []`.
- Move (`git mv`): `src/nthdigit.rs` → `crates/pi-digits/src/nthdigit.rs`; `src/nthdigit2.rs` → `crates/pi-digits/src/nthdigit2.rs`; `src/mem_profile.rs` → `crates/pi-digits/src/mem_profile.rs`; `tests/nthdigit.rs` → `crates/pi-digits/tests/nthdigit.rs`; `tests/nthdigit2.rs` → `crates/pi-digits/tests/nthdigit2.rs`.
- Create: `crates/pi-digits/src/lib.rs`
- Modify: `src/lib.rs`. Replace `pub mod mem_profile; pub mod nthdigit; pub mod nthdigit2;` with `pub use pi_digits::{mem_profile, nthdigit, nthdigit2};`.
- Modify: the moved sources and tests. `crate::pslq::digits_to_bits` doesn't exist in the new crate, so add `pub fn digits_to_bits(d: u32) -> u32 { (d as f64 * std::f64::consts::LOG2_10).ceil() as u32 }` to `crates/pi-digits/src/lib.rs` and point the moved code at `crate::digits_to_bits`. In the moved tests, replace `pihunt::` with `pi_digits::`.
- The CLI tests inside the moved test files use `env!("CARGO_BIN_EXE_pihunt")`, which only exists in the `pihunt` package. Move those CLI test functions into a new `tests/cli_digits.rs` in the root package.

**Interfaces:**
- Produces: crate `pi_digits` with modules `nthdigit`, `nthdigit2`, `mem_profile` (feature-gated as before), and fn `digits_to_bits(u32) -> u32`. `pihunt::{nthdigit, nthdigit2}` paths keep working through the re-export.

- [ ] **Step 1: Baseline.** Run `cargo test --release -j 4 2>&1 | grep "test result"`. Record the passed and ignored totals; expect 102 passed, 8 ignored. Record `pihunt digit` output at positions 1, 762 (`--count 8`), 10000, 100000.
- [ ] **Step 2: Create the crate skeleton.** Write `crates/pi-digits/src/lib.rs`:

```rust
//! Low-memory decimal digit extraction for π: Gourdon's Theorem 1 (`nthdigit`) and the
//! reconstructed Theorem 2 (`nthdigit2`). Shared by the `pihunt` CLI and the site's WASM demo.

#[cfg(feature = "mem-profile")]
pub mod mem_profile;
pub mod nthdigit;
pub mod nthdigit2;

/// Decimal digits → bits, rounded up.
pub fn digits_to_bits(digits: u32) -> u32 {
    (digits as f64 * std::f64::consts::LOG2_10).ceil() as u32
}
```

  Keep whatever `cfg` gating `mem_profile` had in the old `src/lib.rs`. Check with `git show HEAD:src/lib.rs`.
- [ ] **Step 3: Move the files** with `git mv` as listed, fix the paths, and add the workspace config.
- [ ] **Step 4: Run the full suite.** `cargo test --release -j 4 --workspace 2>&1 | grep "test result" | awk '{p+=$4; f+=$6; i+=$8} END {print p, f, i}'`. Expected: the same totals as Step 1, and 0 failed. Run `cargo clippy --workspace --all-targets -j 4`: zero warnings. Re-run the four `pihunt digit` commands and diff against Step 1: identical.
- [ ] **Step 5: Commit.** Message: `refactor: move digit extraction into crates/pi-digits workspace crate`.

---

### Task 2: The `Big` seam (gmp backend) and a bignum-free Theorem 1

**Files:**
- Create: `crates/pi-digits/src/bignum.rs`
- Modify: `crates/pi-digits/src/lib.rs` (add `pub mod bignum;` and the feature guard)
- Modify: `crates/pi-digits/src/nthdigit2.rs`. Replace every non-test `rug::Integer` with `Big`: the functions `bs`, `advance`, the product tree in the ART, and the `PadicBinom` binomial coefficients at the two `Integer::from(k).binomial(i)` sites.
- Modify: `crates/pi-digits/src/nthdigit.rs`. Replace `error_units`' `Integer` arithmetic with exact u128 arithmetic.
- Test: `crates/pi-digits/tests/bignum.rs`

**Interfaces:**
- Produces, in `pi_digits::bignum`: `#[derive(Clone, Debug, PartialEq, Eq)] pub struct Big(..)` with:
  - `Big::from_u64(u64) -> Big`, `Big::zero()`, `Big::one()`
  - `mul(&self, &Big) -> Big`, `add(&self, &Big) -> Big`, `sub(&self, &Big) -> Big` (the caller guarantees `self >= other`)
  - `rem(&self, &Big) -> Big`, `rem_u64(&self, u64) -> u64`, `div_u64_exact(&self, u64) -> Big` (panics in debug if not exact), `mul_u64(&self, u64) -> Big`
  - `to_u64(&self) -> Option<u64>`, `is_zero(&self) -> bool`, `bits(&self) -> u64`
  - `pow_u64(base: u64, exp: u32) -> Big`, `binomial(n: u64, k: u32) -> Big`, `to_decimal_string(&self) -> String`, `shl(&self, u32) -> Big`
  - Everything else in the crate must use only these.

- [ ] **Step 1: Write failing tests** in `crates/pi-digits/tests/bignum.rs`:

```rust
use pi_digits::bignum::Big;

fn lcg(seed: &mut u64) -> u64 {
    *seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
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
        let (a, b, c) = (rand_big(&mut s, words), rand_big(&mut s, words), rand_big(&mut s, words / 2 + 1));
        assert_eq!(a.mul(&b).add(&a.mul(&c)), a.mul(&b.add(&c)), "distributive, words={words}");
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
    assert_eq!(Big::pow_u64(10, 3).to_u64(), Some(1000));
    assert_eq!(Big::binomial(10, 3).to_u64(), Some(120));
    assert_eq!(Big::binomial(100_000_000, 4).to_decimal_string(), "416666641666667083333325000000");
    assert_eq!(Big::pow_u64(2, 70).to_decimal_string(), "1180591620717411303424");
    assert!(Big::zero().is_zero());
    assert_eq!(Big::one().mul_u64(7), Big::from_u64(7));
}
```

  Also add a direct cross-check against GMP. `rug` stays a dev-dependency under both backends (Task 4), so under `pure` this compares dashu against GMP:

```rust
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
    }
}
```

  The expected `binomial(10^8, 4)` string is the value used by `PadicBinom` at realistic N. Before finalising the test, check it independently: `python3 -c "import math;print(math.comb(10**8,4))"`.
- [ ] **Step 2: Run it** with `cargo test --release -j 4 -p pi-digits --test bignum`. Expected: compile error, `bignum` doesn't exist.
- [ ] **Step 3: Implement `bignum.rs` for gmp.** Guard the crate root: `#[cfg(all(feature = "gmp", feature = "pure"))] compile_error!("enable exactly one of `gmp` or `pure`");` and `#[cfg(not(any(feature = "gmp", feature = "pure")))] compile_error!(...)`. The gmp implementation wraps `rug::Integer`:

```rust
//! The only big-integer surface the digit algorithms use. Two backends: `gmp` (rug) for the
//! CLI, `pure` (dashu-int) for WASM. Keep this list minimal: every method here must exist, with
//! identical semantics, in both backends, and `tests/bignum.rs` runs against both.

#[cfg(feature = "gmp")]
mod imp {
    use rug::{Integer, ops::Pow};

    #[derive(Clone, Debug, PartialEq, Eq)]
    pub struct Big(Integer);

    impl Big {
        pub fn from_u64(v: u64) -> Big { Big(Integer::from(v)) }
        pub fn zero() -> Big { Big(Integer::new()) }
        pub fn one() -> Big { Big(Integer::from(1)) }
        pub fn mul(&self, o: &Big) -> Big { Big(Integer::from(&self.0 * &o.0)) }
        pub fn mul_u64(&self, v: u64) -> Big { Big(Integer::from(&self.0 * v)) }
        pub fn add(&self, o: &Big) -> Big { Big(Integer::from(&self.0 + &o.0)) }
        pub fn sub(&self, o: &Big) -> Big { Big(Integer::from(&self.0 - &o.0)) }
        pub fn rem(&self, m: &Big) -> Big { Big(Integer::from(&self.0 % &m.0)) }
        pub fn rem_u64(&self, m: u64) -> u64 { Integer::from(&self.0 % m).to_u64().expect("rem < m") }
        pub fn div_u64_exact(&self, d: u64) -> Big {
            debug_assert!(self.0.is_divisible_u(d as u32) || Integer::from(&self.0 % d) == 0);
            Big(Integer::from(&self.0 / d))
        }
        pub fn to_u64(&self) -> Option<u64> { self.0.to_u64() }
        pub fn is_zero(&self) -> bool { self.0 == 0 }
        pub fn bits(&self) -> u64 { self.0.significant_bits() as u64 }
        pub fn shl(&self, k: u32) -> Big { Big(Integer::from(&self.0 << k)) }
        pub fn pow_u64(base: u64, exp: u32) -> Big { Big(Integer::from(base).pow(exp)) }
        pub fn binomial(n: u64, k: u32) -> Big { Big(Integer::from(n).binomial(k)) }
        pub fn to_decimal_string(&self) -> String { self.0.to_string() }
    }
}

pub use imp::Big;
```

  Remove the `is_divisible_u` part of the debug assertion if it doesn't type-check for large `d`; the `%` check alone is enough.
- [ ] **Step 4: Port `nthdigit2.rs`** to `Big`, mechanically: `Integer::from(&a * &b)` → `a.mul(&b)`; `x % q` → `x.rem(q)`; `Integer::from(&x % &q_big).to_u64().unwrap()` → `x.rem_u64(q)`; `Integer::from(k).binomial(i) * Integer::from(p_pow)` then `% modulus` → `Big::binomial(k, i).mul_u64(p_pow).rem_u64(modulus)`. Tests inside `#[cfg(test)] mod tests` may keep `rug` because `rug` stays a dev-dependency (Task 4 makes that explicit).
- [ ] **Step 5: Make `error_units` bignum-free** in `nthdigit.rs`. It must equal ⌈2^128 / 10^n0⌉ + terms (saturating). For `n0 ≥ 1`, `10^n0` never divides `2^128`, so ⌈2^128/d⌉ = ⌊(2^128 − 1)/d⌋ + 1:

```rust
fn error_units(n0: u32, terms: u64) -> u128 {
    assert!((1..=38).contains(&n0), "n0 out of the u128-certifiable range");
    let d = 10u128.pow(n0);
    (u128::MAX / d + 1).saturating_add(terms as u128)
}
```

  The existing unit test `error_bound_includes_rounding_per_term` must still pass. Add one assertion to it comparing `error_units(20, 0)` against the old rug computation, done inline in the test with `rug`.
- [ ] **Step 6: Run** `cargo test --release -j 4 --workspace` (same totals as Task 1 plus the 2 new bignum tests; 0 failed) and clippy (zero warnings). Diff the four `pihunt digit` outputs: identical.
- [ ] **Step 7: Commit.** Message: `refactor(pi-digits): route all bignum arithmetic through a Big seam; Theorem 1 bignum-free`.

---

### Task 3: `pi_ref`, an exact small-position fallback on `Big`

**Files:**
- Create: `crates/pi-digits/src/pi_ref.rs`
- Modify: `crates/pi-digits/src/nthdigit.rs`. `digits_via_mpfr` becomes a thin wrapper that calls `pi_ref::digits`; rename it `digits_fallback` and update `nthdigit2.rs`'s call sites.
- Test: `crates/pi-digits/tests/pi_ref.rs`

**Interfaces:**
- Produces: `pi_digits::pi_ref::digits(n: u64, count: usize) -> String`, the digits at positions `n+1..=n+count`. It is exact: guard digits plus a boundary check, recomputing with more guard digits if needed.

- [ ] **Step 1: Failing test**, comparing against MPFR (rug is a dev-dependency):

```rust
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
```

- [ ] **Step 2: Run it.** Expected: compile error, `pi_ref` missing.
- [ ] **Step 3: Implement it** with Machin, π = 16·arctan(1/5) − 4·arctan(1/239), in fixed point. Let `P = n + count + guard` decimal digits and `S = 10^P`. Compute `arctan(1/x)·S` as Σ (−1)^k · S / ((2k+1)·x^(2k+1)) using integer division via `div_u64_exact`-free truncating division. Add `pub fn div_u64(&self, d: u64) -> Big` (truncating) to the seam in both backends, and add a test line for it in `tests/bignum.rs`. Keep the running term `t = S/x^(2k+1)` and divide by `x²` each step. The truncation error is ≤ 1 unit per term, so the total error is ≤ number of terms (≪ 10^guard). The digits are certified if adding or subtracting the error bound doesn't change digits `n+1..=n+count`; otherwise, retry with guard + 20. Take digits from `to_decimal_string()` of `floor(π·S)`, left-padded to account for the leading `3`.
- [ ] **Step 4: Swap the fallback.** Point `nthdigit::digits` and `nthdigit2::digits` at `pi_ref::digits`, delete the MPFR fallback body, and run `cargo test --release -j 4 --workspace`. All pass, including the existing small-n MPFR comparison tests.
- [ ] **Step 5: Commit.** Message: `feat(pi-digits): exact integer-Machin fallback for small positions (replaces runtime MPFR)`.

---

### Task 4: The `pure` backend (dashu-int) and feature-gated rayon

**Files:**
- Modify: `crates/pi-digits/Cargo.toml`. `rug` becomes `optional = true`, enabled by feature `gmp = ["dep:rug"]`. Add `dashu-int = { version = "<latest 0.4.x>", optional = true }` with `pure = ["dep:dashu-int"]`. `rayon` becomes optional, `parallel = ["dep:rayon"]`. Add `[dev-dependencies] rug = { same version/features }`, so MPFR references work in tests under both backends.
- Modify: `crates/pi-digits/src/bignum.rs` (add the `#[cfg(feature = "pure")] mod imp` with the same method list)
- Modify: `crates/pi-digits/src/nthdigit.rs`, `nthdigit2.rs`. Wrap every rayon use in a small helper so the sequential path is identical:

```rust
/// `par_iter` when the `parallel` feature is on, a plain `iter` otherwise (WASM has no threads).
#[cfg(feature = "parallel")]
macro_rules! maybe_par_iter { ($e:expr) => { rayon::iter::IntoParallelRefIterator::par_iter($e) } }
#[cfg(not(feature = "parallel"))]
macro_rules! maybe_par_iter { ($e:expr) => { ($e).iter() } }
```

  Add equivalents as needed for `into_par_iter` and `rayon::join` (sequential join: `(a(), b())`). The `fold`/`reduce` pairs used with rayon need a sequential form: for plain iterators, rewrite `.fold(init, f).reduce(init2, g)` as `.fold(init(), f)`, keeping a single code path per call site behind `#[cfg]`.
- Test: the existing suites, run under both feature sets.

**Interfaces:**
- Consumes: the `Big` method list from Task 2 plus `div_u64` from Task 3.
- Produces: `pi-digits` builds with `--no-default-features --features pure` and has no rug or rayon in its non-dev dependency graph (checked in Step 5).

- [ ] **Step 1: Failing build.** Run `cargo test --release -j 4 -p pi-digits --no-default-features --features pure`. Expected: compile errors (no pure `imp`, rayon unconditional).
- [ ] **Step 2: Implement the dashu `imp`.** Every method maps to a `dashu_int::UBig` operation (`UBig` because all values are non-negative). Use `UBig::from(v)`, `&a * &b`, `&a % &b`, `a.bit_len()`, `UBig::from(base).pow(exp as usize)`, `to_string()` for decimal, and `u64::try_from(&x).ok()`. Build `binomial` as a product with exact division: `acc = acc * (n - i) / (i + 1)` for `i` in `0..k`.
- [ ] **Step 3: Gate rayon** with the helpers above, one call site at a time. Keep `cargo test --release -j 4 -p pi-digits` (default features) green after each file.
- [ ] **Step 4: Run all suites under pure.** `cargo test --release -j 4 -p pi-digits --no-default-features --features pure`: all non-ignored tests pass, including `tests/bignum.rs`, which now exercises dashu. Then, individually and under the cap, run each ignored large-position test with `--no-default-features --features pure -- --ignored --exact <name>`. They must pass (same digits).
- [ ] **Step 5: Dependency check.** `cargo tree -p pi-digits --no-default-features --features pure -e normal | grep -E "rug|rayon|gmp-mpfr"` must print nothing.
- [ ] **Step 6: Record the pure-vs-gmp speed.** One back-to-back run each at position 10⁵ for Theorem 2 with the default mem, single-threaded (`RAYON_NUM_THREADS=1` for gmp; pure has no rayon). Put the times in the commit message.
- [ ] **Step 7: Commit.** Message: `feat(pi-digits): pure-Rust dashu backend and optional rayon; suites green under both backends`.

---

### Task 5: The `site/wasm` crate and a Node smoke test

**Files:**
- Modify: root `Cargo.toml`. Add `"site/wasm"` to workspace members.
- Create: `site/wasm/Cargo.toml`:

```toml
[package]
name = "pi-digits-wasm"
version = "0.1.0"
edition = "2024"
publish = false

[lib]
crate-type = ["cdylib"]

[dependencies]
pi-digits = { path = "../../crates/pi-digits", default-features = false, features = ["pure"] }
wasm-bindgen = "0.2"

[profile.release]
# profile settings must live in the workspace root; see Step 2
```

  Cargo only honours `[profile.*]` at the workspace root, so put `[profile.release] lto = true, codegen-units = 1, opt-level = 3` in the root `Cargo.toml` if it's not there already. Check that this doesn't change CLI behaviour; it can only make it faster.
- Create: `site/wasm/src/lib.rs`
- Create: `site/scripts/smoke-wasm.mjs`
- Create: `site/scripts/gen-pi.sh`

**Interfaces:**
- Produces (JS, via wasm-bindgen `--target web` in `site/dist/wasm/`): `digits_thm1(pos: number, count: number): string`, `digits_thm2(pos: number, count: number): string`, `wasm_memory_bytes(): number`. They throw a JS `Error` for `pos < 1`, a non-integer `pos`, or `pos > 200000`.
- Produces: `site/data/pi-200k.txt`, exactly 200 000 ASCII digits of positions 1..=200 000 (no "3.", no newline).

- [ ] **Step 1: Install the tooling once.** `rustup target add wasm32-unknown-unknown`; `cargo install wasm-pack --locked`; install binaryen's `wasm-opt` with `cargo install wasm-opt --locked` if `which wasm-opt` finds nothing.
- [ ] **Step 2: Write the wrapper**, `site/wasm/src/lib.rs`:

```rust
use wasm_bindgen::prelude::*;

const MAX_POS: f64 = 200_000.0;

fn check_pos(pos: f64) -> Result<u64, JsError> {
    if !(pos.is_finite() && pos.fract() == 0.0 && pos >= 1.0 && pos <= MAX_POS) {
        return Err(JsError::new("position must be an integer in 1..=200000"));
    }
    Ok(pos as u64)
}

/// Digits at 1-based positions `pos..pos+count` via Gourdon's Theorem 1.
#[wasm_bindgen]
pub fn digits_thm1(pos: f64, count: u32) -> Result<String, JsError> {
    let p = check_pos(pos)?;
    Ok(pi_digits::nthdigit::digits(p - 1, count as usize))
}

/// Same, via the reconstructed Theorem 2 at its default memory budget.
#[wasm_bindgen]
pub fn digits_thm2(pos: f64, count: u32) -> Result<String, JsError> {
    let p = check_pos(pos)?;
    let n = p - 1;
    Ok(pi_digits::nthdigit2::digits(n, count as usize, pi_digits::nthdigit2::default_mem_bits(n.max(1))))
}

/// Current size of this module's linear memory, in bytes (it only ever grows).
#[wasm_bindgen]
pub fn wasm_memory_bytes() -> f64 {
    (core::arch::wasm32::memory_size(0) * 65536) as f64
}
```

- [ ] **Step 3: Generate the reference digits.** `site/scripts/gen-pi.sh` builds the release CLI, then runs `pihunt digit 1 --count 200000 --method thm1 > site/data/pi-200k.txt`. That uses the certified chunked path, which is exact. Don't commit the output if it's > 200 KB unexpectedly; it should be exactly 200 000 bytes. Verify it with `head -c 5` → `14159`, and check the Feynman point by printing bytes 762–767, expecting `999999`. Check this in to git; it's small and deterministic.
- [ ] **Step 4: Smoke test**, `site/scripts/smoke-wasm.mjs`. Node built-ins only; it loads the `--target web` output:

```js
import { readFile } from 'node:fs/promises';
import init, { digits_thm1, digits_thm2 } from '../dist/wasm/pi_digits_wasm.js';

const wasm = await readFile(new URL('../dist/wasm/pi_digits_wasm_bg.wasm', import.meta.url));
await init({ module_or_path: wasm });
const ref = (await readFile(new URL('../data/pi-200k.txt', import.meta.url), 'utf8')).trim();
let failed = 0;
for (const pos of [1, 5, 762, 1999, 2000, 2001, 10000, 50000]) {
  const want = ref.slice(pos - 1, pos - 1 + 10);
  for (const [name, f] of [['thm1', digits_thm1], ['thm2', digits_thm2]]) {
    const got = f(pos, 10);
    if (got !== want) { failed++; console.error(`FAIL ${name} pos=${pos}: got ${got} want ${want}`); }
  }
}
for (const bad of [0, -1, 1.5, 200001, NaN]) {
  try { digits_thm1(bad, 10); failed++; console.error(`FAIL: pos=${bad} accepted`); } catch { /* expected */ }
}
if (failed) { console.error(`${failed} failure(s)`); process.exit(1); }
console.log('wasm smoke test: all positions match pi-200k.txt; bad positions rejected');
```

  Positions 1999/2000/2001 straddle the fallback boundary (Review Focus 2).
- [ ] **Step 5: Build and run.** `wasm-pack build site/wasm --target web --release --out-dir ../dist/wasm --out-name pi_digits_wasm`, then `wasm-opt -O3` in place on the `.wasm` (skip it if wasm-pack already ran wasm-opt), then `node site/scripts/smoke-wasm.mjs`. Expected: `wasm smoke test: all positions match …`. Record the `.wasm` size.
- [ ] **Step 6: Gitignore** `site/dist/`, then commit `site/wasm`, `site/scripts`, `site/data/pi-200k.txt` and `Cargo.toml`. Message: `feat(site): wasm-bindgen wrapper over pi-digits (pure backend) with Node smoke test`.

---

### Task 6: Demo front end (race and stream) with Web Workers

**Files:**
- Create: `site/src/worker.js`, `site/src/demo.js`, `site/src/demo.test.mjs`
- Modify: `site/src/index.html` (demo markup; created minimally here and filled with content in Task 7)

**Interfaces:**
- Consumes: the WASM exports from Task 5, at `./wasm/pi_digits_wasm.js` relative to the page.
- Worker protocol: page → worker `{run, kind: 'race'|'stream', method: 'thm1'|'thm2', pos, count}`. Worker → page:
  - `{run, type: 'ready'}`
  - `{run, type: 'result', method, pos, digits, ms, memBytes}`
  - `{run, type: 'block', pos, digits, ms, memBytes}` (stream)
  - `{run, type: 'error', message}`
- Produces: from `demo.js`, pure functions (unit-testable in Node) `parsePosition(str) -> {ok: true, pos} | {ok: false, message}`, `checkAgainstReference(ref, pos, digits) -> 'match'|'mismatch'|'unchecked'`, and `lanesAgree(a, b) -> boolean`, plus the DOM wiring.

- [ ] **Step 1: Failing unit tests**, `site/src/demo.test.mjs` (`node --test`):

```js
import test from 'node:test';
import assert from 'node:assert/strict';
import { parsePosition, checkAgainstReference, lanesAgree } from './demo.js';

test('parsePosition accepts plain integers in range', () => {
  assert.deepEqual(parsePosition('1'), { ok: true, pos: 1 });
  assert.deepEqual(parsePosition(' 100000 '), { ok: true, pos: 100000 });
  assert.deepEqual(parsePosition('200000'), { ok: true, pos: 200000 });
  assert.deepEqual(parsePosition('10,000'), { ok: true, pos: 10000 });
});

test('parsePosition rejects everything else with a message', () => {
  for (const bad of ['', '0', '-5', '1.5', '1e5', 'abc', '200001', '99999999999']) {
    const r = parsePosition(bad);
    assert.equal(r.ok, false, `accepted ${JSON.stringify(bad)}`);
    assert.match(r.message, /\S/);
  }
});

test('checkAgainstReference', () => {
  const ref = '1415926535';
  assert.equal(checkAgainstReference(ref, 1, '14159'), 'match');
  assert.equal(checkAgainstReference(ref, 1, '14158'), 'mismatch');
  assert.equal(checkAgainstReference(ref, 8, '53599'), 'unchecked'); // runs past the reference
  assert.equal(checkAgainstReference(null, 1, '14159'), 'unchecked');
});

test('lanesAgree', () => {
  assert.equal(lanesAgree('0123456789', '0123456789'), true);
  assert.equal(lanesAgree('0123456789', '0123456780'), false);
});
```

  `1e5` is rejected on purpose: only plain digits with optional thousands commas are accepted, so the input stays unambiguous.
- [ ] **Step 2: Run** `node --test site/src/`. Expected: FAIL (module missing).
- [ ] **Step 3: Implement the pure functions** at the top of `demo.js`, exported. `parsePosition` strips spaces and commas, requires `/^\d+$/`, and checks 1..=200000. The rejection messages are: "Enter a whole number from 1 to 200,000." (format) and "Positions above 200,000 are beyond the browser demo; use the CLI." (range).
- [ ] **Step 4: Worker**, `site/src/worker.js`, a module worker:

```js
import init, { digits_thm1, digits_thm2, wasm_memory_bytes } from './wasm/pi_digits_wasm.js';

const ready = init();
self.onmessage = async ({ data }) => {
  const { run, kind, method, pos, count } = data;
  try {
    await ready;
    const f = method === 'thm2' ? digits_thm2 : digits_thm1;
    if (kind === 'race') {
      const t0 = performance.now();
      const digits = f(pos, count);
      self.postMessage({ run, type: 'result', method, pos, digits, ms: performance.now() - t0, memBytes: wasm_memory_bytes() });
    } else {
      // Stream: each block is computed from scratch; nothing is carried between blocks.
      for (let p = pos; p + count - 1 <= 200000; p += count) {
        const t0 = performance.now();
        const digits = digits_thm1(p, count);
        self.postMessage({ run, type: 'block', pos: p, digits, ms: performance.now() - t0, memBytes: wasm_memory_bytes() });
      }
    }
  } catch (e) {
    self.postMessage({ run, type: 'error', message: String(e?.message ?? e) });
  }
};
```

- [ ] **Step 5: DOM wiring** in `demo.js`:
  - One `let currentRun = 0`. Every Race, Cancel or Stream action does `currentRun++` and terminates live workers. Message handlers drop any message whose `run !== currentRun` (Review Focus 3).
  - Race creates two workers (`new Worker('./worker.js', { type: 'module' })`) and posts both messages in the same tick.
  - Each lane has a live timer (`requestAnimationFrame`), then shows digits, ms and memory in KiB.
  - When both are done, show "✓ lanes agree" or a prominent red "✗ Theorem 1 and Theorem 2 disagree — please report this" (Review Focus 4). Also show the reference check (✓ / ✗ / not checked).
  - Load `./data/pi-200k.txt` lazily on first use.
  - Above 100 000, show the soft warning "This may take a minute or more in the browser." but still allow the run.
  - Test hook: if `location.hash === '#test-mismatch'`, corrupt lane 2's digits before comparing, so the disagreement UI can be checked by hand.
  - Stream panel: start position input (same `parsePosition`), Start/Stop, append blocks in monospace groups of 10, and show the memory readout after each block.
- [ ] **Step 6: Run** `node --test site/src/` (PASS). Then serve locally: `node site/scripts/serve.mjs site/dist`, a 20-line static server built on `node:http` with correct `application/wasm` and `text/javascript` MIME types; create it in this task. Load the page, race at position 1000 and at 10 000, and cancel mid-race and re-race. Confirm no stale output, that the `#test-mismatch` hook shows the red banner, and that stream shows flat memory for ≥ 20 blocks. Use the Chrome browser tool if available; otherwise report which checks were done and how.
- [ ] **Step 7: Commit.** Message: `feat(site): race and stream demo with run-id guarded workers`.

---

### Task 7: Write-up content, charts, typography

**Files:**
- Modify: `site/src/index.html` (full content)
- Create: `site/src/style.css`, `site/data/benchmarks.json`, `site/scripts/charts.mjs`, `site/vendor/` (KaTeX dist and fonts)

**Interfaces:**
- Consumes: `docs/nthdigit-theorem2.md`, `docs/nthdigit.md`, `docs/findings-bbp-hunt.md` as the sources of truth.
- Produces: `site/data/benchmarks.json`:

```json
{
  "machine": "AMD Ryzen 7 5700G, 6 cores",
  "date": "2026-09-23",
  "rows": [
    { "pos": 10000,    "thm1_s": 0.034, "thm1_mib": 5.2, "thm2_s": 0.016, "thm2_mib": 5.7,  "pidec_s": 3.13 },
    { "pos": 100000,   "thm1_s": 1.51,  "thm1_mib": 5.1, "thm2_s": 0.217, "thm2_mib": 7.3,  "pidec_s": 185.1 },
    { "pos": 1000000,  "thm1_s": 113.4, "thm1_mib": 5.0, "thm2_s": 4.17,  "thm2_mib": 14.1, "pidec_s": 15869 },
    { "pos": 3000000,  "thm1_s": null,  "thm1_mib": null, "thm2_s": 19.3, "thm2_mib": 22.7, "pidec_s": null },
    { "pos": 10000000, "thm1_s": 7107,  "thm1_mib": 4.2, "thm2_s": 111.0, "thm2_mib": 69.2, "pidec_s": null }
  ],
  "memory_history_1e7_mib": { "first_port": 812.1, "streamed": 447, "memo_fix": 69.2 }
}
```

  Every number is copied from the tables in `docs/nthdigit.md` (the Theorem 2 headline table, the Theorem 1 table, and Gourdon's pidec column). Re-check each one against the doc while writing; if the doc disagrees with this block, the doc wins. Fix the JSON and note it in the commit message. The pidec times are Gourdon's own, on a Pentium III 900 MHz, and must be labelled as such everywhere they appear.

- [ ] **Step 1: Vendor the assets.**
  - Download the KaTeX release (`katex.min.css`, `katex.min.js`, `contrib/auto-render.min.js`, fonts/) from the official GitHub release tarball into `site/vendor/katex/`.
  - Download Inter (variable, `InterVariable.woff2` from the official rsms/inter release) and JetBrains Mono (`JetBrainsMono[wght].woff2` from its release) into `site/vendor/fonts/`, along with their OFL license files.
  - Record the versions in `site/vendor/VERSIONS.txt`.
- [ ] **Step 2: Chart generator**, `site/scripts/charts.mjs`. Node built-ins only. It reads `benchmarks.json`, writes two SVG strings, and replaces the `<!-- CHART:time -->` and `<!-- CHART:memory -->` markers in the HTML it's given, writing the result to `site/dist/index.html`.
  - Log-log axes with decade gridlines.
  - Series: Theorem 1, Theorem 2, and pidec 2003 (dashed, labelled "(Pentium III, 2003)"). Null points are skipped, not zeroed.
  - Direct labels at the line ends, no legend box.
  - Colours via CSS custom properties (`var(--series-1)` etc.) so both themes work.
  - `role="img"` with an `aria-label` summarising the chart.
  - Load the `dataviz` skill before writing this step and follow its colour and mark guidance.
- [ ] **Step 3: Content.** Write `index.html` following spec §3, sections 1–13 in order.
  - Maths in KaTeX delimiters (`\( … \)`, `\[ … \]`), rendered by auto-render on load.
  - Derivations in `<details>` blocks.
  - Every claim tagged with a small `proved` / `measured` / `conjecture` badge.
  - The "errors we caught" box lists: the missing fixed-point rounding term in the certification bound (found in review, fixed); the doc position labels mixing CLI and library conventions (both 10⁷ strings MPFR-verified); a verification script that used `int()` on gmpy2 floats, which rounds instead of truncating.
  - The literature table is copied from `docs/nthdigit-theorem2.md` §1, with its URLs.
  - Reproduce section commands: `cargo build --release`, `pihunt digit 1000000` (Theorem 2 default), `pihunt digit 1000000 --method thm1`, `cargo test --release -j 4`.
  - A `<noscript>` block in the demo area explains that the demo needs JavaScript and WebAssembly and links to the benchmark table (Review Focus 5).
  - Load the `frontend-design` skill before writing `style.css`. Direction: clean modern sans (Inter), generous whitespace, a narrow reading column (~70ch) with the demo allowed wider, monospace digits, restrained accent colour, and light/dark via `prefers-color-scheme` plus a toggle. No horizontal scroll at 360 px.
- [ ] **Step 4: Fact check.** For every number in `index.html`, `grep` it in `docs/` or `benchmarks.json`. List any number that can't be traced and fix it or remove it. Include the result, "N numbers checked, all traced", in the commit message.
- [ ] **Step 5: Commit.** Message: `feat(site): write-up content, generated charts, vendored KaTeX and fonts`.

---

### Task 8: `build.sh`, `_headers`, site checks, final verification

**Files:**
- Create: `site/build.sh`, `site/src/_headers`, `site/scripts/check-site.mjs`, `site/README.md`

**Interfaces:**
- Consumes: everything above.
- Produces: `site/dist/`, a complete deployable site. `site/README.md` gives the one-line build and the one-line deploy.

- [ ] **Step 1: `check-site.mjs`** (Node built-ins). Over `site/dist/index.html`:
  - Every `href="#x"` has a matching `id="x"`.
  - Every local `src`/`href` file exists in `dist/`.
  - There are no `http(s)://` references in `<script src>`, `<link href>` or CSS `url()`. External links in `<a href>` are allowed.
  - `_headers` exists and contains `application/wasm` and `Content-Security-Policy`.
  - A `<noscript>` element exists inside the demo section.
  - It exits 1 listing every failure.

  Test it by running it against a deliberately broken copy (delete one anchor target and add one CDN `<script>`): it must report both.
- [ ] **Step 2: `_headers`:**

```
/*
  Content-Security-Policy: default-src 'self'; script-src 'self' 'wasm-unsafe-eval'; style-src 'self' 'unsafe-inline'; img-src 'self' data:; font-src 'self'; worker-src 'self'; connect-src 'self'
  X-Content-Type-Options: nosniff
  Referrer-Policy: strict-origin-when-cross-origin
/wasm/*.wasm
  Content-Type: application/wasm
  Cache-Control: public, max-age=31536000, immutable
/vendor/*
  Cache-Control: public, max-age=31536000, immutable
```

  `'unsafe-inline'` for styles is only there for KaTeX's inline styles. Leave a comment in the README saying so.
- [ ] **Step 3: `build.sh`** (`set -euo pipefail`):
  1. Tooling check with clear install hints.
  2. Clean `site/dist`.
  3. The Task 5 wasm-pack and wasm-opt build into `dist/wasm`.
  4. `gen-pi.sh` if `data/pi-200k.txt` is missing.
  5. `node scripts/smoke-wasm.mjs`.
  6. `node --test src/`.
  7. Copy `src/*` (except tests), `vendor/`, `data/` into `dist/`.
  8. `node scripts/charts.mjs src/index.html > dist/index.html`.
  9. `node scripts/check-site.mjs`.
  10. Print the dist size (`du -sh`) and the wasm size.

  Stop on the first failure.
- [ ] **Step 4: Full verification:**
  - `site/build.sh` from clean: all green.
  - The workspace default suite (`cargo test --release -j 4 --workspace`) passes, plus the pure-backend suite (`cargo test --release -j 4 -p pi-digits --no-default-features --features pure`). Clippy is clean under both.
  - Serve `dist/` locally and do the manual browser pass. Both themes; 360 px width with no horizontal scroll; race at 1 000 and 50 000 with times recorded; cancel/re-race; the `#test-mismatch` banner; stream for 20 blocks; KaTeX renders; charts render in both themes; JS disabled shows the paper and the `<noscript>` message.
  - Report exactly what was observed, with the browser-measured race times.
  - The dist size is under ~1.5 MB (`du -sh site/dist`).
- [ ] **Step 5: README and commit.** `site/README.md` covers the build (`site/build.sh`), deploy (`npx wrangler pages deploy site/dist --project-name <your-project>`), and a note that Pages' git-build image has no Rust toolchain, so upload the prebuilt `dist/`. Message: `feat(site): build script, Pages headers, site checks`.
