# pihunt — design spec

**Date:** 2026-09-22
**Status:** approved design, pending implementation plan

## Goal

Search for new BBP-type digit-extraction formulas for π, especially in base 10 (and powers of 10), using integer relation detection (PSLQ). Built in Rust.

Two stages:

- **Stage 1 — explorer:** correct, well-tested classic PSLQ; basis/reduction engine; TOML-driven batches; JSONL results log. Must rediscover known formulas.
- **Stage 2 — serious search:** resume, fast multi-level PSLQ, precision escalation, exhaustive grids, exclusion-bound reports.

Expected outcome is most likely a null result. The deliverable in that case is a set of precise exclusion claims: "no relation of shape X with max coefficient ≤ C exists (in the reduced basis)".

## Background

A BBP-type formula has the shape

    π = Σ_{k≥0} (1/b^k) · Σ_{s} Σ_{j=1..m} a_{j,s} / (mk + j)^s

with base `b`, period `m`, degrees `s`, integer coefficients `a_{j,s}`. BBP itself is `b=16, m=8, s=1`, coefficients `(4, 0, 0, -2, -1, -1, 0, 0)`, scaled by 1.

Borwein, Borwein & Galway (2004) exclude degree-1 Machin-type BBP formulas for π in non-binary bases. This search targets the ground they did not cover: higher degrees, mixed degrees, larger periods, bases 10^k, and relations that also involve other constants.

## Architecture

Single binary crate `pihunt` at the repo root.

```
src/
├── main.rs         # CLI: run / verify / plan
├── config.rs       # TOML → Batch
├── plan.rs         # Batch → Vec<Job> (grid = cartesian, sample = Sobol/LHS)
├── basis.rs        # Job → Vec<Float> (π, series columns, extras) at N bits
├── reduce.rs       # drop linearly dependent basis columns before the π search
├── pslq/
│   ├── mod.rs      # trait RelationFinder
│   └── classic.rs  # textbook PSLQ, full MPFR precision
├── classify.rs     # Outcome → Hit | Junk | Suspicious | Spurious | Excluded | Inconclusive
├── verify.rs       # residual check against 2× precision columns
├── known.rs        # literature formulas, for known/NEW tagging
├── job.rs          # one job end to end → Record
└── log.rs          # append-only JSONL writer
```

Data flow:

    batch.toml → plan → Vec<Job> → (rayon) → basis → reduce → pslq → classify → verify (Hits only) → log

Boundaries:

- `RelationFinder` is a trait so multi-level PSLQ can be added in stage 2 and diffed against `classic`.
- `basis` knows nothing about PSLQ; it is tested independently against known constant values.
- Jobs are pure and independent: same input → same output.
- Every `Outcome` carries the precision used and the norm bound reached.

## Basis and reduction

### Columns

For a job `(b, m, s ∈ [s_lo, s_hi], extras)`:

- Column 0 is always **π**.
- **Series columns** `S(j,s) = Σ_{k≥0} 1 / (b^k (mk+j)^s)` for `j = 1..m`, each `s` in range.
- **Extras** from a fixed menu: `pi2` (π²), `log2`, `log3`, `log5`, `catalan`, `zeta3`.

Column count `n = 1 + m·(s_hi − s_lo + 1) + |extras|`.

### Precision

- Config gives coefficient bound `C`.
- Auto precision (decimal digits): `ceil(n · log10(C) · 1.5) + 50`. Overridable with an integer ≥ 60. (Prototyping showed factor 1.25 lets spurious ~10⁷-coefficient relations through at n ≈ 45; 1.5 excludes cleanly.)
- Every job builds its columns once at P digits (for PSLQ) and once at 2P digits (for verification).
- Series evaluated with guard bits `ceil(log2(terms)) + 32` above working precision, then rounded.
- Terms per series ≈ `precision_bits / log2(b)` plus margin for the `(mk+j)^s` factor.

### Reduction

Removes rational linear dependencies among non-π columns so every relation found in the main search must involve π.

1. Run PSLQ on all columns except π.
2. If a relation is found, it must have max |coeff| ≤ C **and** verify at 2P precision; otherwise the job ends `Inconclusive` (never drop a column on an unverified relation — prototyping showed unchecked reduction at n ≈ 60 dropping genuine columns on 10⁶–10⁸-coefficient garbage). A trusted relation is recorded as a basis relation (column dropped + coefficients) and the highest-index column with nonzero coefficient is dropped. Basis relations are tagged `known` if they match the known-formula table, otherwise untagged (base-10 log identities are common and uninteresting, so they don't scream NEW).
3. Repeat until PSLQ reports no relation (Excluded).
4. Prepend π and run the main search on the reduced basis.

Dropping a dependent column preserves the Q-span, so any π relation over the original basis still exists over the reduced basis. Exclusion bounds are stated for the **reduced** basis; the log records exactly which columns were dropped and why.

If a reduction step ends inconclusive (precision exhausted / iteration cap), the job ends `Inconclusive` without running the main search.

## PSLQ

Classic Ferguson–Bailey PSLQ in `rug::Float` at job precision.

- Normalise `x`; build `H` (n × n−1, lower trapezoidal), `A = B = I`; initial Hermite reduction.
- Iteration: choose `r` maximising `γ^r · |H_rr|`; swap rows `r, r+1` of `y, H, A, B`; corner fix-up when `r < n−1`; Hermite reduce.
- `γ` default `1.16`, configurable; must be `> sqrt(4/3)` (validated at config load).
- After each iteration, the bound `1 / max_j |H_jj|` is a lower bound on the norm of any integer relation. It is monotone non-decreasing.

### Termination

| Outcome | Trigger | Meaning |
|---|---|---|
| `Relation` | `min |y_i| < 10^-(digits − 30)` | candidate relation = column `i` of `B` |
| `Excluded` | bound `> C · sqrt(n)` | no relation with max \|coefficient\| ≤ C exists |
| `PrecisionExhausted` | max \|A_ij\| exceeds `10^(digits − 30)` | inconclusive; report bound reached |
| `IterationCap` | `max_iterations` reached | inconclusive; report bound reached |

### Classification

For a `Relation`:

1. Normalise: divide by gcd, sign so `a_π > 0`.
2. `a_π = 0` → `Junk` (should be impossible after reduction; logged as a bug signal).
3. `max |a| > C` → `Suspicious`.
4. Otherwise verify: rebuild the basis at 2× precision; residual `|a · x|` must be `< 10^-(2·digits − 60)`. Pass → `Hit`, fail → `Spurious`.
5. `Hit`s are compared against a table of known formulas (BBP base-2^k family, Bailey's base-64 π² formula, known log 2 formulas) and tagged `known` or `NEW`.

`Excluded` is the only null outcome that constitutes an exclusion claim. `PrecisionExhausted` and `IterationCap` are logged as `inconclusive`.

## Batch config

```toml
name    = "base10-scout"
output  = "results/base10-scout.jsonl"
threads = 6                      # default: all cores

[defaults]
coeff_bound      = 10000
precision_digits = "auto"        # or integer
gamma            = 1.16
max_iterations   = 1_000_000
max_columns      = 80            # larger jobs are logged as skipped

[search]
mode    = "grid"                 # "grid" | "sample"
bases   = [10, 100, 1000]
periods = { from = 2, to = 24 }  # inclusive
degrees = [[1, 1], [1, 2], [2, 2]]
extras  = [[], ["log2", "log5"], ["pi2", "log2", "log3", "log5", "catalan", "zeta3"]]

[sample]                         # only used when mode = "sample"
method = "sobol"                 # "sobol" | "lhs"
count  = 200
seed   = 42
```

- **grid:** cartesian product of `bases × periods × degrees × extras`.
- **sample:** one unit-cube dimension per axis (base, period, degree range, extras set); each coordinate floored onto the discrete choices; duplicates dropped; seeded and reproducible. Sample mode is for scouting only — it makes no coverage claim. `count` ≤ 65536 and `seed` is a u32 (Sobol generator limits).

Config validation errors (bad γ, unknown extra, empty axis, `from > to`) fail before any job runs.

## Job identity

`job_id` = blake3 hash (hex, truncated to 16 bytes) of the canonical job parameters: base, period, degree range, **sorted** extras, coeff_bound, resolved precision, γ, finder name, `ALGO_VERSION`.

- Same job → same ID across batches (enables resume/dedupe in stage 2).
- Bump `ALGO_VERSION` whenever math affecting outcomes changes.

## Results log

Append-only JSONL, one line per finished job, written by a single writer thread fed by a channel, flushed per line.

```json
{"job_id":"a3f9c2…","batch":"base10-scout","pihunt_version":"0.1.0",
 "started":"2026-09-22T19:40:00Z","elapsed_ms":84211,
 "params":{"base":10,"period":12,"degrees":[1,2],"extras":["log2","log5"],
           "coeff_bound":10000,"precision_digits":412,"gamma":1.16,
           "finder":"classic","algo_version":1},
 "columns":["pi","S(j=1,s=1)","…","log2","log5"],
 "dropped":[{"column":"S(j=12,s=1)","relation":["1","-3","…"],"tag":"known"}],
 "outcome":"excluded","bound":"1.2e5","iterations":48123,
 "relation":null,"verify":null,"tag":null}
```

- `outcome` ∈ `hit | junk | suspicious | spurious | excluded | inconclusive | skipped`.
- `params` also records `max_iterations`; `verify` is `{"passed": bool, "residual_log10": f64}`; a `note` field says why a job was inconclusive, skipped or junk.
- Integer coefficients and bounds serialised as **strings** (arbitrary size, no silent truncation).
- End of batch: terminal summary with counts per outcome and every `Hit` listed; `NEW` hits highlighted.

## CLI

- `pihunt run <batch.toml>` — execute a batch.
- `pihunt verify <results.jsonl>` — re-verify every `Hit` at 2× precision.
- `pihunt plan <batch.toml>` — dry run: job count, largest `n`, precision range.

## Testing

Unit / integration tests, easiest first:

- **Basis:** `b=16, m=8, s=1` columns match reference values; extras match MPFR constants; precision derivation monotone in `n` and `C`.
- **Reduction:** basis with a planted duplicate (`log2` plus a series equal to `log 2`) drops exactly that column.
- **PSLQ toys:** `[√2, √8] → (2, −1)`; `[log2, log3, log6] → (1, 1, −1)`.
- **Planted relations:** random reals plus one column equal to a small integer combination of the others; recovered exactly across many seeds.
- **No false relations:** random reals over 1000 seeds → `Excluded` or inconclusive, never `Relation`.
- **Monotone bound:** exclusion bound never decreases across iterations.
- **BBP rediscovery:** `b=16, m=8, s=1` → BBP, tagged `known`.
- **π² rediscovery:** `b=64, m=6, s=2` with `pi2` in extras → reduction records Bailey's π² formula as a `basis_relation` dropping `pi2`, tagged `known`. (It has no π term, so it can never be a Hit.)
- **Config:** invalid configs rejected; grid expansion count correct; sample mode deterministic for a fixed seed.
- **Job ID:** stable across runs; insensitive to extras ordering; sensitive to every other parameter.

## Stage 1 done criteria

1. All tests pass, including 1000-seed no-false-relation test.
2. BBP rediscovered end-to-end from a batch TOML, tagged `known`.
3. Bailey's base-64 π² formula rediscovered as a `basis_relation`, tagged `known`.
4. A ~200-job base-10 sample-mode scout batch completes unattended with a valid log.
5. `pihunt verify` re-confirms all Hits in that log.
6. Timing baseline recorded: seconds per job vs `n`.

## Stage 2 roadmap (out of scope for stage 1)

1. Resume: skip job IDs already present in the output log.
2. Multi-level PSLQ as a second `RelationFinder`; must match `classic` outcome class on the full stage-1 corpus before use.
3. Precision escalation: auto-rerun inconclusive jobs at 2× precision.
4. Exhaustive grids over scout-flagged regions; exclusion report generator.
5. Digit extractor for any `NEW` formula.

## Environment

- Rust via mise (`mise.toml`, `rust = "stable"`).
- `CARGO_TARGET_DIR` set by mise to `~/.cache/pihunt/target` — the repo lives under `~/sync`, which is synced.
- System GMP/MPFR (present) used by `rug`.
- Crates: `rug`, `rayon`, `serde`, `serde_json`, `toml`, `clap`, `blake3`, `sobol_burley`, `jiff`. LHS is hand-rolled.
- `results/` is committed (scientific record); `batches/` holds batch configs.

## Stage 2 — as built (2026-09-22)

Built directly from the roadmap above (no separate spec/plan, at the user's request).

- **Multilevel PSLQ** (`src/pslq/multilevel.rs`): two-level, Bailey's pslqm2. f64 inner iterations accumulate an exact integer transform (entries < 2^52; an iteration that would overflow is rolled back from a snapshot). At each sync the transform is applied to the full-precision state, H is re-triangularised with Givens LQ, then fully Hermite-reduced. Classic's internals became the shared `pslq::State`, which alone performs termination checks — the f64 loop never makes claims. The inner loop also stops once its bound estimate reaches the exclusion threshold; without that it ran past exclusions into precision-floor relations (caught by the 1000-seed random test).
- **Validation:** the whole PSLQ test suite runs against both finders, and `tests/equivalence.rs` requires identical verdicts (outcome, relation, tag, dropped columns) on every job of the known and scout batches (200 jobs). ~7× faster than classic; see `docs/timing-baseline.md`.
- **Config:** `defaults.finder = "classic" | "multilevel"` (default classic), `defaults.escalate` (default 2, max 4).
- **Resume + escalation** (`src/runner.rs`): each job runs as a chain of attempts at 1×, 2×, 4×… digits while inconclusive. Attempts whose job ID is already in the output log are reused, so re-running a batch resumes it. Records gain `escalated_from` (serde default, so stage-1 logs still load).
- **Report:** `pihunt report <logs...>` prints markdown: strongest exclusion per shape, hits (NEW ones in their own section), deduplicated basis relations, unresolved shapes.
- **Deferred:** digit extractor (no NEW hit exists); Householder LQ / three-level PSLQ (next speedup).

## Stage 3: precision rule (2026-09-22/23)

`auto_digits`'s flat factor 1.5 (measured at n ≈ 45) was too thin at n ≈ 100: jobs came back
`Inconclusive` (reduction hit a precision-floor relation) or `Suspicious` (main search did the
same), both correctly rejected but expensive to discover — the runner's fix is to escalate and
redo the whole job from scratch. `auto_digits` now grows the factor past n ≈ 36 based on direct
measurement instead of guessing further. Full measurement table, the fitted rule, and a
before/after run of the reported n = 100 example are in `docs/precision-rule.md`.
