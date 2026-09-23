# Findings: the base-10 BBP hunt (stages 1–3)

**Status: closed (2026-09-23).** No base-10 BBP-type formula for π was found. We don't expect one
in this family, for the structural reason explained below. The tooling (`pihunt`) works and is
kept. Search effort has moved on to decimal digit extraction (see the end of this doc).

## What was searched

Shapes π = Σ_k b^(−k) Σ_{j,s} a_{j,s}/(mk+j)^s, optionally with extra constants (log 2, log 3,
log 5, π², Catalan, ζ(3)), searched with PSLQ. Every candidate relation was verified at 2×
precision. Each null result carries a proved coefficient bound in the reduced basis.

| log | shapes | bases | periods | degrees | C | result |
|---|---|---|---|---|---|---|
| `results/known.jsonl` | 24 | 16, 64 | 6–8 | 1, 2 | 10³ | BBP and Bailey's π² formula rediscovered (sanity check) |
| `results/base10-scout.jsonl` | 176 | 10, 100, 1000 | 2–16 | 1..2 | 10³ | 175 excluded, 1 inconclusive |
| `results/base10-grid.jsonl` | 513 | 10, 100, 1000 | 2–20 | 1..1, 1..2, 2..2 | 10⁴ | 513 excluded |
| `results/base10-wide.jsonl` | 1404 | 10, 100, 1000 | 2–40 | 1..1, 2..2, 1..2, 1..3 | 10⁵ | 1306 excluded, 26 suspicious, 72 skipped (>100 columns) |

Reports: `docs/exclusions-base10-grid.md`, `docs/exclusions-base10-wide.md`. The 26 `suspicious`
jobs in the wide sweep predate the fix that escalates them (`8f6bdb1`). They were not re-run,
because the argument below makes that moot.

The only relations found among the non-π columns were the expected family
m·S(b,m,m,1) = b·log(b/(b−1)) scaled, e.g. 3·S(10,6,6,1) = 5·log(10/9).

## Why nothing was found (and why more search won't help)

1. **Degree 1 is already ruled out.** Borwein, Borwein & Galway (2004) exclude degree-1
   ("Machin-type") BBP formulas for π in bases that aren't powers of 2. Base 10 is covered.
2. **Higher degrees can't help, conditionally on a standard conjecture.** A degree-s column
   S(b,m,j,s) is a Q-linear combination of polylogarithm values Li_s(α) at algebraic α. These
   are periods of *weight s*. π is a period of weight 1. The standard expectation (part of the
   Grothendieck period conjecture / motivic picture, unproved) is that Q-linear relations among
   such periods respect weight. If so, in any relation π = (degree-1 part) + (degree-2 part) +
   …, the higher-degree parts must vanish on their own. π would then need a degree-1 formula,
   which item 1 rules out.
3. **So the degree-2/3 exclusions in these logs are consistent with that conjecture rather than
   discoveries waiting to happen.** Finding a mixed-degree base-10 formula for π would *refute*
   a widely believed conjecture. That's possible in principle, not something to spend CPU on.

Caveats, stated plainly: item 2 is conditional. We haven't checked whether BBG's theorem covers
*every* degree-1 BBP shape or only the Machin-type subfamily. Extras such as π² or ζ(3) don't
change the argument, since they are pure weight-2/3 periods.

## Engineering results worth keeping

- **Multilevel PSLQ** (`src/pslq/multilevel.rs`): exact f64 transform, full-precision claims
  only. It is about 7× faster than classic at n ≈ 40, and 3.2 s → 1.0 s at n = 67 after the
  stage-3 work (Householder LQ, snapshot-free inner loop, transposed B, parallel sync). It is
  equivalent to classic on the 200-job corpus. See `docs/timing-baseline.md`.
- **Precision rule** (`docs/precision-rule.md`): the digits factor grows past n ≈ 36. It was
  measured on 66 cases.
- **Soundness machinery**: 2× precision verification of every relation, reduction that never
  trusts an unverified relation, and escalation on inconclusive/suspicious/spurious outcomes.
  Each of these caught a real failure mode during development.
- **Operations**: resume, escalation, sharding (`--shard k/N`), and markdown reports.

## Parked branches (not merged)

- `parked/subsume-pruning`: containment identities (period divisor, base power),
  opt-in pruning, and implied exclusion bounds in the report. Tested (23 tests). Parked because
  it only makes the closed search cheaper.
- `parked/multilevel-dd`: double-double inner loop (`multilevel-dd`). Validated
  equivalent. **Negative result**: half the syncs, but 2–6× slower at every n, because
  double-double square roots and divides dominate the corner rotation.

## Next direction

Attack the cost of *existing* decimal digit extraction (Plouffe/Bellard/Gourdon, ~O(n²) time,
tiny memory) with algorithmic tricks instead of new formulas. The candidate idea is
baby-step/giant-step evaluation of the modular products in the extraction sum (Strassen;
Bostan–Gaudry–Schost), aiming at roughly O(n^1.5) time with O(√n) memory. Step one is a
literature check.
