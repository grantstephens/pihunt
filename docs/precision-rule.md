# Precision rule

`auto_digits(n, C)` picks the working precision (decimal digits) for a job with `n` basis
columns and coefficient bound `C`:

    digits = ceil(n * log10(C) * f(n)) + 50

## The problem with a flat factor

The original factor `f = 1.5` was measured in prototyping at n ≈ 45 (1.25 let spurious
10^7-size relations through there; 1.5 excluded cleanly). It holds up fine through the
sample-scout and first exhaustive-grid regime (n ≲ 40), but thins out badly at n ≈ 100: e.g.
base 10, period 31, degrees 1..3, extras `[pi2, log2, log3, log5, catalan, zeta3]` (n = 100),
C = 10^5 → 800 digits → `MultilevelPslq` came back `Inconclusive`. The *reduction* step (PSLQ
over the non-π columns) ran out of precision with the exclusion bound at only 6.3e4 (needs
`C·√n ≈ 1e6`) and hit a precision-floor "relation" with coefficients ~2.7e10, which the checks
correctly rejected. The runner's escalation then reruns the whole job from scratch at 2×
digits — correct, but expensive to discover by trial at large n.

A second, distinct failure mode showed up in a finished exhaustive sweep at the old flat rule
(C = 10^5): 26 jobs (e.g. base 10, period 20-21, degrees 1..3, and period 30-32, degrees 1..2,
at 530-583 digits) ended `Suspicious` instead of `Inconclusive` — the *main* search (not
reduction) found a "relation" with coefficients far past `coeff_bound`, another precision-floor
near-miss that the checks correctly reject. This is just as much a sign of too few digits as
`Inconclusive`; the measurement below treats both as "not yet resolved" and only accepts
`Excluded` or `Hit` as a genuine answer.

## Measurement

`tests/precision_rule.rs` (`cargo test --release --test precision_rule -- --ignored
--nocapture`) sweeps 22 shapes (bases 10/100/1000, n from 7 to 100 — periods, degree ranges and
extras sets modelled on `batches/*.toml`, including the specific period/degree combinations
flagged `Suspicious` above) × C ∈ {10^3, 10^4, 10^5}. For each (shape, C) it tries factors
1.25, 1.5, 1.75, 2.0, 2.5, 3.0, 3.5 in increasing order against `run_job` with
`MultilevelPslq`, and reports the smallest one whose outcome is `Excluded` or `Hit` (not
`Inconclusive`, `Suspicious`, or `Spurious`).

Measured on an idle 6-core Ryzen 7 5700G (a prior run under a 6-core-busy sweep produced the
same minimal factors, just noisier — one case swung to over 600s from thread contention with
no other change; seconds below should be read as rough, not exact, per `docs/timing-baseline.md`).

| base | period | degrees | #extras | n | C | minimal f | outcome | seconds |
|---|---|---|---|---|---|---|---|---|
| 10 | 6 | 1..1 | 0 | 7 | 1e3 | 1.25 | Excluded | 0.00 |
| 10 | 6 | 1..1 | 0 | 7 | 1e4 | 1.25 | Excluded | 0.00 |
| 10 | 6 | 1..1 | 0 | 7 | 1e5 | 1.25 | Excluded | 0.00 |
| 10 | 10 | 1..1 | 2 | 13 | 1e3 | 1.25 | Excluded | 0.00 |
| 10 | 10 | 1..1 | 2 | 13 | 1e4 | 1.25 | Excluded | 0.01 |
| 10 | 10 | 1..1 | 2 | 13 | 1e5 | 1.25 | Excluded | 0.01 |
| 100 | 10 | 1..1 | 3 | 14 | 1e3 | 1.25 | Excluded | 0.01 |
| 100 | 10 | 1..1 | 3 | 14 | 1e4 | 1.25 | Excluded | 0.01 |
| 100 | 10 | 1..1 | 3 | 14 | 1e5 | 1.25 | Excluded | 0.01 |
| 1000 | 14 | 1..1 | 0 | 15 | 1e3 | 1.25 | Excluded | 0.01 |
| 1000 | 14 | 1..1 | 0 | 15 | 1e4 | 1.25 | Excluded | 0.01 |
| 1000 | 14 | 1..1 | 0 | 15 | 1e5 | 1.25 | Excluded | 0.01 |
| 10 | 8 | 1..2 | 0 | 17 | 1e3 | 1.25 | Excluded | 0.01 |
| 10 | 8 | 1..2 | 0 | 17 | 1e4 | 1.25 | Excluded | 0.02 |
| 10 | 8 | 1..2 | 0 | 17 | 1e5 | 1.25 | Excluded | 0.02 |
| 10 | 12 | 1..1 | 6 | 19 | 1e3 | 1.25 | Excluded | 0.02 |
| 10 | 12 | 1..1 | 6 | 19 | 1e4 | 1.25 | Excluded | 0.02 |
| 10 | 12 | 1..1 | 6 | 19 | 1e5 | 1.25 | Excluded | 0.03 |
| 10 | 12 | 1..2 | 2 | 27 | 1e3 | 1.25 | Excluded | 0.08 |
| 10 | 12 | 1..2 | 2 | 27 | 1e4 | 1.25 | Excluded | 0.11 |
| 10 | 12 | 1..2 | 2 | 27 | 1e5 | 1.25 | Excluded | 0.14 |
| 100 | 16 | 1..2 | 3 | 36 | 1e3 | 1.5 | Excluded | 0.26 |
| 100 | 16 | 1..2 | 3 | 36 | 1e4 | 1.5 | Excluded | 0.37 |
| 100 | 16 | 1..2 | 3 | 36 | 1e5 | 1.5 | Excluded | 0.58 |
| 10 | 16 | 1..2 | 6 | 39 | 1e3 | 1.5 | Excluded | 0.52 |
| 10 | 16 | 1..2 | 6 | 39 | 1e4 | 1.5 | Excluded | 0.64 |
| 10 | 16 | 1..2 | 6 | 39 | 1e5 | 1.25 | Excluded | 0.79 |
| 10 | 20 | 1..2 | 6 | 47 | 1e3 | 1.75 | Excluded | 1.08 |
| 1000 | 20 | 1..2 | 6 | 47 | 1e3 | 1.75 | Excluded | 1.35 |
| 10 | 20 | 1..2 | 6 | 47 | 1e4 | 1.5 | Excluded | 1.25 |
| 1000 | 20 | 1..2 | 6 | 47 | 1e4 | 1.5 | Excluded | 1.14 |
| 10 | 20 | 1..2 | 6 | 47 | 1e5 | 1.5 | Excluded | 1.69 |
| 1000 | 20 | 1..2 | 6 | 47 | 1e5 | 1.5 | Excluded | 1.50 |
| 10 | 16 | 1..3 | 2 | 51 | 1e3 | 1.75 | Excluded | 1.31 |
| 10 | 16 | 1..3 | 2 | 51 | 1e4 | 1.75 | Excluded | 1.92 |
| 10 | 16 | 1..3 | 2 | 51 | 1e5 | 1.5 | Excluded | 2.23 |
| 100 | 24 | 1..2 | 6 | 55 | 1e3 | 1.75 | Excluded | 2.40 |
| 100 | 24 | 1..2 | 6 | 55 | 1e4 | 1.5 | Excluded | 3.09 |
| 100 | 24 | 1..2 | 6 | 55 | 1e5 | 1.5 | Excluded | 3.93 |
| 10 | 20 | 1..3 | 6 | 67 | 1e3 | 2.0 | Excluded | 5.22 |
| 10 | 30 | 1..2 | 6 | 67 | 1e3 | 2.0 | Excluded | 5.97 |
| 10 | 20 | 1..3 | 6 | 67 | 1e4 | 1.75 | Excluded | 6.21 |
| 10 | 30 | 1..2 | 6 | 67 | 1e4 | 1.75 | Excluded | 6.40 |
| 10 | 20 | 1..3 | 6 | 67 | 1e5 | 1.75 | Excluded | 8.76 |
| 10 | 30 | 1..2 | 6 | 67 | 1e5 | 1.75 | Excluded | 10.84 |
| 10 | 31 | 1..2 | 6 | 69 | 1e3 | 2.0 | Excluded | 6.74 |
| 10 | 31 | 1..2 | 6 | 69 | 1e4 | 1.75 | Excluded | 9.48 |
| 10 | 31 | 1..2 | 6 | 69 | 1e5 | 1.75 | Excluded | 11.51 |
| 10 | 21 | 1..3 | 6 | 70 | 1e3 | 2.0 | Excluded | 7.01 |
| 10 | 21 | 1..3 | 6 | 70 | 1e4 | 1.75 | Excluded | 8.79 |
| 10 | 21 | 1..3 | 6 | 70 | 1e5 | 1.75 | Excluded | 105.63 |
| 10 | 32 | 1..2 | 6 | 71 | 1e3 | 2.0 | Excluded | 7.02 |
| 10 | 32 | 1..2 | 6 | 71 | 1e4 | 1.75 | Excluded | 9.95 |
| 10 | 32 | 1..2 | 6 | 71 | 1e5 | 1.75 | Excluded | 11.92 |
| 10 | 24 | 1..3 | 6 | 79 | 1e3 | 2.5 | Excluded | 14.40 |
| 10 | 24 | 1..3 | 6 | 79 | 1e4 | 2.0 | Excluded | 20.63 |
| 10 | 24 | 1..3 | 6 | 79 | 1e5 | 1.75 | Excluded | 24.41 |
| 10 | 28 | 1..3 | 6 | 91 | 1e3 | 2.5 | Excluded | 28.12 |
| 1000 | 28 | 1..3 | 6 | 91 | 1e3 | 2.5 | Excluded | 21.02 |
| 10 | 28 | 1..3 | 6 | 91 | 1e4 | 2.0 | Excluded | 148.43 |
| 1000 | 28 | 1..3 | 6 | 91 | 1e4 | 2.0 | Excluded | 33.13 |
| 10 | 28 | 1..3 | 6 | 91 | 1e5 | 1.75 | Excluded | 40.81 |
| 1000 | 28 | 1..3 | 6 | 91 | 1e5 | 1.75 | Excluded | 40.35 |
| 10 | 31 | 1..3 | 6 | 100 | 1e3 | 2.5 | Excluded | 34.23 |
| 10 | 31 | 1..3 | 6 | 100 | 1e4 | 2.0 | Excluded | 35.29 |
| 10 | 31 | 1..3 | 6 | 100 | 1e5 | 2.0 | Excluded | 35.57 |

(One row's 148.43s and one 105.63s figure are outliers relative to their neighbours —
almost certainly scheduling noise from running ~66 cases in parallel with `rayon` on 6
cores, not a real cost cliff; every neighbouring (n, C) pair with the same minimal `f`
lands 3-5× faster.)

Two things stand out:

1. **Required `f` climbs with `n`, roughly**: 1.25 up to n ≈ 27, 1.5 by n ≈ 36-39, 1.75 by
   n ≈ 47-55, 2.0 by n ≈ 67-71, 2.5 by n ≈ 79-100 (worst case over the tested C's).
2. **Required `f` falls as `C` grows**, for fixed n (e.g. n = 91: f = 2.5 at C = 1e3, 2.0 at
   C = 1e4, 1.75 at C = 1e5). This makes sense: digits already scale with `log10(C)`, so a
   bigger `C` buys more absolute digits at the same `f`; the shortfall that shows up at large n
   is a roughly fixed amount of "PSLQ working precision" that doesn't scale with `log10(C)`.
   Using the worst-case-over-C factor for a given n is therefore safe for every C in the
   tested range (1e3-1e5) — it can only ask for more digits than a larger C strictly needs.

## The new rule

    f(n) = 1.5 + 0.025 * max(0, n - 36)
    digits(n, C) = ceil(n * log10(C) * f(n)) + 50

This is the simplest curve of the shape the task suggested (`f(n) = 1.5 + c·max(0, n - n0)`)
that sits at or above every measured point above, fit to the (36, 1.5) and (79, 2.5) corners
of the worst-case-over-C curve. Checked against the whole table:

| n | required f (worst case over C) | rule f(n) | margin |
|---|---|---|---|
| ≤27 | 1.25 | 1.5 (clamped) | +0.25 |
| 36 | 1.5 | 1.5 | 0 (exact) |
| 39 | 1.5 | 1.575 | +0.075 |
| 47 | 1.75 | 1.775 | +0.025 |
| 51-55 | 1.75 | 1.875-1.975 | +0.125-0.225 |
| 67-71 | 2.0 | 2.275-2.375 | +0.275-0.375 |
| 79 | 2.5 | 2.575 | +0.075 |
| 91 | 2.5 | 2.875 | +0.3 |
| 100 | 2.5 | 3.1 | +0.6 |

Below n = 36 the rule is byte-for-byte the original flat 1.5 (`auto_digits(9, 1000)` is still
91, unchanged), so no small-n job gets more digits than before. At n = 100, C = 1e5 it asks
for 1600 digits, roughly double the old rule's 800 — but that's paid once, up front, instead
of discovered by an `Inconclusive` round-trip through the whole job (including a full
reduction pass) followed by an escalation rerun from scratch.

## Sanity check: the reported n = 100 example

Base 10, period 31, degrees 1..3, extras `[pi2, log2, log3, log5, catalan, zeta3]` (n = 100),
C = 10^5:

Measured directly (`tests/precision_rule.rs::sanity_check_n100_example`, idle 6-core machine,
`MultilevelPslq`, `max_iterations = 300_000`):

| | old rule (flat 1.5) | new rule |
|---|---|---|
| digits | 800 | 1600 |
| outcome | `Inconclusive` — "reduction found a relation that failed checks (max \|coeff\| 27141133321)" | `Excluded` |
| iterations | 129228 (reduction only, then gives up) | 155868 (full job) |
| seconds | 14.68 (then a full-cost rerun at 1600 digits follows, ≈ 62s more) | 62.35 |

The old rule's 800 digits reproduces the reported failure exactly: the reduction step's
"relation" has a ~2.7e10-magnitude coefficient, just as described, and the job ends
`Inconclusive` after 129228 iterations. Under the runner's escalation that's followed by a
full rerun at 1600 digits — coincidentally the same digit count the new rule picks directly —
so this particular case's total cost (~14.7s wasted + ~62s to actually resolve ≈ 77s) drops to
just the ~62s the new rule needs on its first and only attempt. The real payoff is structural,
not this one case's wall-clock: every large-n job in a batch that would have needed a first,
doomed attempt at the old digit count now skips straight to one that resolves.

## What this doesn't fix

Escalation stays as the safety net for anything the measurement above didn't cover (larger n,
smaller C, different extras mixes, or the `Suspicious`/`Spurious` precision-floor artifacts
noted above generally — those are being made to trigger escalation too, independently of this
change). This rule only front-loads the *common* shortfall so most large-n jobs settle on the
first attempt instead of the second.
