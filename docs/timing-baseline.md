# Timing baseline

AMD Ryzen 7 5700G · `coeff_bound = 1000`, auto precision · one job per row (no rayon).
Base 10, shapes from `tests/timing.rs`. Seconds include basis build, reduction and the main search.
Regenerate with `cargo test --release --test timing -- --ignored --nocapture`.

## Stage 2 (2026-09-22) — classic vs multilevel

| n | digits | finder | outcome | iterations | seconds |
|---|---|---|---|---|---|
| 9 | 91 | classic | Excluded | 431 | 0.01 |
| 9 | 91 | multilevel | Excluded | 433 | 0.00 |
| 19 | 136 | classic | Excluded | 2297 | 0.18 |
| 19 | 136 | multilevel | Excluded | 2299 | 0.03 |
| 27 | 172 | classic | Excluded | 5726 | 0.89 |
| 27 | 172 | multilevel | Excluded | 5856 | 0.14 |
| 39 | 226 | classic | Excluded | 12506 | 5.17 |
| 39 | 226 | multilevel | Excluded | 12179 | 0.72 |
| 47 | 262 | classic | Inconclusive | 17916 | 5.60 |
| 47 | 262 | multilevel | Inconclusive | 18016 | 0.79 |
| 55 | 298 | multilevel | Inconclusive | 24880 | 2.06 |
| 67 | 352 | multilevel | Inconclusive | 35570 | 3.2 |

Multilevel is ~7× faster than classic at n = 27–47. Classic rows above n = 47 are skipped
(minutes to hours). The n = 47+ rows are inconclusive at 1× digits (reduction relations
fail checks); escalation resolves those in real runs.

Where multilevel's time goes at n = 67 (30 syncs): f64 inner loop ~1.0 s, sync ~2.2 s, of
which the Givens LQ is ~1.4 s — already allocation-free and near MPFR's multiply speed at
~1170 bits. Next speedups are algorithmic: Householder LQ (~2×), or a three-level variant
that keeps H at medium precision.

## Stage 1 (2026-09-22, `28f4996`) — classic only

| n | digits | outcome | iterations | seconds |
|---|---|---|---|---|
| 9 | 91 | Excluded | 431 | 0.01 |
| 19 | 136 | Excluded | 2297 | 0.18 |
| 27 | 172 | Excluded | 5726 | 0.92 |
| 39 | 226 | Excluded | 12506 | 5.11 |
| 47 | 262 | Inconclusive | 17916 | 5.57 |
