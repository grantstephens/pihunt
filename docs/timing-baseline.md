# Timing baseline — classic PSLQ

2026-09-22 · AMD Ryzen 7 5700G · pihunt `28f4996` · single job per row (no rayon), `coeff_bound = 1000`, auto precision.

Base 10, shapes from `tests/timing.rs`. Seconds include basis build, reduction and the main search.

| n | digits | outcome | iterations | seconds |
|---|---|---|---|---|
| 9 | 91 | Excluded | 431 | 0.01 |
| 19 | 136 | Excluded | 2297 | 0.18 |
| 27 | 172 | Excluded | 5726 | 0.92 |
| 39 | 226 | Excluded | 12506 | 5.11 |
| 47 | 262 | Inconclusive | 17916 | 5.57 |

Scaling is roughly n⁴–n⁵. This is the reference the stage-2 multi-level PSLQ gets measured against.
The n = 47 row ends inconclusive: reduction found a relation that failed the coefficient/2P checks.
