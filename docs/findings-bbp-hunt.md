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

## Bellard's formula as a rediscovery check

`batches/bellard.toml` runs Bellard's 1997 base-2^10 formula through the same pipeline. The
`(-1)^k / 2^(10k)` weight needs a negative base, which `Shape` doesn't support, so pairs of
terms are folded together: base 2^20, period 40, 41 columns, `C = 10^7`, 517 digits. PSLQ finds
it in about 38 s (classic finder, one core) and the hit verifies at 2× precision. `src/known.rs`
holds two forms of it, the directly folded one and the one PSLQ returns after the reduction
drops dependent columns; they differ by zero relations, and `tests/known.rs` checks both
numerically. This is a larger and differently-shaped sanity check than BBP at b=16, m=8.

## A different territory: base-2^k, and the first real use of Latin hypercube sampling

Stages 1–3 above are specifically a **base-10** closure (the `## Why nothing was found` argument
leans on item 1, which only rules out non-power-of-2 bases) — it says nothing about base 2^k,
which is exactly where the Bellard rediscovery above and the original 1996 BBP formula
(`results/known.jsonl`, base 16 = 2^4) both actually live. No batch before this one ever set
`bases` to anything but `[10, 100, 1000]`, so this is genuinely untouched territory, not a rerun
of the closed one.

`batches/base2-exotic.toml`: 17 powers of 2 (2^5 through 2^24, skipping 2^4/2^6/2^20 — already
covered by `known.toml`/`bellard.toml`) plus 7 small non-power-of-2 bases (3, 5, 6, 7, 9, 11, 13)
included at the time "for curiosity," periods 2–48, degrees 1..1 and 1..2, the same extras sets as
the base-10 sweeps. `mode = "sample"` / `method = "lhs"` — `src/plan.rs`'s Latin hypercube sampler
(stratified permutation per axis + jitter, implemented alongside Sobol when grid/sobol/lhs
planning was added) had never actually been exercised by any batch until this one; every prior
sampled sweep used Sobol.

**Correction:** those 7 non-power-of-2 bases were not actually open territory. Borwein, Borwein &
Galway (2004) — "Finding and Excluding b-ary Machin-Type Individual Digit Formulae," *Canad. J.
Math.* 56:897–925 — prove there are no degree-1 Machin-type BBP formulas for π in *any* base that
isn't a power of 2, not just base 10; base 3 is one of their own worked examples. That's the exact
same closure as `## Why nothing was found` item 1 above, just not yet pointed at these specific
base values. The 0 hits found for bases 3/5/6/7/9/11/13 below confirm theory rather than being a
new negative result, and the degree-2 cases at those bases fall under the same conditional
weight-conservation argument (item 2) base 10 already relies on — nothing in that argument is
base-10-specific. `batches/base2-focused.toml` (below) drops them and spends the freed sample
budget on more powers of 2, the only part of this space theory doesn't already foreclose.

Result (`results/base2-exotic.jsonl`, `docs/exclusions-base2-exotic.md`): 150 shapes sampled,
**147 excluded** (proved coefficient bound, no relation), **0 hits** (no π relation, new or
otherwise), **3 unresolved** — one skipped (103 columns, over the 100 cap) and two genuinely
inconclusive even after 3x precision escalation (base 2^24, periods 8 and 16, full extras set;
PSLQ's own reduction found a candidate relation that failed verification both times, the same
"spurious near-relation" category the base-10 wide sweep's suspicious jobs were). The relations
PSLQ *did* find among non-π columns are the same expected family as the base-10 searches —
logarithmic identities among the series and log2/log3/log5, nothing involving π.

Read this the way the base-10 result should *not* be read: as one 150-point sample of an
effectively unbounded (base, period, degree, extras) space, not a closure. There's no
Borwein-Borwein-Galway-style theorem backing a "nothing's here" claim for base 2^k the way there
is for base 10 — quite the opposite, two formulas are already known to live here. A null result
from one sample says only that *this* sample didn't find a third one.

## base2-focused: powers of 2 only, wider periods — a "NEW" hit that wasn't

`batches/base2-focused.toml`: same 17 powers of 2 as `base2-exotic.toml` (no non-power-of-2 bases
— see the correction above), periods widened 2–64 (was 2–48), LHS count raised 150 → 220, new
seed. 219 shapes planned, 190 run (29 skipped, >100-column cap — periods up to 64 push column
counts past it more often than periods up to 48 did).

Result (`results/base2-focused.jsonl`, `docs/exclusions-base2-focused.md`): **188 excluded**, **1
hit** (tagged `derived`, not `NEW` — see below), **1 inconclusive** (base 2^24, period 12, same
spurious-near-relation category as before — PSLQ's reduction found a candidate that failed
verification), **29 skipped**.

The hit looked like a real find at first: `b=256 m=32 s=1..1 x=-`, verified at 2× precision with
residual ~10⁻⁷⁹³ (`8·pi − 64·S(2) + 32·S(8) + 16·S(10) + 16·S(12) − 4·S(18) + 2·S(24) + 1·S(26) +
1·S(28) = 0`, job `026d35c0e5a09b503922e695c708045a`). Checked it by hand before reporting it as
anything: every j in the relation is even, and `S(b,m,2j',1) = ½·S(b,m/2,j',1)`, which collapses
it to base 256, period 16, j' ∈ {1,4,5,6,9,12,13,14} — the classic 1996 BBP residues `{1,4,5,6}`
plus the same set shifted by the half-period. That's exactly what splitting BBP's series by parity
of k produces (`S(b,m,j,1) = S(b²,2m,j,1) + (1/b)·S(b²,2m,m+j,1)`), then period-doubled again.
Verified numerically (mpmath, 300 digits): the PSLQ relation and "classic BBP, radix-split then
period-doubled" both evaluate to the same value, residual ~10⁻³⁰⁰ either way. **It's the 1996 BBP
formula, not a new one.**

This generalizes beyond this one case: any base `b0^t` (any known formula's base raised to any
power) at period `m0·t·d` (its period times `t` times any integer `d`) inherits a mechanical,
non-independent corollary this same way — composing the base-power radix split with a period
multiply. `known::is_derived` (`src/known.rs`, tests in `tests/known.rs` and `tests/classify.rs`)
now checks exactly this before tagging a hit `NEW`: it reconstructs the predicted descendant
relation for every entry in the known-formula table and every `(t, d)` the candidate's shape
admits, and compares it to the actual relation up to the usual gcd/sign normalisation. It caught
this one automatically — `report` now shows `0 NEW` for this batch, tag `derived`, and the stored
record carries a note identifying the source formula and `(t, d)`. It's a small slice of the
parked `subsume-pruning` branch's idea (base-power and period-divisor containment identities),
built because this specific case needed it, not a replacement for that branch's broader pruning.

Net effect: still **zero independent new relations** in base-2^k territory after two LHS sweeps
(300 shapes total). Of the 17 bases swept, 5 (256, 4096, 65536, 262144, 16777216 — powers of 16
and/or 64) are reachable this way from BBP or Bailey's π² formula and will keep producing `derived`
hits rather than useful information at the right period; the other 12 aren't reachable from any
known formula by this mechanism, so a hit there would still mean something.

## base2-clean and base2-degree1: closeness, a regression that found its own limit, and a prior that didn't pan out

Two follow-ups, both negative but worth recording precisely.

**`batches/base2-clean.toml`** dropped the 5 derived-corollary bases and added a new per-job
field, `Outcome::closeness` (`src/pslq/mod.rs`, `src/pslq/state.rs`): the smallest `|y_i|` in
PSLQ's state at termination, divided by the detection threshold `eps`. Always ≥ 1 for a
non-relation; the closer to 1, the nearer that shape came to collapsing onto a relation,
independent of whether it happened to clear the coefficient bound. Pure logging addition, no
effect on exclusion/hit logic (tests in `tests/known.rs` adjacent — see `tests/classify.rs`'s new
case — plus `#[serde(default)]` so old logs without the field still parse). 260 LHS-sampled shapes
over the 12 clean bases, periods 2–64, degrees `[1,1]`/`[1,2]`, same three extras sets as before:
217 excluded, 0 hits, 37 skipped, 0 inconclusive.

Fitting `log10(closeness) ~ columns + precision_digits + log2(base) + period` (ordinary least
squares, 169 excluded jobs with both degree ranges) got R² = 0.997: `closeness` is almost entirely
explained by how big the search was, not by anything π-structural. That's a real result, not a
failure to find one — it means `closeness` alone isn't a usable discovery signal, and ruling that
out cheaply (one regression) beat building the complex-step/gradient machinery that was briefly
considered for the same question and would have hit the same wall (the actual target —
relation-existence at a specific integer shape — isn't a smooth function of base/period at all;
no continuous relaxation of the search makes that fact differentiable in a meaningful way).

**`batches/base2-degree1.toml`** acted on the one lever that *did* hold up: the weight-conservation
conjecture already backing this doc's base-10 closure doesn't depend on the base being 10, so 5 of
the 6 (degree, extras) combinations every batch above had been splitting budget across are
conjecturally reducible to the pure degree-1 case — confirmed by what they kept finding (log
identities, nothing π-related). Pinning `degrees=[[1,1]]`, `extras=[[]]` shrank the space enough
to grid-search exhaustively instead of sampling: all 12 clean bases, every period 2–64, 756 jobs,
no skips (max 65 columns, well under the 100 cap). **756/756 excluded, 0 hits.**

That grid also let us test the one structural prior in hand — BBP (base 2⁴, period 8) and Bellard
(base 2²⁰, period 40) both have period = 2× the base's exponent of 2 — against real data instead
of two examples. It didn't hold: periods within ±4 of `2·log2(base)` have a *higher* mean
dimension-controlled residual than periods more than 20 away (+0.61 vs −0.61), i.e. if anything
slightly farther from a relation, not closer. Two data points was too thin to extrapolate from,
and now there's a reason better than intuition to say so.

Net: the entire open slice of base-2^k territory up to period 64 — every base not already a
corollary of a known formula, pure degree-1, no extras — is now exhaustively checked and
confidently excluded. Extending past period 64 or past base 2²⁴ is the only way to keep searching
here, and nothing above suggests where in that extension to look first.

## multibase: does π need two base towers at once?

Every search above, and both real formulas, fix one base per relation. BBP and Bellard are each
just "π = a Q-linear combination of several weight-1 periods" at a *single* base — nothing
requires those periods to share a base, and mixing columns from two different clean bases is
still weight-1 + weight-1, so the weight-conservation conjecture that ruled out degree/extras
mixing doesn't apply. Untried, and not a stretch.

`examples/multibase.rs`: standalone (not retrofitted into `Shape`/the batch-TOML machinery, which
is one-base-per-shape throughout — reuses `basis::series` and the multilevel PSLQ engine
directly). All `C(12,2) = 66` pairs of the clean bases (deliberately excluding any base already
known to carry a formula or a `known::is_derived` descendant of one — pairing a known base with
anything trivially "succeeds" by rediscovering it with zero coefficients on the other side),
shared period swept 2–24, degree 1, no extras. 1518 jobs.

Result: **1472 excluded, 0 π relations, 46 flagged `basis_relation_no_pi`**. Checked those 46
before trusting them: all land on 7 of the 66 pairs, and every one of those 7 has base exponents
sharing a common factor — e.g. 32 = 2⁵ and 32768 = 2¹⁵ (gcd 5, and 32768 = 32³ exactly); 512 = 2⁹
and 2097152 = 2²¹ (gcd 3, both powers of 8). Same radix-split mechanism `known::is_derived` already
catches, just between two clean bases instead of a clean base and a known formula — pure algebra
from the base list containing a few numbers that are powers of each other, nothing π-related.
Verified one numerically (mpmath, 150 digits, residual ~10⁻¹⁴⁸) rather than taking the pattern on
faith. The other 59 pairs are clean of this and simply excluded.

Net: no evidence π needs two base-2^k towers at once, at least in this slice. Another real,
previously-untried direction, run cleanly, came back empty.

**Arity 3** (`ARITY = 3` in the same file, `results/multibase3.jsonl`): the user's next question
was "what about 7 bases at once?" — declined, and worth recording why. Columns = 1 + arity·m, and
`auto_digits`'s precision factor grows with column count past n ≈ 36; at arity 7, m = 24 that's
169 columns needing ~5,700 digits, well past n ≈ 100 (the range `docs/precision-rule.md` actually
measured). Past validated territory, PSLQ gets *more* likely to turn up a small-coefficient
relation from lattice-reduction noise rather than genuine structure — more free columns, more
chances for a spurious near-miss — and the weight-conservation argument that cleared multi-base
combos in the first place doesn't favour more bases over fewer, so higher arity isn't a more
promising bet, just a less trustworthy one. Arity 3 (`C(12,3) = 220` triples × periods 2–24,
5060 jobs, ≤ 73 columns) is the sane next step: still comfortably inside validated territory.

Result: **4610 excluded, 0 π relations, 450 `basis_relation_no_pi`, 0 inconclusive** (43 min, vs
1:45 for arity 2 — the cost already climbing as expected). All 450 non-π hits are explained by
the triple containing at least one of the 9 base-pairs (out of the 12 clean bases' `C(12,2)=66`)
whose exponents share a common factor — the same radix-split artifact the arity-2 run surfaced,
just inherited into any triple that contains one of those pairs. Checked exhaustively, not
sampled: zero of the 450 are unexplained by that mechanism.

## Status: base-2^k search (2026-10-08)

Summary of everything above, for anyone who doesn't want to read all of it: **no new BBP-type
formula for π was found anywhere in base-2^k territory**, across six independent, theory-motivated
searches, each checked for soundness before being trusted:

| search | shapes | result |
|---|---|---|
| `base2-exotic` (LHS, 17 powers of 2 + 7 non-power-of-2) | 150 | 0 hits, 3 unresolved |
| `base2-focused` (LHS, same 17 bases, wider periods) | 219 | 1 hit — proved to be the 1996 BBP formula in disguise, not new |
| `base2-clean` (LHS, 12 bases not derived from a known formula) | 260 | 0 hits — `closeness` regression (R²=0.997) showed this metric is dimension-driven, not structural |
| `base2-degree1` (exhaustive grid, 12 clean bases × periods 2–64, degree-1 only) | 756 | 0 hits — the entire theoretically-open slice up to period 64, checked completely |
| `multibase` arity 2 (all pairs of the 12 clean bases × periods 2–24) | 1518 | 0 π hits, 46 explained non-π artifacts |
| `multibase` arity 3 (all triples × periods 2–24) | 5060 | 0 π hits, 450 explained non-π artifacts |

**What's actually been ruled out**: every degree-1, no-extras shape at any of the 12 bases that
aren't a mechanical corollary of BBP/Bailey/Bellard, for every period up to 64; every pairwise and
triple-wise combination of those same 12 bases, for every shared period up to 24; every
(degree, extras) combination theory doesn't already predict is futile (confirmed, not assumed —
the base2-clean/degree1 split is why we trust that reduction).

**What's not ruled out, and why we're not chasing it further right now**: periods past 64,
base exponents past 24, and multibase arity past 3. Each is a real extension, not closed by
theory — but each costs more than the last for no principled reason to expect a different answer,
and arity > 3 specifically risks precision past where PSLQ's soundness margin is validated
(`docs/precision-rule.md`). Diminishing returns, not a wall.

**Engineering kept from this phase** (beyond the base-10 phase's list above):
- `known::is_derived` (`src/known.rs`) — catches base-power + period-multiple corollaries of any
  known formula before tagging a hit `NEW`. Caught a real false positive during this work
  (`base2-focused`'s base-256 hit); not hypothetical.
- `Outcome::closeness` (`src/pslq/mod.rs`, `src/pslq/state.rs`) — `min|y_i|/eps` at termination,
  logged for every job. Real signal, just not an independent one once dimension/precision are
  controlled for — that's a tested, reusable fact now, not a guess.
- `examples/multibase.rs` — generic arity-N multi-base search, standalone from `Shape` on purpose.

## Next direction

Attack the cost of *existing* decimal digit extraction (Plouffe/Bellard/Gourdon, ~O(n²) time,
tiny memory) with algorithmic tricks instead of new formulas. The candidate idea is
baby-step/giant-step evaluation of the modular products in the extraction sum (Strassen;
Bostan–Gaudry–Schost), aiming at roughly O(n^1.5) time with O(√n) memory. Step one is a
literature check.
