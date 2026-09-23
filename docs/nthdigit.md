# N-th decimal digit extraction

**Status: implemented, benchmarked, verified against MPFR through 10⁶.** Follows on from
`docs/findings-bbp-hunt.md`'s "next direction" (attacking the cost of existing decimal
digit extraction rather than hunting for new formulas). Unrelated to the PSLQ code
elsewhere in this crate.

## How it works

Xavier Gourdon (2003), Theorem 1. Start from `π/4 = Σ (-1)^k/(2k+1)` and accelerate it with
the Cohen–Villegas–Zagier process using `P(x) = x^M(1-x)^N`, giving

```
S = Σ_{k=0}^{(M+1)N-1} (-1)^k 4/(2k+1)  −  Σ_{k=0}^{N-1} (-1)^k 4 s_k / (2^N (2MN+2k+1))
s_k = Σ_{j=0}^{k} C(N,j)
|S − π| ≤ π/(2eM)^N
```

Choosing `M = 2⌈n/ln³n⌉` (even, ≥ 4) and `N = ⌈(n+n0+1)ln10/ln(2eM)⌉` (even) makes `N ≤ n+2`,
which is the key trick: it makes `4·10^n` and `5^(N-2)·10^(n-N+2)` *integers*, so `frac(10^n
π)` reduces to a sum of ordinary `(integer mod small-integer)/small-integer` fractions —
Gourdon's Proposition 1. Every one of those small integers fits in a `u64`, so the whole
computation is machine-word modular arithmetic, no big-number library required.

The only non-trivial piece is `s_k mod m` (a sum of up to `N` binomial coefficients, mod an
`m` that can be composite and share prime factors with `j ≤ k`, which would otherwise make
`j` non-invertible mod `m`). Algorithm 2 handles this by tracking, for each prime factor
`p_i ≤ k` of `m`, a running *exact* integer `R_i = p_i^{e_i}` (Gourdon proves `e_i ≥ 0`
always, so `R_i` never needs a modular inverse) alongside the usual running products `A`,
`B`, `C` mod `m`, recovering `s_k mod m = C·B⁻¹ mod m` at the end. The complement identity
`Σ_{j≤k} C(N,j) = 2^N − Σ_{j≤N-k-1} C(N,j)` halves the work for `k > N/2`.

## Finding `m`'s prime factors: segmented sieve, not trial division

Algorithm 2 needs the primes `≤ k` dividing each `m_k = 2MN + 2k + 1` before it can run its
`O(k)` binomial loop. The original version found them by trial-dividing `m_k` against
`2..=k`, one `m_k` at a time — correct, and only `O(log m)` memory, but `O(k)` time per `k`,
i.e. `O(N²)` total, and (per the original profiling this task started from) assumed to be
the dominant cost of the whole C-sum.

It's replaced with a **segmented sieve** (`factor_segment` in `src/nthdigit.rs`): sieve the
small primes `≤ √(max m_k)` once per call (`primes_up_to`, plain sieve of Eratosthenes — this
bound is `O(√(MN))`, not `O(N)`, so it's cheap even though it's the one place in this module
whose memory actually scales with `n`), then chunk `k ∈ [0, N)` into segments (size tuned to
give rayon `~32` segments per thread — see `SEGMENTS_PER_THREAD`, and the "segment size
matters more than you'd think" note below) and, for each segment, locate every small prime's
hits via modular arithmetic (`m_k ≡ 0 (mod p) ⟺ k ≡ -c0·inv2(p) (mod p)`, `inv2(p) =
(p+1)/2` for odd `p` — no extended-gcd needed) instead of dividing. Whatever's left of `m_k`
after dividing out every small prime is `1` or a single prime `> √(max m_k)`, included if
it's `≤ k`.

The strict-memory trial-division path is kept (`c_sum_trial_division`, `primes_le_k_dividing`,
`sum_binomials_mod`), selected by the `USE_SIEVED_FACTORING` const in `src/nthdigit.rs`. It
exists for exactly the case Gourdon's theorem is about: if the sieve's extra memory (see
below) is ever actually a problem, flip that flag back to `false` and get the original
`O(log² n)`-memory, slower algorithm. It isn't a problem at any `n` measured here.

**A segmenting trap worth naming:** the first version of this fixed a single `SEGMENT_SIZE =
2^15` for every `N`. At n=10⁴..10⁶, `N` is only in the thousands to hundreds-of-thousands, so
that collapsed the whole C-sum into 1-2 rayon tasks — i.e. *no* parallelism at all — and made
n=10⁵ four to five times *slower* than trial division, before the per-`k` factoring savings
had any chance to matter. Fixed by sizing segments dynamically off `N` and the thread count
(`segment_size_for`) so there are always roughly `32×` as many segments as threads, letting
rayon's work-stealing even out Algorithm 2's per-`k` cost skew (the last segment costs ~`N`
times more than the first).

## Montgomery multiplication in the hot loop

Profiling (see below) found the segmented sieve alone barely moved the needle: prime-factoring
was never actually the dominant cost, the `O(k)` binomial loop itself was, at roughly 25:1
even before any changes. That loop's every step is a handful of `mulmod`s against the *same*
modulus `m` for the whole `Σ_{j=0}^{k}` run, which is exactly what Montgomery multiplication
is for: converting `m`'s operands to "Montgomery form" once and multiplying via `redc`
(a couple of `u64` multiplications and a shift) instead of `mulmod`'s `u128` division for
every single multiply. `Montgomery` in `src/nthdigit.rs` implements this (valid for the odd
`m` every modulus here always is), and `s_k_mod_with_primes` runs its whole accumulation in
Montgomery form, decoding only once at the very end.

One thing that didn't work on the first attempt, left in as a cautionary comment in the code:
re-encoding every prime's running `R_i` into Montgomery form on *every* `j`, regardless of
whether `R_i` actually changed that step, made the "optimisation" measurably *slower* than
plain `mulmod` (extra `redc` calls cost more than the division they were meant to avoid).
Caching each `R_i`'s Montgomery form and only refreshing it when `R_i` itself changes (rare —
most `j` share no factor with most `p_i`) turned it into a real, measured win instead.

## Profiling: where the time actually goes

Instrumented with a throwaway `nthdigit-profile` Cargo feature (`cargo build --features
nthdigit-profile`; adds two `AtomicU64` timers around `factor_segment` and the per-segment
binomial loop, printed by `pihunt digit`'s stderr line — not built by default, zero cost and
zero clippy warnings either way). At n=10⁵ and n=10⁶, with the *original* trial-division
factoring:

| n | factoring (summed across threads) | binomial loop (summed across threads) | loop:factor ratio |
|---|---|---|---|
| 10⁵ | 330.7 ms | 8413.7 ms | ~25× |
| 10⁶ | 32.1 s | 864.7 s | ~27× |

This is the headline finding of this round of work: **trial division was never the dominant
cost** — the `O(k)` binomial loop's `mulmod`-heavy inner steps were, by roughly 25-27× at
both scales tested, essentially flat across a 10× range in `n`. (Both stages are `O(N²)`
overall — trial division's per-step cost is one cheap division-based check, the binomial
loop's is a handful of `u128`-division `mulmod`s per tracked prime factor — so the constant
factor difference, not the asymptotics, explains the gap.) With the sieve in place, factoring
drops to a few milliseconds regardless of `n` (dominated by the one-time `primes_up_to` sieve,
not the O(1)-per-prime segment scans), so **the sieve alone is only worth its ~4% share of
the original total** — the real win had to come from the loop itself, which is why Montgomery
multiplication (above) is the change that actually moves the wall-clock number.

This directly contradicts this task's starting assumption that trial division was reported as
the dominant cost. It's possible that framing came from a coarser measurement (e.g. attributing
the whole per-`k` "factor, then loop" block to "factoring"); either way, the numbers above are
from this session's own instrumentation and are what's reflected in the timings table below.

## Memory bound

Every term of both sums is computed independently from `n`, `M`, `N`, `k` alone, using only
fixed-size `u64`/`u128` scratch. With the segmented sieve (the default —
`USE_SIEVED_FACTORING = true`), that's no longer quite the strict `O(log² n)` Gourdon's
theorem promises: `primes_up_to` holds every prime `≤ √(2MN)`, i.e. `O(π(√(2MN)))` extra
`u64`s (roughly `O(√n) / log n`), plus small `O(segment size)` per-thread scratch during
`factor_segment`. Concretely: a few hundred primes (~2 KiB) at n=10⁵, ~2,300 primes (~18 KiB)
at n=10⁶, an estimated ~13,000 primes (~104 KiB) at n=10⁷, and (extrapolating, not measured)
low hundreds of thousands of primes (a few MiB) by n=10⁹ — small in absolute terms, but a
real, honest departure from the strict bound, and it's the one place in this module whose
memory scales with `n` at all. The strict-memory trial-division path
(`USE_SIEVED_FACTORING = false`) is kept exactly for the case where that departure matters;
see above. Measured peak RSS below stays flat at ~5-6 MiB from n = 10⁴ through n = 10⁶ either
way — this extra allocation is still dwarfed by MPFR/GMP's own baseline footprint at every
`n` tested.

## Fixed-point accumulation and its error budget

Each term `(x mod m)/m` is stored as a `u128` fixed-point fraction `f/2^128`, computed by
exact 128-by-64 long division (`floor(x·2^128/m)`), so its representation error is `< 2^-128`
(one ulp). Terms accumulate via `wrapping_add`/`wrapping_sub`, which is exact arithmetic
modulo 1 on this encoding (we only ever want `frac(B − C)`), so accumulation order doesn't
matter and introduces no additional error. Total error after `(M+1)N + N` terms is
`< ((M+1)N + N)·2^-128`. That is tiny, but not negligible once `n0` gets large: by n ≈ 10⁹
there are ~2^45 terms (~10^-25), so it caps how many digits one evaluation can certify.
The certification bound (`error_units`) therefore adds one ulp per term to Gourdon's
truncation bound `π/(2eM)^N < 10^-(n+n0)` (see `tests/nthdigit.rs::params_are_sane`), `n0`
is capped at 24, and requests longer than 16 digits are split into independent chunks.
(An earlier version omitted the rounding term; that was only unsound for n0 ≳ 28, i.e. for
long `--count` requests or a pathological boundary retry, and is fixed.)

## Digit-boundary safety

`digits(n, count)` never returns a digit it can't certify. It picks `n0 = count + guard`
(guard starts at 4), computes `x = frac_10n_pi(n, n0)`, and re-extracts the requested
`count` digits from `x`, `x + err` and `x − err` (`err` = the `10^-n0` truncation bound plus
one ulp per term, in fixed-point units). If all three agree, the digits are certified and returned; if not (the true value
sits too close to a run of `9`s or `0`s for this guard to resolve), the guard doubles and
the whole computation retries. For `n` below 2000, or too small relative to `n0` for
Gourdon's method to apply (`n < 4·n0`), it falls back to computing π directly with MPFR
(`rug::Float`, `Constant::Pi`), which is cheap there.

**What's uncertain / not battle-tested:** the boundary-safety retry loop has never actually
been *exercised* in these runs — every real position tested landed comfortably clear of a
rounding boundary on the first try (guard = 4), including the Feynman point's run of six
9s (which only affects digits *at* that run, not a boundary the computation itself sits on).
So the retry path is covered by reasoning about the error bounds, not by an observed retry.
It would be worth constructing a deliberately adversarial `n` (a position where the true
digits are `...49999999...` or `...50000000...` right at the `count`-digit cut) to confirm
the doubling loop actually engages and still returns the right answer.

## Timings

AMD Ryzen 7 5700G, 6 cores (shared with other work during this session — see the contention
note below), `cargo build --release`. `pihunt digit <n> --count 10`, wall time from the CLI's
own `Instant` (stderr), peak RSS from `/proc/self/status` `VmHWM`. "Old" = trial division +
plain `mulmod` (this doc's previous numbers, and `USE_SIEVED_FACTORING = false` with the
Montgomery hot loop reverted); "new" = segmented sieve + Montgomery multiplication (this
session's changes, both back-to-back on the same binary build environment moments apart).
Gourdon's own `pidec` figures (Pentium III 900 MHz, single core) alongside for reference —
not an apples-to-apples comparison (different CPU, ~25 years apart, different
implementation), but it's the only published baseline for this algorithm.

| position n | old time | new time | speedup | new peak RSS | Gourdon `pidec` (P3 900MHz) |
|---|---|---|---|---|---|
| 10³ | 0 ms *(MPFR fallback)* | 0 ms *(MPFR fallback)* | — | ~5 MiB | — |
| 10⁴ | ~27-33 ms | ~25-29 ms | ~1.0× (noise-dominated at this size) | ~5 MiB | 3.13 s |
| 10⁵ | ~1.50-1.76 s | ~1.15-1.18 s | **~1.3-1.5×** | ~5.2 MiB | 185.1 s |
| 10⁶ | ~100.7-101.4 s | ~92.0-93.9 s | **~1.08-1.10×** | ~5.2-5.9 MiB | 15,869 s (4h 24m) |
| 10⁷ | 7107 s (1h 58m), old measurement — see note | not re-run this session — see note | est. ~1.1× | — | — (not in Gourdon's table; nearest is 4×10⁶ → 168,191 s) |

**Be honest about the size of this win:** it's a modest ~1.1-1.5×, not the large speedup the
spec's framing (trial division as "the dominant cost") anticipated. The "Profiling" section
above explains why: trial division was already only ~4% of the total, so replacing it with a
sieve can save at most ~4%; the real cost was always the `O(k)` binomial loop's `mulmod`
calls, and Montgomery multiplication (which *is* worth having, per that same profiling) is
what accounts for most of the observed improvement. A ~1.3-1.5× win at n=10⁵ shrinking to
~1.1× by n=10⁶ is consistent with Montgomery's fixed per-multiplication saving being
increasingly dominated by contention noise and other per-`k` overhead (memory traffic,
`strip_power` branches) as `k` grows — not a sign anything regresses; both configurations
scale the same `O(N²)`-ish way.

At n = 10³, `digits()` intentionally takes the MPFR-fallback path (`n < SMALL_N_THRESHOLD =
2000`): Gourdon's method has no advantage that small — his own table starts at 5000 — so
there's no real "Gourdon-algorithm" timing to report at 10³ here (the fallback is instant).
10⁴, 10⁵ and 10⁶ all go through the real Algorithm 1/2 path.

**Shared-machine contention, stated plainly:** this machine ran another agent's benchmark
(a different n-th-digit method, into a separate `target-thm2rs` build, unrelated to this
work) concurrently for part of this session, at one point pinning ~4.8 of the 6 cores and
visibly distorting single-run timings (a sieved+Montgomery n=10⁵ run read 2.85 s during that
window vs. a clean ~1.15-1.18 s once it finished). The ranges above are from runs taken
either before that job started or after it exited (confirmed via `ps`), and the n=10⁶
old-vs-new pair specifically was measured back-to-back (old immediately after new, same
load conditions) to cancel out any remaining drift; even so, background load on this box
fluctuated across the session (`uptime` load average seen anywhere from ~2 to ~8 on 6 cores)
and every number here should be read as "this range, on a noisy shared box," not a
single precise figure.

**10⁷ note:** the spec asked for this row only "if it now runs in <~40 min". The old
algorithm took ~1h 58m at 10⁷; given the ~1.1× measured speedup at the *nearby* n=10⁶ scale
(the regime where trial division's share of the total is largest, so where the sieve+
Montgomery combination has the most to gain), a realistic estimate for the new algorithm is
still roughly ~1h 45m — nowhere near the 40-minute budget, and re-running the old algorithm's
~2h number as a sanity check wasn't worth the wall-clock time on top of everything else this
session already measured back-to-back. The 10⁷ digits from the *old* run (`7259151336`) are
already verified against MPFR (see below); nothing about the sieve or Montgomery correctness
depends on `n`, and both are covered by dedicated unit tests plus the same MPFR checks
at 10⁵/10⁶, so a fresh 10⁷ run would only add a timing data point, not a correctness one.

The flat RSS confirms the memory claim above stays practically true through 10⁶ even with
the sieve's small, honestly-non-zero, n-dependent addition: it's still dwarfed by MPFR/GMP's
own baseline footprint and rayon's thread pool at every `n` tested so far.

## MPFR verification

- n = 10⁵: `cargo test --release --test nthdigit -- --ignored digit_at_1e5_matches_mpfr` — **passes**, digits `6412600243` — unchanged after the sieve/Montgomery rewrite, re-run and reconfirmed this session.
- n = 10⁶: (then) `cargo test --release --test nthdigit -- --ignored digit_at_1e6_matches_mpfr` — **passed**; that test has since been capped at 2·10⁵ (`digit_at_2e5_matches_mpfr`, see "Test sizes and memory" below), the verified 10⁶ digits are recorded here — digits `1309275628` — unchanged, re-run and reconfirmed this session (both 10⁵ and 10⁶ ignored tests together: 75.17 s, dominated by MPFR's own reference computation, not the algorithm under test).
- n = 10⁷: digits `7259151336` — **verified** (with the *old* algorithm; see the 10⁷ timings note above for why this wasn't re-run) against MPFR (gmpy2 `const_pi` at 3.3·10⁷ bits, 10 s), positions 10⁷…10⁷+9.
- n = 1, 762 (Feynman point), and everywhere in `[0, 20000)` (200 sequential + 200 random positions): checked in `digits_match_mpfr_reference`, part of the default `cargo test` run — passes with the new code, and noticeably faster than before (this test alone dropped from ~15-32 s to ~5 s wall time across the whole `cargo test --release --test nthdigit` run, consistent with the timings above once you factor in that most of its 400 evaluations are at small `n` where the win is smaller).
- `digit 1 --count 5` → `14159` and `digit 762 --count 8` → `99999983` were additionally cross-checked against an independent from-scratch Python (`decimal`/Machin and `decimal`/Chudnovsky) π computation before any automated test was written, and n = 2000/10000 against a from-scratch Chudnovsky reference — see the commit message for `src/nthdigit.rs`.

## Theorem 2 (Rust)

**Status: implemented, benchmarked, verified against MPFR through 10⁷.** `src/nthdigit2.rs`
ports the reconstruction in `docs/nthdigit-theorem2.md` (the chunked accumulating remainder
tree / binary-splitting algorithm, `research/thm2/thm2.py`'s Python+gmpy2 prototype) to
production Rust with `rug::Integer` + `rayon`. `pihunt digit <n> --method thm2 --mem <bits>`
(or `--method thm1`, the default — see `src/main.rs`'s `Method` doc comment for why thm1 stays
the default despite losing on speed).

All timings below: same machine as the Theorem-1 table above (Ryzen 7 5700G, 6 cores),
`cargo build --release`, wall time and peak RSS (`VmHWM`) from the CLI's own instrumentation.
Theorem 1's 10⁴/10⁵/10⁶ rows were **re-run back-to-back** with Theorem 2 in this session for a
fair comparison (the machine was shared with another agent's benchmarking run throughout, per
this task's brief — both algorithms felt the same contention). The 10⁷ Theorem-1 comparison
uses the existing measurement from the table above (a fresh 10⁷ Theorem-1 run takes ~2 hours,
out of budget for this session); everything else here is freshly measured.

### Headline: default memory (`mem_bits ≈ 4·√n·log₂10`, i.e. ≈ 4√n decimal digits — doc §6.1's `m ∝ √n` case)

| n | mem_bits | Thm2 time | Thm2 peak RSS | Thm1 time | Thm1 peak RSS | speedup |
|---:|---:|---:|---:|---:|---:|---:|
| 10⁴ | 1 329 | 9 ms | 6.8 MiB | 34 ms | 5.2 MiB | 3.8× |
| 10⁵ | 4 202 | 170 ms | 21.5 MiB | 1.51 s | 5.1 MiB | 8.9× |
| 10⁶ | 13 288 | 3.56 s | 139.6 MiB | 113.4 s | 5.0 MiB | 31.9× |
| 10⁷ | 42 020 | 131.5 s (2m 11s) | 812.1 MiB | 7107 s (1h 58m, from the table above; not re-run) | 4.2 MiB | ~54× |

Digits: 10⁴ → `8566722796`, 10⁵ → `6412600243`, 10⁶ → `1309275628`, 10⁷ → `7259151336` — all
**identical to Theorem 1's output at the same position** (see `tests/nthdigit2.rs` for this
checked automatically at many n/mem_bits combinations) and all **independently verified against
MPFR** (see below).

The speedup grows with `n`, as the doc predicts (`Thm2/Thm1 ∝ 1/(mem_bits · polylog)` roughly,
and `mem_bits` itself grows with `n` in this "default" row): 3.8× at 10⁴ up to ~54× at 10⁷,
close to the reconstruction doc's own Python-prototype-vs-C-Theorem-1 ratios in §6.1 (0.4× at
10⁴ rising to 9.6× at 1.28·10⁶) but *larger* here, because a native Rust ART leaf costs far less
than a Python one relative to Theorem 1's now-also-native inner loop — exactly the "should be
much faster per leaf" the doc's §7 caveats anticipated.

### Other `mem_bits` values (doc §6.2's fixed-n, varying-m experiment)

| n | mem_bits | Thm2 time | Thm2 peak RSS |
|---:|---:|---:|---:|
| 10⁴ | 256 (tiny) | 14 ms | 6.3 MiB |
| 10⁴ | 1 329 (default) | 9 ms | 6.8 MiB |
| 10⁴ | 8 192 (large) | 16 ms | 8.0 MiB |
| 10⁵ | 1 024 (small) | 346 ms | 18.9 MiB |
| 10⁵ | 4 202 (default) | 170 ms | 21.5 MiB |
| 10⁵ | 16 384 (large) | 149 ms | 26.3 MiB |
| 10⁶ | 4 096 (small) | 8.46 s | 116.0 MiB |
| 10⁶ | 13 288 (default) | 3.56 s | 139.6 MiB |
| 10⁶ | 65 536 (large) | 2.78 s | 181.4 MiB |

Same qualitative shape as the prototype's §6.2 table: going from a small to a default `mem_bits`
helps a lot (ART cost drops close to the predicted `1/m`), but pushing well past `4√n` gives
diminishing returns (the ART stops being the bottleneck; the `O(N)`-ish bookkeeping this
implementation doesn't stream — see below — and the p-adic/Lucas passes start to dominate, just
as the doc says). At 10⁴ mem_bits barely matters at all: the whole computation is small enough
that fixed overheads (thread pool spin-up, the `O(N)` factor sieve) swamp the ART's own cost in
either direction.

### Scaling exponent vs the n^1.5·polylog prediction

Least-squares log-log slope across the four "headline" n values above (10⁴, 10⁵, 10⁶, 10⁷,
`mem_bits ∝ √n`): **Theorem 2 ≈ 1.38**, matching the reconstruction doc's own measured range
(1.26–1.43 across its various counters, §6.1) and consistent with the predicted local slope of
`n^1.5/log²(n/m)` (≈1.4-1.5 over this range, since the `log²(n/m)` denominator grows slowly and
eats a bit of the naive 1.5 exponent). Theorem 1 over the same three re-run points (10⁴–10⁶)
comes out at **≈1.76** here, close to its historical `n²·loglog n/log²n` slope (~1.8, both in
this doc's own earlier table and the reconstruction doc's C-baseline measurement). Four points
is a thin fit — this is a sanity check that the scaling is in the right ballpark, not a precise
exponent measurement.

### MPFR verification

- n = 10⁵, 10⁶: `cargo test --release --test nthdigit2 -- --ignored` — **passes**, `6412600243`
  / `1309275628`, checked against `rug`/MPFR the same way as the Theorem-1 table above.
- n = 10⁶, 10⁷: independently cross-checked with **gmpy2** (not `rug`/MPFR — a different library
  binding, in `research/thm2`'s own venv) via `gmpy2.const_pi()` at precision `⌈(pos+30)·log₂10⌉
  + 16` bits and slicing its `gmpy2.digits(pi, 10)` mantissa string at `[pos:pos+10]`: both match
  exactly (`1309275628`, `7259151336`), 0.7 s and 11.8 s respectively to compute the reference.
- n = 10⁴, and every position `tests/nthdigit2.rs::digits_match_mpfr_reference_across_mem_bits`
  and `digits_match_theorem1_across_positions_and_mem_bits` cover (positions 0..2000 by steps of
  37, the Feynman point, ~100 random positions up to 20000, and ~45 random positions up to
  ~42000 across three `mem_bits` regimes) — part of the default `cargo test` run, all green.

### Memory: what's actually `O(mem_bits)` here (read this before trusting the RSS numbers)

The peak-RSS numbers above are real measurements, and they *do* show the ART's own working set
scaling with `mem_bits` rather than `n` (10⁶ at mem_bits=65536 uses *more* memory than mem_bits
=4096, correctly). But **this implementation's total peak RSS is not `O(mem_bits)`** the way
Theorem 1's is `O(log² n)` — it grows with `n` too (6.8 MiB at 10⁴ up to 812 MiB at 10⁷), because
the `m_k` factorisation table and the sorted ART item list are held in full (`O(N)` words)
before any chunk runs, matching a caveat the reconstruction doc states about its own prototype
(§7: "bookkeeping memory in the prototype is `O(N)`, not `O(m)`... the prototype does not [stream
it]"). This port carries the same caveat forward rather than fixing it — see the module docs in
`src/nthdigit2.rs` for exactly what's `O(mem_bits)` (the product tree + recurrence state per ART
chunk, and the p-adic tables, built and dropped one prime at a time) versus what's `O(N)` (the
factor table and item list). Fixing this for real means factoring `m_k` chunk-by-chunk with a
segmented sieve restricted to each chunk's numeric window instead of sieving `[0, N)` up front
(doc §4.4's last paragraph) — routine, but out of scope for this session.

Practically: Theorem 2 is unambiguously the faster algorithm from `n ≈ 10⁴` upward on this
machine, by a growing margin, and its digits check out against both Theorem 1 and two
independent MPFR bindings through `n = 10⁷`. Its memory story is *not* yet the `O(mem_bits)`
headline the theorem promises — plan for `O(N)`-ish RSS (hundreds of MiB by `n = 10⁷`) until the
streaming item-generation described in doc §4.4 gets implemented.

## New tests added for this round

- `sieved_factorisation_matches_trial_division` (`src/nthdigit.rs`): `factor_segment`'s output
  against `primes_le_k_dividing` (used as ground truth via an unbounded trial division) for
  every `k` across four `(M, N)` pairs and three segment sizes (including ones that don't
  evenly divide `N`, to exercise the last-short-segment case).
- `sieved_factorisation_handles_repeated_factor_and_large_cofactor`: a hand-picked
  `m = 3² · 5 · 9973` (a repeated small factor plus a prime cofactor well above `√m`, so it's
  never in `small_primes`) — the case the task specifically asked to cover.
- `sieved_algorithm2_matches_trial_division_algorithm2`: end-to-end, `sum_binomials_mod`
  (trial division) vs `sum_binomials_mod_sieved` (sieve) for every `k < N` across three
  `(M, N)` pairs, including the `k > N/2` complement-trick branch.
- `montgomery_mul_matches_mulmod`: `Montgomery`'s `encode`/`decode`/`mul` against plain
  `mulmod`, fixed edge cases plus 5,000 random odd moduli and operands.

All four are part of the default `cargo test --release` run (no `--ignored` needed).

## Test sizes and memory

The suite has to run safely on an 8 GiB machine, including with `--include-ignored`.
Measured peak RSS per test (release build, run one at a time; 2026-09-23):

| test | peak RSS | time |
|---|---:|---:|
| every non-ignored test | ≤ 33 MiB | ≤ 6 s |
| `nthdigit2::digit_at_1e6_matches_mpfr` (ignored) | 70 MiB | 3 s |
| `nthdigit::digit_at_1e6_matches_mpfr` (ignored, **now capped at 2·10⁵**) | 11 MiB | 74 s |
| `nthdigit2::digit_at_1e7_matches_mpfr` (ignored, **removed**) | 447 MiB | 110 s |
| `precision_rule::measure_precision_rule` (ignored) | 79 MiB | 265 s |

No single test comes near the limit. The OOM came from running heavy work concurrently:
the ignored large-position checks ran in parallel within one test binary, each using every
core. They're now serialized by a lock. The 10⁷ check is gone from the suite; its digits
(`7259151336`) are MPFR-verified above. To re-check it by hand:
`pihunt digit 10000001 --method thm2` (~2 min, ~450 MiB).

