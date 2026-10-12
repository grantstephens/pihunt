# N-th decimal digit extraction

## Position convention (read this first)

- **CLI:** `pihunt digit <pos>` is 1-based: `pihunt digit 1 --count 5` → `14159` (3.**1**4159…).
- **Library:** `nthdigit::digits(n, count)` / `nthdigit2::digits(n, count, mem)` return positions
  `n+1 ..= n+count`, so `pihunt digit <pos>` calls them with `n = pos − 1`.
- Every digit string quoted in this doc is a **CLI** position unless it says otherwise:

| `pihunt digit` position | digits | position + 1 |
|---:|---|---|
| 10⁴ | `8566722796` | `5667227966` |
| 10⁵ | `6412600243` | `4126002437` |
| 10⁶ | `1309275628` | `3092756283` |
| 3·10⁶ | `3697067915` | `6970679155` |
| 10⁷ | `7259151336` | `2591513361` |

(10⁴–3·10⁶ regenerated with the CLI on 2026-09-23, Theorem 1 and 2 agreeing; 10⁷ both columns
MPFR-verified.)

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

**Update (2026-10-08): the *other* hot loop never got this treatment, and by Theorem 2's era it
was the dominant one.** `b_sum` — the `B` part of the `B - C` split, shared verbatim by Theorem 1
and Theorem 2 — computed `10^n mod m` by plain binary exponentiation through `powmod`/`mulmod`'s
`u128 % u128` the whole time. That was never the bottleneck *for Theorem 1*, where `C`
(`s_k_mod_with_primes`, Montgomery'd above) dominates at ~25:1. Theorem 2 replaces `C` with the
much faster chunked remainder tree (`src/nthdigit2.rs`), which left `b_sum`'s un-Montgomery'd
`powmod` as the new dominant cost without anyone having profiled `--method thm2` specifically to
notice — `perf record` on a 1-thread `n=3e6` run found **49.6% of all CPU time in `__umodti3`**
(the compiler's software 128-bit division routine — x86-64 has no hardware instruction for a full
128÷128 divide) plus another chunk in `b_sum`'s own closure calling it, roughly 70% combined.

Fix: `powmod_mont` (`src/nthdigit.rs`), the same binary-exponentiation algorithm as `powmod` but
through `Montgomery::mul`/`redc` instead of `mulmod`, tested against `powmod` for agreement
across random `(m, base, exp)` including `m=1`. `b_sum` builds a fresh `Montgomery` per `k`
(`m = 2k+1` changes every iteration, so there's no running state to reuse the way
`s_k_mod_with_primes` does) and calls `powmod_mont` in place of `powmod`; the one remaining
`mulmod` per `k` (`4 % m` times the result) was left as plain `mulmod` — one division per `k`
instead of up to ~46 (one per `mulmod` call inside the old `powmod`'s binary exponentiation) is
already most of the win, and Montgomery-ing a single multiply wasn't worth the extra
encode/decode round trip. Verified unchanged output (`digit 3000000`'s digits identical,
`cargo test`'s full suite plus the slow MPFR cross-checks at 10⁵/10⁶/2×10⁵ all pass) before
trusting the speedup.

Measured at `n = 3×10⁶`, same machine, back-to-back before/after on the same build: 1 thread
**103.3 s → 81.7-86.1 s** (≈ 20-26%), default 6 threads **44.1 s → 18.0 s** (≈ 2.4×, bigger than
the 1-thread gain since the division work this removed didn't parallelise any better than the
multiplication work replacing it, so cutting its *total* cost matters more once there are more
threads contending for the same memory bandwidth). Re-profiling confirmed `__umodti3` is gone
from the hot-symbol list entirely; the new top cost is `Montgomery::redc` itself (real
multiplication work, not division), plus GMP's own `__gmpn_mul_basecase_zen` (the ART's bignum
work, unchanged, just a bigger *share* of a smaller total). Benefits Theorem 1 too, since
`b_sum` is shared code — not separately re-measured here, but the mechanism doesn't care which
theorem called it.

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

- position 10⁵: `cargo test --release --test nthdigit -- --ignored digit_at_1e5_matches_mpfr` — **passes** (the test checks library offset n = 10⁵, i.e. `4126002437`; the CLI string at position 10⁵ is `6412600243`) — unchanged after the sieve/Montgomery rewrite, re-run and reconfirmed this session.
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

**2026-09-23 update — a memory-composition profiling pass** (`--features mem-profile`, a
counting global allocator plus an explicit per-structure breakdown; see `src/mem_profile.rs` and
`src/nthdigit2.rs`'s "Memory" module docs for the full writeup) found that the streamed
implementation's `O(N)` bookkeeping wasn't actually `cofactor` (which the streaming commit's docs
blamed) but two other things: `PadicBinom`'s query memo tables leaking across an entire prime's
worth of queries instead of one query's recursion tree (the dominant term by far — ~40 of ~58 MiB
at `n=1e6`), and the small-prime Lucas/p-adic need lists genuinely costing `O(N loglog sqrt(max
m_k))` words (not the `O(pi(sqrt(max m_k)))` the docs claimed) because they stored every matching
`(k,t)` pair instead of the compact arithmetic-progression range each one actually is. Both fixed
(`SideRanges` + `PadicBinom::clear_query_memo`) — see below for the three-stage before/after/now
numbers this leaves `cofactor` as the (now much smaller) dominant named term.

### Headline: default memory (`mem_bits ≈ 4·√n·log₂10`, i.e. ≈ 4√n decimal digits — doc §6.1's `m ∝ √n` case)

Thm1 columns are unchanged from before this session (Theorem 1 wasn't touched). Thm2 columns are
freshly measured after the 2026-09-23 memory fix; the "before" sub-columns show the same `n`'s
peak RSS at each of the two prior stages, where known, for context.

| n | mem_bits | Thm2 time | Thm2 peak RSS (now) | Thm2 peak RSS (streamed, pre-fix) | Thm2 peak RSS (pre-streaming) | Thm1 time | Thm1 peak RSS | speedup |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10⁴ | 1 329 | 16 ms | 5.7 MiB | — | 6.8 MiB | 34 ms | 5.2 MiB | 2.1× |
| 10⁵ | 4 202 | 217 ms | 7.3 MiB | — | 21.5 MiB | 1.51 s | 5.1 MiB | 7.0× |
| 10⁶ | 13 288 | 4.17 s | 14.1 MiB | 67.4 MiB (measured this session) | 139.6 MiB | 113.4 s | 5.0 MiB | 27.2× |
| 3·10⁶ | 23 015 | 19.3 s | 22.7 MiB | 62.2 MiB (measured this session) | — (not previously benchmarked) | not re-run | — | — |
| 10⁷ | 42 020 | 111.0 s (1m 51s) | 69.2 MiB | 447 MiB (task brief's figure) | 812.1 MiB | 7107 s (1h 58m, from the table above; not re-run) | 4.2 MiB | ~64× |

"Pre-streaming" is this table's old numbers (predates the chunked-item-generation streaming
commits this session started from); "streamed, pre-fix" is what the streaming commits left
behind (447 MiB at 10⁷ per this session's task brief; the 10⁴/10⁵/3·10⁶ streamed-pre-fix figures
were never separately benchmarked, hence the dashes — only 10⁶/10⁷ were measured directly this
session before applying the fix). The 10⁴/10⁵ Thm1 speedup numbers dropped slightly from the
previous table's 3.8×/8.9× — not a regression, just this session's fresh timing run landing on a
noisier point on a shared machine (see the RSS columns: Thm2 itself got *smaller*, not slower, at
every `n`); the headline scaling trend (speedup growing with `n`) is unchanged.

Digits (see "Position convention" at the top of this doc — CLI positions, i.e. `pihunt digit
10^k` = positions 10^k … 10^k+9): 10⁴ → `8566722796`, 10⁵ → `6412600243`, 10⁶ → `1309275628`,
3·10⁶ → `3697067915`, 10⁷ → `7259151336` — all **identical to Theorem 1's output at the same
position** (see `tests/nthdigit2.rs` for this checked automatically at many n/mem_bits
combinations) and all **independently verified against MPFR** (see below).

**On an earlier "off-by-one correction":** one session briefly relabelled the 10⁷ row as
`2591513361`. That value is correct too, but it is `pihunt digit 10000001` (library call
`digits(10_000_000, 10)`), i.e. one position later than every other row. Nothing was ever wrong
in the code; only labels mixed the CLI (1-based position) and library (offset `n`, digits start at
`n+1`) conventions. Both 10⁷ strings are MPFR-verified:
`digit 10000000` → `7259151336`, `digit 10000001` → `2591513361`.

The speedup grows with `n`, as the doc predicts (`Thm2/Thm1 ∝ 1/(mem_bits · polylog)` roughly,
and `mem_bits` itself grows with `n` in this "default" row): 3.8× at 10⁴ up to ~54× at 10⁷,
close to the reconstruction doc's own Python-prototype-vs-C-Theorem-1 ratios in §6.1 (0.4× at
10⁴ rising to 9.6× at 1.28·10⁶) but *larger* here, because a native Rust ART leaf costs far less
than a Python one relative to Theorem 1's now-also-native inner loop — exactly the "should be
much faster per leaf" the doc's §7 caveats anticipated.

**2026-09-25 update — `cofactor` bounded too.** The 2026-09-23 fix above left `cofactor` (a `Vec`
that held every leftover-large-prime Lucas need for the whole run) as the dominant *named* term at
large `n` — still genuinely `O(N)`, not `O(mem_bits)` (see `docs/nthdigit-theorem2.md` §7's
"Update" and "Update 2"). `CofactorResolver` (`src/nthdigit2.rs`) now resolves these needs in
`O(mem_bits)`-sized batches as they're found, instead of collecting the whole run's worth before
resolving any — see that module's "Memory" docs for the mechanism. Measured peak RSS, same machine
as above, `--release`, before/after run back-to-back so machine load isn't a confound:

| n | mem_bits | Thm2 time (before) | Thm2 peak RSS (before) | Thm2 time (after) | Thm2 peak RSS (after) |
|---:|---:|---:|---:|---:|---:|
| 10⁶ | 13 288 | 4.58 s | 13.9 MiB | 4.74 s | 10.8 MiB |
| 3·10⁶ | 23 015 | 19.5 s | 23.4 MiB | 20.0 s | 15.2 MiB |
| 10⁷ | 42 020 | 119.0 s | 68.9 MiB | 118.3-120.5 s (2 runs) | 19.8-24.0 MiB (2 runs) |

(These "before" figures were re-measured this session rather than reused from the 2026-09-23 row
above, specifically to control for machine load: an early comparison against the *old*,
differently-timed 111.0 s/69.2 MiB figures looked like a ~9% regression, which a same-session,
same-build baseline showed was just this machine currently running ~7% slower across the board,
not a real cost of the fix.) Time is unchanged within run-to-run noise at every size; memory drops
25-70%, growing with `n` since `cofactor` was a larger fraction of the total the bigger `n` got.
10⁴/10⁵ aren't shown: `cofactor` was never a measurable fraction of their (already small) peak RSS.

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

- positions 10⁵, 10⁶: `cargo test --release --test nthdigit2 -- --ignored` — **passes** (library offsets n = 10⁵, 10⁶; CLI strings at those positions `6412600243` / `1309275628`), checked against `rug`/MPFR the same way as the Theorem-1 table above.
- n = 10⁶, 10⁷: independently cross-checked with **gmpy2** (not `rug`/MPFR — a different library
  binding, in `research/thm2`'s own venv) via `gmpy2.const_pi()` at precision `⌈(pos+30)·log₂10⌉
  + 16` bits and slicing its `gmpy2.digits(pi, 10)` mantissa string at `[pos:pos+10]`: both match
  exactly (`1309275628`, `7259151336`), 0.7 s and 11.8 s respectively to compute the reference.
- n = 10⁴, and every position `tests/nthdigit2.rs::digits_match_mpfr_reference_across_mem_bits`
  and `digits_match_theorem1_across_positions_and_mem_bits` cover (positions 0..2000 by steps of
  37, the Feynman point, ~100 random positions up to 20000, and ~45 random positions up to
  ~42000 across three `mem_bits` regimes) — part of the default `cargo test` run, all green.

### Comparison against Bellard's existing n-th-digit implementations

Every other comparison in this doc is Theorem 1 vs Theorem 2 vs each other, plus one external,
unreproducible baseline (`pidec` on a 2003 Pentium III). `research/bellard/` adds a real,
independently-authored O(n²)-time algorithm to the mix: Fabrice Bellard's 1997 `pi.c`/`pi1.c`
(<https://bellard.org/pi/pi_n2/pi_n2.html> is the article; `pi.c` is a refinement of Plouffe's
1996 any-base n-th-digit algorithm down to O(n²) time, O(log n) memory, and `pi1.c` is a ~4×
faster variant using Gosper's series for the same quantity — not a different algorithm). Both
re-measured on this doc's own machine (AMD Ryzen 7 5700G, 6 cores, `gcc -O3 -march=native`) rather
than reusing `research/bellard/results.txt`'s original cloud-VM numbers, which weren't comparable
to anything else in this table:

| n | Theorem 2 | Theorem 1 | `pi1.c` (Gosper, 1997) | `pi.c` (1997) | `pidec` (2003, different HW) |
|---:|---:|---:|---:|---:|---:|
| 10⁴ | 16 ms | 24 ms | 2.10 s | 7.87 s | 3.13 s |
| 2·10⁴ | 49 ms | 73 ms | 8.03 s | 30.15 s | — |

All four local methods agree on digits at both positions (`8566722796` at 10⁴, `8203856539` at
2·10⁴) — an independent correctness cross-check for free, since `pi.c`/`pi1.c` are a completely
separate author's implementation of a completely separate algorithm. `pi1.c`'s 1997-vs-1997
speedup over `pi.c` here (3.7-3.75×) matches `research/bellard/README.md`'s own figure (~4.4×)
closely enough to be the same effect, not a fluke.

O(n²) means there's no point extending this table rightward: `pi.c` crosses 5 minutes by n=5·10⁴
(`research/bellard/results.txt`), and Theorem 2's whole premise is being asymptotically faster
than exactly this family of algorithm — the gap only widens from here, as the rest of this doc's
tables already show against Theorem 1 (the same O(n²) class). `research/bellard/rust/` separately
ports both `pi.c` and `pi1.c` to Rust — each with a naive (line-for-line) variant and a
Montgomery/batched-inversion/threaded `fast` one — as its own case study in how much of
Bellard-vs-Theorem-2 is "1997 C" vs "O(n²) algorithm". It doesn't change the conclusion above,
since every variant stays O(n²), but the numbers (same machine, same session;
`research/bellard/README.md` has the full writeup) are a good demonstration of how much headroom
"just engineering" leaves on the table for an unchanged O(n²) algorithm:

| | naive, 1 thread | fast, 1 thread | fast, 6 threads |
|---|---:|---:|---:|
| `pi.c` (C) | 7.87 s | — | — |
| `pi.c` (Rust) | 3.07 s | 1.29 s | 0.24 s |
| `pi1.c` (C) | 2.10 s | — | — |
| `pi1.c` (Rust) | 0.91 s | 0.23 s | 0.06 s |

(n=10⁴ throughout; the two `fast` columns reflect a follow-up optimization below, so the 6-thread
column is tighter than a first pass would show, and the sub-100ms entries there carry more
run-to-run noise than the rest.) Reading it two ways: `pi1`'s *naive* Rust port alone is already
~2.3× faster than `pi.c`'s original C — a bigger language/compiler gap than `pi.c`'s own naive
port shows (~1.4×), plausibly because `pi1.c`'s heavier per-k inner loop (four factor-stripping
terms instead of two, plus a costlier modular inverse) gives whatever's driving that gap more
work to amortise over. And fully optimized, `pi1` fast on 6 threads is **~53× faster than the
original `pi1.c`** (n=2·10⁴, the more reliable scale: 8.03 s → 0.152 s), and ~5.6× faster than
`pi.c`'s own fast variant on the same cores — Gosper's algorithmic edge and the engineering
compound rather than trade off. `pi1`'s `fast` variant has one necessary exception: the single
prime `a == 2` can't use Montgomery arithmetic (it needs an odd modulus, and `av = 2^vmax` there
never is), so it falls back to the naive per-k computation just for that one prime out of
thousands — unmeasurable cost, not a correctness compromise. Every variant cross-checked against
its C program's real output (not just known reference digits) at five positions spanning four
orders of magnitude.

**A follow-up optimization pass, and one rejected.** `Mont::mul`'s final conditional subtract
(`if u >= n { u - n } else { u }`) was a plain, data-dependent, ~50/50 branch — and, measured,
*not* something LLVM was already compiling branchlessly. Replacing it with a wrapping-subtract +
sign-mask (`reduce_once` in `rust/src/main.rs`) is where the `fast`-column improvement above over
the first pass comes from: a real ~12% (`pi.c`) to ~20-25% (`pi1`, bigger because its heavier
per-k loop calls `Mont::mul` more often) speedup, confirmed with repeated A/B timing, not a
one-off. The same trick applied to the batch accumulator's equivalent branch
(`s += t; if s >= av { s -= av }` in `flush!`) was tried too, on the same reasoning — and measured
a small, *reproducible regression* (~1.5-2%) instead, so it was reverted. Not every data-dependent
branch is a free win just because it looks like the one that was.

The bigger idea on the list, SIMD across primes (four different primes in lockstep via AVX2,
vectorising the Montgomery-multiply chain that runs on every `k`), was also tried —
`--formula pi1-simd4`, a separate formula so it can't touch the `pi`/`pi1` numbers above.
Carefully verified (the core 4-lane multiply unit-tested against scalar REDC on 20,000 random
inputs; the full algorithm unit-tested per-lane against `prime_frac_fast_gosper`, not just
checked end-to-end) and then measured: it's consistently **~7-8% slower** than `pi1` fast, not
faster, across every size and thread count tried. Kept in the tree as a documented negative
result rather than quietly dropped. Root-caused afterward with `perf` (installed partway through
this investigation specifically to stop guessing): the vector accumulators were being spilled to
the stack and reloaded every single iteration rather than kept in registers — one reload alone
was ~25% of the hot function's own samples, and static disassembly confirmed it isn't register
scarcity (only 5 of 16 YMM registers were even in use). Four targeted fixes were tried — isolating
the update into its own call, shrinking the per-lane buffers, and two escalating `asm!`
register-pin hints (one-sided, then symmetric) — and all four produced the *identical* spill,
down to the exact instruction. A real fix would need the whole update rewritten as one opaque
hand-written assembly block rather than a hint around existing codegen — scoped out as a
materially larger, correctness-critical undertaking for an uncertain slice of an already-modest
gap. A third, structurally different angle — vectorising the batched-inverse extraction in
`flush_lane` instead of the main loop, which has no loop-carried dependency and a single shared
modulus (neither property the main loop's spill needed) — was also tried: measured ~29% *slower*,
traced to the cost of gathering `pa[be[i]]` for four data-dependent indices rather than the
vectorisation itself (confirmed by isolating a scalar-but-still-two-pass version, which matched
the original baseline exactly). See `research/bellard/README.md`'s "Update" sections for the full
measurements across all three angles.

Two other things linked from bellard.org/pi are **not** in this table on purpose: the "exotic"
1997 PSLQ-discovered BBP-type series for π (already captured, for a different reason, as a
PSLQ-hunt rediscovery check — see `docs/findings-bbp-hunt.md`'s "Bellard's formula as a
rediscovery check") and Gourevitch & Guillera's 2003 binomial-sum paper on proving such formulas.
Neither is an n-th-digit *algorithm* to benchmark — one's a closed-form value for π itself (useful
for binary splitting, not single-digit extraction), the other's a proof technique. The 2009
"2700 billion digits" record on the same page ships no source and was already out of scope before
this round (`research/bellard/README.md`'s own table says so).

`pi_bin.c` (n-th *binary*/hex digit, a different base entirely — not comparable to the decimal
table above) is our own implementation of bellard.org/pi's `pi_bin` article, since the page itself
ships no code for it. Re-measured here too: 0.29 s at 10⁶ hex digits, 3.30 s at 10⁷ — both
cross-checked against known hex digits (`6C65E52C` at the 10⁶th, `7AF5863E` at the 10⁷th),
consistent with `research/bellard/results.txt`'s original cloud-VM figures scaled by roughly the
same ~3.6-3.7× this machine runs faster across the board (not independently re-verified for
O(d log d) scaling beyond what `results.txt` already measured).

### Memory: what's actually `O(mem_bits)` here (read this before trusting the RSS numbers)

**This section is now current as of the 2026-09-23 memory-profiling session** (superseding the
history below it, kept for the record). Two rounds of work happened before this table's current
numbers:

1. **Pre-streaming** (812 MiB at `n=10⁷`): the original port held a full `O(N)` factorisation
   table and a full sorted ART item list before any chunk ran, matching a caveat the
   reconstruction doc states about its own Python prototype (§7).
2. **Streaming** (447 MiB at `n=10⁷`, per this session's task brief — not separately re-derived
   here, see `src/nthdigit2.rs`'s git history): `stream_needs_and_chunks` removed the global
   factor table and item list, generating Main-item chunks on the fly. Its own docs claimed the
   *only* remaining `O(N)` piece was the cofactor `Vec` (doc §7's caveat) — **this claim was
   wrong**, see below.
3. **This session's fix** (69.2 MiB at `n=10⁷`): profiling with `--features mem-profile` (a
   counting global allocator + an explicit per-structure breakdown, `src/mem_profile.rs` /
   `nthdigit2::log_needs_sizes`) at `n=1e6` and `n=3e6` found the streaming step's actual
   dominant costs were **not** the cofactor `Vec`:
   - `PadicBinom`'s per-query memo tables (`memo_s`/`memo_c`) were leaking across an entire
     prime's worth of queries instead of clearing between them — for small primes with many
     queries (e.g. `p=3` at `n=1e6`: 34 557 queries), this alone measured ~40 of the ~58 MiB
     peak. Fixed: clear both maps before each top-level query
     (`PadicBinom::clear_query_memo`) — cheap and provably correctness-neutral.
   - `lucas_small`/`padic_small` stored an explicit `Vec<(k,t[,e])>` per key: `O(N loglog
     sqrt(max m_k))` words total (Mertens' third theorem), not the `O(pi(sqrt(max m_k)))` the
     docs claimed — measured 12.9+2.4 MiB at `n=1e6`, 36.7+6.1 MiB at `n=3e6`. Fixed: `SideRanges`
     — a fixed prime's matching `k` form one arithmetic progression per mirror side, so only the
     first/last `t` needs storing; consumers regenerate `(k,t)` pairs on the fly.

   `cofactor` itself — the piece the streaming docs blamed — genuinely is `O(N)` (most of its
   primes are singleton occurrences that don't compact into a range the way a shared small prime
   does) but is *small* in absolute terms at the `n` this runs at: ~1.5 MiB at `n=1e6`, ~6 MiB at
   `n=3e6`. It's now the dominant *named* remaining term, but no longer the dominant *actual*
   one — the rest of peak RSS at these sizes is `rug::Integer`/rayon/allocator overhead this
   session didn't itemise further (see `src/nthdigit2.rs`'s module docs for the full writeup).

**Bound achieved:** `O(mem_bits · log(mem_bits) · threads + pi(sqrt(max m_k)) + chunks +
cofactor_count)`. `cofactor_count` is the one term that scales with `N` rather than `mem_bits` —
still not the `O(mem_bits)` the theorem promises in the strict asymptotic sense, but empirically
tens of MiB through `n = 10⁷` rather than hundreds, and every other term this session identified
really is bounded independent of `N`. Fixing `cofactor_count` for real would need either an
external (disk-backed) sort of the cofactor needs by target, or a smarter per-prime
classification that doesn't require discovering the cofactor's value before knowing whether it's
"small" (doc §4.4's last paragraph gestures at a segmented-sieve-over-`q` approach) — future
work, not implemented here; see `src/nthdigit2.rs`'s module docs for the detailed argument for
why it resists the same compaction the small primes got.

Practically: Theorem 2 is unambiguously the faster algorithm from `n ≈ 10⁴` upward on this
machine, by a growing margin, its digits check out against both Theorem 1 and two independent
MPFR bindings through `n = 10⁷` (see the correction note above the headline table), and as of
this session its memory footprint is genuinely small in absolute terms — tens of MiB, not
hundreds, through `n = 10⁷` — even though the strict `O(mem_bits)` bound the theorem promises
isn't fully achieved (see `cofactor_count` above). See "Should the CLI default switch to thm2?"
below for what this changes about the `--method` default.

#### Should the CLI default switch to thm2?

Not changed this session — but the case for it is now much stronger than when `thm1` was picked
as the default (see `src/main.rs`'s `Method` doc comment for the original reasoning: thm2's
"still `O(N)`-ish bookkeeping" was the stated blocker). That blocker is now mostly gone: thm2 is
faster at every `n` tested from `10⁴` up (2.1×-64× in the table above) *and* smaller in absolute
peak RSS than it was (69 MiB vs Theorem 1's 4.2 MiB at `n=10⁷` — Theorem 1 is still smaller in
absolute terms, being genuinely `O(log² n)`, but "69 MiB to be 64× faster" is a trade nearly any
caller would take, and it no longer risks hundreds of MiB to hundreds of MiB it did before this
session). Two things still argue for leaving the default alone rather than flipping it in this
same commit: (1) thm2's memory is still not the *asymptotic* `O(mem_bits)` the theorem promises
(`cofactor_count` above) — safe in practice through `n=10⁷` but not proven safe at, say, `n=10⁹`
without further work; (2) switching a default is a separate decision from a memory-bugfix
session, and deserves its own explicit call rather than riding along here. Recommendation for a
follow-up: switch the default to `thm2`, keeping `thm1` available via `--method thm1` for
anyone who wants the strictly-`O(log² n)` guarantee.

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
| `nthdigit2::digit_at_1e6_matches_mpfr` (ignored) | ~15 MiB (2026-09-23 memory fix; was 70 MiB) | 4.5 s |
| `nthdigit::digit_at_1e6_matches_mpfr` (ignored, **now capped at 2·10⁵**) | 11 MiB | 74 s |
| `nthdigit2::digit_at_1e7_matches_mpfr` (ignored, **removed**) | was 447 MiB pre-fix; ~70 MiB now (not re-added as a test — see below) | 110 s |
| `precision_rule::measure_precision_rule` (ignored) | 79 MiB | 265 s |

No single test comes near the limit. The OOM came from running heavy work concurrently:
the ignored large-position checks ran in parallel within one test binary, each using every
core. They're now serialized by a lock. The 10⁷ check is gone from the suite; its digits are
MPFR-verified above. To re-check by hand: `pihunt digit 10000000` → `7259151336` (~2 min,
~70 MiB; Theorem 2 is the default method).

