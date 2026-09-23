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

## Memory bound

Every term of both sums is computed independently from `n`, `M`, `N`, `k` alone, using only
fixed-size `u64`/`u128` scratch (a handful of primes ≤ k dividing `m`, so `O(log m) =
O(log n)` words). Nothing is ever sized by `N` or `n` — no sieve, no per-term vector, no
table of binomial coefficients. So per-thread memory is `O(log² n)` (Gourdon's theorem),
and total memory is that times the (fixed, small) thread count. Measured peak RSS below
confirms this: it's flat at ~5 MiB from n = 10³ through n = 10⁶.

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

AMD Ryzen 7 5700G, 6 cores, `cargo build --release`. `pihunt digit <n> --count 10`, wall
time from the CLI's own `Instant` (stderr), peak RSS from `/proc/self/status` `VmHWM`.
Gourdon's own `pidec` figures (Pentium III 900 MHz, single core) alongside for reference —
not an apples-to-apples comparison (different CPU, ~25 years apart, different
implementation), but it's the only published baseline for this algorithm.

| position n | our time | our peak RSS | Gourdon `pidec` (P3 900MHz) |
|---|---|---|---|
| 10³ | 0 ms *(MPFR fallback — see below)* | 4.9 MiB | — |
| 10⁴ | 36 ms | 4.9 MiB | 3.13 s |
| 10⁵ | 2.15 s | 5.0 MiB | 185.1 s |
| 10⁶ | 116.9 s (1m 57s) | 4.9 MiB | 15,869 s (4h 24m) |
| 10⁷ | 7107 s (1h 58m) — see note | 4.2 MiB | — (not in Gourdon's table; nearest is 4×10⁶ → 168,191 s) |

At n = 10³, `digits()` intentionally takes the MPFR-fallback path (`n < SMALL_N_THRESHOLD =
2000`): Gourdon's method has no advantage that small — his own table starts at 5000 — so
there's no real "Gourdon-algorithm" timing to report at 10³ here (the fallback is instant).
10⁴, 10⁵ and 10⁶ all go through the real Algorithm 1/2 path.

**10⁷ note, stated plainly:** the spec asked for this row only "if it finishes within ~30
min". It didn't — it took ~1h 58m, about 60× the n=10⁶ time (roughly in line with the
n=10⁵→10⁶ ratio of ~54×, so the observed scaling is consistent, just past the time budget).
It's included above anyway because it did finish and the number is real; its digits were
later verified against MPFR (see below). Take the 10⁷ row as "it works and the memory
stays flat," not as a benchmark result to compare against a 30-minute budget.

The very flat RSS confirms the `O(log² n)` memory claim directly: essentially all of the
~5 MiB is MPFR/GMP's own baseline footprint and rayon's thread pool, not anything that
scales with `n`.

## MPFR verification

- n = 10⁵: `cargo test --release --test nthdigit -- --ignored digit_at_1e5_matches_mpfr` — **passes**, digits `6412600243`.
- n = 10⁶: `cargo test --release --test nthdigit -- --ignored digit_at_1e6_matches_mpfr` — **passes**, digits `1309275628`, run twice (once before an unrelated session interruption, once after, both green); the MPFR reference computation itself takes ~100-115 s (computing π to ~1,000,000+ decimal digits of MPFR precision, once, to check against).
- n = 10⁷: digits `7259151336` — **verified** afterwards against MPFR (gmpy2 `const_pi` at 3.3·10⁷ bits, 10 s), positions 10⁷…10⁷+9.
- n = 1, 762 (Feynman point), and everywhere in `[0, 20000)` (200 sequential + 200 random positions): checked in `digits_match_mpfr_reference`, part of the default `cargo test` run.
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
