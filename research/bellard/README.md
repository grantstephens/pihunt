# Bellard `pi.c` / `pi1.c` / `pi_bin` benchmark

`pi.c` is F. Bellard's 1997 O(n²)-time, O(log n)-memory decimal digit extractor
(<https://bellard.org/pi/pi.c>, unmodified). `bench.sh` compiles it and times it against
`target/release/pihunt digit` at the same position. Measured pi.c times are in `results.txt`
(including an original cloud-VM pass and a later re-measurement on the actual
`docs/nthdigit.md` benchmark machine, so the two sets of numbers aren't directly comparable to
each other — each is labelled).

Scaling is ~n² (10⁴ → 14 s, 2·10⁴ → 54 s, 5·10⁴ → 324 s on the original cloud VM), so keep
n ≤ ~5·10⁴. For the full writeup comparing these against Theorem 1/Theorem 2 side by side
(same machine, same session), see `docs/nthdigit.md`'s "Comparison against Bellard's existing
n-th-digit implementations".

## Rust port (`rust/`)

`cd rust && cargo build --release && target/release/bellard-pi 10000 [--formula pi|pi1|pi1-simd4] [--variant naive|fast] [--threads T]`

* `--formula pi` (default), `naive`: a line-for-line port of `pi.c` (u64 `%`, a Euclid inverse
  per k): ~1.4× faster than the C build, just from the compiler/integer division.
* `--formula pi` `fast` (default variant) keeps Bellard's algorithm and O(log n)-ish memory but
  uses Montgomery u32 multiplication, batched inversion (one Euclid per 512 terms), a sieve for
  the primes, and threads over the primes. About 4× faster single-threaded, ~12× on 4 cores
  (n=2·10⁴: 41 s → 3.1 s).
* `--formula pi1`, `naive`: a line-for-line port of `pi1.c` (Gosper's series) — four running
  terms and a special `a == 2` branch (ported as its own `inv_mod2`, not simplified away) instead
  of `pi.c`'s two terms. Measured **~2.2-2.3× faster than the C build** single-threaded (n=10⁴:
  2.10 s → 0.91 s; n=2·10⁴: 8.03 s → 3.63 s) — a bigger language/compiler gap than `pi.c`'s own
  naive port shows, consistent with `pi1.c`'s heavier per-k inner loop (four factor-stripping
  terms and a costlier modular inverse) giving the same per-iteration overhead more to amortise.
  (A first attempt at this measurement looked like Rust running *slower* than C — that was a
  stale pre-built binary sitting in this directory's own `target/` from before the `--formula`
  flag existed, silently running `pi.c`'s naive path instead and still printing correct-looking
  digits, since both formulas compute the same digits of π. Worth remembering: a surprising
  benchmark result is a reason to check what actually ran, not to report it.)
* `--formula pi1`, `fast` (default): same Montgomery/batched-inversion/threaded treatment as
  `pi.c`'s `fast`, adapted for four running terms instead of two (`bk[i]` holds `25k-3`, the
  per-term scalar Gosper's formula multiplies by, instead of `pi.c`'s bare `k`). One exception:
  `a == 2` can't go through Montgomery arithmetic at all (it needs an odd modulus; `av = 2^vmax`
  for that branch never is), so that single prime — one out of thousands up to `3*big_n` — falls
  back to `prime_frac_naive_gosper` instead; negligible cost. Figures below include the
  branchless-`Mont::mul` optimisation from the Update section further down (it affects both
  `fast` variants, since they share `Mont::mul`). ~4.3× faster than `pi1`'s own naive port
  single-threaded at n=2·10⁴ (3.63 s → 0.85 s) — more than `pi.c` fast gets over its own naive
  (~2.7×), plausibly because `pi1`'s heavier per-k loop calls `Mont::mul` more often per
  contributing term, giving it more to amortise. All the way down: **`pi1` fast, 6 cores, vs the
  original `pi1.c` C program is ~53×** at n=2·10⁴ (8.03 s → 0.152 s) — and **~5.6× faster than
  `pi.c`'s own fast variant** on the same cores (0.848 s → 0.152 s), i.e. Gosper's algorithmic
  edge and the Montgomery/threading engineering compound rather than trade off.

All four variants are cross-checked against their respective C program's actual output (not just
the known reference digits) at n = 1, 100, 1000, 10⁴, 2·10⁴.

Still O(n²); `pihunt digit` (Gourdon-style) remains far faster at these sizes.

### Update (2026-10-07, same session): one real micro-optimization found, one tried and rejected

Investigated further optimizations in both `fast` variants. No `perf` available in this
environment, so this was reasoning-from-structure plus direct A/B timing (several repeats each,
to separate signal from run-to-run noise) rather than profiler-guided.

* **`Mont::mul`'s conditional subtract, made branchless: a real, substantial win.** The Montgomery
  reduction's last step was a plain `if u >= n { u - n } else { u }` — data-dependent, roughly
  50/50, and (measured) *not* something LLVM was already turning into a branchless sequence.
  Replacing it with `reduce_once` (wrapping-subtract + sign-mask, no branch) measured, at n=2·10⁴,
  single-threaded: `pi.c` fast 5.43 s → 4.77 s (**12% faster**), `pi1` fast 1.06 s → 0.85 s
  (**19-20% faster** — bigger, since `pi1`'s heavier per-k loop calls `Mont::mul` more often per
  contributing term). Same shape on 6 cores: `pi.c` fast 0.961 s → 0.848 s (12%), `pi1` fast
  0.202 s → 0.152 s (25%). This one change benefits every `fast` run, since both formulas share
  `Mont::mul`.
* **The same trick on the batch accumulator (`s += t; if s >= av { s -= av }` in `flush!`): tried,
  measured, rejected.** Same pattern, same reasoning for why it should help — but it didn't.
  Several repeated A/B runs showed a small, reproducible *regression* (~1.5-2%) from switching
  this specific call site to branchless, not an improvement. Kept as plain branches in both
  `fast` variants; `reduce_once` is used only in `Mont::mul`, where it actually measured faster.
  (See `reduce_once`'s doc comment in `rust/src/main.rs` for the detail — not every data-dependent
  branch is a free win just because it looks like the one that was.)

The `fast`-variant bullets above already reflect these numbers. Sub-100ms multi-threaded timings
(e.g. n=10⁴ on 6 cores) carry more run-to-run noise than the rest — thread startup and OS
scheduling jitter start to matter at that scale — which is why the ratios above lean on n=2·10⁴.

Things considered and not pursued, with reasons: wider `u64` Montgomery arithmetic (would need
double-width intermediates for the multiply, no evidence the current `u32` modulus range is
actually the constraint); retuning `BLK` (the batch size already amortises the one-inversion-per-
block cost well below the per-term `Mont::mul` cost, so this isn't where the time goes). Neither
was measured, just reasoned about — flagged here as open directions, not closed ones. Unlike the
third item on the original list here, SIMD across primes, which *was* tried: see below.

### Update (2026-10-07, same session): tried the SIMD redesign. It doesn't help.

`--formula pi1-simd4` (`simd4` module in `rust/src/main.rs`) processes four *different* odd
primes in lockstep via AVX2, vectorising the num/den Montgomery-multiply chain that runs on every
`k` regardless of whether `k` contributes to the sum — the one part of `prime_frac_fast_gosper`
with no per-lane-divergent control flow, and (reasoned, not yet measured at the time) the
dominant cost. `a == 2` and any 0-3 leftover primes after grouping into fours still go through
the existing, already-validated scalar paths; this is a distinct `--formula`, not a `pi1` variant,
specifically so the `pi`/`pi1` numbers everywhere else in this README stay exactly reproducible.

**What's actually new here**: `montmul4`, a 4-lane Montgomery multiply (`a[i]*b[i]/R mod n[i]`,
each lane its own modulus — REDC has no cross-lane dependency, so this part generalises for
free) via `_mm256_mul_epu32`/`_mm256_mullo_epi32`, with the same branchless conditional-subtract
trick as `reduce_once`. Factor-stripping (`DIVN`) and the batched-inverse extraction stay scalar,
one lane at a time, reusing the per-prime bookkeeping shape the scalar `fast` variants already
use — only the uniform, always-executes part is vectorised.

**Verified carefully, given how much new `unsafe` this is**, before trusting any performance
number: `montmul4` unit-tested against the scalar REDC on 20,000 random `(a,b,n)` triples;
`frac4` (the full 4-lane algorithm) unit-tested *per lane* against `prime_frac_fast_gosper` called
individually on the same primes — stronger than an end-to-end digit check, since it can't hide
two compensating bugs the way matching final digits could; then the full CLI cross-checked
against `pi1` at six positions, both 1 and 6 threads. All green.

**Measured result: `pi1-simd4` is consistently ~7-8% *slower* than `pi1` fast**, not faster — at
n=2·10⁴ (1 thread: 0.85 s → 0.90-0.92 s; 6 cores: 0.152 s → 0.16-0.17 s) and at n=5·10⁴ (1 thread:
4.79 s → 5.17 s). Reproducible across repeated runs and both thread counts, not noise.

Two follow-up attempts, both measured and neither enough to close the gap: (1) the first version
round-tripped `num`/`den` through scalar form every single `k` even on iterations where no lane
needed it (unpack, then immediately re-pack the same values next iteration) — fixed to carry them
as vectors across the whole loop, unpacking only when a lane actually has a contribution to push;
measured ~4% faster, still behind. (2) `pack()` used `_mm256_set_epi64x`, which likely lowers to a
multi-instruction insert sequence; switched to `_mm256_cvtepu32_epi64` (a single zero-extend
instruction, the operation AVX2 actually has for "4x u32 -> 4x u64") — measured no further change,
meaning `pack` wasn't the bottleneck either.

### Update (2026-10-07, later same session): root-caused with `perf`, two targeted fixes tried, neither helped

`perf` wasn't available in the sandbox this was first built in; it is now (installed via pacman,
with the user's explicit go-ahead first, since that's a system-level change beyond the repo).
That turned "best-effort explanation" into an actual measurement.

`perf record`/`report`/`annotate` on `pi1-simd4` at n=2·10⁴: `frac4_avx2` (the vectorised k-loop)
is ~70-72% of total cycles, `flush_lane` (the scalar batched-inverse extraction, believed cheap —
wrongly) is ~28-33%. Inside `frac4_avx2`, one single instruction —
`vpmuludq 0x1d0(%rsp),%ymm0,%ymm0`, reloading `numv` from its stack slot to start the next `k`'s
multiply — accounted for **~25-28% of that function's own samples on its own**. `numv`/`denv`
were not being kept in YMM registers across the loop at all: stored back to the stack at the end
of every iteration (`vmovdqu %ymm0,0x1c0(%rsp)`) and immediately reloaded from the same address a
few instructions later to begin the next one — a textbook store-to-load-forwarding stall, paid on
every single `k`.

Two targeted fixes, both informed by this measurement, both tried, neither closed the gap:

1. **Isolate the update into its own `#[inline(never)]` call**, on the theory that competing
   with all the surrounding scalar bookkeeping (four lanes' worth of `kq1..kq4`, `v`, and the
   `ba`/`bb`/`bk`/`be` buffers) for registers was forcing the spill. Measured: *slower*
   (~0.90 s → ~0.95-1.03 s at n=2·10⁴) — the real function-call boundary (and whatever the
   calling convention does with a 2x `__m256i` return) cost more than the isolation saved.
   Reverted.
2. **Shrink the per-lane batch buffers** (a new, smaller `BLK` just for this module, 64 instead
   of 512) on the theory that a ~26 KB stack frame (four lanes' worth of those buffers) was
   itself why the allocator wouldn't commit registers to anything for the whole loop's duration.
   Measured: no real change (~0.89-0.93 s, within noise of the baseline) — and re-profiling
   confirmed the *exact same* spill instruction, same ~28% share, completely unmoved. Reverted.

### Update (2026-10-07, later still): pinpointed further with static disassembly; tried the assembly route; still unresolved

Two more things worth separating out, since "it's a register-allocation issue" was still a bit
vague: *is it actually register scarcity*, and *can a small, targeted assembly fix pin it down*?

**Static disassembly first** (`objdump -d`, no profiler needed for this part): `frac4_avx2` uses
only **5 of the 16 available YMM registers** (`ymm0`-`ymm4`) in the entire function, and both
`montmul4` calls are confirmed inlined (no `call` instruction for them anywhere in the hot loop —
the only `call`s left are `log@GLIBC` in one-time setup code, `malloc`/`grow_one` for the `pa`
table `Vec`s, and an unreachable divide-by-zero panic path). So this is **not** register
exhaustion — there were 11 completely unused YMM registers the whole time. The decision to spill
`numv`/`denv` is a genuine choice, not a resource constraint, and the spill is specifically
*loop-carried*: the reload at the top of one `k` reads back exactly what the previous `k` wrote,
confirmed by reading the addresses off by hand (`0x1c0(%rsp)` written once per iteration, read
back once per iteration, nothing else touches it in between).

**Then the assembly route** (the user's own suggestion, tried before concluding further). Two
escalating experiments, both using `std::arch::asm!` to force `numv`/`denv` into specific,
otherwise-unused registers (`ymm5`/`ymm6`) at precise program points, *without* touching any of
the real arithmetic — a bare register-pin hint, to test cheaply whether that alone would be
enough before committing to rewriting the whole REDC chain in raw asm:

1. **Pin at the write end only** (right after the 4th `mo4.mul` call, forcing the result into
   `ymm5`/`ymm6`). The hint was honored *exactly where asked* — `vblendvpd` (the REDC's final
   conditional-move step) now writes directly into `ymm5`/`ymm6` instead of wherever it landed
   before — but the very next instructions store them straight back to the *same* `0x1a0`/`0x1c0`
   stack slots anyway. The compiler respected the local hint and then immediately overrode its
   effect for the loop as a whole.
2. **Pin at both ends** (read, at the top of the update, and write, at the bottom), symmetric.
   Identical outcome: same spill instruction, same stack addresses, same ~28% cost, and timing
   unchanged (if anything, marginally worse from the extra `asm!` overhead: ~0.93 s vs ~0.90-0.92 s
   baseline).

**So: four independent, targeted interventions — function isolation, smaller buffers, a one-sided
register pin, a symmetric register pin — produced the *identical* spill, down to the exact
instruction, every time.** That consistency is itself the finding: this isn't register pressure,
isn't stack-frame size, and isn't fixable by hinting at the use sites. The loop-carried value's
storage location is being decided somewhere upstream of instruction selection — plausibly during
LLVM's loop-vectorization or SSA/PHI-node construction for this specific, heavily-branching loop
body (DIVN's nested `if`/`while`/`loop` across four lanes) — a decision local `asm!` hints at
individual program points simply don't reach.

A fix that actually works would need to remove LLVM's opportunity to make that decision at all:
the *entire* per-`k` update (not just a hint around existing codegen, the whole arithmetic chain)
written as one opaque `asm!` block, so `numv`/`denv` are never, at any point, a normal
LLVM-visible SSA value with a storage decision to make — they'd exist only as register operands
Rust hands to one black-box instruction sequence. That's a real rewrite of Montgomery REDC in
raw AVX2 assembly (proper operand/clobber declarations, scratch registers for `t`/`m`/`mn`/`u`,
and — critically — full re-validation from scratch, since a hand-written REDC that's subtly
wrong would fail silently rather than loudly). Scoped out of this round: a correctness-critical,
from-scratch assembly implementation is a materially different, larger undertaking than the
hint-based experiments above, for an uncertain amount of an already-modest (~8%) gap. Flagged
precisely, not attempted, for whoever picks this up next with more time to spend on it.

Kept in the tree as a correct, tested, but *not recommended* code path — a real answer to "does
this redesign help" (no, measured, here's why it was expected to, and here's exactly where the
gap lives down to the instruction) rather than a guess either way.

### Update (2026-10-07, later still): a third, structurally different angle — also didn't pay off

The main k-loop wasn't the only candidate: `flush_lane` (the batched-inverse extraction) is
~28-33% of total time and was completely untouched by any of the above. Unlike `numv`/`denv`,
its final extraction step (`ba[i]*ib[i]*bk[i]*pa[be[i]]` for each contributing item `i`) has no
dependency between different `i` once `ib[i]` is known — not loop-carried, same modulus for the
whole batch. Both properties the main loop's register-spill issue needed and didn't have, so this
looked like a genuinely different, more promising target, not just another attempt at the same
thing.

Implemented: split the existing single reverse pass (which computed `ib[i]` and immediately
consumed it) into two — an unchanged sequential pass computing and storing every `ib[i]`, then a
second pass processing groups of 4 independent items via `montmul4` (same modulus broadcast
across all 4 lanes, reusing the already-validated primitive directly).

**Measured: ~29% slower** (0.91 s → 1.15 s at n=2·10⁴), a clear regression, not noise.
Diagnosed by testing the two changes separately: reverting to a **scalar** extraction but keeping
the two-pass split (computing and storing `ib[]` first, same as the vectorised version, just not
vectorising the second pass) measured identical to the original single-pass baseline (~0.91 s) —
so the split itself costs nothing. The entire regression was the vectorisation. Most likely
culprit: unlike the main loop's `t_num1` etc. (already plain, contiguous `[u32; 4]` arrays coming
straight out of scalar bookkeeping), this extraction needs `pa[be[i]]` for four different,
data-dependent indices — a genuine gather, not a contiguous load — and assembling that into a
vector plus unpacking the result back out costs more than the three scalar multiplies it replaces.

Reverted cleanly (confirmed identical before/after via `git diff --stat`). Three substantially
different angles tried this session — cross-prime vectorisation of the main loop, four targeted
fixes for its register spill (including the assembly route), and same-modulus batched extraction
— and none beat the scalar `fast` variant. That's a reasonably thorough search of the approaches
that looked promising on paper; closing the remaining ~8% would need either hand-written assembly
for the whole per-`k` update (scoped out above as a materially larger, correctness-critical
undertaking) or a genuine hardware gather instruction for the `pa[be[i]]` lookup specifically,
neither attempted here.

## Other algorithms on bellard.org/pi

| file | what | status |
|---|---|---|
| `pi.c` | decimal digit, O(n²), O(log n) memory (Jan 1997) | Bellard's source, unmodified |
| `pi1.c` | same, using Gosper's series (Feb 1997) | Bellard's source, unmodified; ~4.4× faster than `pi.c` (n=5·10⁴: 74 s vs 324 s), same digits |
| `pi_bin.c` | n-th *binary* digit via the `pi_bin` article's formula | **our implementation** (the page ships no code); O(d log d), 10⁷ hex digits in 12 s |
| 2009 "2700 billion digits" | y-cruncher-class Chudnovsky run | no source published, not benchmarkable |

`pi_bin` answers a different question (base 2/16, not base 10), so `bench.sh` does not run it
against `pihunt`; time it directly with `gcc -O3 -o pi_bin pi_bin.c -lm && ./pi_bin 1000000`.
`bench.sh` now runs `pi.c` and `pi1.c` side by side and flags any digit mismatch between them.
