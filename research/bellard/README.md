# Bellard `pi.c` benchmark

`pi.c` is F. Bellard's 1997 O(n²)-time, O(log n)-memory decimal digit extractor
(<https://bellard.org/pi/pi.c>, unmodified). `bench.sh` compiles it and times it against
`target/release/pihunt digit` at the same position. Measured pi.c times are in `results.txt`.

Scaling is ~n² (10⁴ → 14 s, 2·10⁴ → 54 s, 5·10⁴ → 324 s), so keep n ≤ ~5·10⁴.
For `pihunt` timings at the same sizes see `docs/nthdigit.md` (≈25 ms at 10⁴, ≈1.15 s at 10⁵).

## Rust port (`rust/`)

`cd rust && cargo build --release && target/release/bellard-pi 10000 [--variant naive|fast] [--threads T]`

* `naive` is a line-for-line port (u64 `%`, a Euclid inverse per k): ~1.4× faster than the C
  build, just from the compiler/integer division.
* `fast` (default) keeps Bellard's algorithm and O(log n)-ish memory but uses Montgomery u32
  multiplication, batched inversion (one Euclid per 512 terms), a sieve for the primes, and
  threads over the primes. About 4× faster single-threaded, ~12× on 4 cores (n=2·10⁴: 41 s → 3.1 s).

Still O(n²); `pihunt digit` (Gourdon-style) remains far faster at these sizes.
