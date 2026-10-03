# Bellard `pi.c` / `pi1.c` / `pi_bin` benchmark

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
