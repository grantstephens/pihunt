# Bellard `pi.c` benchmark

`pi.c` is F. Bellard's 1997 O(n²)-time, O(log n)-memory decimal digit extractor
(<https://bellard.org/pi/pi.c>, unmodified). `bench.sh` compiles it and times it against
`target/release/pihunt digit` at the same position. Measured pi.c times are in `results.txt`.

Scaling is ~n² (10⁴ → 14 s, 2·10⁴ → 54 s, 5·10⁴ → 324 s), so keep n ≤ ~5·10⁴.
For `pihunt` timings at the same sizes see `docs/nthdigit.md` (≈25 ms at 10⁴, ≈1.15 s at 10⁵).
