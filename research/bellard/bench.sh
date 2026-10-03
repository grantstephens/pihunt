#!/usr/bin/env bash
# Benchmark Bellard's 1997 pi.c (O(n^2), O(log n) memory) against `pihunt digit`.
# usage: research/bellard/bench.sh [n ...]      (default: 1000 10000 20000)
# pi.c's position n is the same as `pihunt digit n` (1-based; both print 9-10 digits from n).
set -euo pipefail
here=$(cd "$(dirname "$0")" && pwd)
root=$(cd "$here/../.." && pwd)
exe=${TMPDIR:-/tmp}/pi_bellard
gcc -O3 -march=native -o "$exe" "$here/pi.c" -lm
pihunt=$root/target/release/pihunt
[ -x "$pihunt" ] || echo "note: $pihunt not built; skipping pihunt column" >&2
now() { date +%s.%N; }
printf '%-8s %-12s %-12s %-12s %s\n' n bellard_s pihunt_s bellard_dig pihunt_dig
for n in "${@:-1000 10000 20000}"; do
  for n in $n; do
    t0=$(now); b=$("$exe" "$n" | awk '{print $NF}'); t1=$(now)
    p=-; tp=-
    if [ -x "$pihunt" ]; then
      t2=$(now); p=$("$pihunt" digit "$n" --count 9 2>/dev/null); t3=$(now)
      tp=$(echo "$t3 - $t2" | bc)
    fi
    printf '%-8s %-12s %-12s %-12s %s\n' "$n" "$(echo "$t1 - $t0" | bc)" "$tp" "$b" "$p"
  done
done
