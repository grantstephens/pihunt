#!/usr/bin/env bash
# Generates site/data/pi-200k.txt: the first 200,000 decimal digits of pi, cross-checked against
# an independent MPFR computation. See crates/pi-digits/examples/gen_pi.rs for the method and why
# this doesn't just shell out to the `pihunt digit` CLI.
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$repo_root"

cargo run --release -p pi-digits --example gen_pi

out="$repo_root/site/data/pi-200k.txt"
size=$(wc -c < "$out")
if [ "$size" -ne 200000 ]; then
    echo "error: $out is $size bytes, expected exactly 200000" >&2
    exit 1
fi

echo "head: $(head -c 5 "$out")"
echo "Feynman point (bytes 762-767): $(tail -c +762 "$out" | head -c 6)"
