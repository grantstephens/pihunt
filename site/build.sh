#!/usr/bin/env bash
# Builds the pihunt digits site into site/dist: a ready-to-deploy static site (paper write-up +
# in-browser WASM demo). Run from anywhere; every path is resolved relative to this script's own
# location. See site/README.md for the one-line build and deploy commands.
#
# Env: PATH must include ~/.cargo/bin (cargo, wasm-pack). CARGO_TARGET_DIR for both the wasm
# crate and (only if site/data/pi-200k.txt is missing) the root workspace are set below to
# directories outside the repo -- this worktree lives under ~/sync, which is synced, and must
# never receive build output. Override with WASM_CARGO_TARGET_DIR / ROOT_CARGO_TARGET_DIR if
# needed.
set -euo pipefail

site_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
dist_dir="$site_dir/dist"

export PATH="$HOME/.cargo/bin:$PATH"
wasm_target_dir="${WASM_CARGO_TARGET_DIR:-$HOME/.cache/pihunt/target-wasm}"
root_target_dir="${ROOT_CARGO_TARGET_DIR:-$HOME/.cache/pihunt/target-site}"

log() { printf '==> %s\n' "$*"; }
die() { printf 'error: %s\n' "$*" >&2; exit 1; }

# --- Step 1: tooling check, with install hints --------------------------------------------------
log "checking tooling"
command -v cargo >/dev/null 2>&1 \
  || die "cargo not found on PATH. Install Rust: https://rustup.rs"
command -v wasm-pack >/dev/null 2>&1 \
  || die "wasm-pack not found on PATH. Install: cargo install wasm-pack"
command -v node >/dev/null 2>&1 \
  || die "node not found on PATH. Install Node.js >= 18 (needed for 'node --test'): https://nodejs.org"
node_major="$(node -e 'process.stdout.write(process.versions.node.split(".")[0])')"
if [ "$node_major" -lt 18 ]; then
  die "node $(node --version) is too old; need Node >= 18 for 'node --test'"
fi

# --- Step 2: clean dist ---------------------------------------------------------------------
log "cleaning $dist_dir"
rm -rf "$dist_dir"
mkdir -p "$dist_dir"

# --- Step 3: wasm-pack build (wasm-pack already runs wasm-opt for a --release build; do not
# run wasm-opt again separately) -----------------------------------------------------------------
log "building site/wasm with wasm-pack (--release, target web)"
(
  cd "$site_dir/wasm"
  CARGO_TARGET_DIR="$wasm_target_dir" wasm-pack build --target web --release \
    --out-dir ../dist/wasm --out-name pi_digits_wasm
)

# --- Step 4: generate site/data/pi-200k.txt if it's missing --------------------------------------
if [ ! -f "$site_dir/data/pi-200k.txt" ]; then
  log "site/data/pi-200k.txt missing; generating it (gen-pi.sh)"
  CARGO_TARGET_DIR="$root_target_dir" "$site_dir/scripts/gen-pi.sh"
else
  log "site/data/pi-200k.txt already present; skipping generation"
fi

# --- Step 5: wasm smoke test ---------------------------------------------------------------------
log "running wasm smoke test"
node "$site_dir/scripts/smoke-wasm.mjs"

# --- Step 6: JS unit tests (pure functions in site/src, e.g. demo.test.mjs) ----------------------
log "running node --test over site/src"
node --test "$site_dir/src"

# --- Step 7: copy src/* (excluding dev-only files), vendor/, data/ into dist --------------------
log "copying static assets into dist"
for f in "$site_dir"/src/*; do
  base="$(basename "$f")"
  case "$base" in
    demo.test.mjs|package.json) continue ;;
  esac
  cp -R "$f" "$dist_dir/"
done
cp -R "$site_dir/vendor" "$dist_dir/vendor"
cp -R "$site_dir/data" "$dist_dir/data"

# --- Step 8: splice the benchmark charts into index.html -----------------------------------------
log "rendering charts into dist/index.html"
node "$site_dir/scripts/charts.mjs" "$site_dir/src/index.html" > "$dist_dir/index.html"

# --- Step 9: static site checks -------------------------------------------------------------------
log "running check-site.mjs"
node "$site_dir/scripts/check-site.mjs" "$dist_dir"

# --- Step 10: report sizes ---------------------------------------------------------------------
log "build complete"
du -sh "$dist_dir"
du -sh "$dist_dir"/wasm/*.wasm
