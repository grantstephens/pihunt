# Digits site: design spec

**Date:** 2026-09-23
**Status:** design approved in conversation, pending written-spec review

## Goal

Publish the Theorem 2 work as a static site the user deploys to Cloudflare Pages. It has two parts:

1. A paper-style write-up of the reconstruction and verified implementation of Xavier Gourdon's unpublished 2003 "Theorem 2": decimal digits of π at position n using O(m) memory in about O(n²/m · polylog) time.
2. A live in-browser demo that races Theorem 1 against Theorem 2, both compiled to WebAssembly from the same code the CLI uses, and streams digits in constant memory.

**Audience:** readers who know what BBP is (math/CS).
**Success:** a visitor can understand the claim and see it work. The in-browser ✓ checks and the reproduction commands let them verify it themselves.

## Decisions (from the brainstorm)

| question | decision |
|---|---|
| scope | Paper-style Theorem 2 write-up plus live demo |
| demo engine | Both theorems, raced side by side |
| location | `site/` in this repo (monorepo), so the demo runs the exact algorithm code |
| code sharing | Approach A: extract a `pi-digits` crate with switchable bignum backends |
| typography | Clean modern sans (Inter, Helvetica/system fallback); monospace digits (JetBrains Mono); both self-hosted |
| deploy | User runs `wrangler pages deploy site/dist` (direct upload); no git-integrated Pages build |

## 1. Crate split and bignum backends

### Layout

```
Cargo.toml              # becomes a workspace; pihunt stays the root package
crates/pi-digits/
  src/lib.rs
  src/bignum.rs         # the single backend seam
  src/nthdigit.rs       # Theorem 1 (moved from src/)
  src/nthdigit2.rs      # Theorem 2 (moved from src/)
  src/pi_ref.rs         # small-position fallback, integer Machin, backend-agnostic
  src/mem_profile.rs    # moved from src/, behind its existing feature
  tests/                # tests/nthdigit.rs and tests/nthdigit2.rs move here
site/wasm/              # wasm-bindgen wrapper crate (section 2)
```

`pihunt` depends on `pi-digits` with features `gmp` and `parallel`. It re-exports `nthdigit` and `nthdigit2`, so `pihunt::nthdigit…` paths and the CLI behave exactly as today. The PSLQ/BBP code stays where it is.

### Features on `pi-digits`

- `gmp` (default): `rug`, i.e. today's behaviour.
- `pure`: a pure-Rust bignum. `dashu-int` is preferred; `num-bigint` is the fallback if `dashu-int` blocks.
- `parallel` (default): rayon. Turned off for WASM.
- Exactly one of `gmp` or `pure` must be enabled, otherwise `compile_error!`.

### The seam

`bignum.rs` defines a newtype `Big` exposing only the operations the algorithms use: construction from and conversion to `u64`, multiply, add, subtract, remainder by `Big` and by `u64`, power, and comparison. There is one implementation per backend.

The certification path is not behind the seam. That covers the u128 fixed-point accumulation, `error_units`, and the digit-boundary retry loop. It is plain integer code and must not change.

### Small-position fallback

`digits_via_mpfr` (used below position 2000, or when n0 exceeds the certifiable cap) is replaced in both backends by `pi_ref`: π computed with integer Machin arithmetic on `Big`, at `n + count + guard` digits. It is tested against MPFR. The certification semantics are unchanged: this path is exact.

## 2. WASM module and demo

### `site/wasm`

A `wasm-bindgen` `cdylib` depending on `pi-digits` with `default-features = false, features = ["pure"]`. It exports:

- `digits_thm1(pos: f64, count: u32) -> String` (1-based position, like the CLI)
- `digits_thm2(pos: f64, count: u32) -> String`, using `default_mem_bits`
- `wasm_memory_bytes() -> f64`

The build is `wasm-pack build --target web --release` with LTO, followed by `wasm-opt -O3`.

### Race panel

- The user enters a position. Presets: 1 000, 10 000, 50 000, 100 000. Soft warning above 10⁵; hard cap at 2·10⁵.
- **Race** starts two Web Workers at the same moment, one per theorem, each running single-threaded WASM.
- Each lane shows the algorithm, a live elapsed timer, then the 10 digits, final time and WASM memory.
- **Cancel** terminates both workers.
- **Checks:**
  - The two lanes must agree; a disagreement is flagged prominently.
  - Results up to position 200 000 are checked against a shipped `pi-200k.txt` and marked ✓ or ✗.

### Stream panel

- Start position p. One worker computes successive 10-digit blocks, each from scratch, using Theorem 1.
- Blocks are appended live, with a flat memory readout. A Stop button ends it.
- This is the original "keep churning, don't store previous digits" goal.

### Failure handling

- No WASM or Worker support: a plain message linking to the CLI and the benchmark tables.
- A worker error: shown in its lane, with no hang.

### Expected performance

These are estimates only. Real numbers are measured during the build and the page reports those. At 10⁵: Theorem 1 roughly 10–20 s, Theorem 2 a few seconds.

## 3. Page content and structure

One long page:

1. Hero: title (e.g. *"Decimal digits of π in O(√n) memory: Gourdon's unpublished 2003 Theorem 2, reconstructed"*), a two-sentence summary, then the race.
2. Race panel, then stream panel.
3. Background: BBP gives hex digits; decimal is hard; the timeline Plouffe 1996 (n³) → Bellard 1997 (n²) → Gourdon 2003 Theorem 1 (n²/log²n, O(log²n) memory).
4. The claim: Theorem 2 as stated; "details will be added soon"; the literature table.
5. The reconstruction: CVZ-accelerated series, (P, T, D) recurrence, partial-fraction split into main / Lucas / p-adic (the new p-adic Lucas recursion), chunked remainder tree. Key equations inline; derivations in collapsible blocks; every claim tagged **proved / measured / conjecture**.
6. Complexity: the bound, and log²(n/m) as the fingerprint matching Gourdon's claim.
7. Implementation: u128 fixed-point certification with every rounding counted; the memory story 812 → 447 → 69 MiB at 10⁷; the remaining O(N) cofactor term.
8. Benchmarks: two log-log charts (time vs position, memory vs position; Theorem 1, Theorem 2, 2003 pidec) plus the full table.
9. Verification: MPFR agreement at 10⁵/10⁶/10⁷ with the position-convention table; the two-backend test runs; the browser ✓ checks; an "errors we caught" box (the certification hole, the off-by-one labels, the rounding-not-truncating verification script).
10. Limits and caveats: the cofactor term; "reconstructed and verified", not a formal proof; that it is Gourdon's own method is a conjecture.
11. Reproduce it: exact `cargo` and `pihunt digit` commands.
12. Postscript: one paragraph on the base-10 BBP hunt, its null result and the weight argument.
13. References.

**Sources of truth:** `docs/nthdigit-theorem2.md`, `docs/nthdigit.md`, `docs/findings-bbp-hunt.md`. Every number on the page must trace to those or to `site/data/benchmarks.json`.

**Presentation:**
- Inter (variable woff2) for body and headings, with a Helvetica/system-sans fallback stack.
- JetBrains Mono for digits and code.
- KaTeX vendored for the maths.
- Light and dark themes; phone-width layout with no horizontal scroll.
- Charts are inline SVG generated at build time from `benchmarks.json`.
- No third-party requests at runtime.

## 4. Build and deploy

```
site/
  src/        index.html, style.css, app.js, worker.js
  data/       benchmarks.json, pi-200k.txt (generated if missing)
  wasm/       wasm-bindgen crate
  vendor/     KaTeX, fonts
  scripts/    charts.mjs (no npm deps), gen-pi.sh, check-site.mjs
  build.sh    produces site/dist/
  dist/       gitignored
```

`site/build.sh` does the following, in order:
1. `wasm-pack` build plus `wasm-opt`.
2. Generate `pi-200k.txt` if missing, using the CLI's exact MPFR path.
3. WASM smoke test (section 5).
4. Render charts into the HTML.
5. Copy assets.
6. Write `_headers`: long cache for hashed assets, `application/wasm`, CSP `default-src 'self'` with the minimal additions KaTeX/WASM need. No COOP/COEP.
7. Run `check-site.mjs`.

**One-time tooling**, installed through the user's rustup/cargo, not system-wide: the `wasm32-unknown-unknown` target, `wasm-pack`, and binaryen's `wasm-opt` if missing.

**Deploy:** `npx wrangler pages deploy site/dist --project-name <name>`.

**Size budget:** about 1 MB total. WASM ~200–400 KB. The WASM and the π digits load lazily when the demo is first used.

## 5. Testing and verification

- **Refactor is behaviour-preserving:** with the default `gmp` backend, all existing suites pass unchanged (102 tests), clippy is clean, and the CLI prints the same digits at 10⁴, 10⁵, 10⁶.
- **Pure backend parity:** the full `pi-digits` suite runs with `--no-default-features --features pure`, still against MPFR (dev-dependency). A dedicated test compares `Big` operations across backends on random inputs, including multi-thousand-bit values. The ignored 10⁶ check runs once per backend.
- **WASM smoke test:** Node loads the built module and checks both theorems at positions 1, 762, 10⁴ and 5·10⁴ against `pi-200k.txt`. It runs inside `build.sh`.
- **Site checks:** `check-site.mjs` verifies internal links, anchors and `_headers`. A manual pass in a real browser covers both themes and phone width, via a local server and the Chrome tool, reporting what was observed.
- **Memory safety rules** (this machine has 7.7 GiB): heavy runs happen one at a time under `ulimit -v`; never `--include-ignored` across the workspace; `cargo -j 4`.

## Out of scope

- Multi-threaded WASM (SharedArrayBuffer, COOP/COEP).
- Browser positions above 2·10⁵.
- Git-integrated Pages build.
- Analytics.
- Porting the PSLQ/BBP machinery to the site.
