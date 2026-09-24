# pihunt digits site

A static site: the Theorem 2 write-up plus an in-browser WebAssembly demo that races it against
the classical (Theorem 1) method, built from the same Rust code the `pihunt` CLI uses.

## Build

```sh
site/build.sh
```

Builds `site/wasm` with `wasm-pack` (release, `--target web`; wasm-pack already runs `wasm-opt`
on the output, so the script doesn't run it again), generates `site/data/pi-200k.txt` if it's
missing, runs the wasm smoke test and the JS unit tests, copies everything into `site/dist`
(excluding the dev-only `site/src/demo.test.mjs` and `site/src/package.json`), splices the
benchmark charts into `dist/index.html`, and runs `site/scripts/check-site.mjs` over the result.
Requires `cargo`, `wasm-pack`, and Node.js >= 18 on `PATH`; the script checks for these up front
and prints an install hint if any is missing. Never builds into the repo tree itself — set
`WASM_CARGO_TARGET_DIR` / `ROOT_CARGO_TARGET_DIR` to override the default `CARGO_TARGET_DIR`s
under `~/.cache/pihunt/`.

`site/dist` is the complete, self-contained deployable output.

## Serve locally

```sh
node site/scripts/serve.mjs
```

Serves `site/dist` (falling back to `site/src` and `site/` itself for files dist doesn't have,
e.g. while iterating without a full rebuild) at `http://localhost:8787/`. Useful for a manual
browser pass (light/dark theme, narrow viewport, the race/cancel/mismatch-banner/stream demo,
KaTeX rendering, chart rendering, and the `<noscript>` fallback with JS disabled) before deploying.

## Deploy

```sh
npx wrangler pages deploy site/dist --project-name <your-project>
```

Cloudflare Pages' git-integrated build image has no Rust toolchain (no `cargo`, no `wasm-pack`),
so this site is **not** built by Pages from source — build locally with `site/build.sh` first and
upload the prebuilt `site/dist` directory directly, as above.

## Notes

- `site/src/_headers` (copied into `dist/_headers` by the build) sets the Cloudflare Pages
  response headers, including the Content-Security-Policy. `style-src` includes `'unsafe-inline'`
  solely because KaTeX renders math by setting inline `style="..."` attributes (and individual
  `element.style.*` properties) on the spans it generates — there is no other inline styling on
  this page, and no inline `<script>` anywhere (`script-src` has no `'unsafe-inline'`).
  `'wasm-unsafe-eval'` is required for `WebAssembly.instantiate`/`instantiateStreaming`.
- `site/scripts/check-site.mjs` re-validates the build output (anchors, local asset existence, no
  external `http(s)://` script/link/CSS references, `_headers` present and correct, the demo
  section's `<noscript>` fallback). It runs as the last step of `build.sh` and can also be run
  standalone: `node site/scripts/check-site.mjs [dist-dir]`.
