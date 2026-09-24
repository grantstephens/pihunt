// Node smoke test for the site/wasm build. Node built-ins only; it loads the `--target web`
// output produced by `wasm-pack build site/wasm --target web --release --out-dir ../dist/wasm`.
import { readFile } from 'node:fs/promises';
import init, { digits_thm1, digits_thm2 } from '../dist/wasm/pi_digits_wasm.js';

const wasm = await readFile(new URL('../dist/wasm/pi_digits_wasm_bg.wasm', import.meta.url));
await init({ module_or_path: wasm });
const ref = (await readFile(new URL('../data/pi-200k.txt', import.meta.url), 'utf8')).trim();
let failed = 0;
for (const pos of [1, 5, 762, 1999, 2000, 2001, 10000, 50000]) {
  const want = ref.slice(pos - 1, pos - 1 + 10);
  for (const [name, f] of [['thm1', digits_thm1], ['thm2', digits_thm2]]) {
    const got = f(pos, 10);
    if (got !== want) { failed++; console.error(`FAIL ${name} pos=${pos}: got ${got} want ${want}`); }
  }
}
for (const bad of [0, -1, 1.5, 200001, NaN]) {
  try { digits_thm1(bad, 10); failed++; console.error(`FAIL: pos=${bad} accepted`); } catch { /* expected */ }
}
if (failed) { console.error(`${failed} failure(s)`); process.exit(1); }
console.log('wasm smoke test: all positions match pi-200k.txt; bad positions rejected');
