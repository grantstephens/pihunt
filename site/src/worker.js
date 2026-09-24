import init, { digits_thm1, digits_thm2, wasm_memory_bytes } from './wasm/pi_digits_wasm.js';

const ready = init();
self.onmessage = async ({ data }) => {
  const { run, kind, method, pos, count } = data;
  try {
    await ready;
    const f = method === 'thm2' ? digits_thm2 : digits_thm1;
    if (kind === 'race') {
      const t0 = performance.now();
      const digits = f(pos, count);
      self.postMessage({ run, type: 'result', method, pos, digits, ms: performance.now() - t0, memBytes: wasm_memory_bytes() });
    } else {
      // Stream: each block is computed from scratch; nothing is carried between blocks.
      let p = pos;
      for (; p + count - 1 <= 200000; p += count) {
        const t0 = performance.now();
        const digits = digits_thm1(p, count);
        self.postMessage({ run, type: 'block', pos: p, digits, ms: performance.now() - t0, memBytes: wasm_memory_bytes() });
      }
      // The next block would overrun the 200,000-digit browser demo range: stop with a
      // visible message instead of just falling silent.
      self.postMessage({ run, type: 'end', message: 'reached the end of the browser demo range (position 200,000)' });
    }
  } catch (e) {
    self.postMessage({ run, type: 'error', message: String(e?.message ?? e) });
  }
};
