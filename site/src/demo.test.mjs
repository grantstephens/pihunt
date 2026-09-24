import test from 'node:test';
import assert from 'node:assert/strict';
import { parsePosition, checkAgainstReference, lanesAgree, browserSupportsDemo } from './demo.js';

test('parsePosition accepts plain integers in range', () => {
  assert.deepEqual(parsePosition('1'), { ok: true, pos: 1 });
  assert.deepEqual(parsePosition(' 100000 '), { ok: true, pos: 100000 });
  assert.deepEqual(parsePosition('200000'), { ok: true, pos: 200000 });
  assert.deepEqual(parsePosition('10,000'), { ok: true, pos: 10000 });
});

test('parsePosition rejects everything else with a message', () => {
  for (const bad of ['', '0', '-5', '1.5', '1e5', 'abc', '200001', '99999999999']) {
    const r = parsePosition(bad);
    assert.equal(r.ok, false, `accepted ${JSON.stringify(bad)}`);
    assert.match(r.message, /\S/);
  }
});

test('checkAgainstReference', () => {
  const ref = '1415926535';
  assert.equal(checkAgainstReference(ref, 1, '14159'), 'match');
  assert.equal(checkAgainstReference(ref, 1, '14158'), 'mismatch');
  assert.equal(checkAgainstReference(ref, 8, '53599'), 'unchecked'); // runs past the reference
  assert.equal(checkAgainstReference(null, 1, '14159'), 'unchecked');
});

test('lanesAgree', () => {
  assert.equal(lanesAgree('0123456789', '0123456789'), true);
  assert.equal(lanesAgree('0123456789', '0123456780'), false);
});

test('browserSupportsDemo', () => {
  assert.equal(browserSupportsDemo({ WebAssembly: {}, Worker: function () {} }), true);
  assert.equal(browserSupportsDemo({ Worker: function () {} }), false); // no WebAssembly
  assert.equal(browserSupportsDemo({ WebAssembly: {} }), false); // no Worker
  assert.equal(browserSupportsDemo({}), false);
  // Plain Node (no DOM globals) must not accidentally report support.
  assert.equal(browserSupportsDemo(globalThis), false);
});
