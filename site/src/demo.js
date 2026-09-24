// Demo front end: race (Theorem 1 vs Theorem 2) and stream panels.
//
// This module is imported by `node --test` (see demo.test.mjs) as well as by
// the page (`<script type="module" src="./demo.js">`). The pure functions
// below must not touch `document`/`window` at import time; all DOM wiring is
// guarded behind `typeof document !== 'undefined'` so this file stays
// importable in plain Node.

// --- Pure functions (unit-tested in Node) -----------------------------------

const POSITION_CAP = 200000;
const POSITION_FORMAT_MESSAGE = 'Enter a whole number from 1 to 200,000.';
const POSITION_RANGE_MESSAGE = 'Positions above 200,000 are beyond the browser demo; use the CLI.';

/**
 * Parse a user-entered position string. Accepts plain digits with optional
 * surrounding whitespace and thousands commas; rejects everything else
 * (including scientific notation, decimals, and negative numbers) so the
 * input stays unambiguous.
 * @param {string} input
 * @returns {{ok: true, pos: number} | {ok: false, message: string}}
 */
export function parsePosition(input) {
  const cleaned = String(input).replace(/[\s,]/g, '');
  if (!/^\d+$/.test(cleaned)) {
    return { ok: false, message: POSITION_FORMAT_MESSAGE };
  }
  const pos = Number(cleaned);
  if (pos < 1) {
    return { ok: false, message: POSITION_FORMAT_MESSAGE };
  }
  if (pos > POSITION_CAP) {
    return { ok: false, message: POSITION_RANGE_MESSAGE };
  }
  return { ok: true, pos };
}

/**
 * Compare computed digits against the known reference (site/data/pi-200k.txt,
 * positions 1..=200000). Returns 'unchecked' when there is no reference text,
 * or when the requested range runs past what the reference covers.
 * @param {string | null | undefined} ref
 * @param {number} pos 1-based
 * @param {string} digits
 * @returns {'match' | 'mismatch' | 'unchecked'}
 */
export function checkAgainstReference(ref, pos, digits) {
  if (ref == null) return 'unchecked';
  const start = pos - 1;
  const end = start + digits.length;
  if (start < 0 || end > ref.length) return 'unchecked';
  return ref.slice(start, end) === digits ? 'match' : 'mismatch';
}

/**
 * Whether two computed digit strings (Theorem 1 vs Theorem 2 lanes) agree.
 * @param {string} a
 * @param {string} b
 * @returns {boolean}
 */
export function lanesAgree(a, b) {
  return a === b;
}

// --- Internal helpers --------------------------------------------------------

function groupDigits(str, size = 10) {
  const groups = [];
  for (let i = 0; i < str.length; i += size) groups.push(str.slice(i, i + size));
  return groups.join(' ');
}

function formatMs(ms) {
  return `${ms.toFixed(1)} ms`;
}

function formatKiB(bytes) {
  return `${(bytes / 1024).toFixed(1)} KiB`;
}

// Flips the last digit so the corrupted string is guaranteed to differ.
// Used only behind the #test-mismatch hook, to exercise the disagreement UI.
function corruptDigits(digits) {
  if (!digits) return digits;
  const last = digits[digits.length - 1];
  const replacement = last === '0' ? '1' : '0';
  return digits.slice(0, -1) + replacement;
}

const RACE_DIGIT_COUNT = 20;
const STREAM_BLOCK_SIZE = 1000;
const SOFT_WARNING_THRESHOLD = 100000;
const SOFT_WARNING_MESSAGE = 'This may take a minute or more in the browser.';
const AGREE_MESSAGE = '✓ lanes agree';
const MISMATCH_MESSAGE = '✗ Theorem 1 and Theorem 2 disagree — please report this';

// --- DOM wiring ---------------------------------------------------------------
// Only runs in a browser; never touches `document` at import time so this
// file stays importable by `node --test`.

if (typeof document !== 'undefined') {
  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', initDemo);
  } else {
    initDemo();
  }
}

function initDemo() {
  let currentRun = 0;
  let liveWorkers = [];
  let referenceTextPromise = null;

  function loadReferenceText() {
    if (!referenceTextPromise) {
      referenceTextPromise = fetch('./data/pi-200k.txt')
        .then((r) => (r.ok ? r.text() : null))
        .then((t) => (typeof t === 'string' ? t.trim() : null))
        .catch(() => null);
    }
    return referenceTextPromise;
  }

  // Every Race, Cancel or Stream action bumps the run id and terminates any
  // live workers, so stale messages (checked via `data.run !== currentRun`)
  // are dropped by the handlers below.
  function bumpRun() {
    currentRun += 1;
    for (const w of liveWorkers) {
      try {
        w.terminate();
      } catch {
        /* already gone */
      }
    }
    liveWorkers = [];
    return currentRun;
  }

  function startLaneTimer(lane, runId) {
    lane._timerStopped = false;
    const t0 = performance.now();
    function tick() {
      if (lane._timerStopped || runId !== currentRun) return;
      lane.timerEl.textContent = formatMs(performance.now() - t0);
      lane._raf = requestAnimationFrame(tick);
    }
    tick();
  }

  function stopLaneTimer(lane) {
    lane._timerStopped = true;
    if (lane._raf) cancelAnimationFrame(lane._raf);
  }

  initRacePanel();
  initStreamPanel();

  // --- Race panel ---

  function initRacePanel() {
    const posInput = document.getElementById('race-pos');
    const startBtn = document.getElementById('race-start');
    const cancelBtn = document.getElementById('race-cancel');
    const errorEl = document.getElementById('race-error');
    const warningEl = document.getElementById('race-warning');
    const agreementEl = document.getElementById('lanes-agreement');
    const referenceEl = document.getElementById('reference-check');
    const presetBtns = document.querySelectorAll('[data-race-preset]');

    const lanes = {
      thm1: {
        digitsEl: document.getElementById('thm1-digits'),
        msEl: document.getElementById('thm1-ms'),
        memEl: document.getElementById('thm1-mem'),
        timerEl: document.getElementById('thm1-timer'),
      },
      thm2: {
        digitsEl: document.getElementById('thm2-digits'),
        msEl: document.getElementById('thm2-ms'),
        memEl: document.getElementById('thm2-mem'),
        timerEl: document.getElementById('thm2-timer'),
      },
    };

    // { runId, pos, results: { thm1?, thm2? }, referenceText? }
    let raceState = null;

    presetBtns.forEach((btn) => {
      btn.addEventListener('click', () => {
        posInput.value = btn.getAttribute('data-race-preset');
      });
    });

    cancelBtn.addEventListener('click', () => {
      bumpRun();
      for (const lane of Object.values(lanes)) stopLaneTimer(lane);
    });

    startBtn.addEventListener('click', () => {
      startRace();
    });

    function setError(msg) {
      if (msg) {
        errorEl.textContent = msg;
        errorEl.classList.remove('hidden');
      } else {
        errorEl.textContent = '';
        errorEl.classList.add('hidden');
      }
    }

    function setWarning(show) {
      warningEl.textContent = SOFT_WARNING_MESSAGE;
      warningEl.classList.toggle('hidden', !show);
    }

    function startRace() {
      const parsed = parsePosition(posInput.value);
      if (!parsed.ok) {
        setError(parsed.message);
        return;
      }
      setError(null);
      setWarning(parsed.pos > SOFT_WARNING_THRESHOLD);

      const runId = bumpRun();
      agreementEl.textContent = '';
      agreementEl.className = '';
      referenceEl.textContent = '';
      referenceEl.className = '';

      raceState = { runId, pos: parsed.pos, results: {} };

      for (const [method, lane] of Object.entries(lanes)) {
        lane.digitsEl.textContent = '';
        lane.msEl.textContent = '–';
        lane.memEl.textContent = '–';
        startLaneTimer(lane, runId);

        const worker = new Worker('./worker.js', { type: 'module' });
        liveWorkers.push(worker);
        worker.onmessage = ({ data }) => {
          if (data.run !== currentRun) return;
          handleRaceMessage(method, lane, data);
        };
        worker.onerror = (e) => {
          if (runId !== currentRun) return;
          handleRaceMessage(method, lane, { type: 'error', message: e.message });
        };
        worker.postMessage({ run: runId, kind: 'race', method, pos: parsed.pos, count: RACE_DIGIT_COUNT });
      }

      loadReferenceText().then((text) => {
        if (runId !== currentRun) return;
        raceState.referenceText = text;
        maybeUpdateReferenceCheck();
      });
    }

    function handleRaceMessage(method, lane, data) {
      stopLaneTimer(lane);
      if (data.type === 'error') {
        lane.digitsEl.textContent = `error: ${data.message}`;
        lane.msEl.textContent = '–';
        lane.memEl.textContent = '–';
        raceState.results[method] = { error: data.message };
      } else {
        lane.digitsEl.textContent = groupDigits(data.digits);
        lane.msEl.textContent = formatMs(data.ms);
        lane.memEl.textContent = formatKiB(data.memBytes);
        raceState.results[method] = { digits: data.digits, ms: data.ms, memBytes: data.memBytes };
      }
      maybeShowAgreement();
      maybeUpdateReferenceCheck();
    }

    function maybeShowAgreement() {
      const { thm1, thm2 } = raceState.results;
      if (!thm1 || !thm2) return;
      if (thm1.error || thm2.error) {
        agreementEl.textContent = `⚠ error: ${thm1.error ?? thm2.error}`;
        agreementEl.className = 'error';
        return;
      }
      let compareDigits2 = thm2.digits;
      if (typeof location !== 'undefined' && location.hash === '#test-mismatch') {
        compareDigits2 = corruptDigits(compareDigits2);
      }
      const agree = lanesAgree(thm1.digits, compareDigits2);
      agreementEl.textContent = agree ? AGREE_MESSAGE : MISMATCH_MESSAGE;
      agreementEl.className = agree ? 'match' : 'mismatch';
    }

    function maybeUpdateReferenceCheck() {
      if (!raceState || !('referenceText' in raceState)) return;
      const thm1Result = raceState.results.thm1;
      if (!thm1Result || thm1Result.error) return;
      const status = checkAgainstReference(raceState.referenceText, raceState.pos, thm1Result.digits);
      const label = {
        match: '✓ matches known digits',
        mismatch: '✗ MISMATCH vs known digits',
        unchecked: 'reference: not checked',
      }[status];
      referenceEl.textContent = label;
      referenceEl.className = status === 'mismatch' ? 'mismatch' : status === 'match' ? 'match' : '';
    }
  }

  // --- Stream panel ---

  function initStreamPanel() {
    const posInput = document.getElementById('stream-pos');
    const startBtn = document.getElementById('stream-start');
    const stopBtn = document.getElementById('stream-stop');
    const errorEl = document.getElementById('stream-error');
    const outputEl = document.getElementById('stream-output');
    const memEl = document.getElementById('stream-mem');

    function setError(msg) {
      if (msg) {
        errorEl.textContent = msg;
        errorEl.classList.remove('hidden');
      } else {
        errorEl.textContent = '';
        errorEl.classList.add('hidden');
      }
    }

    stopBtn.addEventListener('click', () => {
      bumpRun();
    });

    startBtn.addEventListener('click', () => {
      const parsed = parsePosition(posInput.value);
      if (!parsed.ok) {
        setError(parsed.message);
        return;
      }
      setError(null);

      const runId = bumpRun();
      outputEl.textContent = '';
      memEl.textContent = '–';

      const worker = new Worker('./worker.js', { type: 'module' });
      liveWorkers.push(worker);
      worker.onmessage = ({ data }) => {
        if (data.run !== currentRun) return;
        if (data.type === 'block') {
          const block = document.createElement('span');
          block.className = 'stream-block';
          block.textContent = groupDigits(data.digits);
          outputEl.appendChild(block);
          memEl.textContent = formatKiB(data.memBytes);
        } else if (data.type === 'error') {
          setError(data.message);
        }
      };
      worker.onerror = (e) => {
        if (runId !== currentRun) return;
        setError(e.message);
      };
      worker.postMessage({ run: runId, kind: 'stream', method: 'thm1', pos: parsed.pos, count: STREAM_BLOCK_SIZE });
    });
  }
}
