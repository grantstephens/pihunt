#!/usr/bin/env node
// Renders the two benchmark charts (time vs position, memory vs position) as
// inline SVG and splices them into an HTML file at build time.
//
// Usage: node scripts/charts.mjs <input-html-path> > <output-html-path>
//
// Node built-ins only, no dependencies. Reads site/data/benchmarks.json
// (relative to this script), replaces the `<!-- CHART:time -->` and
// `<!-- CHART:memory -->` markers in the given HTML file, and writes the
// result to stdout.

import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const dataPath = path.join(here, '..', 'data', 'benchmarks.json');

// --- Colour: fixed categorical slots (see dataviz skill's palette.md),
// referenced as CSS custom properties so both themes resolve correctly.
// slot 1 (blue) = Theorem 1, slot 2 (orange) = Theorem 2, slot 3 (aqua) =
// Gourdon's 2003 pidec. Validated: light/dark categorical, 3-slot all-pairs
// (validate_palette.js), worst adjacent CVD deltaE 9.2 light / 9.4 dark.
const SERIES = {
  thm1: { varName: '--series-1', label: 'Theorem 1' },
  thm2: { varName: '--series-2', label: 'Theorem 2' },
  pidec: { varName: '--series-3', label: 'pidec 2003 (Pentium III, 2003)' },
};

const WIDTH = 640;
const HEIGHT = 400;
const MARGIN = { top: 24, right: 132, bottom: 44, left: 64 };
const PLOT_W = WIDTH - MARGIN.left - MARGIN.right;
const PLOT_H = HEIGHT - MARGIN.top - MARGIN.bottom;

function escapeXml(s) {
  return String(s)
    .replace(/&/g, '&amp;')
    .replace(/</g, '&lt;')
    .replace(/>/g, '&gt;')
    .replace(/"/g, '&quot;');
}

// Smallest power of 10 <= v, and largest power of 10 >= v.
function decadeFloor(v) {
  return Math.pow(10, Math.floor(Math.log10(v)));
}
function decadeCeil(v) {
  return Math.pow(10, Math.ceil(Math.log10(v)));
}

function log10Scale(value, domainMin, domainMax, rangeMin, rangeMax) {
  const t = (Math.log10(value) - Math.log10(domainMin)) / (Math.log10(domainMax) - Math.log10(domainMin));
  return rangeMin + t * (rangeMax - rangeMin);
}

function decadeTicks(min, max) {
  const lo = Math.round(Math.log10(min));
  const hi = Math.round(Math.log10(max));
  const ticks = [];
  for (let e = lo; e <= hi; e++) ticks.push(Math.pow(10, e));
  return ticks;
}

function formatPos(v) {
  return v.toLocaleString('en-US');
}

function formatSeconds(v) {
  if (v >= 1) return `${v >= 1000 ? v.toLocaleString('en-US') : v} s`;
  return `${(v * 1000).toLocaleString('en-US')} ms`;
}

function formatMiB(v) {
  return `${v} MiB`;
}

/**
 * Build one log-log line chart as an SVG string.
 * @param {object} opts
 * @param {{pos:number, values: Record<string, number|null>}[]} opts.points row data
 * @param {string[]} opts.seriesKeys which SERIES keys to draw, in draw order
 * @param {string} opts.yLabel axis title
 * @param {(v:number)=>string} opts.formatY tick/end-label formatter for y values
 * @param {string} opts.ariaLabel full aria-label summarising the chart
 * @param {Set<string>} opts.dashed series keys to draw dashed
 */
function renderChart({ points, seriesKeys, yLabel, formatY, ariaLabel, dashed }) {
  const xs = points.map((p) => p.pos);
  const allY = [];
  for (const p of points) {
    for (const k of seriesKeys) {
      const v = p.values[k];
      if (v != null) allY.push(v);
    }
  }
  const xMin = decadeFloor(Math.min(...xs));
  const xMax = decadeCeil(Math.max(...xs));
  const yMin = decadeFloor(Math.min(...allY));
  const yMax = decadeCeil(Math.max(...allY));

  const xTicks = decadeTicks(xMin, xMax);
  const yTicks = decadeTicks(yMin, yMax);

  const X = (pos) => MARGIN.left + log10Scale(pos, xMin, xMax, 0, PLOT_W);
  const Y = (v) => MARGIN.top + PLOT_H - log10Scale(v, yMin, yMax, 0, PLOT_H);

  const parts = [];

  // Decade gridlines (recessive, hairline) -- vertical (x) then horizontal (y).
  for (const t of xTicks) {
    const x = X(t);
    parts.push(
      `<line class="chart-grid" x1="${x.toFixed(1)}" y1="${MARGIN.top}" x2="${x.toFixed(1)}" y2="${(MARGIN.top + PLOT_H).toFixed(1)}" />`
    );
    parts.push(
      `<text class="chart-tick chart-tick-x" x="${x.toFixed(1)}" y="${(MARGIN.top + PLOT_H + 18).toFixed(1)}" text-anchor="middle">${escapeXml(formatPos(t))}</text>`
    );
  }
  for (const t of yTicks) {
    const y = Y(t);
    parts.push(
      `<line class="chart-grid" x1="${MARGIN.left}" y1="${y.toFixed(1)}" x2="${(MARGIN.left + PLOT_W).toFixed(1)}" y2="${y.toFixed(1)}" />`
    );
    parts.push(
      `<text class="chart-tick chart-tick-y" x="${(MARGIN.left - 8).toFixed(1)}" y="${(y + 4).toFixed(1)}" text-anchor="end">${escapeXml(formatY(t))}</text>`
    );
  }

  // Axis baselines.
  parts.push(
    `<line class="chart-axis" x1="${MARGIN.left}" y1="${(MARGIN.top + PLOT_H).toFixed(1)}" x2="${(MARGIN.left + PLOT_W).toFixed(1)}" y2="${(MARGIN.top + PLOT_H).toFixed(1)}" />`
  );
  parts.push(
    `<line class="chart-axis" x1="${MARGIN.left}" y1="${MARGIN.top}" x2="${MARGIN.left}" y2="${(MARGIN.top + PLOT_H).toFixed(1)}" />`
  );

  // Axis titles.
  parts.push(
    `<text class="chart-axis-title" x="${(MARGIN.left + PLOT_W / 2).toFixed(1)}" y="${HEIGHT - 6}" text-anchor="middle">Digit position n</text>`
  );
  parts.push(
    `<text class="chart-axis-title" x="16" y="${(MARGIN.top + PLOT_H / 2).toFixed(1)}" text-anchor="middle" transform="rotate(-90 16 ${(MARGIN.top + PLOT_H / 2).toFixed(1)})">${escapeXml(yLabel)}</text>`
  );

  // Series lines, markers, and direct end-labels. Null points are skipped
  // (break the line) rather than zeroed.
  for (const key of seriesKeys) {
    const meta = SERIES[key];
    const segPoints = points
      .map((p) => ({ pos: p.pos, v: p.values[key] }))
      .filter((p) => p.v != null);
    if (segPoints.length === 0) continue;

    const pathD = segPoints
      .map((p, i) => `${i === 0 ? 'M' : 'L'}${X(p.pos).toFixed(1)},${Y(p.v).toFixed(1)}`)
      .join(' ');
    const dashAttr = dashed.has(key) ? ' stroke-dasharray="6 5"' : '';
    parts.push(
      `<path class="chart-line" d="${pathD}" fill="none" stroke="var(${meta.varName})" stroke-width="2"${dashAttr} />`
    );
    for (const p of segPoints) {
      parts.push(
        `<circle class="chart-marker" cx="${X(p.pos).toFixed(1)}" cy="${Y(p.v).toFixed(1)}" r="4" fill="var(${meta.varName})" stroke="var(--chart-surface)" stroke-width="2" />`
      );
    }
    // Direct label at the line's last plotted point.
    const last = segPoints[segPoints.length - 1];
    const lx = X(last.pos) + 8;
    const ly = Y(last.v);
    parts.push(
      `<text class="chart-series-label" x="${lx.toFixed(1)}" y="${(ly + 4).toFixed(1)}" fill="var(${meta.varName})">${escapeXml(meta.label)}</text>`
    );
  }

  return (
    `<svg class="chart-svg" viewBox="0 0 ${WIDTH} ${HEIGHT}" xmlns="http://www.w3.org/2000/svg" ` +
    `role="img" aria-label="${escapeXml(ariaLabel)}">` +
    parts.join('') +
    `</svg>`
  );
}

function buildCharts(bench) {
  const points = bench.rows.map((r) => ({
    pos: r.pos,
    values: { thm1: r.thm1_s, thm2: r.thm2_s, pidec: r.pidec_s },
  }));
  const memPoints = bench.rows.map((r) => ({
    pos: r.pos,
    values: { thm1: r.thm1_mib, thm2: r.thm2_mib },
  }));

  // Derived from `bench.rows` rather than hardcoded, so these aria-label sentences can't go
  // stale the way they did when benchmarks.json was updated but these string literals weren't.
  const firstRow = bench.rows[0];
  const lastRow = bench.rows[bench.rows.length - 1];

  const timeChart = renderChart({
    points,
    seriesKeys: ['pidec', 'thm1', 'thm2'],
    yLabel: 'Wall-clock time',
    formatY: formatSeconds,
    dashed: new Set(['pidec']),
    ariaLabel:
      'Log-log chart of wall-clock time versus digit position, for Theorem 1, Theorem 2, and Gourdon’s 2003 pidec (Pentium III, 2003). Theorem 2 is fastest at every measured position from 10,000 upward, and the gap widens with position: at ' +
      formatPos(lastRow.pos) +
      ', Theorem 1 takes ' +
      formatSeconds(lastRow.thm1_s) +
      ' and Theorem 2 takes ' +
      formatSeconds(lastRow.thm2_s) +
      '. pidec has no data point past 1,000,000.',
  });

  const memChart = renderChart({
    points: memPoints,
    seriesKeys: ['thm1', 'thm2'],
    yLabel: 'Peak memory (MiB)',
    formatY: formatMiB,
    dashed: new Set(),
    ariaLabel:
      'Log-log chart of peak memory versus digit position, for Theorem 1 and Theorem 2. Theorem 1 stays flat at roughly 5 MiB across all positions; Theorem 2 grows with position, from ' +
      formatMiB(firstRow.thm2_mib) +
      ' at ' +
      formatPos(firstRow.pos) +
      ' to ' +
      formatMiB(lastRow.thm2_mib) +
      ' at ' +
      formatPos(lastRow.pos) +
      ', trading bounded memory for a large speed advantage.',
  });

  return { timeChart, memChart };
}

function main() {
  const inputPath = process.argv[2];
  if (!inputPath) {
    process.stderr.write('usage: node scripts/charts.mjs <input-html-path> > <output-html-path>\n');
    process.exit(1);
  }
  const html = readFileSync(inputPath, 'utf8');
  const bench = JSON.parse(readFileSync(dataPath, 'utf8'));
  const { timeChart, memChart } = buildCharts(bench);

  let out = html;
  if (!out.includes('<!-- CHART:time -->')) {
    process.stderr.write('charts.mjs: <!-- CHART:time --> marker not found\n');
    process.exit(1);
  }
  if (!out.includes('<!-- CHART:memory -->')) {
    process.stderr.write('charts.mjs: <!-- CHART:memory --> marker not found\n');
    process.exit(1);
  }
  out = out.replace('<!-- CHART:time -->', timeChart);
  out = out.replace('<!-- CHART:memory -->', memChart);

  process.stdout.write(out);
}

main();
