#!/usr/bin/env node
// Static consistency checks over a built site/dist. Node built-ins only.
// Usage: node scripts/check-site.mjs [dist-dir]   (default: ../dist relative to this script)
//
// Checks, all over dist/index.html unless noted:
//   1. every href="#x" has a matching id="x"
//   2. every local script src / link href file exists under dist/
//   3. no http(s):// references in <script src>, <link href>, or linked CSS url() (external
//      <a href> links are fine and are not checked for existence)
//   4. dist/_headers exists and contains "application/wasm" and "Content-Security-Policy"
//   5. a <noscript> element exists inside the demo section (aria-labelledby="demo-heading")
//
// Exits 1 and lists every failure if any check fails; exits 0 and prints a summary otherwise.

import { readFileSync, existsSync } from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';

const here = path.dirname(fileURLToPath(import.meta.url));
const distDir = path.resolve(process.argv[2] ?? path.join(here, '..', 'dist'));
const indexPath = path.join(distDir, 'index.html');

const failures = [];
const fail = (msg) => failures.push(msg);

if (!existsSync(indexPath)) {
  console.error(`check-site: ${indexPath} does not exist`);
  process.exit(1);
}

const html = readFileSync(indexPath, 'utf8');

function isLocalRef(value) {
  return (
    !!value &&
    !/^([a-z][a-z0-9+.-]*:)?\/\//i.test(value) &&
    !value.startsWith('#') &&
    !value.startsWith('mailto:') &&
    !value.startsWith('data:')
  );
}

function stripFragmentQuery(value) {
  return value.split('#')[0].split('?')[0];
}

function localPathFor(baseDir, value) {
  return path.join(baseDir, stripFragmentQuery(value).replace(/^\.\//, ''));
}

// --- Check 1: href="#x" anchors resolve to a matching id="x" ---------------------------------
const ids = new Set();
for (const m of html.matchAll(/\bid="([^"]+)"/g)) ids.add(m[1]);

for (const m of html.matchAll(/\bhref="#([^"]+)"/g)) {
  const target = m[1];
  if (!ids.has(target)) {
    fail(`href="#${target}" has no matching id="${target}" in index.html`);
  }
}

// --- Checks 2 & 3: local src/href existence + no http(s) in script/link, gather CSS to scan ---
const cssFilesToScan = [];

for (const tagMatch of html.matchAll(/<(script|link|img)\b[^>]*>/gi)) {
  const tag = tagMatch[0];
  const tagName = tagMatch[1].toLowerCase();
  const attrName = tagName === 'link' ? 'href' : 'src';
  const attrMatch = tag.match(new RegExp(`\\b${attrName}="([^"]*)"`, 'i'));
  if (!attrMatch) continue;
  const value = attrMatch[1];

  if (!isLocalRef(value)) {
    if (/^https?:\/\//i.test(value)) {
      fail(`<${tagName} ${attrName}="${value}"> is an external http(s) reference (must be local)`);
    }
    continue;
  }

  const filePath = localPathFor(distDir, value);
  if (!existsSync(filePath)) {
    fail(`<${tagName} ${attrName}="${value}"> -> missing file dist/${path.relative(distDir, filePath)}`);
    continue;
  }
  if (tagName === 'link' && stripFragmentQuery(value).endsWith('.css')) {
    cssFilesToScan.push(filePath);
  }
}

// <a href> to local (non-#, non-external) targets: allowed to be relative, but if present must
// point at a real file. External http(s) <a href> links are explicitly allowed and untouched.
for (const m of html.matchAll(/<a\b[^>]*\bhref="([^"]*)"[^>]*>/gi)) {
  const value = m[1];
  if (!isLocalRef(value)) continue; // external or "#..." — fine
  const filePath = localPathFor(distDir, value);
  if (!existsSync(filePath)) {
    fail(`<a href="${value}"> -> missing file dist/${path.relative(distDir, filePath)}`);
  }
}

// --- Check 3 (continued): CSS url() references, for every stylesheet linked from index.html --
// Per the spec this only checks for http(s):// leakage, not existence: a stylesheet is allowed
// to list format() fallbacks (e.g. woff/ttf after woff2) for fonts that were deliberately not
// vendored, since browsers skip fetching a url() whose format() they don't need.
for (const cssPath of cssFilesToScan) {
  let css;
  try {
    css = readFileSync(cssPath, 'utf8');
  } catch (err) {
    fail(`could not read ${path.relative(distDir, cssPath)}: ${err.message}`);
    continue;
  }
  const cssRel = path.relative(distDir, cssPath);
  for (const m of css.matchAll(/url\(\s*(['"]?)([^'")]+)\1\s*\)/g)) {
    const value = m[2];
    if (/^https?:\/\//i.test(value)) {
      fail(`${cssRel}: url(${value}) is an external http(s) reference`);
    }
  }
}

// --- Check 4: _headers exists and has the required directives --------------------------------
const headersPath = path.join(distDir, '_headers');
if (!existsSync(headersPath)) {
  fail('dist/_headers does not exist');
} else {
  const headers = readFileSync(headersPath, 'utf8');
  if (!headers.includes('application/wasm')) fail('dist/_headers is missing "application/wasm"');
  if (!headers.includes('Content-Security-Policy')) fail('dist/_headers is missing "Content-Security-Policy"');
}

// --- Check 5: a <noscript> element exists inside the demo section ----------------------------
// Robustly extract the section with aria-labelledby="demo-heading" by tag-depth scanning, since
// it contains nested <section> elements (race-panel, stream-panel) that a naive regex would stop
// at early.
function extractBalancedSection(source, openTagStart) {
  let pos = openTagStart;
  let depth = 0;
  while (pos < source.length) {
    const nextOpen = source.indexOf('<section', pos);
    const nextClose = source.indexOf('</section>', pos);
    if (nextClose === -1) return null;
    if (nextOpen !== -1 && nextOpen < nextClose) {
      depth++;
      pos = nextOpen + '<section'.length;
    } else {
      depth--;
      pos = nextClose + '</section>'.length;
      if (depth === 0) return source.slice(openTagStart, pos);
    }
  }
  return null;
}

const demoHeadingIdx = html.indexOf('aria-labelledby="demo-heading"');
if (demoHeadingIdx === -1) {
  fail('no section with aria-labelledby="demo-heading" found (cannot locate the demo section)');
} else {
  const sectionStart = html.lastIndexOf('<section', demoHeadingIdx);
  const demoSection = sectionStart === -1 ? null : extractBalancedSection(html, sectionStart);
  if (!demoSection) {
    fail('could not extract the demo section (unbalanced <section> tags)');
  } else if (!/<noscript\b/i.test(demoSection)) {
    fail('no <noscript> element found inside the demo section');
  }
}

// --- Report ------------------------------------------------------------------------------------
if (failures.length > 0) {
  console.error(`check-site: ${failures.length} failure(s) in ${path.relative(process.cwd(), distDir) || '.'}:`);
  for (const f of failures) console.error(`  - ${f}`);
  process.exit(1);
}

console.log(`check-site: all checks passed (${path.relative(process.cwd(), distDir) || '.'})`);
