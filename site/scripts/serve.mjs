#!/usr/bin/env node
// Tiny static file server for local demo preview. Node built-ins only.
//
// Serves, in order of preference: site/dist, then site/src, then site itself
// -- so /wasm/* comes from dist (the wasm-pack build output), pages/JS from
// src, and /data/* from the repo's site/data. Usage: node serve.mjs [port]
import { createServer } from 'node:http';
import { readFile, stat } from 'node:fs/promises';
import { fileURLToPath } from 'node:url';
import path from 'node:path';

const here = path.dirname(fileURLToPath(import.meta.url));
const siteRoot = path.resolve(here, '..');
const roots = [path.join(siteRoot, 'dist'), path.join(siteRoot, 'src'), siteRoot];

const MIME_TYPES = {
  '.html': 'text/html',
  '.js': 'text/javascript',
  '.mjs': 'text/javascript',
  '.css': 'text/css',
  '.wasm': 'application/wasm',
  '.txt': 'text/plain',
  '.json': 'application/json',
};

function contentTypeFor(filePath) {
  return MIME_TYPES[path.extname(filePath).toLowerCase()] ?? 'application/octet-stream';
}

// Resolves a request path against each root in turn, refusing anything that
// would escape that root (no directory traversal outside the roots).
async function resolveFile(rawUrl) {
  let pathname = decodeURIComponent(rawUrl.split('?')[0]);
  if (pathname === '/') pathname = '/index.html';

  const normalized = path.normalize(pathname);
  if (normalized.split(path.sep).includes('..')) return null;

  for (const root of roots) {
    const candidate = path.join(root, normalized);
    const relative = path.relative(root, candidate);
    if (relative.startsWith('..') || path.isAbsolute(relative)) continue; // escaped root
    try {
      const st = await stat(candidate);
      if (st.isFile()) return candidate;
    } catch {
      // not found under this root; try the next one
    }
  }
  return null;
}

const port = Number(process.argv[2]) || 8787;

const server = createServer(async (req, res) => {
  try {
    const file = await resolveFile(req.url ?? '/');
    if (!file) {
      res.writeHead(404, { 'Content-Type': 'text/plain' });
      res.end('Not found');
      return;
    }
    const body = await readFile(file);
    res.writeHead(200, { 'Content-Type': contentTypeFor(file) });
    res.end(body);
  } catch {
    res.writeHead(500, { 'Content-Type': 'text/plain' });
    res.end('Internal error');
  }
});

server.listen(port, () => {
  console.log(`serving ${roots.join(', ')} at http://localhost:${port}/`);
});
