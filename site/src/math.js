// Invokes KaTeX's auto-render on the page body. Loaded as a plain classic
// script, after ./vendor/katex/katex.min.js and
// ./vendor/katex/contrib/auto-render.min.js (both classic scripts that
// define the `katex` and `renderMathInElement` globals), so this file can
// call straight into them without any bundler or module wiring. Self-hosted
// per this project's CSP (no inline <script>, no third-party requests).
//
// Wide display equations are handled in CSS alone (style.css's
// `.katex-display { overflow-x: auto; ... }`): they scroll horizontally
// inside their own box at full, readable size rather than being shrunk to
// fit. Nothing here needs to run after render.
(function () {
  function run() {
    if (typeof renderMathInElement !== 'function') return;
    renderMathInElement(document.body, {
      delimiters: [
        { left: '\\[', right: '\\]', display: true },
        { left: '\\(', right: '\\)', display: false },
      ],
      throwOnError: false,
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', run);
  } else {
    run();
  }
})();
