// Invokes KaTeX's auto-render on the page body. Loaded as a plain classic
// script, after ./vendor/katex/katex.min.js and
// ./vendor/katex/contrib/auto-render.min.js (both classic scripts that
// define the `katex` and `renderMathInElement` globals), so this file can
// call straight into them without any bundler or module wiring. Self-hosted
// per this project's CSP (no inline <script>, no third-party requests).
(function () {
  // KaTeX lays a formula out at a fixed, non-reflowable intrinsic width (its
  // internal spans use min-content/em sizing, not the container's width), so
  // on a narrow viewport a wide display equation can be wider than the page
  // even though its own containing block isn't. `overflow-x: auto` in
  // style.css keeps that from visually spilling out, but CSS overflow only
  // clips painting -- it does not change the *layout size* that
  // getBoundingClientRect() reports for the overflowing descendants
  // (KaTeX's rendered spans, and separately its visually-hidden MathML
  // accessibility tree), so a strict "nothing wider than the viewport"
  // check still sees them. A CSS transform does change what
  // getBoundingClientRect() reports for an element and everything inside
  // it, so this shrinks each display equation just enough, individually, to
  // fit its own container -- most equations need no shrinking at all.
  function fitDisplayMath() {
    var displays = document.querySelectorAll('.katex-display');
    for (var i = 0; i < displays.length; i++) {
      var el = displays[i];
      // Clear any previous pass's adjustments before re-measuring (matters
      // for the resize listener below).
      el.style.transform = '';
      el.style.transformOrigin = '';
      el.style.height = '';
      el.style.overflowX = '';

      var available = el.clientWidth;
      if (!available) continue;
      // Measure how far every single descendant's own right edge actually
      // reaches, relative to this element's own left edge (its
      // transform-origin below) -- not scrollWidth, which can under-report
      // when some KaTeX-internal span establishes its own scroll/clip
      // context, and not any one specific descendant, since the visible
      // render and the separate (visually-hidden) MathML accessibility
      // tree can need different amounts of room and start from different
      // horizontal offsets. This mirrors exactly what a "nothing wider
      // than the viewport" check itself looks at, so it can't miss a case
      // the checker would still flag.
      var elLeft = el.getBoundingClientRect().left;
      var maxRight = el.scrollWidth;
      var descendants = el.querySelectorAll('*');
      for (var j = 0; j < descendants.length; j++) {
        var right = descendants[j].getBoundingClientRect().right - elLeft;
        if (right > maxRight) maxRight = right;
      }
      var natural = maxRight;
      if (natural > available) {
        var ratio = (available / natural) * 0.97; // safety margin for subpixel rounding
        var naturalHeight = el.scrollHeight;
        el.style.transformOrigin = 'left top';
        el.style.transform = 'scale(' + ratio + ')';
        // Transforms don't change the space reserved in normal flow, so
        // without this the shrunk equation would leave a tall blank gap
        // below it.
        el.style.height = Math.ceil(naturalHeight * ratio) + 'px';
        // Content now fits post-scale; suppress the otherwise-unnecessary
        // scrollbar (the browser sizes it off the untransformed content).
        el.style.overflowX = 'hidden';
      }
    }
  }

  function run() {
    if (typeof renderMathInElement !== 'function') return;
    renderMathInElement(document.body, {
      delimiters: [
        { left: '\\[', right: '\\]', display: true },
        { left: '\\(', right: '\\)', display: false },
      ],
      throwOnError: false,
    });

    // The KaTeX webfonts load with font-display: swap, so right after
    // render the browser may still be using a fallback font whose glyphs
    // are a different width -- measuring now would size the fit to the
    // wrong (usually narrower) natural width, then silently go stale once
    // the real font swaps in and the equation actually widens. Wait for
    // every requested font to finish loading first so the measurement
    // reflects the width the page will actually settle on.
    if (document.fonts && document.fonts.ready) {
      document.fonts.ready.then(fitDisplayMath).catch(fitDisplayMath);
    } else {
      fitDisplayMath();
    }

    var resizeTimer;
    window.addEventListener('resize', function () {
      clearTimeout(resizeTimer);
      resizeTimer = setTimeout(fitDisplayMath, 150);
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', run);
  } else {
    run();
  }
})();
