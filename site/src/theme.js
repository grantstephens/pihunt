// Light/dark theme toggle. Self-hosted, loaded as a plain classic script (no
// inline <script>, no inline event handlers, per this project's CSP).
//
// Default behaviour follows the OS/browser `prefers-color-scheme`, via CSS
// alone (see style.css). This file only handles the explicit override: a
// `data-theme="light"` or `data-theme="dark"` attribute on <html>, set by
// clicking the toggle and remembered in localStorage for next visit.
(function () {
  const STORAGE_KEY = 'pihunt-theme';
  const root = document.documentElement;

  function apply(theme) {
    if (theme === 'light' || theme === 'dark') {
      root.setAttribute('data-theme', theme);
    } else {
      root.removeAttribute('data-theme');
    }
  }

  function readStored() {
    try {
      return localStorage.getItem(STORAGE_KEY);
    } catch {
      return null;
    }
  }

  function writeStored(theme) {
    try {
      if (theme) localStorage.setItem(STORAGE_KEY, theme);
      else localStorage.removeItem(STORAGE_KEY);
    } catch {
      /* private mode / storage disabled: toggle still works for this load */
    }
  }

  function prefersDark() {
    return typeof matchMedia === 'function' && matchMedia('(prefers-color-scheme: dark)').matches;
  }

  function currentEffectiveTheme() {
    const explicit = root.getAttribute('data-theme');
    if (explicit) return explicit;
    return prefersDark() ? 'dark' : 'light';
  }

  function updateButton(button) {
    const effective = currentEffectiveTheme();
    const isDark = effective === 'dark';
    button.setAttribute('aria-pressed', String(isDark));
    const label = button.querySelector('.theme-toggle-label');
    if (label) label.textContent = isDark ? 'Light mode' : 'Dark mode';
  }

  // Apply any stored preference immediately (before DOMContentLoaded would
  // be ideal to avoid a flash, but this file loads at the end of <body> per
  // the CSP's no-inline-script rule, so a brief flash on first paint is the
  // accepted trade-off).
  apply(readStored());

  function init() {
    const button = document.getElementById('theme-toggle');
    if (!button) return;
    updateButton(button);
    button.addEventListener('click', () => {
      const next = currentEffectiveTheme() === 'dark' ? 'light' : 'dark';
      apply(next);
      writeStored(next);
      updateButton(button);
    });
  }

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init);
  } else {
    init();
  }
})();
