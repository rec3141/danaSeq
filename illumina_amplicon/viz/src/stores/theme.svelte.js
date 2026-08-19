/**
 * Light and dark, for the page and for the plots.
 *
 * The page's own colours come from CSS variables that `.dark` on <html> swaps,
 * so Tailwind classes need no help. Plotly, phylocanvas and the hand-drawn
 * heatmap take colours as values rather than classes, and have to be handed the
 * palette; `chrome` is that palette, and reading it inside the effect that draws
 * a plot is enough to make the plot redraw when the theme changes.
 *
 * Only chrome lives here — backgrounds, axes, gridlines, unhighlighted marks.
 * The colours that carry data (taxon palettes, cluster colours, the selection
 * highlight) mean the same thing in both themes and are not themed.
 */

const STORAGE_KEY = 'danaseq-viz-theme';

const CHROME = {
  dark: {
    plotBg: 'rgba(2, 6, 15, 1)',
    canvasBg: '#0f172a',
    axis: '#94a3b8',
    faint: '#64748b',
    grid: '#1e293b',
    cell: '#334155',      // a heatmap cell with nothing to say
    masked: '#0f172a',    // a cell filtered out
    rule: '#475569',      // dendrogram branches
    node: [148, 163, 184, 255],
    stroke: [71, 85, 105, 255],
    label: [100, 116, 139, 255],
    highlight: [220, 225, 235, 180],
    collapsed: '#1e293b',
  },
  light: {
    plotBg: 'rgba(255, 255, 255, 1)',
    canvasBg: '#ffffff',
    axis: '#475569',
    faint: '#64748b',
    grid: '#e2e8f0',
    cell: '#cbd5e1',
    masked: '#f8fafc',
    rule: '#94a3b8',
    node: [71, 85, 105, 255],
    stroke: [148, 163, 184, 255],
    label: [100, 116, 139, 255],
    highlight: [30, 41, 59, 180],
    collapsed: '#cbd5e1',
  },
};

function initial() {
  try {
    const saved = localStorage.getItem(STORAGE_KEY);
    if (saved === 'light' || saved === 'dark') return saved;
  } catch { /* private browsing, or no storage at all */ }
  return window.matchMedia?.('(prefers-color-scheme: light)').matches ? 'light' : 'dark';
}

export const theme = $state({ mode: initial() });

/** The palette for whatever is drawing itself rather than wearing a class. */
export function chrome() {
  return CHROME[theme.mode] ?? CHROME.dark;
}

export function setTheme(mode) {
  theme.mode = mode === 'light' ? 'light' : 'dark';
  document.documentElement.classList.toggle('dark', theme.mode === 'dark');
  try { localStorage.setItem(STORAGE_KEY, theme.mode); } catch { /* nothing to do */ }
}

export function toggleTheme() {
  setTheme(theme.mode === 'dark' ? 'light' : 'dark');
}

// The inline script in index.html sets the class before first paint; this keeps
// the class and the store agreeing when the module loads.
setTheme(theme.mode);
