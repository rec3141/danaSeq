/**
 * Persisted /heatmap UI settings — rank / column order / label / value mode /
 * top-N / k-clusters. Lifted out of HeatmapView's $state so the user's choices
 * survive switching tabs and reload (localStorage). Defaults match the prior
 * inline initial values.
 */
import { writable } from 'svelte/store';

const KEY = 'microscape.heatmapSettings';

const DEFAULTS = {
  rankIdx: 4,        // 0..5 — Phylum/Class/Order/Family/Genus/Species, default Genus
  orderIdx: 0,       // 0..(orderModes-1), 0 = default sort
  labelIdx: 0,       // 0..(labelModes-1)
  rowOrderIdx: 0,    // 0=relabund, 1=taxonomic, 2=ward
  valueMode: 'pct',  // 'pct' | 'log'
  topN: 30,          // top-N taxa shown
  kClusters: 4,      // Ward k-cut
};

function initial() {
  try {
    const v = localStorage.getItem(KEY);
    if (!v) return { ...DEFAULTS };
    const parsed = JSON.parse(v);
    return { ...DEFAULTS, ...parsed };
  } catch (_) { return { ...DEFAULTS }; }
}

function persist(state) {
  try { localStorage.setItem(KEY, JSON.stringify(state)); } catch (_) {}
}

function makeHeatmapSettings() {
  const store = writable(initial());
  const { subscribe, update } = store;
  return {
    subscribe,
    /** Patch one or more keys. */
    patch(diff) {
      update(s => {
        const next = { ...s, ...diff };
        persist(next);
        return next;
      });
    },
    reset() {
      const d = { ...DEFAULTS };
      persist(d);
      update(() => d);
    },
  };
}

export const heatmapSettings = makeHeatmapSettings();
