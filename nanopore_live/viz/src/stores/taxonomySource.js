/**
 * Global toggle selecting which taxonomy the viz displays — Kraken (NCBI) or
 * GTDB (from the per-read sendsketch classification). Persisted to
 * localStorage so a user's choice survives reload. Default: 'gtdb'.
 *
 * Consumed by stores/data.js (derived sampleTaxonomy / taxonomySunburst) and
 * stores/taxonomy.js (derived rankOrder) — flipping this value
 * reactively rebuilds every taxonomy-derived view (Map rings, Heatmap,
 * Sunburst).
 */
import { writable } from 'svelte/store';

const KEY = 'microscape.taxonomySource';
const VALID = new Set(['kraken', 'gtdb']);

function initial() {
  try {
    const v = localStorage.getItem(KEY);
    if (v && VALID.has(v)) return v;
  } catch (_) { /* localStorage unavailable (SSR / private mode) */ }
  return 'gtdb';
}

function makeTaxonomySource() {
  const store = writable(initial());
  const { subscribe, set, update } = store;
  return {
    subscribe,
    set(v) {
      if (!VALID.has(v)) return;
      set(v);
      try { localStorage.setItem(KEY, v); } catch (_) {}
    },
    update,
  };
}

export const taxonomySource = makeTaxonomySource();
