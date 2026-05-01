import { writable, derived } from 'svelte/store';

// Set of sample keys ("flowcell:barcode")
export const cartItems = writable(new Set());

// When true, all views filter to only cart items
export const cartActive = writable(false);

export function addToCart(sampleId) {
  cartItems.update(s => { s.add(sampleId); return new Set(s); });
}

export function removeFromCart(sampleId) {
  cartItems.update(s => { s.delete(sampleId); return new Set(s); });
}

export function toggleCart(sampleId) {
  cartItems.update(s => {
    if (s.has(sampleId)) s.delete(sampleId);
    else s.add(sampleId);
    return new Set(s);
  });
}

export function clearCart() {
  cartItems.set(new Set());
}

export const cartCount = derived(cartItems, $c => $c.size);

// ---- Read selections ----
// Map<selectionName, { readIds: Set<readId>, samples: Map<sampleId, count> }>
// Per-sample counts are pre-aggregated at add time so /map can iterate samples
// rather than the (much larger) read list on every render.
export const readSelections = writable(new Map());

// Which selection is currently being visualized on /map. null = none active.
// Default is set to the most-recently-added selection by addReadSelection().
export const activeReadSelection = writable(null);

let _selectionCounter = 0;

export function addReadSelection(readIds, allReads) {
  if (!readIds || !allReads) return null;
  const ids = readIds instanceof Set ? readIds : new Set(readIds);
  if (ids.size === 0) return null;

  // Pre-aggregate per-sample counts by iterating the full read list once.
  const samples = new Map();
  for (const r of allReads) {
    if (ids.has(r.id)) {
      samples.set(r.sample, (samples.get(r.sample) || 0) + 1);
    }
  }

  _selectionCounter += 1;
  const name = `Selection ${_selectionCounter}`;
  readSelections.update(m => {
    const next = new Map(m);
    next.set(name, { readIds: ids, samples });
    return next;
  });
  // Auto-activate the newly added selection.
  activeReadSelection.set(name);
  return name;
}

export function removeReadSelection(name) {
  readSelections.update(m => {
    const next = new Map(m);
    next.delete(name);
    return next;
  });
  activeReadSelection.update(active => {
    if (active === name) {
      // Pick any remaining selection, or null.
      let fallback = null;
      readSelections.update(m => {
        for (const k of m.keys()) { fallback = k; break; }
        return m;
      });
      return fallback;
    }
    return active;
  });
}

export function clearReadSelections() {
  readSelections.set(new Map());
  activeReadSelection.set(null);
}

export const readSelectionCount = derived(readSelections, $r => $r.size);
