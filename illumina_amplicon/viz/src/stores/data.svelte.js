/**
 * Reactive data stores for microscape visualization.
 *
 * Svelte 5 module-level $state can't be exported AND reassigned.
 * Solution: wrap all state in a single exported object whose
 * properties are mutated (not reassigned).
 */

export const store = $state({
  // Core data
  samples: [],
  asvs: [],
  counts: { data: [], samples: [], asvs: [] },
  network: { edges: [] },
  taxonomy: {},
  treeNewick: '',      // phylogeny (NJ) tree
  provenance: null,    // per-step read counts (data/provenance.json)
  heatmapNewick: '',   // Ward clustering tree (from heatmap data)
  heatmap: null,
  aggCounts: null,     // Pre-aggregated counts per taxonomy level

  // Selection
  selectedSample: null,
  selectedAsv: null,

  // UI
  loading: true,
  error: null,
});

/** Group colors as RGBA for regl-scatterplot */
export const GROUP_COLORS = {
  prokaryote: [0.3, 0.5, 1.0, 0.8],
  eukaryote: [1.0, 0.3, 0.3, 0.8],
  chloroplast: [0.2, 0.85, 0.4, 0.8],
  mitochondria: [0.2, 0.9, 0.9, 0.8],
  unknown: [0.6, 0.6, 0.6, 0.5],
};

/** Group colors as hex for UI */
export const GROUP_HEX = {
  prokaryote: '#4d80ff',
  eukaryote: '#ff4d4d',
  chloroplast: '#33d966',
  mitochondria: '#33e6e6',
  unknown: '#999999',
};

// ── Taxonomy coloring ─────────────────────────────────────────────────────

/** Generate N perceptually spaced colors using golden-angle HSL */
function generatePalette(n) {
  const colors = [];
  const golden = 137.508;
  const startHue = 120;  // start at green
  for (let i = 0; i < n; i++) {
    const h = (startHue + i * golden) % 360;
    const s = 55 + (i % 3) * 15;  // 55-85% saturation
    const l = 45 + (i % 4) * 8;   // 45-69% lightness
    colors.push(`hsl(${h}, ${s}%, ${l}%)`);
  }
  return colors;
}

/**
 * Build a color map: taxon name → hex color for the top N taxa at a level.
 * Returns { colorMap: {name: hex}, ranked: [{name, count, color}] }
 */
/**
 * Predicate over ASVs for the current taxonomy filter: (asvId, lineage?) → bool.
 * Returns null when there is no filter (caller should not filter at all).
 *
 * `filterLevel` is set when the filter came from a drill-down click or an
 * autocomplete pick, so the filter names a taxon at that exact rank and is
 * matched against that rank alone. Matching it as a substring of the whole
 * lineage instead makes "SAR" (the eukaryotic supergroup) also select SAR11,
 * SAR86 and SAR116, which sit on a different branch entirely. A hand-typed
 * filter has no level and stays a free regex over the lineage and ASV id.
 */
export function makeTaxonMatcher(taxonFilter = '', filterLevel = '') {
  if (!taxonFilter) return null;
  const db = Object.keys(store.taxonomy)[0];
  const levels = db ? (store.taxonomy[db]?.levels || []) : [];
  const assignments = db ? (store.taxonomy[db]?.assignments || {}) : {};
  const lower = taxonFilter.toLowerCase();

  const levelIdx = filterLevel ? levels.indexOf(filterLevel) : -1;
  if (levelIdx >= 0) {
    return (asvId) => (assignments[asvId]?.[levelIdx] || '').toLowerCase() === lower;
  }

  let re;
  try { re = new RegExp(taxonFilter, 'i'); } catch { re = null; }
  return (asvId, lineage) => {
    const fullTax = lineage ?? (assignments[asvId]?.filter(Boolean).join(';') || '');
    const id = asvId ?? '';
    return re
      ? (re.test(fullTax) || re.test(id))
      : (fullTax.toLowerCase().includes(lower) || id.toLowerCase().includes(lower));
  };
}

/**
 * Taxon names at `level` whose ASVs match `taxonFilter`.
 *
 * Needed because a drill-down filter names an ANCESTOR ("Bacteria") while the
 * pre-aggregated counts are keyed by a DESCENDANT ("Pseudomonadota"). Matching
 * the filter against the child name drops every point. Returns null when there
 * is no filter (caller should not filter at all).
 */
export function taxaMatchingFilter(level, taxonFilter = '', filterLevel = '') {
  if (!taxonFilter || !level || level === '_asv' || level === 'group') return null;
  const db = Object.keys(store.taxonomy)[0];
  if (!db || !store.taxonomy[db]) return null;

  const levels = store.taxonomy[db].levels || [];
  const assignments = store.taxonomy[db].assignments || {};
  const levelIdx = levels.indexOf(level);
  if (levelIdx < 0) return null;

  const matches = makeTaxonMatcher(taxonFilter, filterLevel);
  const out = new Set();
  for (const asvId in assignments) {
    const vals = assignments[asvId] || [];
    if (vals[levelIdx] && matches(asvId)) out.add(vals[levelIdx]);
  }
  return out;
}


// Memo for buildTaxColorMap. It is called many times with the same arguments in one
// render cycle, and each call walks every ASV assignment.
let _taxColorCache = { key: '', result: { colorMap: {}, ranked: [] } };

export function buildTaxColorMap(level, taxonFilter = '', filterLevel = '') {
  const cacheKey = `${level}|${taxonFilter}|${filterLevel}`;
  if (_taxColorCache.key === cacheKey) return _taxColorCache.result;

  const db = Object.keys(store.taxonomy)[0];
  if (!db || !store.taxonomy[db]) return { colorMap: {}, ranked: [] };

  const levels = store.taxonomy[db].levels || [];
  const assignments = store.taxonomy[db].assignments || {};

  // Build a set of ASV IDs that pass the taxonomy filter
  let filteredAsvIds = null;
  const matches = makeTaxonMatcher(taxonFilter, filterLevel);
  if (matches) {
    filteredAsvIds = new Set();
    for (const asvId in assignments) {
      if (matches(asvId)) filteredAsvIds.add(asvId);
    }
  }

  // Special case: color by individual ASV ID
  if (level === '_asv') {
    const asvIds = store.asvs
      .map(a => a.id)
      .filter(id => id && (!filteredAsvIds || filteredAsvIds.has(id)));
    const palette = generatePalette(asvIds.length);
    const colorMap = {};
    // Map, not .find() per element — this ran O(n^2) over every ASV.
    const asvLookup = new Map(store.asvs.map(a => [a.id, a]));
    const ranked = asvIds.map((id, i) => {
      colorMap[id] = palette[i];
      return { name: id, count: asvLookup.get(id)?.total_reads ?? 0, color: palette[i] };
    });
    ranked.sort((a, b) => b.count - a.count);
    const result = { colorMap, ranked };
    _taxColorCache = { key: cacheKey, result };
    return result;
  }

  const levelIdx = levels.indexOf(level);
  if (levelIdx < 0) return { colorMap: {}, ranked: [] };

  // Count ASVs per taxon at this level (only filtered ASVs)
  const counts = {};
  for (const asvId in assignments) {
    if (filteredAsvIds && !filteredAsvIds.has(asvId)) continue;
    const val = assignments[asvId]?.[levelIdx] || 'unclassified';
    counts[val] = (counts[val] || 0) + 1;
  }

  // Rank by count, assign a unique color to every taxon
  const ranked = Object.entries(counts)
    .map(([name, count]) => ({ name, count }))
    .sort((a, b) => b.count - a.count);

  const palette = generatePalette(ranked.length);
  const colorMap = {};
  for (let i = 0; i < ranked.length; i++) {
    ranked[i].color = palette[i];
    colorMap[ranked[i].name] = palette[i];
  }

  const result = { colorMap, ranked };
  _taxColorCache = { key: cacheKey, result };
  return result;
}

/**
 * Determine the effective color level: if the taxon filter matches a taxon
 * at the current level, drill down to the next level to show diversity within.
 */
export function getEffectiveColorLevel(colorByLevel, taxonFilter, filterLevel = '') {
  if (!taxonFilter || colorByLevel === 'group') return colorByLevel;

  const db = Object.keys(store.taxonomy)[0];
  if (!db || !store.taxonomy[db]) return colorByLevel;

  const levels = store.taxonomy[db].levels || [];
  const assignments = store.taxonomy[db].assignments || {};
  const levelIdx = levels.indexOf(colorByLevel);
  if (levelIdx < 0) return colorByLevel;

  // Check if the filter matches a taxon name at the current level
  const taxaAtLevel = new Set();
  for (const asvId in assignments) {
    const val = assignments[asvId]?.[levelIdx];
    if (val) taxaAtLevel.add(val);
  }

  // If filter exactly matches one taxon at this level, drill down
  let matchCount = 0;
  // First try exact string match (handles special chars like parentheses)
  for (const t of taxaAtLevel) {
    if (t.toLowerCase() === taxonFilter.toLowerCase()) matchCount++;
  }
  // If no exact match, try as regex — but only for a hand-typed filter. A picked
  // taxon is an exact name at filterLevel; matching it loosely here re-opens the
  // cross-branch bleed the matcher exists to prevent.
  if (matchCount === 0 && !filterLevel) {
    try {
      const escaped = taxonFilter.replace(/[.*+?^${}()|[\]\\]/g, '\\$&');
      const re = new RegExp(`^${escaped}$`, 'i');
      for (const t of taxaAtLevel) {
        if (re.test(t)) matchCount++;
      }
    } catch {
      return colorByLevel;
    }
  }

  if (matchCount === 1) {
    if (levelIdx < levels.length - 1) {
      // Check if the next level down actually has data for the filtered set
      const nextLevel = levels[levelIdx + 1];
      const nextLevelIdx = levelIdx + 1;
      const matches = makeTaxonMatcher(taxonFilter, filterLevel);
      let hasNextLevelData = false;
      for (const asvId in assignments) {
        if (assignments[asvId]?.[nextLevelIdx] && matches(asvId)) {
          hasNextLevelData = true;
          break;
        }
      }
      if (hasNextLevelData) return nextLevel;
    }
    // At deepest level or next level has no data — color by individual ASV
    return '_asv';
  }

  return colorByLevel;
}

/**
 * Find which taxonomy level a taxon name belongs to.
 * Returns the level name (e.g. 'Phylum') or null if not found.
 */
export function findTaxonLevel(taxonName) {
  const db = Object.keys(store.taxonomy)[0];
  if (!db || !store.taxonomy[db] || !taxonName) return null;

  const levels = store.taxonomy[db].levels || [];
  const assignments = store.taxonomy[db].assignments || {};
  const lower = taxonName.toLowerCase();

  for (let levelIdx = 0; levelIdx < levels.length; levelIdx++) {
    for (const asvId in assignments) {
      const val = assignments[asvId]?.[levelIdx];
      if (val && val.toLowerCase() === lower) {
        return levels[levelIdx];
      }
    }
  }
  return null;
}

/**
 * Get the hex color for an ASV given a taxonomy level and color map.
 */
export function getAsvColor(asvId, level, colorMap) {
  if (level === '_asv') {
    return colorMap[asvId] || '#475569';
  }

  const db = Object.keys(store.taxonomy)[0];
  if (!db || !store.taxonomy[db]) return '#475569';

  const levels = store.taxonomy[db].levels || [];
  const assignments = store.taxonomy[db].assignments || {};
  const levelIdx = levels.indexOf(level);
  if (levelIdx < 0) return '#475569';

  const val = assignments[asvId]?.[levelIdx] || 'unclassified';
  return colorMap[val] || '#475569';
}

/** Convert hex to regl-scatterplot RGBA [0-1] */
export function hexToRegl(hex) {
  if (!hex || hex[0] !== '#') return [0.4, 0.4, 0.4, 0.5];
  const r = parseInt(hex.slice(1, 3), 16) / 255;
  const g = parseInt(hex.slice(3, 5), 16) / 255;
  const b = parseInt(hex.slice(5, 7), 16) / 255;
  return [r, g, b, 0.8];
}

/**
 * Get color for a sample or ASV based on cluster assignment.
 * mode: 'sampleCluster' or 'asvCluster'
 */
export function getClusterColor(id, mode, k) {
  const heatmap = store.heatmap;
  if (!heatmap) return '#475569';

  const clusters = mode === 'sampleCluster'
    ? heatmap.sampleClusters?.[String(k)]
    : heatmap.asvClusters?.[String(k)];

  if (!clusters || !(id in clusters)) return '#475569';

  const cid = clusters[id];
  const hue = ((cid - 1) * 137.508) % 360;
  return `hsl(${hue}, 70%, 55%)`;
}

/** Convert hex to phylocanvas RGBA [0-255] */
export function hexToRgba255(hex) {
  if (!hex || hex[0] !== '#') return [71, 85, 105, 255];
  const r = parseInt(hex.slice(1, 3), 16);
  const g = parseInt(hex.slice(3, 5), 16);
  const b = parseInt(hex.slice(5, 7), 16);
  return [r, g, b, 255];
}

// ── Data loading ────────────────────────────────────────────────────────────

async function fetchJson(url) {
  // Try .gz first
  try {
    const gzRes = await fetch(url + '.gz');
    if (gzRes.ok) {
      const buf = await gzRes.arrayBuffer();
      const bytes = new Uint8Array(buf);
      if (bytes[0] === 0x1f && bytes[1] === 0x8b) {
        const ds = new DecompressionStream('gzip');
        const reader = new Blob([buf]).stream().pipeThrough(ds).getReader();
        const chunks = [];
        while (true) {
          const { done, value } = await reader.read();
          if (done) break;
          chunks.push(value);
        }
        const combined = new Uint8Array(chunks.reduce((a, c) => a + c.length, 0));
        let offset = 0;
        for (const c of chunks) { combined.set(c, offset); offset += c.length; }
        return JSON.parse(new TextDecoder().decode(combined));
      }
      return JSON.parse(new TextDecoder().decode(buf));
    }
  } catch (_) { /* fall through */ }

  const res = await fetch(url);
  if (!res.ok) throw new Error(`${res.status} ${url}`);
  return res.json();
}

export async function loadData() {
  store.loading = true;
  store.error = null;

  try {
    const [samples, asvs, counts, network, taxonomy, treeNewick, provenance] = await Promise.all([
      fetchJson('./data/samples.json').catch(() => []),
      fetchJson('./data/asvs.json').catch(() => []),
      fetchJson('./data/counts.json').catch(() => ({ data: [], samples: [], asvs: [] })),
      fetchJson('./data/network.json').catch(() => ({ edges: [] })),
      fetchJson('./data/taxonomy.json').catch(() => ({})),
      fetch('./data/tree.nwk').then(r => r.ok ? r.text() : '').catch(() => ''),
      fetchJson('./data/provenance.json').catch(() => null),
    ]);

    store.samples = samples;
    store.asvs = asvs;
    _taxColorCache = { key: '', result: { colorMap: {}, ranked: [] } };  // memo depends on these
    store.counts = counts;
    store.network = network;
    store.taxonomy = taxonomy;
    store.treeNewick = treeNewick.trim();
    store.provenance = provenance || null;

    // Load heatmap data (async, non-blocking)
    fetch('./data/heatmap.json.gz')
      .then(r => r.ok ? r.json() : null)
      .then(d => {
        if (d) {
          store.heatmap = d;
          if (d.asvWardNewick) store.heatmapNewick = d.asvWardNewick;
        }
      })
      .catch(() => {});

    // Load pre-aggregated counts
    fetch('./data/aggregated_counts.json.gz')
      .then(r => r.ok ? r.json() : null)
      .then(d => { if (d) store.aggCounts = d; })
      .catch(() => {});
  } catch (e) {
    store.error = e.message;
    console.error('[microscape] Data load failed:', e);
  }

  store.loading = false;
}

// ── Helpers ─────────────────────────────────────────────────────────────────

export function countsBySample() {
  const map = new Map();
  const data = store.counts?.data;
  const sampleIds = store.counts?.samples;
  if (!data || !sampleIds) return map;
  for (const row of data) {
    const [si, ai, count, prop] = row;
    const sampleId = sampleIds[si];
    if (!map.has(sampleId)) map.set(sampleId, []);
    map.get(sampleId).push({ asv_idx: ai, count, proportion: prop });
  }
  return map;
}

export function countsByAsv() {
  const map = new Map();
  const data = store.counts?.data;
  if (!data) return map;
  for (const row of data) {
    const [si, ai, count, prop] = row;
    if (!map.has(ai)) map.set(ai, []);
    map.get(ai).push({ sample_idx: si, count, proportion: prop });
  }
  return map;
}
