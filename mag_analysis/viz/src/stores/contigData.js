// Per-contig data: contig explorer rows (+ t-SNE/UMAP), per-sample depths and
// gene features.
//
// Runs written by the current preprocessor store these chunked and columnar
// (see preprocess/viz_chunks.py): a manifest per data set plus part files, with
// contigs ordered by length (longest first). Only the chunks needed for the
// requested minimum contig length are fetched; lowering the minimum fetches
// more. Older runs have single files (contig_explorer.json, contig_tsne.json,
// contig_umap.json, contig_sample_depths.json, genes.json); when a manifest is
// absent those are loaded whole, as before.
//
// Rows keep the legacy shape as far as consumers can tell: row.id, row.length,
// row.kaiju_phylum, ... Chunked rows read their fields through prototype
// getters over typed-array columns, so a million rows do not each carry forty
// properties. Every row also has _row, its position in the full length-sorted
// order, which indexes per-sample depth columns.

import { writable, get } from 'svelte/store';
import { fetchJSON } from './fetchJson.js';

const DATA = 'data/';

// The Contigs view starts at this minimum length; shorter contigs load on demand.
export const DEFAULT_MIN_CONTIG_LENGTH = 10000;

export const contigExplorer = writable(null);
export const sampleDepths = writable(null);
export const contigLoadState = writable({ loading: false, done: 0, total: 0, error: null });
export const geneSearch = writable({
  index: null, version: 0, loading: false, shardsDone: 0, shardsNeeded: 0,
  shardsTotal: 0, minLength: null, legacy: false,
});

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

async function tryManifest(name, format) {
  try {
    const m = await fetchJSON(DATA + name);
    return m && m.format === format ? m : null;
  } catch (e) {
    return null;
  }
}

async function mapLimit(items, limit, fn) {
  let next = 0;
  const workers = Array.from({ length: Math.min(limit, items.length) }, async () => {
    while (next < items.length) {
      const item = items[next++];
      await fn(item);
    }
  });
  await Promise.all(workers);
}

function decodeIds(encoded, enc) {
  if (!enc) return encoded;
  const out = new Array(encoded.length);
  const { prefix, width } = enc;
  let n = 0;
  for (let i = 0; i < encoded.length; i++) {
    n += encoded[i];
    out[i] = prefix + String(n).padStart(width, '0');
  }
  return out;
}

// Sort order shared by every chunked data set: length descending, then id.
function compareKey(lenA, idA, lenB, idB) {
  if (lenA !== lenB) return lenB - lenA;
  return idA < idB ? -1 : idA > idB ? 1 : 0;
}

// Serialize loads so overlapping requests never fetch a chunk twice.
let queue = Promise.resolve();
function enqueue(fn) {
  const p = queue.then(fn);
  queue = p.catch(() => {});
  return p;
}

// ---------------------------------------------------------------------------
// Contig explorer
// ---------------------------------------------------------------------------

let explorerManifest;          // undefined: not probed; null: legacy layout
let explorerProbe = null;
const chunkRows = [];          // decoded rows per chunk index
let loadedChunks = 0;          // chunks 0..loadedChunks-1 are loaded
let legacyLoaded = false;

function probeExplorer() {
  if (!explorerProbe) {
    explorerProbe = tryManifest('contig_explorer_manifest.json', 'danaseq-contig-explorer')
      .then(m => { explorerManifest = m; return m; });
  }
  return explorerProbe;
}

export function decodeExplorerChunk(man, body) {
  const vocab = man.vocab || {};
  const cols = body.columns;
  const props = {};
  for (const name of man.column_order) {
    const col = cols[name];
    if (!col) continue;
    if (vocab[name]) {
      const voc = vocab[name];
      const arr = voc.length <= 65536 ? Uint16Array.from(col) : Uint32Array.from(col);
      props[name] = { get() { return voc[arr[this._i]]; }, enumerable: true };
    } else {
      const arr = new Float64Array(col.length);
      for (let i = 0; i < col.length; i++) arr[i] = col[i] == null ? NaN : col[i];
      props[name] = {
        get() { const v = arr[this._i]; return v === v ? v : null; },
        enumerable: true,
      };
    }
  }
  function Row(id, i, row) { this.id = id; this._i = i; this._row = row; }
  Object.defineProperties(Row.prototype, props);
  const ids = decodeIds(body.ids, man.id_encoding);
  const rows = new Array(ids.length);
  for (let i = 0; i < ids.length; i++) rows[i] = new Row(ids[i], i, body.offset + i);
  return rows;
}

// Number of leading chunks holding every contig of at least minLength.
function chunksNeeded(man, minLength) {
  const chunks = man.chunks;
  let need = 0;
  while (need < chunks.length && chunks[need].max_length >= minLength) need++;
  return Math.max(need, Math.min(1, chunks.length));
}

function publishExplorer() {
  const man = explorerManifest;
  const contigs = [].concat(...chunkRows.slice(0, loadedChunks));
  const total = man.chunks.length;
  // Every contig at least this long is loaded (ties can straddle a chunk edge).
  const nextMax = loadedChunks < total ? man.chunks[loadedChunks].max_length : 0;
  contigExplorer.set({
    contigs,
    chunked: true,
    complete: loadedChunks >= total,
    n_total: man.n_contigs,
    loaded_min_length: loadedChunks >= total ? 0 : nextMax + 1,
    has_tsne: !!man.has_tsne,
    has_umap: !!man.has_umap,
    pca_variance_explained: man.pca_variance_explained,
    transform: man.transform,
    extents: man.extents,
    vocab: man.vocab,
    bin_min_length: man.bin_min_length,
  });
}

async function ensureExplorerChunks(need) {
  const man = explorerManifest;
  if (need <= loadedChunks) return;
  const todo = [];
  for (let k = loadedChunks; k < need; k++) if (!chunkRows[k]) todo.push(k);
  let done = loadedChunks + (need - loadedChunks - todo.length);
  contigLoadState.set({ loading: true, done, total: need, error: null });
  try {
    await mapLimit(todo, 4, async (k) => {
      const body = await fetchJSON(DATA + man.chunks[k].file);
      chunkRows[k] = decodeExplorerChunk(man, body);
      done++;
      contigLoadState.set({ loading: true, done, total: need, error: null });
    });
  } catch (e) {
    console.error('Failed to load contig chunk:', e);
    contigLoadState.set({ loading: false, done, total: need, error: e.message });
  }
  // Publish the contiguous loaded prefix, even after a partial failure.
  let k = loadedChunks;
  while (k < need && chunkRows[k]) k++;
  if (k > loadedChunks) {
    loadedChunks = k;
    publishExplorer();
  }
  if (get(contigLoadState).loading) {
    contigLoadState.set({ loading: false, done: loadedChunks, total: need, error: null });
  }
  if (depthsWanted) await ensureDepthChunks();
}

async function loadLegacyExplorer() {
  if (legacyLoaded) return;
  legacyLoaded = true;
  contigLoadState.set({ loading: true, done: 0, total: 1, error: null });
  try {
    const [data, tsneEmb, umapEmb] = await Promise.all([
      fetchJSON(DATA + 'contig_explorer.json'),
      fetchJSON(DATA + 'contig_tsne.json').catch(() => null),
      fetchJSON(DATA + 'contig_umap.json').catch(() => null),
    ]);
    // Merge embeddings into contig records (each file is {contig_id: [x, y]})
    let nTsne = 0, nUmap = 0;
    const contigs = data?.contigs || [];
    contigs.forEach((c, i) => {
      c._row = i;
      if (tsneEmb?.[c.id]) { c.tsne_x = tsneEmb[c.id][0]; c.tsne_y = tsneEmb[c.id][1]; nTsne++; }
      if (umapEmb?.[c.id]) { c.umap_x = umapEmb[c.id][0]; c.umap_y = umapEmb[c.id][1]; nUmap++; }
    });
    contigExplorer.set({
      ...data,
      contigs,
      chunked: false,
      complete: true,
      n_total: contigs.length,
      loaded_min_length: 0,
      has_tsne: nTsne > 0,
      has_umap: nUmap > 0,
    });
    contigLoadState.set({ loading: false, done: 1, total: 1, error: null });
  } catch (e) {
    console.error('Failed to load contig explorer:', e);
    contigLoadState.set({ loading: false, done: 0, total: 1, error: e.message });
  }
}

// Load contigs of at least minLength (all contigs for legacy runs, which have
// a single file). Resolves once they are in the contigExplorer store.
export function loadContigExplorer(minLength = DEFAULT_MIN_CONTIG_LENGTH) {
  return enqueue(async () => {
    const man = await probeExplorer();
    if (!man) return loadLegacyExplorer();
    return ensureExplorerChunks(chunksNeeded(man, minLength));
  });
}

// Load enough contigs to include every member of the named bin (any binner).
export function ensureContigsForBin(binName) {
  return enqueue(async () => {
    const man = await probeExplorer();
    if (!man || !binName) return;
    let minLen = Infinity;
    for (const mins of Object.values(man.bin_min_length || {})) {
      if (mins[binName] != null && mins[binName] < minLen) minLen = mins[binName];
    }
    if (minLen < Infinity) return ensureExplorerChunks(chunksNeeded(man, minLen));
  });
}

// ---------------------------------------------------------------------------
// Per-sample depths
// ---------------------------------------------------------------------------

let depthsWanted = false;
let depthManifest;             // undefined: not probed; null: absent
const depthChunks = [];        // per chunk: [{ rows: Int32Array, vals: Float32Array } per sample]

function depthsEmpty() {
  return { samples: [], counts: {}, maxDepth: 0, depthOf: () => 0 };
}

function publishDepths() {
  const dm = depthManifest;
  const samples = dm.samples;
  const counts = {};
  const loaded = [];
  for (let k = 0; k < depthChunks.length; k++) if (depthChunks[k]) loaded.push(k);
  samples.forEach((name, s) => {
    let n = 0;
    for (const k of loaded) n += dm.chunks[k].nonzero[s];
    counts[name] = n;
  });
  const nTotal = explorerManifest.n_contigs;
  const cache = new Map();
  function column(s) {
    let col = cache.get(s);
    if (!col) {
      col = new Float32Array(nTotal);
      for (const k of loaded) {
        const { rows, vals } = depthChunks[k][s];
        for (let j = 0; j < rows.length; j++) col[rows[j]] = vals[j];
      }
      cache.clear();  // keep one sample column (nTotal floats) at a time
      cache.set(s, col);
    }
    return col;
  }
  sampleDepths.set({
    samples,
    maxDepth: dm.maxDepth,
    counts,
    depthOf: (c, s) => (s >= 0 && c._row != null ? column(s)[c._row] : 0),
  });
}

async function ensureDepthChunks() {
  if (depthManifest === undefined) {
    depthManifest = await tryManifest('contig_sample_depths_manifest.json', 'danaseq-contig-sample-depths');
  }
  const dm = depthManifest;
  if (!dm) {
    if (get(sampleDepths) === null) sampleDepths.set(depthsEmpty());
    return;
  }
  const todo = [];
  for (let k = 0; k < loadedChunks; k++) if (!depthChunks[k]) todo.push(k);
  if (!todo.length && get(sampleDepths) !== null) return;
  await mapLimit(todo, 4, async (k) => {
    try {
      const body = await fetchJSON(DATA + dm.chunks[k].file);
      depthChunks[k] = body.idx.map((idx, s) => {
        const rows = new Int32Array(idx.length);
        let r = body.offset;
        for (let j = 0; j < idx.length; j++) { r += idx[j]; rows[j] = r; }
        return { rows, vals: Float32Array.from(body.val[s]) };
      });
    } catch (e) {
      console.warn('Per-sample depth chunk not available:', e.message);
    }
  });
  publishDepths();
}

async function loadLegacySampleDepths() {
  try {
    const data = await fetchJSON(DATA + 'contig_sample_depths.json');
    const samples = data.samples || [];
    const depths = data.depths || {};
    const counts = {};
    const nz = new Array(samples.length).fill(0);
    for (const vals of Object.values(depths)) {
      for (let s = 0; s < vals.length; s++) if (vals[s] > 0) nz[s]++;
    }
    samples.forEach((name, s) => { counts[name] = nz[s]; });
    sampleDepths.set({
      samples,
      maxDepth: data.maxDepth,
      counts,
      depthOf: (c, s) => depths[c.id]?.[s] ?? 0,
    });
  } catch (e) {
    console.warn('Per-sample depth data not available:', e.message);
    sampleDepths.set(depthsEmpty());
  }
}

// Load per-sample depths for the loaded contigs; contigs loaded later get
// theirs automatically. The store value is
//   { samples, maxDepth, counts: {sample: nonzero contigs}, depthOf(row, sampleIndex) }.
export function loadSampleDepths() {
  if (depthsWanted) return queue;
  depthsWanted = true;
  return enqueue(async () => {
    const man = await probeExplorer();
    if (!man) return loadLegacySampleDepths();
    return ensureDepthChunks();
  });
}

// ---------------------------------------------------------------------------
// Genes
// ---------------------------------------------------------------------------

let genesManifest;             // undefined: not probed; null: legacy genes.json
let genesProbe = null;
let legacyGenesPromise = null;
const shardCache = new Map();  // shard index -> promise of {contig_id: [feature]}
const SHARD_CACHE_MAX = 3;

function probeGenes() {
  if (!genesProbe) {
    genesProbe = tryManifest('genes_manifest.json', 'danaseq-genes')
      .then(m => { genesManifest = m; return m; });
  }
  return genesProbe;
}

function loadLegacyGenes() {
  if (!legacyGenesPromise) {
    legacyGenesPromise = fetchJSON(DATA + 'genes.json').catch(e => {
      console.warn('Gene data not available:', e.message);
      return {};
    });
  }
  return legacyGenesPromise;
}

function fetchShard(k) {
  return fetchJSON(DATA + genesManifest.shards[k].file);
}

// Shards viewed in the detail panel stay cached (a few, least recently used).
function getShard(k) {
  let p = shardCache.get(k);
  if (p) {
    shardCache.delete(k);
    shardCache.set(k, p);
    return p;
  }
  p = fetchShard(k).catch(e => { shardCache.delete(k); throw e; });
  shardCache.set(k, p);
  while (shardCache.size > SHARD_CACHE_MAX) shardCache.delete(shardCache.keys().next().value);
  return p;
}

// Index of the shard holding the contig with this (length, id), or -1.
function findShard(man, length, id) {
  const shards = man.shards;
  let lo = 0, hi = shards.length;
  while (lo < hi) {
    const mid = (lo + hi) >> 1;
    const [l, i] = shards[mid].last;
    if (compareKey(l, i, length, id) < 0) lo = mid + 1;
    else hi = mid;
  }
  if (lo >= shards.length) return -1;
  const [fl, fi] = shards[lo].first;
  return compareKey(fl, fi, length, id) <= 0 ? lo : -1;
}

// Gene features of one contig. Takes an explorer row (or anything with id and
// length); a bare id works only for legacy runs, whose genes.json has no shards.
export async function loadContigGenes(contig) {
  const id = typeof contig === 'string' ? contig : contig?.id;
  if (!id) return null;
  const man = await probeGenes();
  if (!man) return (await loadLegacyGenes())[id] || null;
  const length = typeof contig === 'string' ? null : contig.length;
  if (length == null) return null;
  const k = findShard(man, length, id);
  if (k < 0) return null;
  try {
    return (await getShard(k))[id] || null;
  } catch (e) {
    console.warn('Gene shard not available:', e.message);
    return null;
  }
}

// Search index: Map<contigId, lowercased gene names + products joined by \0>
function indexGenes(index, genesByContig) {
  for (const [contigId, genes] of Object.entries(genesByContig)) {
    const parts = [];
    for (const gene of genes) {
      if (gene.g) parts.push(gene.g);
      if (gene.p) parts.push(gene.p);
    }
    if (parts.length) index.set(contigId, parts.join('\0').toLowerCase());
  }
}

const searchIndex = new Map();
const searchedShards = new Set();
let legacySearchBuilt = false;

function setSearch(patch) {
  const cur = get(geneSearch);
  geneSearch.set({ ...cur, ...patch, index: searchIndex, version: cur.version + 1 });
}

// Build the gene search index for the contigs currently loaded in the explorer:
// the gene shards covering those contigs are fetched one at a time and indexed
// as they arrive. Call again after more contigs load to extend it. Legacy runs
// index their whole genes.json.
export function ensureGeneSearchIndex() {
  return enqueue(async () => {
    const man = await probeGenes();
    if (!man) {
      if (legacySearchBuilt) return;
      setSearch({ loading: true, legacy: true });
      indexGenes(searchIndex, await loadLegacyGenes());
      legacySearchBuilt = true;
      setSearch({ loading: false, shardsDone: 1, shardsNeeded: 1, shardsTotal: 1, minLength: 0 });
      return;
    }
    const exp = get(contigExplorer);
    const shards = man.shards;
    let need = shards.length;
    let minLength = 0;
    if (exp?.chunked && !exp.complete && exp.contigs.length) {
      const last = exp.contigs[exp.contigs.length - 1];
      need = 0;
      while (need < shards.length &&
             compareKey(shards[need].first[0], shards[need].first[1], last.length, last.id) <= 0) need++;
      minLength = exp.loaded_min_length;
    }
    const todo = [];
    for (let k = 0; k < need; k++) if (!searchedShards.has(k)) todo.push(k);
    setSearch({ loading: todo.length > 0, shardsDone: need - todo.length, shardsNeeded: need,
                shardsTotal: shards.length, minLength });
    for (const k of todo) {
      try {
        // Not via the detail-panel cache: indexed shards are dropped after use.
        indexGenes(searchIndex, shardCache.has(k) ? await shardCache.get(k) : await fetchShard(k));
      } catch (e) {
        console.warn('Gene shard not available:', e.message);
      }
      searchedShards.add(k);
      setSearch({ shardsDone: get(geneSearch).shardsDone + 1 });
    }
    setSearch({ loading: false });
  });
}
