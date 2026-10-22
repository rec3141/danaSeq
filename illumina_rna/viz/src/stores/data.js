import { writable, get } from 'svelte/store';
import { fetchJSON } from './fetchJson.js';

// index.json lists the references; each has tables written by
// bin/summarize_expression.py: samples, bins, functions, top_genes, and
// (with a contig embedding) contig_map + one contig_expr layer per sample.
export const index     = writable(null);
export const refName   = writable(null);
export const samples   = writable(null);  // { samples, rows[], n_genes, n_contigs, genes_detected[] }
export const bins      = writable(null);  // { samples, bins: [{bin, length, n_genes, completeness, contamination, classification, reads[], tpm[]}] }
export const funcs     = writable(null);  // { samples, kinds: { ko|cog|cazy: [{id, tpm[]}] } }
export const topGenes  = writable(null);  // { samples, genes: [{id, contig, product, gene, bin, ko, reads[], tpm[]}] }
export const contigMap = writable(null);  // { samples, id[], x[], y[], length[], bin[], bins[], total_reads[], ... }

export const loading = writable(true);
export const error   = writable(null);

// Files are listed with their .json.gz name; fetchJSON wants the base name
// and tries <name>.gz first.
const base = f => `./data/${f.replace(/\.gz$/, '')}`;

export async function loadIndex() {
  loading.set(true);
  error.set(null);
  try {
    const idx = await fetchJSON('./data/index.json');
    index.set(idx);
    const first = idx.references?.[0]?.name;
    if (first) await loadReference(first);
  } catch (e) {
    error.set(e.message);
    console.error('Data load error:', e);
  } finally {
    loading.set(false);
  }
}

let current = null;
export async function loadReference(name) {
  const ref = get(index)?.references?.find(r => r.name === name);
  if (!ref || current === name) return;
  current = name;
  loading.set(true);
  try {
    const load = k => ref.files[k] ? fetchJSON(base(ref.files[k])) : Promise.resolve(null);
    const [s, b, f, g] = await Promise.all([load('samples'), load('bins'), load('functions'), load('top_genes')]);
    samples.set(s); bins.set(b); funcs.set(f); topGenes.set(g);
    contigMap.set(null);
    layerCache.clear();
    refName.set(name);
  } catch (e) {
    error.set(e.message);
  } finally {
    loading.set(false);
  }
}

// ── Contig map: loaded when the Explorer opens; per-sample layers on demand ──
const layerCache = new Map();
export const contigMapState = writable({ loading: false, error: null, available: true });

export async function loadContigMap() {
  const ref = get(index)?.references?.find(r => r.name === get(refName));
  if (!ref?.files?.contig_map) { contigMapState.set({ loading: false, error: null, available: false }); return; }
  if (get(contigMap)) return;
  contigMapState.set({ loading: true, error: null, available: true });
  try {
    const m = await fetchJSON(base(ref.files.contig_map));
    contigMap.set(m);
    contigMapState.set({ loading: false, error: null, available: true });
  } catch (e) {
    contigMapState.set({ loading: false, error: e.message, available: true });
  }
}

// Reads per map contig for sample j, as a Float32Array aligned to the map.
export async function loadContigLayer(map, j) {
  if (!layerCache.has(j)) {
    const raw = await fetchJSON(base(map.files[j]));
    const reads = new Float32Array(map.n);
    let pos = 0;
    raw.idx.forEach((d, k) => { pos += d; reads[pos] = raw.reads[k]; });  // delta-decoded
    layerCache.set(j, reads);
  }
  return layerCache.get(j);
}

// ── Shared helpers ───────────────────────────────────────────────────────────
export const mean = a => (a.length ? a.reduce((x, y) => x + y, 0) / a.length : 0);
export const median = a => {
  const v = a.filter(x => x != null && !Number.isNaN(x)).sort((x, y) => x - y);
  return v.length ? v[Math.floor(v.length / 2)] : null;
};
export const fmtInt = n => (n == null ? '-' : Math.round(n).toLocaleString());
export const fmtPct = (n, d = 1) => (n == null ? '-' : `${n.toFixed(d)}%`);
export const fmtTpm = n => (n == null ? '-' : n >= 100 ? Math.round(n).toLocaleString() : n.toFixed(1));
export function fmtBp(v) {
  if (v == null) return '-';
  if (v >= 1e9) return (v / 1e9).toFixed(2) + ' Gb';
  if (v >= 1e6) return (v / 1e6).toFixed(1) + ' Mb';
  if (v >= 1e3) return (v / 1e3).toFixed(1) + ' Kb';
  return Math.round(v) + ' bp';
}

// Same tiers as the MAG dashboard
export function qualityTier(b) {
  if (b?.completeness == null) return null;
  if (b.completeness >= 90 && b.contamination < 5) return 'HQ';
  if (b.completeness >= 50 && b.contamination < 10) return 'MQ';
  return 'LQ';
}

// GTDB lineage helpers ("d__Bacteria;p__Pseudomonadota;...")
const RANK_PREFIX = { d: 'domain', p: 'phylum', c: 'class', o: 'order', f: 'family', g: 'genus', s: 'species' };
export function gtdbRanks(c) {
  const out = {};
  for (const part of (c || '').split(';')) {
    const [p, v] = part.split('__');
    if (RANK_PREFIX[p] && v) out[RANK_PREFIX[p]] = v;
  }
  return out;
}
export const lastRank = c => {
  const r = gtdbRanks(c);
  return r.species || r.genus || r.family || r.order || r.class || r.phylum || r.domain || null;
};

export function downloadTsv(filename, header, rows) {
  const body = rows.map(r => r.map(v => (v == null ? '' : v)).join('\t')).join('\n');
  const blob = new Blob([header.join('\t') + '\n' + body + '\n'], { type: 'text/tab-separated-values' });
  const url = URL.createObjectURL(blob);
  const a = document.createElement('a');
  a.href = url; a.download = filename; a.click();
  URL.revokeObjectURL(url);
}
