import { writable } from 'svelte/store';

// index.json lists the references; each has four tables written by
// bin/summarize_expression.py: samples, bins, functions, top_genes.
export const index    = writable(null);
export const refName  = writable(null);
export const samples  = writable(null);  // { samples, rows[], n_genes, n_contigs, genes_detected[] }
export const bins     = writable(null);  // { samples, bins: [{bin, length, n_genes, completeness, contamination, classification, reads[], tpm[]}] }
export const funcs    = writable(null);  // { samples, kinds: { ko|cog|cazy: [{id, tpm[]}] } }
export const topGenes = writable(null);  // { samples, genes: [{id, contig, product, gene, bin, ko, reads[], tpm[]}] }

export const loading = writable(true);
export const error   = writable(null);

async function fetchJSON(url) {
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Failed to load ${url}: ${res.status}`);
  // .json.gz is served either with Content-Encoding: gzip (already inflated
  // by the browser) or as raw bytes; peek at the magic number to tell.
  const buf = new Uint8Array(await res.arrayBuffer());
  if (buf[0] === 0x1f && buf[1] === 0x8b) {
    const ds = new DecompressionStream('gzip');
    const text = await new Response(new Blob([buf]).stream().pipeThrough(ds)).text();
    return JSON.parse(text);
  }
  return JSON.parse(new TextDecoder().decode(buf));
}

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
  let idx;
  index.subscribe(v => (idx = v))();
  const ref = idx?.references?.find(r => r.name === name);
  if (!ref || current === name) return;
  current = name;
  loading.set(true);
  try {
    const get = k => ref.files[k] ? fetchJSON(`./data/${ref.files[k]}`) : Promise.resolve(null);
    const [s, b, f, g] = await Promise.all([get('samples'), get('bins'), get('functions'), get('top_genes')]);
    samples.set(s); bins.set(b); funcs.set(f); topGenes.set(g);
    refName.set(name);
  } catch (e) {
    error.set(e.message);
  } finally {
    loading.set(false);
  }
}

// ── shared helpers ──────────────────────────────────────────────────────────
export const mean = a => (a.length ? a.reduce((x, y) => x + y, 0) / a.length : 0);
export const median = a => {
  const v = a.filter(x => x != null && !Number.isNaN(x)).sort((x, y) => x - y);
  return v.length ? v[Math.floor(v.length / 2)] : null;
};
export const fmtInt = n => (n == null ? '—' : Math.round(n).toLocaleString());
export const fmtPct = n => (n == null ? '—' : `${n.toFixed(1)}%`);
export const fmtTpm = n => (n == null ? '—' : n >= 100 ? Math.round(n).toLocaleString() : n.toFixed(1));

// Sequential blue for magnitude on the dark surface: near-zero recedes into
// the surface, high values step toward the light end.
export const SEQ_BLUE = [
  [0.0, '#0f172a'], [0.12, '#0d366b'], [0.3, '#1c5cab'],
  [0.5, '#2a78d6'], [0.7, '#5598e7'], [0.85, '#86b6ef'], [1.0, '#cde2fb'],
];
