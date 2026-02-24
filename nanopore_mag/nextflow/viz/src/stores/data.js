import { writable, get } from 'svelte/store';

// Individual data stores
export const overview = writable(null);
export const mags = writable(null);
export const checkm2All = writable(null);
export const taxonomySunburst = writable(null);
export const keggHeatmap = writable(null);
export const scgHeatmap = writable(null);
export const coverage = writable(null);
export const mgeSummary = writable(null);
export const mgePerBin = writable(null);
export const eukaryotic = writable(null);
export const contigLengths = writable(null);
// Lazy-loaded (large files)
export const contigExplorer = writable(null);

export const loading = writable(true);
export const error = writable(null);

async function fetchJSON(url) {
  // Try pre-compressed .json.gz first — works with any static server, no server config needed.
  // DecompressionStream is supported in all modern browsers (Chrome 80+, Firefox 113+, Safari 16.4+).
  const gzRes = await fetch(url + '.gz');
  if (gzRes.ok) {
    const ds = new DecompressionStream('gzip');
    const text = await new Response(gzRes.body.pipeThrough(ds)).text();
    return JSON.parse(text);
  }
  // Fall back to plain JSON (older preprocess runs, dev environments)
  const res = await fetch(url);
  if (!res.ok) throw new Error(`Failed to load ${url}: ${res.status}`);
  return res.json();
}

export async function loadAllData() {
  loading.set(true);
  error.set(null);
  try {
    const [
      overviewData,
      magsData,
      sunburstData,
      keggData,
      scgData,
      coverageData,
      mgeSumData,
      mgePerBinData,
      eukData,
      contigLenData,
    ] = await Promise.all([
      fetchJSON('/data/overview.json'),
      fetchJSON('/data/mags.json'),
      fetchJSON('/data/taxonomy_sunburst.json'),
      fetchJSON('/data/kegg_heatmap.json'),
      fetchJSON('/data/scg_heatmap.json').catch(() => null),
      fetchJSON('/data/coverage.json'),
      fetchJSON('/data/mge_summary.json'),
      fetchJSON('/data/mge_per_bin.json'),
      fetchJSON('/data/eukaryotic.json'),
      fetchJSON('/data/contig_lengths.json'),
    ]);

    overview.set(overviewData);
    mags.set(magsData);
    taxonomySunburst.set(sunburstData);
    keggHeatmap.set(keggData);
    if (scgData) scgHeatmap.set(scgData);
    coverage.set(coverageData);
    mgeSummary.set(mgeSumData);
    mgePerBin.set(mgePerBinData);
    eukaryotic.set(eukData);
    contigLengths.set(contigLenData);
  } catch (e) {
    error.set(e.message);
    console.error('Data loading error:', e);
  } finally {
    loading.set(false);
  }
}

// Lazy load checkm2_all (3MB+ with all binner bins)
let checkm2Loading = false;
export async function loadCheckm2All() {
  if (checkm2Loading || get(checkm2All) !== null) return;
  checkm2Loading = true;
  try {
    const data = await fetchJSON('/data/checkm2_all.json');
    checkm2All.set(data);
  } catch (e) {
    console.error('Failed to load checkm2_all:', e);
  } finally {
    checkm2Loading = false;
  }
}

// Lazy load gene features (12MB+ keyed by contig ID)
let genesData = null;
let genesLoading = false;
export async function loadContigGenes(contigId) {
  if (!genesData && !genesLoading) {
    genesLoading = true;
    try {
      genesData = await fetchJSON('/data/genes.json');
    } catch (e) {
      console.warn('Gene data not available:', e.message);
      genesData = {};
    } finally {
      genesLoading = false;
    }
  }
  // Wait if another call is loading
  while (genesLoading) {
    await new Promise(r => setTimeout(r, 50));
  }
  return genesData?.[contigId] || null;
}

// Lazy load contig explorer + embeddings (separate files)
let contigLoading = false;
export async function loadContigExplorer() {
  if (contigLoading || get(contigExplorer) !== null) return;
  contigLoading = true;
  try {
    const [data, embeddings] = await Promise.all([
      fetchJSON('/data/contig_explorer.json'),
      fetchJSON('/data/contig_embeddings.json').catch(() => null),
    ]);
    // Merge embeddings (t-SNE/UMAP) into contig records
    if (embeddings && data?.contigs) {
      let nTsne = 0, nUmap = 0;
      for (const c of data.contigs) {
        const emb = embeddings[c.id];
        if (emb) {
          if (emb.tsne_x !== undefined) { c.tsne_x = emb.tsne_x; c.tsne_y = emb.tsne_y; nTsne++; }
          if (emb.umap_x !== undefined) { c.umap_x = emb.umap_x; c.umap_y = emb.umap_y; nUmap++; }
        }
      }
      data.has_tsne = nTsne > 0;
      data.has_umap = nUmap > 0;
    }
    contigExplorer.set(data);
  } catch (e) {
    console.error('Failed to load contig explorer:', e);
  } finally {
    contigLoading = false;
  }
}
