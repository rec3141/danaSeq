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

// Live pipeline status (polled from pipeline_status.json)
export const pipelineStatus = writable(null);

let statusPollTimer = null;

export async function startStatusPolling(intervalMs = 30000) {
  if (statusPollTimer) return;
  await refreshPipelineStatus();
  statusPollTimer = setInterval(refreshPipelineStatus, intervalMs);
}

export function stopStatusPolling() {
  if (statusPollTimer) {
    clearInterval(statusPollTimer);
    statusPollTimer = null;
  }
}

async function refreshPipelineStatus() {
  try {
    // Cache-bust so vite preview doesn't serve stale data
    const res = await fetch('/data/pipeline_status.json?t=' + Date.now());
    if (res.ok) {
      const data = await res.json();
      pipelineStatus.set(data);
      // Update overview store's process statuses for DAG coloring
      overview.update(curr => curr ? {
        ...curr,
        processes: Object.fromEntries(
          Object.entries(data.processes).map(([k, v]) => [k, v.status])
        ),
        pipeline_total: data.pipeline_total,
        pipeline_completed: data.pipeline_completed,
        pipeline_running: data.pipeline_running,
        pipeline_pending: data.pipeline_pending,
        pipeline_failed: data.pipeline_failed,
        pipeline_skipped: data.pipeline_skipped,
      } : curr);
      if (!data.pipeline_active) stopStatusPolling();
    }
  } catch (e) { /* pipeline_status.json may not exist yet */ }
}

async function fetchJSON(url) {
  // Try pre-compressed .json.gz first — works with any static server, no server config needed.
  // If the server sets Content-Encoding: gzip (e.g. vite preview), the browser auto-decompresses
  // and we can just parse directly. Otherwise we decompress manually via DecompressionStream.
  const gzRes = await fetch(url + '.gz');
  if (gzRes.ok) {
    const ce = gzRes.headers.get('Content-Encoding');
    if (ce && ce.includes('gzip')) {
      // Browser already decompressed transparently
      return gzRes.json();
    }
    // Manual decompression (plain static server, no Content-Encoding header)
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

// Load all genes at once (for search index). Reuses genesData closure.
let genesLoadPromise = null;
export async function loadAllGenes() {
  if (genesData) return genesData;
  if (genesLoadPromise) return genesLoadPromise;
  genesLoadPromise = (async () => {
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
    while (genesLoading) {
      await new Promise(r => setTimeout(r, 50));
    }
    return genesData;
  })();
  return genesLoadPromise;
}

// Build search index: Map<contigId, string> where string = lowercased gene names + products joined by \0
export function buildGeneSearchIndex(allGenes) {
  const index = new Map();
  for (const [contigId, genes] of Object.entries(allGenes)) {
    const parts = [];
    for (const gene of genes) {
      if (gene.g) parts.push(gene.g);
      if (gene.p) parts.push(gene.p);
    }
    if (parts.length) {
      index.set(contigId, parts.join('\0').toLowerCase());
    }
  }
  return index;
}

// Lazy load per-sample depth data (sidecar to contig_explorer)
export const sampleDepths = writable(null);

let sampleDepthsLoading = false;
export async function loadSampleDepths() {
  if (sampleDepthsLoading || get(sampleDepths) !== null) return;
  sampleDepthsLoading = true;
  try {
    const data = await fetchJSON('/data/contig_sample_depths.json');
    sampleDepths.set(data);
  } catch (e) {
    console.warn('Per-sample depth data not available:', e.message);
    sampleDepths.set({ samples: [], depths: {} });
  } finally {
    sampleDepthsLoading = false;
  }
}

// Lazy load contig explorer + embeddings (separate files per method)
let contigLoading = false;
export async function loadContigExplorer() {
  if (contigLoading || get(contigExplorer) !== null) return;
  contigLoading = true;
  try {
    const [data, tsneEmb, umapEmb] = await Promise.all([
      fetchJSON('/data/contig_explorer.json'),
      fetchJSON('/data/contig_tsne.json').catch(() => null),
      fetchJSON('/data/contig_umap.json').catch(() => null),
    ]);
    // Merge embeddings into contig records (each file is {contig_id: [x, y]})
    if (data?.contigs) {
      let nTsne = 0, nUmap = 0;
      for (const c of data.contigs) {
        if (tsneEmb?.[c.id]) { c.tsne_x = tsneEmb[c.id][0]; c.tsne_y = tsneEmb[c.id][1]; nTsne++; }
        if (umapEmb?.[c.id]) { c.umap_x = umapEmb[c.id][0]; c.umap_y = umapEmb[c.id][1]; nUmap++; }
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
