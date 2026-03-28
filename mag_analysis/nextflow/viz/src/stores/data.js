import { writable, get } from 'svelte/store';

// Individual data stores
export const overview = writable(null);
export const mags = writable(null);
export const binQuality = writable(null);
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
export const phyloTree = writable(null);
export const biosynthetic = writable(null);
export const ecosystemServices = writable(null);

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
    const res = await fetch('data/pipeline_status.json?t=' + Date.now());
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
  const gzRes = await fetch(url + '.gz?t=' + Date.now());
  const gzContentType = (gzRes.headers.get('Content-Type') || '').toLowerCase();
  if (gzRes.ok && !gzContentType.includes('text/html')) {
    const ce = (gzRes.headers.get('Content-Encoding') || '').toLowerCase();
    if (ce.includes('gzip') || ce.includes('br') || ce.includes('deflate')) {
      // Browser already decompressed transparently
      return gzRes.json();
    }
    // Sniff first two bytes — if they match the gzip magic number (1f 8b),
    // decompress manually; otherwise the server already decoded for us.
    const buf = await gzRes.arrayBuffer();
    const header = new Uint8Array(buf, 0, 2);
    if (header[0] === 0x1f && header[1] === 0x8b) {
      const ds = new DecompressionStream('gzip');
      const text = await new Response(
        new Blob([buf]).stream().pipeThrough(ds)
      ).text();
      return JSON.parse(text);
    }
    // Not actually gzipped (server decoded it) — parse as plain JSON
    const text = new TextDecoder().decode(buf);
    return JSON.parse(text);
  }
  // Fall back to plain JSON (older preprocess runs, dev environments)
  const res = await fetch(url + '?t=' + Date.now());
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
      fetchJSON('data/overview.json'),
      fetchJSON('data/mags.json').catch(() => null),
      fetchJSON('data/taxonomy_sunburst.json').catch(() => null),
      fetchJSON('data/kegg_heatmap.json').catch(() => null),
      fetchJSON('data/scg_heatmap.json').catch(() => null),
      fetchJSON('data/coverage.json').catch(() => null),
      fetchJSON('data/mge_summary.json').catch(() => null),
      fetchJSON('data/mge_per_bin.json').catch(() => null),
      fetchJSON('data/eukaryotic.json').catch(() => null),
      fetchJSON('data/contig_lengths.json').catch(() => null),
    ]);

    overview.set(overviewData);
    if (magsData) mags.set(magsData);
    if (sunburstData) taxonomySunburst.set(sunburstData);
    if (keggData) keggHeatmap.set(keggData);
    if (scgData) scgHeatmap.set(scgData);
    if (coverageData) coverage.set(coverageData);
    if (mgeSumData) mgeSummary.set(mgeSumData);
    if (mgePerBinData) mgePerBin.set(mgePerBinData);
    if (eukData) eukaryotic.set(eukData);
    if (contigLenData) contigLengths.set(contigLenData);
  } catch (e) {
    error.set(e.message);
    console.error('Data loading error:', e);
  } finally {
    loading.set(false);
  }
}

// Lazy load bin quality data (3MB+ with all binner bins)
let binQualityLoading = false;
export async function loadBinQuality() {
  if (binQualityLoading || get(binQuality) !== null) return;
  binQualityLoading = true;
  try {
    const data = await fetchJSON('data/checkm2_all.json');
    binQuality.set(data);
  } catch (e) {
    console.error('Failed to load checkm2_all:', e);
  } finally {
    binQualityLoading = false;
  }
}

// Lazy load gene features (12MB+ keyed by contig ID)
let genesData = null;
let genesLoading = false;
export async function loadContigGenes(contigId) {
  if (!genesData && !genesLoading) {
    genesLoading = true;
    try {
      genesData = await fetchJSON('data/genes.json');
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
        genesData = await fetchJSON('data/genes.json');
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
    const data = await fetchJSON('data/contig_sample_depths.json');
    sampleDepths.set(data);
  } catch (e) {
    console.warn('Per-sample depth data not available:', e.message);
    sampleDepths.set({ samples: [], depths: {} });
  } finally {
    sampleDepthsLoading = false;
  }
}

// Lazy load biosynthetic data (antiSMASH BGC regions)
let biosyntheticLoading = false;
export async function loadBiosynthetic() {
  if (biosyntheticLoading || get(biosynthetic) !== null) return;
  biosyntheticLoading = true;
  try {
    const data = await fetchJSON('data/biosynthetic.json');
    biosynthetic.set(data);
  } catch (e) {
    console.warn('Biosynthetic data not available:', e.message);
    biosynthetic.set({ n_regions: 0, type_counts: {}, regions: [], per_bin: {} });
  } finally {
    biosyntheticLoading = false;
  }
}

// Lazy load ecosystem services (ECOSSDB)
let esLoading = false;
export async function loadEcosystemServices() {
  if (esLoading || get(ecosystemServices) !== null) return;
  esLoading = true;
  try {
    const data = await fetchJSON('data/ecosystem_services.json');
    // Also try loading SDG data
    try {
      const sdg = await fetchJSON('data/es_sdg.json');
      if (sdg) data.sdg = sdg;
    } catch (e) {
      console.warn('SDG data not available');
    }
    ecosystemServices.set(data);
  } catch (e) {
    console.warn('Ecosystem services data not available:', e.message);
    ecosystemServices.set(null);
  } finally {
    esLoading = false;
  }
}

// Lazy load phylotree (GTDB-Tk phylogenetic classification)
let phyloLoading = false;
export async function loadPhyloTree() {
  if (phyloLoading || get(phyloTree) !== null) return;
  phyloLoading = true;
  try {
    const data = await fetchJSON('data/phylotree.json');
    phyloTree.set(data);
  } catch (e) {
    console.warn('Phylotree data not available:', e.message);
    phyloTree.set({ hierarchy: null, bins: [], newick: {} });
  } finally {
    phyloLoading = false;
  }
}

// Lazy load contig explorer + embeddings (separate files per method)
let contigLoading = false;
export async function loadContigExplorer() {
  if (contigLoading || get(contigExplorer) !== null) return;
  contigLoading = true;
  try {
    const [data, tsneEmb, umapEmb] = await Promise.all([
      fetchJSON('data/contig_explorer.json'),
      fetchJSON('data/contig_tsne.json').catch(() => null),
      fetchJSON('data/contig_umap.json').catch(() => null),
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
