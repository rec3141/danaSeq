import { writable, get } from 'svelte/store';
import { fetchJSON } from './fetchJson.js';

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
// Sample rows (samples.json) and sample metadata keyed by their FLOWCELL:barcodeNN
// id (metadata.json: shipped in the bundle or uploaded by a lab admin).
export const samples = writable(null);
export const sampleMetadata = writable(null);
// Lazy-loaded (large files); per-contig stores live in contigData.js
export const phyloTree = writable(null);
export const biosynthetic = writable(null);
export const ecosystemServices = writable(null);

export const loading = writable(true);
export const error = writable(null);

// Live pipeline status (polled from pipeline_status.json)
export const pipelineStatus = writable(null);

let statusPollTimer = null;
let statusPolling = false;
let statusFailures = 0;

export async function startStatusPolling(intervalMs = 30000) {
  if (statusPolling) return;
  statusPolling = true;
  statusFailures = 0;
  await refreshPipelineStatus();
  if (statusPolling) statusPollTimer = setInterval(refreshPipelineStatus, intervalMs);
}

export function stopStatusPolling() {
  statusPolling = false;
  if (statusPollTimer) {
    clearInterval(statusPollTimer);
    statusPollTimer = null;
  }
}

// Polling stops when the run is finished (pipeline_active false), when the file
// is absent (404: published runs often have none), or after two consecutive
// unusable responses (other HTTP errors, or an HTML fallback page instead of
// JSON). A network error that yields no response at all is treated as
// transient and polling continues.
async function refreshPipelineStatus() {
  let res;
  try {
    // Cache-bust so vite preview doesn't serve stale data
    res = await fetch('data/pipeline_status.json?t=' + Date.now());
  } catch (e) {
    return;
  }
  let data = null;
  if (res.ok) {
    try { data = await res.json(); } catch (e) { data = null; }
  }
  if (!data?.processes) {
    statusFailures++;
    if (res.status === 404 || statusFailures >= 2) stopStatusPolling();
    return;
  }
  statusFailures = 0;
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
      samplesData,
      metadataData,
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
      fetchJSON('data/samples.json').catch(() => null),
      fetchJSON('data/metadata.json').catch(() => null),
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
    if (samplesData) samples.set(samplesData);
    if (metadataData) sampleMetadata.set(metadataData);
  } catch (e) {
    error.set(e.message);
    console.error('Data loading error:', e);
  } finally {
    loading.set(false);
  }
}

// Lazy load bin quality data (3MB+ with all binner bins)
// One shared load, so a caller that needs the data can await one in flight
let binQualityLoading = null;
export function loadBinQuality() {
  if (get(binQuality) !== null) return Promise.resolve();
  binQualityLoading ??= fetchJSON('data/checkm2_all.json')
    .then(data => binQuality.set(data))
    .catch(e => console.error('Failed to load checkm2_all:', e))
    .finally(() => { binQualityLoading = null; });
  return binQualityLoading;
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
let phyloLoading = null;
export function loadPhyloTree() {
  if (get(phyloTree) !== null) return Promise.resolve();
  phyloLoading ??= fetchJSON('data/phylotree.json')
    .then(data => phyloTree.set(data))
    .catch(e => {
      console.warn('Phylotree data not available:', e.message);
      phyloTree.set({ hierarchy: null, bins: [], newick: {} });
    })
    .finally(() => { phyloLoading = null; });
  return phyloLoading;
}

// GTDB-Tk classification as a bin taxonomy source (taxonomy.gtdbtk). Bundles
// preprocessed before it was written there carry it only in phylotree.json's
// bins, so fill it in from there the first time it's asked for.
export async function ensureGtdbTaxonomy() {
  await loadBinQuality();
  const bins = get(binQuality);
  if (!bins?.length || bins.some(b => b.taxonomy?.gtdbtk)) return;
  await loadPhyloTree();
  const lineageByName = new Map((get(phyloTree)?.bins || []).map(b => [b.name, b.lineage]));
  if (!lineageByName.size) return;
  binQuality.set(bins.map(b => {
    const lineage = lineageByName.get(b.name);
    return lineage && Object.keys(lineage).length
      ? { ...b, taxonomy: { ...b.taxonomy, gtdbtk: lineage } }
      : b;
  }));
}

export * from './contigData.js';
