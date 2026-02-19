import { writable, get } from 'svelte/store';

// Individual data stores
export const overview = writable(null);
export const mags = writable(null);
export const checkm2All = writable(null);
export const taxonomySunburst = writable(null);
export const keggHeatmap = writable(null);
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

// Lazy load contig explorer (large file)
let contigLoading = false;
export async function loadContigExplorer() {
  if (contigLoading || get(contigExplorer) !== null) return;
  contigLoading = true;
  try {
    const data = await fetchJSON('/data/contig_explorer.json');
    contigExplorer.set(data);
  } catch (e) {
    console.error('Failed to load contig explorer:', e);
  } finally {
    contigLoading = false;
  }
}
