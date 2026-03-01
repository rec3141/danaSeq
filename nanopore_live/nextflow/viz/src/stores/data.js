import { writable, get } from 'svelte/store';

export const overview = writable(null);
export const samples = writable(null);
export const sampleTsne = writable(null);
export const sampleTaxonomy = writable(null);
export const taxonomySunburst = writable(null);
export const sampleFunction = writable(null);
export const metadata = writable(null);

// Lazy-loaded (large)
export const readExplorer = writable(null);

export const loading = writable(true);
export const error = writable(null);

async function fetchJSON(url) {
  // Try pre-compressed .json.gz first (only trust it if content-type is JSON or gzip)
  try {
    const gzRes = await fetch(url + '.gz');
    if (gzRes.ok) {
      const ct = gzRes.headers.get('Content-Type') || '';
      if (ct.includes('json') || ct.includes('gzip') || ct.includes('octet-stream')) {
        const ce = gzRes.headers.get('Content-Encoding');
        if (ce && ce.includes('gzip')) return gzRes.json();
        const ds = new DecompressionStream('gzip');
        const text = await new Response(gzRes.body.pipeThrough(ds)).text();
        return JSON.parse(text);
      }
      // Content-Type is text/html (Vite SPA fallback) — not a real .gz file
    }
  } catch (_) { /* fall through to plain JSON */ }
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
      samplesData,
      tsneData,
      taxData,
      sunburstData,
      funcData,
      metaData,
    ] = await Promise.all([
      fetchJSON('/data/overview.json'),
      fetchJSON('/data/samples.json'),
      fetchJSON('/data/sample_tsne.json'),
      fetchJSON('/data/sample_taxonomy.json'),
      fetchJSON('/data/taxonomy_sunburst.json'),
      fetchJSON('/data/sample_function.json'),
      fetchJSON('/data/metadata.json').catch(() => null),
    ]);

    overview.set(overviewData);
    samples.set(samplesData);
    sampleTsne.set(tsneData);
    sampleTaxonomy.set(taxData);
    taxonomySunburst.set(sunburstData);
    sampleFunction.set(funcData);
    if (metaData) metadata.set(metaData);
  } catch (e) {
    error.set(e.message);
    console.error('Data loading error:', e);
  } finally {
    loading.set(false);
  }
}

// Lazy-load read explorer data (potentially huge)
let readLoading = false;
export async function loadReadExplorer() {
  if (readLoading || get(readExplorer) !== null) return;
  readLoading = true;
  try {
    const data = await fetchJSON('/data/read_explorer.json');
    readExplorer.set(data);
  } catch (e) {
    console.error('Failed to load read explorer:', e);
  } finally {
    readLoading = false;
  }
}
