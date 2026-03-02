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

/** Normalize a barcode value to "barcodeNN" format.
 *  Accepts: "1", "01", "barcode1", "barcode01" → "barcode01". */
function normalizeBarcode(raw) {
  if (!raw) return '';
  const s = String(raw).trim().toLowerCase();
  const m = s.match(/^(?:barcode)?(\d+)$/);
  if (!m) return '';
  return 'barcode' + m[1].padStart(2, '0');
}

/** Normalize a column header to its canonical name.
 *  Case-insensitive aliases: latitude→lat, longitude/long/lng→lon,
 *  sample_id→flowcell. */
function normalizeColumn(name) {
  const lower = name.trim().toLowerCase();
  const aliases = {
    'latitude': 'lat', 'longitude': 'lon', 'long': 'lon', 'lng': 'lon',
    'sample_id': 'flowcell',
  };
  return aliases[lower] || lower;
}

/** Parse a metadata TSV string and update the metadata store.
 *  Requires flowcell and barcode columns. Canonical sample key: FLOWCELL_barcodeNN.
 *  Validates against known sample IDs before replacing existing data.
 *  Returns { matched, total, error? }. */
export function loadMetadataTsv(text) {
  const lines = text.split('\n').filter(l => l.trim());
  if (lines.length < 2) return { matched: 0, total: 0, error: 'File has no data rows' };

  const rawHeader = lines[0].split('\t');
  const header = rawHeader.map(normalizeColumn);
  const hasFlowcell = header.includes('flowcell');
  const hasBarcode = header.includes('barcode');
  if (!hasFlowcell || !hasBarcode) {
    return { matched: 0, total: 0, error: 'Requires both flowcell (or sample_id) and barcode columns' };
  }

  const result = {};
  for (let i = 1; i < lines.length; i++) {
    const fields = lines[i].split('\t');
    const row = {};
    header.forEach((h, j) => { row[h] = fields[j] ?? ''; });

    const flowcell = (row.flowcell || '').trim();
    const barcode = normalizeBarcode(row.barcode);
    if (!flowcell || !barcode) continue;
    const key = `${flowcell}_${barcode}`;

    const parsed = {};
    for (const [k, v] of Object.entries(row)) {
      if (k === 'flowcell' || k === 'barcode') continue;
      const num = Number(v);
      if (v !== '' && !isNaN(num)) parsed[k] = num;
      else if (v) parsed[k] = v;
    }
    result[key] = parsed;
  }

  const total = Object.keys(result).length;
  if (total === 0) return { matched: 0, total: 0, error: 'No rows parsed' };

  // Validate: check how many keys match known samples
  const knownSamples = get(samples);
  const knownIds = knownSamples ? new Set(knownSamples.map(s => s.id)) : null;
  const matched = knownIds ? Object.keys(result).filter(k => knownIds.has(k)).length : total;

  if (matched === 0) {
    return { matched: 0, total, error: `0/${total} rows matched known samples` };
  }

  metadata.set(result);
  return { matched, total };
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
