import { writable, derived, get } from 'svelte/store';
import { taxonomySource } from './taxonomySource.js';

export const overview = writable(null);
export const samples = writable(null);
export const sampleTsne = writable(null);
// Raw taxonomy payloads from JSON, each a {kraken, gtdb} envelope. Old flat
// builds (no source key) are treated as kraken for back-compat.
export const sampleTaxonomyRaw = writable(null);
export const taxonomySunburstRaw = writable(null);
export const sampleFunction = writable(null);
export const metadata = writable(null);

// Derived views that track the active taxonomy source. Consumers across the
// SPA import these and stay source-agnostic.
function pickSource(raw, source) {
  if (!raw) return null;
  if (raw.kraken || raw.gtdb) return raw[source] ?? raw.kraken ?? raw.gtdb ?? null;
  // Back-compat: pre-GTDB flat-shape JSON — treat as Kraken regardless of source.
  return raw;
}

export const sampleTaxonomy = derived(
  [sampleTaxonomyRaw, taxonomySource],
  ([$raw, $src]) => pickSource($raw, $src),
);

export const taxonomySunburst = derived(
  [taxonomySunburstRaw, taxonomySource],
  ([$raw, $src]) => pickSource($raw, $src),
);

// Lazy-loaded (large)
export const readExplorer = writable(null);

export const loading = writable(true);
export const error = writable(null);

async function fetchJSON(url) {
  // Resolve absolute /data/* paths against Vite's BASE_URL so the SPA works
  // when hosted under a subdirectory (e.g. /LoganSkylar/) as well as at /.
  // In dev / single-site builds BASE_URL is '/' so this is a no-op.
  if (url.startsWith('/data/')) {
    const base = (import.meta.env.BASE_URL || '/').replace(/\/$/, '');
    url = base + url;
  }
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
  const ct = res.headers.get('Content-Type') || '';
  if (ct.includes('html')) return null;  // Vite SPA fallback — file doesn't exist yet
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
      fetchJSON('./data/overview.json').catch(() => null),
      fetchJSON('./data/samples.json').catch(() => null),
      fetchJSON('./data/sample_tsne.json').catch(() => null),
      fetchJSON('./data/sample_taxonomy.json').catch(() => null),
      fetchJSON('./data/taxonomy_sunburst.json').catch(() => null),
      fetchJSON('./data/sample_function.json').catch(() => null),
      fetchJSON('./data/metadata.json').catch(() => null),
    ]);

    if (overviewData) overview.set(overviewData);
    if (samplesData) samples.set(samplesData);
    if (tsneData) sampleTsne.set(tsneData);
    if (taxData) sampleTaxonomyRaw.set(taxData);
    if (sunburstData) taxonomySunburstRaw.set(sunburstData);
    if (funcData) sampleFunction.set(funcData);
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
 *  samp_name/flowcell_barcode/id→sample_id (combined-column path),
 *  lat_lon kept as-is for split handling in row parser. */
function normalizeColumn(name) {
  const lower = name.trim().toLowerCase();
  const aliases = {
    'latitude': 'lat', 'longitude': 'lon', 'long': 'lon', 'lng': 'lon',
    'samp_name': 'sample_id', 'flowcell_barcode': 'sample_id', 'id': 'sample_id',
  };
  return aliases[lower] || lower;
}

/** Split a combined sample identifier into [flowcell, barcode].
 *  Accepts FLOWCELL:barcodeNN, FLOWCELL_barcodeNN, FLOWCELL-barcodeNN.
 *  Last separator wins so flowcells containing underscores still parse. */
function splitCombinedId(raw) {
  if (!raw) return ['', ''];
  const s = String(raw).trim();
  if (!s) return ['', ''];
  // Find the last occurrence of any separator: ':', '_', '-'.
  let lastIdx = -1;
  for (let i = s.length - 1; i >= 0; i--) {
    const c = s[i];
    if (c === ':' || c === '_' || c === '-') { lastIdx = i; break; }
  }
  if (lastIdx < 0) return [s, ''];
  return [s.slice(0, lastIdx), s.slice(lastIdx + 1)];
}

/** Try to parse a single value as a date and return ISO YYYY-MM-DD, or null. */
function tryParseDate(raw, mode) {
  const s = String(raw).trim();
  if (!s) return null;
  let m;
  if (mode === 'yyyymmdd') {
    m = s.match(/^(\d{4})(\d{2})(\d{2})$/);
    if (m) return `${m[1]}-${m[2]}-${m[3]}`;
  } else if (mode === 'iso') {
    m = s.match(/^(\d{4})-(\d{2})-(\d{2})(?:[T ].*)?$/);
    if (m) return `${m[1]}-${m[2]}-${m[3]}`;
  } else if (mode === 'slashed-iso') {
    m = s.match(/^(\d{4})\/(\d{2})\/(\d{2})$/);
    if (m) return `${m[1]}-${m[2]}-${m[3]}`;
  } else if (mode === 'mdy') {
    m = s.match(/^(\d{1,2})\/(\d{1,2})\/(\d{4})$/);
    if (m) {
      const mo = m[1].padStart(2, '0'), d = m[2].padStart(2, '0');
      if (+mo >= 1 && +mo <= 12 && +d >= 1 && +d <= 31) return `${m[3]}-${mo}-${d}`;
    }
  } else if (mode === 'dmy') {
    m = s.match(/^(\d{1,2})\/(\d{1,2})\/(\d{4})$/);
    if (m) {
      const d = m[1].padStart(2, '0'), mo = m[2].padStart(2, '0');
      if (+mo >= 1 && +mo <= 12 && +d >= 1 && +d <= 31) return `${m[3]}-${mo}-${d}`;
    }
  }
  return null;
}

/** Detect whether a column is dominantly a single date format and, if so,
 *  return a function (raw → ISO) to apply across the whole column. ≥80% of
 *  non-empty values must match. Conservative: returns null on mixed/ambiguous. */
function detectDateFormat(values) {
  const nonEmpty = values.filter(v => v != null && String(v).trim() !== '');
  if (nonEmpty.length === 0) return null;
  const modes = ['iso', 'yyyymmdd', 'slashed-iso', 'mdy', 'dmy'];
  const stats = {};
  for (const mode of modes) {
    let hits = 0;
    for (const v of nonEmpty) if (tryParseDate(v, mode)) hits++;
    stats[mode] = hits;
  }
  // Disambiguate MDY vs DMY: prefer DMY only if it disambiguates (any value
  // with month > 12 in MDY interpretation forces DMY).
  const slashCount = nonEmpty.filter(v => /^\d{1,2}\/\d{1,2}\/\d{4}$/.test(String(v).trim())).length;
  if (slashCount > 0) {
    let needsDmy = false;
    for (const v of nonEmpty) {
      const m = String(v).trim().match(/^(\d{1,2})\/(\d{1,2})\/\d{4}$/);
      if (m && +m[1] > 12 && +m[2] <= 12) { needsDmy = true; break; }
    }
    const slashMode = needsDmy ? 'dmy' : 'mdy';
    if (stats[slashMode] / nonEmpty.length >= 0.8) {
      return (v) => tryParseDate(v, slashMode) ?? v;
    }
  }
  // Try unambiguous formats.
  for (const mode of ['iso', 'yyyymmdd', 'slashed-iso']) {
    if (stats[mode] / nonEmpty.length >= 0.8) {
      return (v) => tryParseDate(v, mode) ?? v;
    }
  }
  return null;
}

/** Build a downloadable metadata template TSV string.
 *  Header: flowcell, barcode, samp_name, collection_date, lat, lon, env_medium, depth_m.
 *  One row per known sample (from the samples store), with flowcell/barcode/samp_name
 *  pre-filled from the canonical FLOWCELL:barcodeNN id. Other columns left empty. */
export function downloadMetadataTemplate() {
  const cols = ['flowcell', 'barcode', 'samp_name', 'collection_date', 'lat', 'lon', 'env_medium', 'depth_m'];
  const rows = [cols.join('\t')];
  const samps = get(samples);
  if (samps) {
    for (const s of samps) {
      const id = s.id || '';
      const sepIdx = id.indexOf(':');
      const flowcell = sepIdx >= 0 ? id.slice(0, sepIdx) : id;
      const barcode = sepIdx >= 0 ? id.slice(sepIdx + 1) : '';
      rows.push([flowcell, barcode, id, '', '', '', '', ''].join('\t'));
    }
  }
  return rows.join('\n') + '\n';
}

/** Parse a delimited (CSV/TSV) string into rows of string fields.
 *  Honors quoted fields with embedded delimiters, escaped double-quotes
 *  ("" → "), and CRLF/LF line endings. Empty lines are skipped. */
function parseDelimited(text, delim) {
  const rows = [];
  let row = [];
  let field = '';
  let inQuotes = false;
  let fieldStarted = false;
  const n = text.length;
  for (let i = 0; i < n; i++) {
    const c = text[i];
    if (inQuotes) {
      if (c === '"') {
        if (text[i + 1] === '"') { field += '"'; i++; continue; }
        inQuotes = false;
        continue;
      }
      field += c;
      continue;
    }
    if (c === '"' && !fieldStarted) { inQuotes = true; fieldStarted = true; continue; }
    if (c === delim) { row.push(field); field = ''; fieldStarted = false; continue; }
    if (c === '\r') continue;
    if (c === '\n') {
      row.push(field);
      if (row.some(v => v !== '')) rows.push(row);
      row = []; field = ''; fieldStarted = false;
      continue;
    }
    field += c;
    fieldStarted = true;
  }
  if (field !== '' || row.length > 0) {
    row.push(field);
    if (row.some(v => v !== '')) rows.push(row);
  }
  return rows;
}

/** Heuristically detect tab vs comma delimiter from the first line of text. */
function detectDelimiter(text) {
  let tabs = 0, commas = 0, inQuotes = false;
  const lookahead = Math.min(text.length, 4096);
  for (let i = 0; i < lookahead; i++) {
    const c = text[i];
    if (c === '"') inQuotes = !inQuotes;
    else if (!inQuotes) {
      if (c === '\n') { if (tabs > 0 || commas > 0) break; }
      else if (c === '\t') tabs++;
      else if (c === ',') commas++;
    }
  }
  return tabs > 0 ? '\t' : ',';
}

/** Process pre-parsed metadata rows ([header, ...dataRows] of string arrays)
 *  and update the metadata store.
 *  Accepts either flowcell + barcode columns, or a single combined column
 *  (sample_id / samp_name / flowcell_barcode / id) with values like
 *  FAR84275:barcode00, FAR84275_barcode00, or FAR84275-barcode00.
 *  Canonical sample key: FLOWCELL:barcodeNN.
 *  Validates against known sample IDs before replacing existing data.
 *  Returns { matched, total, error? }. */
function processMetadataRows(rows) {
  if (!rows || rows.length < 2) return { matched: 0, total: 0, error: 'File has no data rows' };

  const rawHeader = rows[0];
  const header = rawHeader.map(normalizeColumn);
  const hasFlowcell = header.includes('flowcell');
  const hasBarcode = header.includes('barcode');
  const hasCombined = header.includes('sample_id');
  const useCombined = !(hasFlowcell && hasBarcode) && hasCombined;
  if (!useCombined && !(hasFlowcell && hasBarcode)) {
    return {
      matched: 0, total: 0,
      error: 'Requires either flowcell + barcode columns, or a combined sample_id/samp_name/flowcell_barcode/id column',
    };
  }

  // First pass: parse rows into a list (preserving column values as raw strings)
  // so we can run a column-wide date-format autodetect before final coercion.
  const rawRows = [];
  for (let i = 1; i < rows.length; i++) {
    const fields = rows[i];
    const row = {};
    header.forEach((h, j) => { row[h] = fields[j] != null ? String(fields[j]) : ''; });

    let flowcell, barcode;
    if (useCombined) {
      const [fc, bc] = splitCombinedId(row.sample_id);
      flowcell = fc.trim();
      barcode = normalizeBarcode(bc);
    } else {
      flowcell = (row.flowcell || '').trim();
      barcode = normalizeBarcode(row.barcode);
    }
    if (!flowcell || !barcode) continue;
    const key = `${flowcell}:${barcode}`;

    // Strip identifier columns from the data payload.
    const data = {};
    for (const [k, v] of Object.entries(row)) {
      if (k === 'flowcell' || k === 'barcode' || k === 'sample_id') continue;
      // Split lat_lon into lat / lon (MIxS spec: "lat lon" or "lat,lon").
      if (k === 'lat_lon') {
        const sv = String(v).trim();
        if (sv) {
          const parts = sv.split(/[\s,]+/).filter(Boolean);
          if (parts.length >= 2) {
            data.lat = parts[0];
            data.lon = parts[1];
          }
        }
        continue;
      }
      data[k] = v;
    }
    rawRows.push({ key, data });
  }

  if (rawRows.length === 0) return { matched: 0, total: 0, error: 'No rows parsed' };

  // Column-wide date autodetect: collect all column names, then for each
  // non-identifier column, sniff for a single dominant date format.
  const allCols = new Set();
  for (const r of rawRows) for (const k of Object.keys(r.data)) allCols.add(k);
  const dateFormatters = {};
  for (const col of allCols) {
    const colValues = rawRows.map(r => r.data[col]).filter(v => v !== undefined);
    const fmt = detectDateFormat(colValues);
    if (fmt) dateFormatters[col] = fmt;
  }

  // Final pass: coerce values (numbers stay numbers; date columns get ISO strings).
  const result = {};
  for (const { key, data } of rawRows) {
    const parsed = {};
    for (const [k, v] of Object.entries(data)) {
      if (v === '' || v == null) continue;
      if (dateFormatters[k]) {
        parsed[k] = dateFormatters[k](v);
        continue;
      }
      const num = Number(v);
      if (!isNaN(num)) parsed[k] = num;
      else parsed[k] = v;
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

/** Parse a metadata TSV/CSV string (delimiter auto-detected) and update the
 *  metadata store. Returns { matched, total, error? }. */
export function loadMetadataTsv(text) {
  const delim = detectDelimiter(text);
  const rows = parseDelimited(text, delim);
  return processMetadataRows(rows);
}

/** Parse an XLSX workbook (first sheet) and update the metadata store.
 *  Takes an ArrayBuffer. Returns { matched, total, error? } via Promise. */
export async function loadMetadataXlsx(arrayBuffer) {
  let XLSX;
  try {
    XLSX = await import('xlsx');
  } catch (e) {
    return { matched: 0, total: 0, error: 'XLSX support failed to load' };
  }
  let wb;
  try {
    wb = XLSX.read(arrayBuffer, { type: 'array', cellDates: true });
  } catch (e) {
    return { matched: 0, total: 0, error: 'Could not read XLSX file (corrupt or unsupported)' };
  }
  const sheetName = wb.SheetNames[0];
  if (!sheetName) return { matched: 0, total: 0, error: 'Workbook contains no sheets' };
  const ws = wb.Sheets[sheetName];
  const rawRows = XLSX.utils.sheet_to_json(ws, { header: 1, defval: '', raw: true, blankrows: false });
  // Stringify cell values: Date → ISO YYYY-MM-DD, everything else via String().
  const rows = rawRows.map(r => r.map(v => {
    if (v == null) return '';
    if (v instanceof Date) {
      const y = v.getUTCFullYear();
      const m = String(v.getUTCMonth() + 1).padStart(2, '0');
      const d = String(v.getUTCDate()).padStart(2, '0');
      return `${y}-${m}-${d}`;
    }
    return String(v);
  }));
  return processMetadataRows(rows);
}

/** Dispatcher: import metadata from a File. Detects type by extension and
 *  routes to the correct parser. Returns { matched, total, error? }. */
export async function loadMetadataFile(file) {
  if (!file) return { matched: 0, total: 0, error: 'No file selected' };
  const name = (file.name || '').toLowerCase();
  const dot = name.lastIndexOf('.');
  const ext = dot >= 0 ? name.slice(dot + 1) : '';
  if (ext === 'xlsx') {
    const buf = await file.arrayBuffer();
    return loadMetadataXlsx(buf);
  }
  if (ext === 'tsv' || ext === 'csv' || ext === 'txt') {
    const text = await file.text();
    return loadMetadataTsv(text);
  }
  return {
    matched: 0, total: 0,
    error: `Unsupported file type${ext ? ` (.${ext})` : ''} — use .tsv, .csv, .txt, or .xlsx`,
  };
}

// Lazy-load read explorer data (potentially huge)
let readLoading = false;
export async function loadReadExplorer() {
  if (readLoading || get(readExplorer) !== null) return;
  readLoading = true;
  try {
    const data = await fetchJSON('./data/read_explorer.json');
    readExplorer.set(data);
  } catch (e) {
    console.error('Failed to load read explorer:', e);
  } finally {
    readLoading = false;
  }
}
