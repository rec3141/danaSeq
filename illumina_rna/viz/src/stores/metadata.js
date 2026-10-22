import { writable, derived, get } from 'svelte/store';
import { samples, fmtPct } from './data.js';
import { fetchJSON } from './fetchJson.js';
import { wardOrder } from '../lib/cluster.js';

// ── Sample metadata ──────────────────────────────────────────────────────────
// A table with one row per sample, from (first wins):
//   1. a file uploaded in this browser (TSV / CSV / XLSX; kept in localStorage
//      until removed),
//   2. data/metadata.json: { "<sample>": { column: value } }, the shape both the
//      pipeline (--metadata) and microscape.app's run-editor upload produce
//      (the server answers that URL with the editor's table when there is one).
// Both become { source, header[], rows[][] } of strings; matching samples and
// typing the columns happen here, so every path behaves the same.

export const metaRaw = writable(null);      // { source, origin: 'bundle'|'upload', header, rows }
let bundled = null;

const storeKey = () => `danaseq-rna-metadata:${location.pathname}`;
function saveUpload(raw) {
  try { raw ? localStorage.setItem(storeKey(), JSON.stringify(raw)) : localStorage.removeItem(storeKey()); } catch { /* storage blocked */ }
}

/** { "<sample>": { column: value } } → { header, rows } with a sample column. */
function keyedToRaw(obj, source) {
  const cols = [];
  for (const row of Object.values(obj)) for (const k of Object.keys(row ?? {})) if (!cols.includes(k)) cols.push(k);
  return { source, origin: 'bundle', header: ['sample', ...cols],
    rows: Object.entries(obj).map(([id, row]) => [id, ...cols.map(c => (row?.[c] == null ? '' : String(row[c])))]) };
}

export async function initMetadata() {
  try {
    const m = await fetchJSON('./data/metadata.json');
    if (m && typeof m === 'object' && Object.keys(m).length) bundled = keyedToRaw(m, 'metadata.json');
  } catch { /* none shipped or uploaded on microscape.app */ }
  let upload = null;
  try { upload = JSON.parse(localStorage.getItem(storeKey()) || 'null'); } catch { upload = null; }
  metaRaw.set(upload?.header ? upload : bundled);
}

export function clearUploadedMetadata() {
  saveUpload(null);
  metaRaw.set(bundled);
}

export const hasBundledMetadata = () => !!bundled;

// Sample names as the pipeline writes them, and loose forms a sheet may use
// (file names, R1 suffixes, case).
const loose = s => String(s).trim().toLowerCase()
  .replace(/\.(fastq|fq)(\.gz)?$/, '').replace(/[_.-](r?[12])$/, '');
const ID_COLUMNS = ['sample', 'sample_id', 'sample_name', 'samp_name', 'sampleid', 'name', 'id', 'library', 'library_id'];

const normCol = h => String(h).trim();
const isId = h => ID_COLUMNS.includes(h.trim().toLowerCase().replace(/\s+/g, '_'));

// ── Column typing ──
const num = v => {
  const s = String(v).trim().replace(/,/g, '');
  return s !== '' && !isNaN(Number(s)) ? Number(s) : null;
};
const DATE_RES = [
  [/^(\d{4})-(\d{1,2})-(\d{1,2})(?:[T ].*)?$/, m => [m[1], m[2], m[3]]],
  [/^(\d{4})\/(\d{1,2})\/(\d{1,2})$/, m => [m[1], m[2], m[3]]],
  [/^(\d{4})(\d{2})(\d{2})$/, m => [m[1], m[2], m[3]]],
];
function isoDate(v, dmy) {
  const s = String(v).trim();
  for (const [re, f] of DATE_RES) {
    const m = s.match(re);
    if (m) { const [y, mo, d] = f(m); return check(y, mo, d); }
  }
  const m = s.match(/^(\d{1,2})[/.](\d{1,2})[/.](\d{4})$/);
  if (m) return dmy ? check(m[3], m[2], m[1]) : check(m[3], m[1], m[2]);
  return null;
}
function check(y, mo, d) {
  if (+mo < 1 || +mo > 12 || +d < 1 || +d > 31) return null;
  return `${y}-${String(mo).padStart(2, '0')}-${String(d).padStart(2, '0')}`;
}

/** Type one column from its non-empty string values: number (≥90% numeric),
 *  date (≥80% one date form; 8-digit numbers count as YYYYMMDD) or category. */
function typeColumn(values) {
  const v = values.filter(x => x != null && String(x).trim() !== '');
  if (!v.length) return null;
  const dmy = v.some(x => { const m = String(x).trim().match(/^(\d{1,2})[/.](\d{1,2})[/.]\d{4}$/); return m && +m[1] > 12; });
  const dates = v.filter(x => isoDate(x, dmy)).length;
  const eightDigit = v.every(x => /^\d{8}$/.test(String(x).trim()));
  if (dates / v.length >= 0.8 && (eightDigit || !v.every(x => num(x) != null))) return { type: 'date', dmy };
  if (v.filter(x => num(x) != null).length / v.length >= 0.9) return { type: 'number' };
  return { type: 'category' };
}

/** Matched, typed metadata for the loaded samples:
 *  { source, origin, idColumn, columns: [{key, type, levels?, min?, max?}],
 *    bySample: {name: {key: value}}, matched, unmatchedRows[], missingSamples[] } */
export const sampleMeta = derived([metaRaw, samples], ([$raw, $samples]) => {
  if (!$raw?.header?.length || !$raw.rows?.length) return null;
  const names = $samples?.samples ?? [];
  const exact = new Map(names.map(n => [n, n]));
  const byLoose = new Map(names.map(n => [loose(n), n]));

  const header = $raw.header.map(normCol);
  let idCol = header.findIndex(isId);
  if (idCol < 0) idCol = 0;

  // Rows → samples; a sample named twice keeps its first row
  const matchedRows = new Map();
  const unmatchedRows = [];
  for (const r of $raw.rows) {
    const id = String(r[idCol] ?? '').trim();
    if (!id) continue;
    const name = exact.get(id) ?? byLoose.get(loose(id));
    if (!name) { unmatchedRows.push(id); continue; }
    if (!matchedRows.has(name)) matchedRows.set(name, r);
  }

  // lat_lon "lat lon" splits into two numeric columns (MIxS)
  const cols = [];
  header.forEach((h, j) => {
    if (j === idCol || !h) return;
    if (/^lat[_ ]?lon$/i.test(h)) {
      cols.push({ key: 'lat', src: j, part: 0 }, { key: 'lon', src: j, part: 1 });
    } else cols.push({ key: h, src: j });
  });
  const cell = (r, c) => {
    const v = r[c.src] == null ? '' : String(r[c.src]).trim();
    return c.part == null ? v : (v.split(/[\s,]+/)[c.part] ?? '');
  };

  const bySample = Object.fromEntries([...matchedRows.keys()].map(n => [n, {}]));
  const columns = [];
  for (const c of cols) {
    const vals = [...matchedRows.values()].map(r => cell(r, c));
    const t = typeColumn(vals);
    if (!t) continue;   // empty for every matched sample
    const col = { key: c.key, type: t.type };
    for (const [n, r] of matchedRows) {
      const s = cell(r, c);
      if (s === '') continue;
      bySample[n][c.key] = t.type === 'number' ? num(s)
        : t.type === 'date' ? (isoDate(/^\d{8}$/.test(s) ? `${s.slice(0, 4)}-${s.slice(4, 6)}-${s.slice(6)}` : s, t.dmy) ?? s)
        : s;
    }
    const present = Object.values(bySample).map(m => m[c.key]).filter(v => v != null);
    if (t.type === 'category') {
      col.levels = [...new Set(present)].sort((a, b) => String(a).localeCompare(String(b), undefined, { numeric: true }));
    } else {
      const nums = present.map(v => (t.type === 'date' ? Date.parse(v) : v)).filter(Number.isFinite);
      col.min = Math.min(...nums); col.max = Math.max(...nums);
    }
    columns.push(col);
  }

  return {
    source: $raw.source, origin: $raw.origin, idColumn: header[idCol], columns, bySample,
    matched: matchedRows.size, unmatchedRows,
    missingSamples: names.filter(n => !matchedRows.has(n)),
  };
});

// ── Loading files ──
function parseDelimited(text, delim) {
  const rows = [];
  let row = [], field = '', inQuotes = false, started = false;
  for (let i = 0; i < text.length; i++) {
    const c = text[i];
    if (inQuotes) {
      if (c === '"') { if (text[i + 1] === '"') { field += '"'; i++; } else inQuotes = false; }
      else field += c;
      continue;
    }
    if (c === '"' && !started) { inQuotes = true; started = true; }
    else if (c === delim) { row.push(field); field = ''; started = false; }
    else if (c === '\r') continue;
    else if (c === '\n') {
      row.push(field);
      if (row.some(v => v !== '')) rows.push(row);
      row = []; field = ''; started = false;
    } else { field += c; started = true; }
  }
  if (field !== '' || row.length) { row.push(field); if (row.some(v => v !== '')) rows.push(row); }
  return rows;
}
const firstLine = t => t.slice(0, t.indexOf('\n') < 0 ? t.length : t.indexOf('\n'));

async function readXlsx(buf) {
  const XLSX = await import('xlsx');
  const wb = XLSX.read(buf, { type: 'array', cellDates: true });
  const ws = wb.Sheets[wb.SheetNames[0]];
  if (!ws) throw new Error('Workbook has no sheets');
  return XLSX.utils.sheet_to_json(ws, { header: 1, defval: '', raw: true, blankrows: false })
    .map(r => r.map(v => (v instanceof Date ? v.toISOString().slice(0, 10) : v == null ? '' : String(v))));
}

/** Load a metadata file picked by the viewer. Returns { matched, total, error? }. */
export async function loadMetadataFile(file) {
  const ext = (file.name.split('.').pop() || '').toLowerCase();
  let rows;
  try {
    if (ext === 'xlsx' || ext === 'xls') rows = await readXlsx(await file.arrayBuffer());
    else if (['tsv', 'csv', 'txt'].includes(ext)) {
      const text = (await file.text()).replace(/^﻿/, '');
      rows = parseDelimited(text, firstLine(text).includes('\t') ? '\t' : ',');
    } else return { error: `Unsupported file type (.${ext}): use .tsv, .csv, .txt or .xlsx` };
  } catch (e) {
    return { error: `Could not read ${file.name}: ${e.message}` };
  }
  if (!rows || rows.length < 2) return { error: 'The file needs a header row and at least one sample row' };
  const raw = { source: file.name, origin: 'upload', header: rows[0].map(String), rows: rows.slice(1) };
  // Only replace the current table if some rows name a known sample
  const prev = get(metaRaw);
  metaRaw.set(raw);
  const m = get(sampleMeta);
  if (!m?.matched) {
    metaRaw.set(prev);
    return { error: `No row names a sample of this run (sample names look like ${get(samples)?.samples?.[0] ?? '?'})` };
  }
  saveUpload(raw);
  return { matched: m.matched, total: m.matched + m.unmatchedRows.length, columns: m.columns.length };
}

/** Template TSV: one row per sample with everything already known about it:
 *  flowcell / lane / barcode parsed from the sample name (FLOWCELL_Lnn_BARCODE)
 *  and the values of the loaded metadata table, plus empty MIxS-style columns. */
export function metadataTemplate() {
  const names = get(samples)?.samples ?? [];
  const m = get(sampleMeta);
  const fromName = n => { const x = n.match(/^(.+)_(L\d+)_([^_]+)$/); return x ? { flowcell: x[1], lane: x[2], barcode: x[3] } : {}; };
  const known = (m?.columns ?? []).map(c => c.key);
  const cols = ['sample', ...new Set([...(names.some(n => fromName(n).flowcell) ? ['flowcell', 'lane', 'barcode'] : []), ...known,
    'collection_date', 'lat', 'lon', 'depth_m', 'env_medium', 'site'])];
  return [cols.join('\t'), ...names.map(n => {
    const v = { ...fromName(n), ...(m?.bySample[n] ?? {}), sample: n };
    return cols.map(c => v[c] ?? '').join('\t');
  })].join('\n') + '\n';
}

export function metadataTsv() {
  const m = get(sampleMeta);
  if (!m) return '';
  const names = get(samples)?.samples ?? [];
  const keys = m.columns.map(c => c.key);
  return [['sample', ...keys].join('\t'),
    ...names.filter(n => m.bySample[n]).map(n => [n, ...keys.map(k => m.bySample[n][k] ?? '')].join('\t'))].join('\n') + '\n';
}

// ── Sample annotations: one choice shared by every view ──────────────────────
// Built-in attributes (sequencing run, read-flow QC) plus each metadata column.
// The chosen one colours the sample t-SNEs, the strip above the heatmaps and
// the composition bars, and can order the samples there.

const runOf = n => n.split('_').slice(0, -1).join('_') || n;

// Categorical hues in fixed order (Tailwind 400 steps on the dark surface);
// levels past the palette fold into one grey "Other".
export const CATEGORY_PALETTE = ['#22d3ee', '#f472b6', '#fbbf24', '#34d399', '#a78bfa', '#fb923c',
  '#60a5fa', '#f87171', '#a3e635', '#e879f9', '#2dd4bf', '#facc15'];
export const OTHER_COLOR = '#64748b';
export const MISSING_COLOR = '#334155';

let annotKey0 = 'run';
try { annotKey0 = localStorage.getItem('danaseq-rna-annot') || 'run'; } catch { /* storage blocked */ }
export const sampleAnnotKey = writable(annotKey0);
sampleAnnotKey.subscribe(k => { try { localStorage.setItem('danaseq-rna-annot', k); } catch { /* storage blocked */ } });

// Viridis stops for numeric annotations (same scale as the heatmaps)
const VIRIDIS = ['#440154', '#482878', '#3e4989', '#31688e', '#26828e', '#1f9e89', '#35b779', '#6ece58', '#b5de2b', '#fde725'];
function viridis(t) {
  const x = Math.max(0, Math.min(1, t)) * (VIRIDIS.length - 1);
  const i = Math.min(VIRIDIS.length - 2, Math.floor(x)), f = x - i;
  const p = h => [1, 3, 5].map(o => parseInt(h.slice(o, o + 2), 16));
  const a = p(VIRIDIS[i]), b = p(VIRIDIS[i + 1]);
  return '#' + a.map((v, k) => Math.round(v + (b[k] - v) * f).toString(16).padStart(2, '0')).join('');
}

/** All annotations: [{ key, label, group, type, value(name), display(name) }]. */
export const sampleAnnotations = derived([samples, sampleMeta], ([$samples, $meta]) => {
  const rows = Object.fromEntries(($samples?.rows ?? []).map(r => [r.sample, r]));
  const out = [
    { key: 'run', label: 'Sequencing run', group: 'Sequencing', type: 'category', value: runOf },
    ...[['pct_rrna_total', '% rRNA'], ['pct_mapped', '% mapped'], ['pct_cds', '% on CDS']]
      .map(([k, label]) => ({ key: k, label, group: 'Read QC', type: 'number', unit: '%', value: n => rows[n]?.[k] ?? null })),
  ];
  for (const c of $meta?.columns ?? []) {
    out.push({ key: `meta:${c.key}`, label: c.key, group: 'Metadata', type: c.type,
      value: n => $meta.bySample[n]?.[c.key] ?? null });
  }
  for (const a of out) {
    a.display = n => { const v = a.value(n); return v == null ? '–' : a.unit ? fmtPct(v) : String(v); };
  }
  return out;
});

/** The chosen annotation with colours for a sample list:
 *  { ...annotation, colorOf(name), legend: { kind: 'category', items: [{label, color, n}] }
 *    | { kind: 'number'|'date', min, max, minLabel, maxLabel }, order(names) } */
export const sampleAnnot = derived([sampleAnnotations, sampleAnnotKey, samples], ([$all, $key, $samples]) => {
  const a = $all.find(x => x.key === $key) ?? $all[0];
  if (!a) return null;
  const names = $samples?.samples ?? [];
  const vals = names.map(a.value);
  let colorOf, legend;
  if (a.type === 'category') {
    const counts = new Map();
    vals.forEach(v => { if (v != null) counts.set(String(v), (counts.get(String(v)) || 0) + 1); });
    const levels = [...counts.keys()].sort((x, y) => x.localeCompare(y, undefined, { numeric: true }));
    const col = new Map(levels.map((l, i) => [l, i < CATEGORY_PALETTE.length ? CATEGORY_PALETTE[i] : OTHER_COLOR]));
    colorOf = n => { const v = a.value(n); return v == null ? MISSING_COLOR : col.get(String(v)) ?? OTHER_COLOR; };
    const items = levels.slice(0, CATEGORY_PALETTE.length).map(l => ({ label: l, color: col.get(l), n: counts.get(l) }));
    const rest = levels.slice(CATEGORY_PALETTE.length);
    if (rest.length) items.push({ label: `Other (${rest.length})`, color: OTHER_COLOR, n: rest.reduce((s, l) => s + counts.get(l), 0) });
    legend = { kind: 'category', items, levels, colorFor: l => col.get(String(l)) ?? OTHER_COLOR };
  } else {
    const toNum = v => (v == null ? null : a.type === 'date' ? Date.parse(v) : v);
    const nums = vals.map(toNum).filter(Number.isFinite);
    const min = Math.min(...nums), max = Math.max(...nums);
    const t = v => (max > min ? (v - min) / (max - min) : 0.5);
    colorOf = n => { const v = toNum(a.value(n)); return Number.isFinite(v) ? viridis(t(v)) : MISSING_COLOR; };
    const lab = v => (a.type === 'date' ? new Date(v).toISOString().slice(0, 10) : a.unit ? fmtPct(v) : +v.toPrecision(4));
    legend = nums.length ? { kind: a.type, min, max, minLabel: lab(min), maxLabel: lab(max), colorAt: viridis } : { kind: a.type };
  }
  const missing = vals.filter(v => v == null).length;
  // Sample order by the annotation (missing last), ties by name
  const order = list => list.map((n, j) => [n, j]).sort(([x], [y]) => {
    const vx = a.value(x), vy = a.value(y);
    if (vx == null || vy == null) return (vx == null) - (vy == null) || x.localeCompare(y);
    const c = typeof vx === 'number' && typeof vy === 'number' ? vx - vy : String(vx).localeCompare(String(vy), undefined, { numeric: true });
    return c || x.localeCompare(y);
  }).map(([, j]) => j);
  return { ...a, colorOf, legend, missing, order, viridis };
});

/** Heatmap strip (D3Heatmap colAnnotation) for samples `names`. */
export const annotStrip = (a, names) =>
  (a ? { label: a.label, colors: names.map(a.colorOf), values: names.map(a.display) } : null);

/** Column order for an entities × samples heatmap: Ward on log1p values
 *  (reusing a row clustering's column order when there is one) or by the annotation. */
export function sampleColumnOrder(mode, annot, names, matrix, clustColOrder = null) {
  if (mode === 'annot' && annot) return annot.order(names);
  if (clustColOrder) return clustColOrder;
  if (!matrix.length || names.length < 3) return names.map((_, j) => j);
  return wardOrder(names.map((_, j) => matrix.map(r => Math.log1p(r[j] ?? 0))));
}
