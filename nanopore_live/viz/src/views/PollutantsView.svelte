<script>
  import LeafletMap from '../components/charts/LeafletMap.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, metadata } from '../stores/data.js';
  import { cartItems, cartActive, toggleCart, addToCart } from '../stores/cart.js';
  import { selectedSample } from '../stores/selection.js';

  // ---------------------------------------------------------------------------
  // Data: ./data/pollutants.json (viz/preprocess/preprocess.py
  // build_pollutant_payload). Shape:
  //   features: [[db, category, pathway, gene], ...]
  //   dbs:      {db: {params, n_samples, n_genes_hit}}
  //   samples:  {sid: {n_proteins: {db: N}, counts: {featureIdx: hits}}}
  // Each protein is counted once per DB (best-scoring target).
  // ---------------------------------------------------------------------------
  let payload = $state(null);
  let loadState = $state('loading'); // loading | ok | missing

  async function fetchMaybeGz(url) {
    try {
      const gz = await fetch(url + '.gz');
      if (gz.ok && (gz.headers.get('Content-Type') || '').match(/json|gzip|octet/)) {
        const ce = gz.headers.get('Content-Encoding') || '';
        if (ce.includes('gzip')) return await gz.json();
        const ds = new DecompressionStream('gzip');
        return JSON.parse(await new Response(gz.body.pipeThrough(ds)).text());
      }
    } catch (_) { /* fall through */ }
    const res = await fetch(url);
    if (!res.ok) throw new Error(`HTTP ${res.status}`);
    if ((res.headers.get('Content-Type') || '').includes('html')) throw new Error('not found');
    return res.json();
  }

  $effect(() => {
    const base = (import.meta.env.BASE_URL || '/').replace(/\/$/, '');
    fetchMaybeGz(`${base}/data/pollutants.json`)
      .then(p => {
        payload = p; loadState = 'ok';
        // Methylmercury is the headline concern: open on hgcA when present.
        if ((p.features || []).some(f => f[0] === 'Hg-MATE' && f[3] === 'hgcA')) selectFeatured('hgcA');
      })
      .catch(e => { console.warn('pollutants.json not available:', e.message); loadState = 'missing'; });
  });

  // Display metadata for known databases. Unknown DBs still render, with
  // their raw name and the search parameters recorded in the payload.
  const DB_INFO = {
    'Hg-MATE': {
      label: 'Hg-MATE',
      blurb: 'Mercury methylation genes hgcA / hgcB (methylmercury production)',
      cite: 'Gionfriddo et al. 2021, Hg-MATE-Db.v1, doi:10.25573/serc.13105370',
      cutoff: 'hmmsearch: hgcA ≥95 bits + cap-helix motif; hgcB ≥80 bits',
      color: '#fb7185',
    },
    'CANT-HYD': {
      label: 'CANT-HYD',
      blurb: 'Hydrocarbon-degradation marker HMMs (aerobic + anaerobic alkane/aromatic)',
      cite: 'Khot et al. 2022, Front. Microbiol. 12:764058',
      cutoff: 'hmmsearch, noise cutoffs (--cut_nc)',
      color: '#fbbf24',
    },
    'Mapping markers': {
      label: 'Fecal markers (read mapping)',
      blurb: 'Human fecal / sewage indicator genomes (crAssphage) mapped with minimap2',
      cite: 'Dutilh et al. 2014; Stachler et al. 2017 (crAssphage as sewage marker)',
      cutoff: 'reads aligned at ≥90% identity; always per million reads',
      color: '#a78bfa',
    },
    'KEGG-xeno': {
      label: 'KEGG xenobiotics',
      blurb: 'KOfam profiles for KEGG xenobiotics-degradation pathways (BTEX, PAH, organohalides, atrazine, nitroaromatics, steroids …)',
      cite: 'Aramaki et al. 2020 (KOfam); KEGG br08901 class 1.11',
      cutoff: 'hmmsearch, kofam_scan per-KO thresholds (full / domain score); core-metabolism KOs excluded',
      color: '#38bdf8',
    },
    'BacMet': {
      label: 'BacMet',
      blurb: 'Experimentally confirmed metal & biocide resistance genes (Hg, As, Cu, Zn/Cd/Co, Cr, Ni, Ag, QACs …)',
      cite: 'Pal et al. 2014, Nucleic Acids Res. 42:D737 (BacMet2 EXP)',
      cutoff: 'DIAMOND blastp, ≥90% identity, ≥30 aa alignment (BacMet-Scan defaults)',
      color: '#94a3b8',
    },
    'PAZy': {
      label: 'PAZy',
      blurb: 'Experimentally validated plastic-active enzymes (PET, PA, PUR, PLA, PBAT, PHA …)',
      cite: 'Buchholz et al. 2022 Proteins; Feuerriegel et al. 2026 NAR',
      cutoff: 'DIAMOND blastp, ≥50% identity, ≥50% query+subject cover, E≤1e-5',
      color: '#f472b6',
    },
    'RDase': {
      label: 'Reductive dehalogenases',
      blurb: 'rdhA — organohalide respiration (chlorinated ethenes, chlorophenols, PCBs)',
      cite: 'NCBIfam TIGR02486 + Pfam PF13486',
      cutoff: 'hmmsearch: TIGR02486 trusted cutoff (rdhA); PF13486 gathering cutoff (rdhA-like)',
      color: '#2dd4bf',
    },
    'CARD': {
      label: 'CARD (ARGs)',
      blurb: 'Antibiotic resistance genes — emerging-contaminant indicator',
      cite: 'Alcock et al. 2023, Nucleic Acids Res. 51:D690 (CARD 4.0.1)',
      cutoff: 'DIAMOND blastp, per-model curated bit-score cutoffs (RGI "Strict")',
      color: '#fb923c',
    },
    'HADEG': {
      label: 'HADEG',
      blurb: 'Aerobic hydrocarbon & plastic degradation and biosurfactant genes',
      cite: 'Rojas-Vargas et al. 2023, Comput. Biol. Chem. 107:107966',
      cutoff: 'DIAMOND blastp, ≥50% identity, ≥50% query+subject cover, E≤1e-5',
      color: '#34d399',
    },
  };
  const DB_FALLBACK_COLORS = ['#a78bfa', '#fb923c', '#38bdf8', '#f472b6', '#2dd4bf', '#818cf8'];
  function dbInfo(db) {
    if (DB_INFO[db]) return DB_INFO[db];
    const i = Object.keys(payload?.dbs || {}).indexOf(db);
    return { label: db, blurb: '', cite: '', cutoff: payload?.dbs?.[db]?.params || '',
             color: DB_FALLBACK_COLORS[Math.max(0, i) % DB_FALLBACK_COLORS.length] };
  }

  // Featured quick-picks: genes worth one click. Each `match` tests a feature
  // [db, category, pathway, gene]; only picks with data are shown.
  const FEATURED = [
    { id: 'hgcA',   label: 'hgcA (Hg methylation)', title: 'HgcA with the conserved cap-helix motif (Hg-MATE) — methylmercury production potential',
      match: f => f[0] === 'Hg-MATE' && f[3] === 'hgcA' },
    { id: 'hgcAB',  label: 'hgcAB same read', title: 'Reads carrying both hgcA and hgcB (counted in reads)',
      match: f => f[0] === 'Hg-MATE' && f[3] === 'hgcAB pair (same read)' },
    { id: 'crass',  label: 'crAssphage (sewage)', title: 'Reads mapping to crAssphage — human fecal / sewage marker (per million reads)',
      match: f => f[0] === 'Mapping markers' && /crass/i.test(f[2] + f[3]) },
    { id: 'alkB',   label: 'alkB / CYP153',     title: 'Medium-chain alkane hydroxylases (CANT-HYD AlkB + CYP153)',
      match: f => f[0] === 'CANT-HYD' && (f[3] === 'alkB' || f[3] === 'CYP153') },
    { id: 'almA',   label: 'almA / ladA',       title: 'Long-chain alkane monooxygenases (CANT-HYD)',
      match: f => f[0] === 'CANT-HYD' && ['almA', 'ladA', 'ladB'].includes(f[3]) },
    { id: 'rhd',    label: 'Aromatic RHDs',     title: 'Aromatic ring-hydroxylating dioxygenases (CANT-HYD NdoB/C, MAH, non-NdoB)',
      match: f => f[0] === 'CANT-HYD' && /^(ndoB|ndoC|MAH|non-ndoB)/.test(f[3]) },
    { id: 'anaer',  label: 'Anaerobic activation', title: 'Fumarate-addition and carboxylation genes (CANT-HYD assA/bssA/nmsA/abcA/ncA)',
      match: f => f[0] === 'CANT-HYD' && f[2].startsWith('Anaerobic') },
    { id: 'plastic', label: 'Plastics', title: 'Plastic-active enzymes (PAZy, plus HADEG plastics)',
      match: f => (f[0] === 'PAZy') || (f[0] === 'HADEG' && f[1] === 'Plastics') },
    { id: 'mer',    label: 'mer (Hg resistance)', title: 'Mercury resistance genes (BacMet: merA/merB/merR/merT …)',
      match: f => f[0] === 'BacMet' && f[2].startsWith('Mercury') },
    { id: 'rdh',    label: 'rdhA (organohalides)', title: 'Reductive dehalogenases — organohalide respiration',
      match: f => f[0] === 'RDase' },
  ];

  // ---- Selection: a node in the db -> category -> pathway -> gene tree -----
  // `sel` = { db, category, pathway, gene } with nulls for "any"; or a
  // featured pick { featured: id }.
  let sel = $state({ db: null, category: null, pathway: null, gene: null, featured: null });
  let expanded = $state(new Set());

  function selectNode(node) {
    sel = { db: node.db ?? null, category: node.category ?? null, pathway: node.pathway ?? null,
            gene: node.gene ?? null, featured: null };
  }
  function selectFeatured(id) {
    sel = { db: null, category: null, pathway: null, gene: null, featured: sel.featured === id ? null : id };
  }
  function toggleExpand(key) {
    const s = new Set(expanded);
    if (s.has(key)) s.delete(key); else s.add(key);
    expanded = s;
  }

  let features = $derived(payload?.features || []);

  // Read-unit features (same-read gene pairs such as hgcAB, and mapping
  // markers such as crAssphage) count reads, not proteins. A same-read pair
  // only enters a selection when picked explicitly (gene level / featured), so
  // it never double counts the underlying hgcA/hgcB proteins; read-only DBs
  // (mapping markers) behave like any other DB.
  let pairSet = $derived(new Set(payload?.read_features || payload?.pair_features || []));
  let readOnlyDbs = $derived.by(() => {
    const byDb = {};
    features.forEach((f, i) => { (byDb[f[0]] ??= []).push(pairSet.has(i)); });
    return new Set(Object.entries(byDb).filter(([, v]) => v.every(Boolean)).map(([k]) => k));
  });
  function featureMatches(f, s, i) {
    if (i != null && pairSet.has(i) && !readOnlyDbs.has(f[0])
        && !(s.featured ? FEATURED.find(x => x.id === s.featured)?.match(f) : s.gene === f[3])) return false;
    if (s.featured) {
      const fp = FEATURED.find(x => x.id === s.featured);
      return fp ? fp.match(f) : false;
    }
    return (!s.db || f[0] === s.db) && (!s.category || f[1] === s.category)
        && (!s.pathway || f[2] === s.pathway) && (!s.gene || f[3] === s.gene);
  }

  // Indices of features in the current selection
  let selFeatureIdx = $derived.by(() => {
    const out = new Set();
    features.forEach((f, i) => { if (featureMatches(f, sel, i)) out.add(i); });
    return out;
  });

  let selLabel = $derived.by(() => {
    if (sel.featured) return FEATURED.find(x => x.id === sel.featured)?.label || sel.featured;
    const parts = [sel.db ? dbInfo(sel.db).label : 'All databases', sel.category, sel.pathway, sel.gene].filter(Boolean);
    return parts.join(' › ');
  });

  // ---- Normalisation -------------------------------------------------------
  // per_mprot: hits per million proteins searched (per DB, summed across DBs)
  // per_mread: hits per million reads; hits: raw protein count.
  const NORMS = [
    { id: 'per_mprot', label: 'per M proteins', title: 'Hits per million annotated proteins searched (per database)' },
    { id: 'per_mread', label: 'per M reads',    title: 'Hits per million reads in the sample' },
    { id: 'hits',      label: 'raw hits',       title: 'Number of reads carrying a gene in the selection' },
  ];
  let normIdx = $state(0);
  let norm = $derived(NORMS[normIdx]);

  let readCount = $derived.by(() => {
    const m = {};
    for (const s of $samples || []) m[s.id] = s.read_count || 0;
    return m;
  });

  // Value of an arbitrary feature-index set in one sample.
  function sampleValue(sid, idxSet, mode) {
    const ps = payload?.samples?.[sid];
    if (!ps) return null;  // not searched
    let hits = 0, rate = 0;
    for (const [k, n] of Object.entries(ps.counts)) {
      const i = +k;
      if (!idxSet.has(i)) continue;
      hits += n;
      if (pairSet.has(i)) {   // read-unit feature: always per million reads
        if (readCount[sid] > 0) rate += n / readCount[sid] * 1e6;
        continue;
      }
      const np = ps.n_proteins[features[i][0]] || 0;
      if (np > 0) rate += n / np * 1e6;
    }
    if (mode === 'hits') return hits;
    if (mode === 'per_mread') return readCount[sid] > 0 ? hits / readCount[sid] * 1e6 : 0;
    return rate;
  }
  function sampleHits(sid, idxSet) { return sampleValue(sid, idxSet, 'hits'); }

  // Proteins searched for the DB(s) in the current selection.
  function sampleProteins(sid) {
    const ps = payload?.samples?.[sid];
    if (!ps) return null;
    const dbsInSel = new Set([...selFeatureIdx].map(i => features[i][0]));
    if (sel.db) dbsInSel.add(sel.db);
    let n = 0;
    for (const [db, v] of Object.entries(ps.n_proteins)) if (!dbsInSel.size || dbsInSel.has(db)) n = Math.max(n, v);
    return n;
  }

  function fmtNorm(v) {
    if (v == null) return '—';
    if (v === 0) return '0';
    if (norm.id === 'hits') return v.toLocaleString();
    if (v >= 100) return v.toFixed(0);
    if (v >= 1) return v.toFixed(1);
    return v.toFixed(2);
  }

  // ---- Tree (db -> category -> pathway -> gene) with cohort totals --------
  // Totals honour the cart filter so the tree reflects the samples in view.
  let cohortIds = $derived.by(() => {
    if (!payload) return [];
    let ids = Object.keys(payload.samples);
    if ($cartActive && $cartItems.size > 0) ids = ids.filter(id => $cartItems.has(id));
    return ids;
  });

  let tree = $derived.by(() => {
    if (!payload) return [];
    const fTot = new Array(features.length).fill(0);
    const fSamp = new Array(features.length).fill(0);
    for (const sid of cohortIds) {
      for (const [k, n] of Object.entries(payload.samples[sid].counts)) {
        fTot[+k] += n;
        if (n > 0) fSamp[+k] += 1;
      }
    }
    const root = new Map();
    const node = (map, key, extra) => {
      if (!map.has(key)) map.set(key, { key, ...extra, total: 0, idx: [], children: new Map() });
      return map.get(key);
    };
    features.forEach((f, i) => {
      const [db, cat, pw, gene] = f;
      const a = node(root, db, { db, level: 'db', name: db });
      const b = node(a.children, cat, { db, category: cat, level: 'category', name: cat });
      const c = node(b.children, pw, { db, category: cat, pathway: pw, level: 'pathway', name: pw });
      const d = node(c.children, gene, { db, category: cat, pathway: pw, gene, level: 'gene', name: gene });
      for (const n of (pairSet.has(i) && !readOnlyDbs.has(db) ? [d] : [a, b, c, d])) { n.total += fTot[i]; n.idx.push(i); }
      d.nSamples = Math.max(d.nSamples || 0, fSamp[i]);
    });
    // Include searched DBs that have no hits at all
    for (const db of Object.keys(payload.dbs || {})) node(root, db, { db, level: 'db', name: db });
    const toArr = (m, path) => [...m.values()]
      .map(n => ({ ...n, path: path + '|' + n.key, children: toArr(n.children, path + '|' + n.key) }))
      .sort((x, y) => y.total - x.total || x.name.localeCompare(y.name));
    return toArr(root, '');
  });

  // Flattened visible rows of the tree
  let treeRows = $derived.by(() => {
    const out = [];
    const walk = (nodes, depth) => {
      for (const n of nodes) {
        out.push({ ...n, depth });
        if (expanded.has(n.path)) walk(n.children, depth + 1);
      }
    };
    walk(tree, 0);
    return out;
  });

  function isSelected(n) {
    return !sel.featured && (sel.db ?? null) === (n.db ?? null)
      && (sel.category ?? null) === (n.category ?? null)
      && (sel.pathway ?? null) === (n.pathway ?? null)
      && (sel.gene ?? null) === (n.gene ?? null);
  }

  let featuredAvail = $derived(FEATURED.filter(fp => features.some(f => fp.match(f))));

  // ---- Map ------------------------------------------------------------------
  let colorMode = $state('value');   // 'value' | 'flowcell' | metadata column
  let sizeBy = $state('value');      // 'value' | 'fixed'
  let sizeScale = $state(1.0);
  let hideUnsearched = $state(true);

  let metaCols = $derived.by(() => {
    if (!$metadata) return [];
    const cols = new Set();
    for (const m of Object.values($metadata)) for (const k of Object.keys(m)) {
      if (!['lat', 'lon'].includes(k)) cols.add(k);
    }
    return [...cols].sort();
  });
  const COLOR_CYCLE_BASE = ['value', 'flowcell'];
  let colorCycle = $derived([...COLOR_CYCLE_BASE, ...metaCols]);
  function cycleColor() {
    const i = colorCycle.indexOf(colorMode);
    colorMode = colorCycle[(i + 1) % colorCycle.length];
  }
  let colorLabel = $derived(colorMode === 'value' ? `Colour: ${norm.label}` : `Colour: ${colorMode.replace(/_/g, ' ')}`);

  let allMarkers = $derived.by(() => {
    if (!$samples || !$metadata || !payload) return [];
    return $samples.map(s => {
      const m = $metadata[s.id];
      if (!m || m.lat == null || m.lon == null) return null;
      const searched = !!payload.samples[s.id];
      if (!searched && hideUnsearched) return null;
      const v = searched ? sampleValue(s.id, selFeatureIdx, norm.id) : null;
      return {
        id: s.id, lat: m.lat, lon: m.lon, flowcell: s.flowcell,
        read_count: s.read_count,
        value: v ?? 0,
        hits: searched ? sampleHits(s.id, selFeatureIdx) : 0,
        searched,
        ...m,
      };
    }).filter(Boolean);
  });

  let markers = $derived.by(() => {
    let f = allMarkers;
    if ($cartActive && $cartItems.size > 0) f = f.filter(m => $cartItems.has(m.id));
    return f;
  });

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8'];
  function viridisColor(t) {
    const r = Math.round(Math.min(255, Math.max(0, 68 + t * 187)));
    const g = Math.round(Math.min(255, Math.max(0, 1 + t * 230)));
    const b = Math.round(Math.min(255, Math.max(0, 84 + (t < 0.5 ? t * 200 : (1 - t) * 200))));
    return `rgb(${r},${g},${b})`;
  }

  let colorBy = $derived(colorMode);
  let markerColorMap = $derived.by(() => {
    const map = {};
    const vals = markers.map(m => m[colorBy]);
    const numeric = vals.filter(v => typeof v === 'number');
    if (colorBy === 'value' || (numeric.length && numeric.length >= vals.filter(v => v != null).length * 0.8)) {
      // log scale for abundance (heavy-tailed); linear for metadata numbers
      const tf = colorBy === 'value' ? (v => Math.log1p(v)) : (v => v);
      const nums = numeric.map(tf);
      const lo = Math.min(...nums), hi = Math.max(...nums), rg = hi - lo || 1;
      for (const m of markers) {
        const v = m[colorBy];
        if (colorBy === 'value' && (!m.searched)) { map[m.id] = '#475569'; continue; }
        if (colorBy === 'value' && v === 0) { map[m.id] = '#334155'; continue; }
        map[m.id] = typeof v === 'number' ? viridisColor((tf(v) - lo) / rg) : '#475569';
      }
      return map;
    }
    const cats = [...new Set(vals.filter(v => v != null).map(String))].sort();
    cats.forEach((c, i) => { map[c] = PALETTE[i % PALETTE.length]; });
    return map;
  });

  let noGeoData = $derived(!$metadata || allMarkers.length === 0);
  let lassoIds = $state(null);
  function addLassoToCart() {
    if (lassoIds) { for (const id of lassoIds) addToCart(id); lassoIds = null; }
  }

  // ---- Detail panel: breakdown of the selection for the clicked sample ------
  let selectedDetail = $derived.by(() => {
    if (!$selectedSample || !$samples) return null;
    const s = $samples.find(x => x.id === $selectedSample);
    if (!s) return null;
    return { ...s, ...($metadata?.[s.id] || {}) };
  });

  // Next tree level below the selection: db -> category -> pathway -> gene
  let breakdownLevel = $derived(sel.featured ? 3 : sel.gene ? 3 : sel.pathway ? 3 : sel.category ? 2 : sel.db ? 1 : 0);

  let selectedBreakdown = $derived.by(() => {
    if (!selectedDetail || !payload) return null;
    const ps = payload.samples[selectedDetail.id];
    if (!ps) return { searched: false };
    const groups = new Map();
    for (const [k, n] of Object.entries(ps.counts)) {
      const i = +k;
      if (!selFeatureIdx.has(i)) continue;
      const f = features[i];
      const key = breakdownLevel === 0 ? f[0] : f[breakdownLevel];
      const np = pairSet.has(i) ? (selectedDetail.read_count || 0) : (ps.n_proteins[f[0]] || 0);
      const g = groups.get(key) || { name: key, db: f[0], hits: 0, rate: 0 };
      g.hits += n;
      if (np > 0) g.rate += n / np * 1e6;
      groups.set(key, g);
    }
    const rows = [...groups.values()].map(g => ({
      ...g,
      value: norm.id === 'hits' ? g.hits
           : norm.id === 'per_mread' ? (selectedDetail.read_count > 0 ? g.hits / selectedDetail.read_count * 1e6 : 0)
           : g.rate,
    })).sort((a, b) => b.value - a.value);
    return { searched: true, rows, n_proteins: ps.n_proteins };
  });

  // ---- Tables -----------------------------------------------------------------
  let sampleRows = $derived.by(() => {
    if (!$samples || !payload) return [];
    let rows = $samples.filter(s => payload.samples[s.id]).map(s => {
      const m = $metadata?.[s.id] || {};
      return {
        id: s.id, station: m.station ?? null, date: m.date ?? m.collection_date ?? null,
        read_count: s.read_count,
        n_proteins: sampleProteins(s.id),
        hits: sampleHits(s.id, selFeatureIdx),
        value: sampleValue(s.id, selFeatureIdx, norm.id),
      };
    });
    if ($cartActive && $cartItems.size > 0) rows = rows.filter(r => $cartItems.has(r.id));
    return rows.sort((a, b) => (b.value ?? 0) - (a.value ?? 0));
  });

  let sampleColumns = $derived([
    { key: 'id', label: 'Sample' },
    { key: 'station', label: 'Station' },
    { key: 'date', label: 'Date' },
    { key: 'read_count', label: 'Reads', render: v => v?.toLocaleString() ?? '-' },
    { key: 'n_proteins', label: 'Proteins searched', render: v => v?.toLocaleString() ?? '-' },
    { key: 'hits', label: 'Hits', render: v => v ? v.toLocaleString() : '0' },
    { key: 'value', label: norm.label, render: v => fmtNorm(v) },
  ]);

  // Gene table: every gene in the current selection, across the cohort
  let geneRows = $derived.by(() => {
    if (!payload) return [];
    const agg = new Map();
    for (const sid of cohortIds) {
      const ps = payload.samples[sid];
      for (const [k, n] of Object.entries(ps.counts)) {
        const i = +k;
        if (!selFeatureIdx.has(i) || n <= 0) continue;
        const f = features[i];
        const key = f.join('|');
        const np = pairSet.has(i) ? (readCount[sid] || 0) : (ps.n_proteins[f[0]] || 0);
        const g = agg.get(key) || { key, db: f[0], category: f[1], pathway: f[2], gene: f[3], hits: 0, n_samples: 0, rate_sum: 0 };
        g.hits += n; g.n_samples += 1;
        if (np > 0) g.rate_sum += n / np * 1e6;
        agg.set(key, g);
      }
    }
    const nCohort = Math.max(1, cohortIds.length);
    return [...agg.values()].map(g => ({ ...g, mean_rate: g.rate_sum / nCohort, prevalence: g.n_samples / nCohort }))
      .sort((a, b) => b.hits - a.hits);
  });

  const geneColumns = [
    { key: 'gene', label: 'Gene' },
    { key: 'db', label: 'DB' },
    { key: 'category', label: 'Category' },
    { key: 'pathway', label: 'Pathway' },
    { key: 'hits', label: 'Total hits', render: v => v.toLocaleString() },
    { key: 'n_samples', label: 'Samples' },
    { key: 'prevalence', label: 'Prevalence', render: v => `${(v * 100).toFixed(0)}%` },
    { key: 'mean_rate', label: 'Mean per M prot (reads*)', render: v => v >= 1 ? v.toFixed(1) : v.toFixed(2) },
  ];

  let totals = $derived.by(() => {
    let hits = 0, withHits = 0;
    for (const r of sampleRows) { hits += r.hits || 0; if (r.hits > 0) withHits++; }
    return { hits, withHits, n: sampleRows.length };
  });

  // ---- Mercury (hgcAB) summary panel ----------------------------------------
  function geneIdx(db, gene) {
    const out = new Set();
    features.forEach((f, i) => { if (f[0] === db && f[3] === gene) out.add(i); });
    return out;
  }
  let mercury = $derived.by(() => {
    if (!payload || !features.some(f => f[0] === 'Hg-MATE')) return null;
    const ia = geneIdx('Hg-MATE', 'hgcA'), ib = geneIdx('Hg-MATE', 'hgcB'), ip = geneIdx('Hg-MATE', 'hgcAB pair (same read)');
    let nSearched = 0, nA = 0, nB = 0, nPair = 0, hitsA = 0;
    const top = [];
    for (const sid of cohortIds) {
      const ps = payload.samples[sid];
      if (ps.n_proteins['Hg-MATE'] == null) continue;
      nSearched++;
      const a = sampleHits(sid, ia), b = sampleHits(sid, ib), p = sampleHits(sid, ip);
      if (a > 0) { nA++; hitsA += a; top.push({ sid, a, rate: sampleValue(sid, ia, 'per_mprot') }); }
      if (b > 0) nB++;
      if (p > 0) nPair++;
    }
    top.sort((x, y) => y.rate - x.rate);
    const taxa = {};
    for (const [sid, t] of Object.entries(payload.carriers?.['Hg-MATE|hgcA'] || {})) {
      if (!cohortIds.includes(sid)) continue;
      for (const [k, n] of Object.entries(t)) taxa[k] = (taxa[k] || 0) + n;
    }
    const topTaxa = Object.entries(taxa).sort((x, y) => y[1] - x[1]).slice(0, 8);
    return { nSearched, nA, nB, nPair, hitsA, top: top.slice(0, 5), topTaxa };
  });

  let selectedCarriers = $derived.by(() => {
    if (!selectedDetail || !payload?.carriers) return [];
    const out = [];
    for (const key of ['Hg-MATE|hgcA', 'Hg-MATE|hgcB']) {
      const t = payload.carriers[key]?.[selectedDetail.id];
      if (t) out.push({ gene: key.split('|')[1], taxa: Object.entries(t).sort((a, b) => b[1] - a[1]) });
    }
    return out;
  });

  const BTN_ON = 'border-cyan-400 bg-cyan-400/10 text-cyan-400';
  const BTN_OFF = 'border-slate-600 text-slate-400 hover:border-slate-500 hover:text-slate-200';
</script>

<div class="space-y-6">
  {#if loadState === 'loading'}
    <div class="text-xs text-slate-500">Loading pollutant marker data…</div>
  {:else if loadState === 'missing'}
    <div class="bg-slate-800/60 border border-slate-700 rounded-md px-3 py-2 text-xs text-slate-400">
      No pollutant marker data for this dataset. Run the pipeline with
      <code class="bg-slate-900 px-1 rounded text-slate-300">--pollutant_dbs</code> (CANT-HYD, HADEG …)
      and an annotator, then re-run the viz preprocess to produce
      <code class="bg-slate-900 px-1 rounded text-slate-300">pollutants.json</code>.
    </div>
  {:else}
    {#if mercury}
      <!-- Methylmercury potential: hgcAB is the key Canadian-lake concern, so it
           gets its own summary above the generic database cards. -->
      <div class="rounded-lg border border-rose-400/40 bg-rose-400/5 p-3 text-xs">
        <div class="flex items-baseline justify-between gap-3 flex-wrap">
          <div>
            <span class="text-rose-300 font-semibold text-sm">Mercury methylation potential (hgcAB)</span>
            <span class="text-slate-400 ml-2">microbes with hgcA + hgcB convert Hg(II) to neurotoxic methylmercury</span>
          </div>
          <div class="flex gap-2">
            <button class="px-2 py-1 rounded-md border transition-colors {sel.featured === 'hgcA' ? BTN_ON : BTN_OFF}" onclick={() => selectFeatured('hgcA')}>Map hgcA</button>
            {#if mercury.nPair > 0}
              <button class="px-2 py-1 rounded-md border transition-colors {sel.featured === 'hgcAB' ? BTN_ON : BTN_OFF}" onclick={() => selectFeatured('hgcAB')}>Map hgcAB reads</button>
            {/if}
          </div>
        </div>
        <div class="flex flex-wrap gap-x-6 gap-y-1 mt-2">
          <span class="text-slate-300">hgcA in <span class="font-mono text-rose-200">{mercury.nA}</span>/{mercury.nSearched} samples</span>
          <span class="text-slate-300">hgcB in <span class="font-mono text-rose-200">{mercury.nB}</span></span>
          <span class="text-slate-300">hgcA+hgcB on one read in <span class="font-mono text-rose-200">{mercury.nPair}</span></span>
          <span class="text-slate-300">hgcA-carrying reads <span class="font-mono text-rose-200">{mercury.hitsA.toLocaleString()}</span></span>
        </div>
        {#if mercury.top.length}
          <div class="mt-2 grid gap-2 md:grid-cols-2">
            <div>
              <div class="text-[10px] uppercase tracking-wide text-slate-500">Highest hgcA (per M proteins)</div>
              {#each mercury.top as t (t.sid)}
                <button class="block font-mono text-slate-300 hover:text-cyan-300" onclick={() => selectedSample.set(t.sid)}>
                  {t.sid} <span class="text-slate-500">{fmtNorm(t.rate)} · {t.a} hits</span>
                </button>
              {/each}
            </div>
            {#if mercury.topTaxa.length}
              <div>
                <div class="text-[10px] uppercase tracking-wide text-slate-500">hgcA carrier taxa (reads, all samples)</div>
                {#each mercury.topTaxa as [k, n] (k)}
                  <div class="text-slate-300">{k} <span class="text-slate-500 font-mono">{n}</span></div>
                {/each}
              </div>
            {/if}
          </div>
        {:else}
          <div class="mt-1 text-slate-500 italic">No hgcA with the cap-helix motif detected in the samples in view.</div>
        {/if}
      </div>
    {/if}

    <!-- Database summary cards -->
    <div class="grid gap-3" style="grid-template-columns: repeat(auto-fit, minmax(16rem, 1fr));">
      {#each Object.entries(payload.dbs) as [db, d] (db)}
        {@const info = dbInfo(db)}
        <button
          class="text-left rounded-lg border p-3 transition-colors bg-slate-800
            {sel.db === db && !sel.category && !sel.featured ? 'border-cyan-400' : 'border-slate-700 hover:border-slate-500'}"
          onclick={() => selectNode({ db })}
        >
          <div class="flex items-baseline justify-between gap-2">
            <span class="font-semibold text-sm" style="color: {info.color}">{info.label}</span>
            <span class="text-[10px] text-slate-500 font-mono">{d.n_samples} samples · {d.n_genes_hit} genes hit</span>
          </div>
          {#if info.blurb}<div class="text-xs text-slate-300 mt-1">{info.blurb}</div>{/if}
          <div class="text-[10px] text-slate-500 mt-1">{info.cutoff}</div>
          {#if info.cite}<div class="text-[10px] text-slate-600 italic">{info.cite}</div>{/if}
        </button>
      {/each}
    </div>

    <!-- Controls -->
    <div class="flex items-center gap-2 flex-wrap text-xs">
      <span class="text-slate-500 uppercase tracking-wide text-[10px]">Selection</span>
      <span class="px-2 py-1 rounded-md border border-amber-400/60 bg-amber-400/10 text-amber-300 font-medium max-w-[28rem] truncate" title={selLabel}>{selLabel}</span>
      {#if sel.db || sel.featured}
        <button class="px-2 py-1 rounded-md border {BTN_OFF}" onclick={() => selectNode({})} title="Select all databases">✕</button>
      {/if}
      {#each featuredAvail as fp (fp.id)}
        <button class="px-2 py-1 rounded-md border transition-colors {sel.featured === fp.id ? BTN_ON : BTN_OFF}"
          onclick={() => selectFeatured(fp.id)} title={fp.title}>{fp.label}</button>
      {/each}

      <span class="text-slate-600 mx-1">|</span>
      <button class="px-3 py-1 rounded-md border transition-colors {BTN_ON}"
        onclick={() => normIdx = (normIdx + 1) % NORMS.length} title={norm.title}>
        {norm.label} &#x25BE;
      </button>
      <button class="px-3 py-1 rounded-md border transition-colors {colorMode === 'value' ? BTN_ON : BTN_OFF}"
        onclick={cycleColor} title={`Marker colour — click to cycle: ${colorCycle.join(' → ')}`}>
        {colorLabel} &#x25BE;
      </button>
      <button class="px-3 py-1 rounded-md border transition-colors {sizeBy === 'value' ? BTN_ON : BTN_OFF}"
        onclick={() => sizeBy = sizeBy === 'value' ? 'fixed' : 'value'} title="Marker size">
        Size: {sizeBy === 'value' ? norm.label : 'fixed'}
      </button>
      <div class="text-slate-400 flex items-center gap-1">
        <div class="single-range relative w-16 h-5 flex items-center">
          <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
          <input type="range" min="0.3" max="10" step="0.1" bind:value={sizeScale} />
        </div>
        <span class="text-slate-500 w-8 font-mono">{sizeScale.toFixed(1)}x</span>
      </div>
      <label class="flex items-center gap-1 text-slate-400" title="Hide samples whose proteins were not searched (no annotation)">
        <input type="checkbox" bind:checked={hideUnsearched} class="rounded bg-slate-800 border-slate-600 text-cyan-400" />
        searched only
      </label>
      <span class="text-slate-500">
        {markers.length} mapped samples
        {#if $cartActive && $cartItems.size > 0}<span class="text-cyan-400">(cart filtered)</span>{/if}
        <span class="text-slate-600 ml-1">shift-drag to select</span>
      </span>
      {#if lassoIds}
        <button class="px-3 py-1 rounded-md border {BTN_ON} hover:bg-cyan-400/20" onclick={addLassoToCart}>
          Add {lassoIds.length} to Cart
        </button>
      {/if}
    </div>

    {#if noGeoData}
      <div class="bg-slate-800/60 border border-slate-700 rounded-md px-3 py-2 text-xs text-slate-400">
        Upload metadata with <code class="bg-slate-900 px-1 rounded text-slate-300">lat</code> + <code class="bg-slate-900 px-1 rounded text-slate-300">lon</code> columns to plot samples on the map; the tables below work without it.
      </div>
    {/if}

    <div class="flex gap-6">
      <!-- Tree: database > category > pathway > gene -->
      <div class="w-72 shrink-0 h-[500px] flex flex-col bg-slate-800 rounded-lg border border-slate-700 p-2 text-xs">
        <div class="flex items-center justify-between px-1 mb-1">
          <span class="text-[10px] uppercase tracking-wide text-slate-500">Database › category › pathway › gene</span>
          <span class="text-[10px] text-slate-500">hits</span>
        </div>
        <ul class="flex-1 overflow-y-auto -mx-1">
          {#each treeRows as n (n.path)}
            <li>
              <div class="flex items-center gap-1 rounded px-1 py-0.5 {isSelected(n) ? 'bg-cyan-400/15' : 'hover:bg-slate-700/50'}"
                   style="padding-left: {0.25 + n.depth * 0.8}rem">
                {#if n.children.length}
                  <button class="w-4 text-slate-500 hover:text-slate-200 shrink-0" onclick={() => toggleExpand(n.path)}
                    title={expanded.has(n.path) ? 'Collapse' : 'Expand'}>{expanded.has(n.path) ? '▾' : '▸'}</button>
                {:else}
                  <span class="w-4 shrink-0"></span>
                {/if}
                <button class="flex-1 min-w-0 text-left truncate
                    {n.level === 'db' ? 'font-semibold' : ''} {n.level === 'gene' ? 'font-mono' : ''}
                    {isSelected(n) ? 'text-cyan-300' : n.total ? 'text-slate-200' : 'text-slate-500'}"
                  style={n.level === 'db' ? `color: ${dbInfo(n.db).color}` : ''}
                  onclick={() => { selectNode(n); if (n.children.length && !expanded.has(n.path)) toggleExpand(n.path); }}
                  title={n.name}>{n.level === 'db' ? dbInfo(n.db).label : n.name}</button>
                <span class="text-[10px] font-mono shrink-0 {n.total ? 'text-slate-400' : 'text-slate-600'}">{n.total.toLocaleString()}</span>
              </div>
            </li>
          {/each}
        </ul>
      </div>

      <div class="flex-1 h-[500px] rounded-lg overflow-hidden border border-slate-700">
        <LeafletMap
          {markers}
          colorMap={markerColorMap}
          {colorBy}
          sizeBy={sizeBy === 'value' ? 'value' : 'fixed'}
          {sizeScale}
          sizeExp={2}
          nudgeMeters={100}
          onMarkerClick={(id) => selectedSample.set(id)}
          onSelect={(ids) => lassoIds = ids}
          exportName={`danaseq_pollutants_${(sel.featured || sel.gene || sel.pathway || sel.category || sel.db || 'all').replace(/[^\w.-]+/g, '_')}`}
        />
      </div>

      <!-- Detail panel -->
      <div class="w-72 shrink-0 bg-slate-800 rounded-lg border border-slate-700 p-4 space-y-3 h-fit max-h-[500px] overflow-y-auto">
        {#if selectedDetail}
          <div class="flex items-center justify-between gap-2">
            <h3 class="text-sm font-semibold text-cyan-400 font-mono truncate">{selectedDetail.id}</h3>
            <button class="text-xs px-2 py-1 rounded border transition-colors shrink-0
                {$cartItems.has(selectedDetail.id) ? 'bg-cyan-400/20 text-cyan-400 border-cyan-400/40' : 'text-slate-400 border-slate-600 hover:text-cyan-400 hover:border-cyan-400/40'}"
              onclick={() => toggleCart(selectedDetail.id)}>
              {$cartItems.has(selectedDetail.id) ? 'In Cart' : '+ Cart'}
            </button>
          </div>
          {#if selectedDetail.station || selectedDetail.date}
            <div class="text-[11px] text-slate-400">{[selectedDetail.station, selectedDetail.date, selectedDetail.depth_m != null ? `${selectedDetail.depth_m} m` : null].filter(Boolean).join(' · ')}</div>
          {/if}
          <div class="border border-amber-400/30 bg-amber-400/5 rounded-md p-2 space-y-1">
            <div class="text-xs text-amber-300 font-medium truncate" title={selLabel}>{selLabel}</div>
            {#if !selectedBreakdown?.searched}
              <div class="text-[10px] text-slate-500 italic">proteins not searched for this sample</div>
            {:else if selectedBreakdown.rows.length === 0}
              <div class="text-[10px] text-slate-500 italic">no hits in this selection</div>
            {:else}
              {@const mx = Math.max(...selectedBreakdown.rows.map(r => r.value), 1e-12)}
              {#each selectedBreakdown.rows.slice(0, 25) as r (r.name + r.db)}
                <div class="text-[11px]">
                  <div class="flex justify-between gap-2">
                    <span class="truncate {breakdownLevel === 3 ? 'font-mono' : ''} text-slate-200" title={r.name}>{breakdownLevel === 0 ? dbInfo(r.name).label : r.name}</span>
                    <span class="font-mono text-slate-400 shrink-0">{fmtNorm(r.value)}</span>
                  </div>
                  <div class="h-1 rounded bg-slate-700 overflow-hidden">
                    <div class="h-1 rounded" style="width: {r.value / mx * 100}%; background: {dbInfo(r.db).color}"></div>
                  </div>
                </div>
              {/each}
              {#if selectedBreakdown.rows.length > 25}
                <div class="text-[10px] text-slate-500">+{selectedBreakdown.rows.length - 25} more</div>
              {/if}
            {/if}
          </div>
          {#each selectedCarriers as c (c.gene)}
            <div class="border border-rose-400/30 bg-rose-400/5 rounded-md p-2">
              <div class="text-xs text-rose-300 font-medium mb-1"><span class="font-mono">{c.gene}</span> carrier reads by taxon</div>
              {#each c.taxa.slice(0, 10) as [k, n] (k)}
                <div class="flex justify-between text-[11px]"><span class="truncate text-slate-200" title={k}>{k}</span><span class="font-mono text-slate-400">{n}</span></div>
              {/each}
            </div>
          {/each}
          <div class="grid grid-cols-2 gap-2 text-xs">
            <div class="text-slate-400">Reads</div><div class="text-slate-200 font-mono">{selectedDetail.read_count?.toLocaleString() ?? '-'}</div>
            {#if selectedBreakdown?.searched}
              {#each Object.entries(selectedBreakdown.n_proteins) as [db, n] (db)}
                <div class="text-slate-400">Proteins ({dbInfo(db).label})</div><div class="text-slate-200 font-mono">{n.toLocaleString()}</div>
              {/each}
            {/if}
            {#if selectedDetail.lat != null}
              <div class="text-slate-400">Lat / Lon</div><div class="text-slate-200 font-mono">{selectedDetail.lat}, {selectedDetail.lon}</div>
            {/if}
          </div>
        {:else}
          <p class="text-slate-500 text-xs">Click a map marker or table row to see the selection broken down for that sample.</p>
          <p class="text-slate-600 text-[11px]">Counts are reads carrying a gene in the selection (each annotated protein's best hit per database; a gene split into several ORFs on one read counts once). Per-M-protein rates divide by the proteins searched; same-read pairs (hgcAB) and mapping markers (crAssphage) are always per million reads.</p>
        {/if}
      </div>
    </div>

    <div class="flex flex-wrap items-baseline gap-x-6 gap-y-1 text-xs px-2">
      <span class="text-slate-400 uppercase tracking-wide text-[10px]">Totals (searched samples):</span>
      <span class="text-slate-300">Samples <span class="text-slate-100 font-mono">{totals.n.toLocaleString()}</span></span>
      <span class="text-slate-300">with hits <span class="text-slate-100 font-mono">{totals.withHits.toLocaleString()}</span></span>
      <span class="text-amber-300">Hits <span class="text-slate-100 font-mono">{totals.hits.toLocaleString()}</span></span>
    </div>

    <div class="grid gap-6 xl:grid-cols-2">
      <div>
        <h3 class="text-xs uppercase tracking-wide text-slate-500 mb-2">Samples — {selLabel}</h3>
        <DataTable
          columns={sampleColumns}
          rows={sampleRows}
          onRowClick={(row) => selectedSample.set(row.id)}
          selectedId={$selectedSample}
          idKey="id"
          maxHeight="360px"
          exportFilename="pollutants_samples"
          actionLabel={(row) => $cartItems.has(row.id) ? 'In Cart' : '+ Cart'}
          actionFn={(row) => toggleCart(row.id)}
          actionStyle={(row) => $cartItems.has(row.id)
            ? 'text-[10px] px-2 py-0.5 rounded border bg-cyan-400/20 text-cyan-400 border-cyan-400/40 transition-colors'
            : 'text-[10px] px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-cyan-400 hover:border-cyan-400/40 transition-colors'}
        />
      </div>
      <div>
        <h3 class="text-xs uppercase tracking-wide text-slate-500 mb-2">Genes — {selLabel}</h3>
        <DataTable
          columns={geneColumns}
          rows={geneRows}
          onRowClick={(row) => selectNode({ db: row.db, category: row.category, pathway: row.pathway, gene: row.gene })}
          idKey="key"
          maxHeight="360px"
          exportFilename="pollutants_genes"
        />
      </div>
    </div>
  {/if}
</div>

<style>
  .single-range input[type="range"] {
    -webkit-appearance: none;
    appearance: none;
    background: transparent;
    position: absolute;
    width: 100%;
    height: 100%;
    margin: 0;
    padding: 0;
    cursor: pointer;
  }
  .single-range input[type="range"]::-webkit-slider-thumb {
    -webkit-appearance: none;
    height: 14px;
    width: 14px;
    border-radius: 50%;
    background: #22d3ee;
    cursor: pointer;
    border: 2px solid #0f172a;
    box-shadow: 0 0 3px rgba(0,0,0,0.4);
  }
  .single-range input[type="range"]::-moz-range-thumb {
    height: 14px;
    width: 14px;
    border-radius: 50%;
    background: #22d3ee;
    cursor: pointer;
    border: 2px solid #0f172a;
    box-shadow: 0 0 3px rgba(0,0,0,0.4);
  }
  .single-range input[type="range"]::-webkit-slider-runnable-track { height: 0; }
  .single-range input[type="range"]::-moz-range-track { height: 0; background: transparent; }
</style>
