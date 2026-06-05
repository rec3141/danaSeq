<script>
  import LeafletMap from '../components/charts/LeafletMap.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, metadata, sampleTaxonomy, readExplorer, loadReadExplorer } from '../stores/data.js';
  import { cartItems, cartActive, toggleCart, addToCart, readSelections, activeReadSelection } from '../stores/cart.js';
  import { selectedSample } from '../stores/selection.js';
  import { taxNav, activeSubTaxa, ancestorInfo, nextRank, prevRank, rankOrder } from '../stores/taxonomy.js';
  import { sampleClusters, sampleClusterColors } from '../stores/clusters.js';
  import TaxonomyDrillNav from '../components/layout/TaxonomyDrillNav.svelte';

  // Per-species AIS payload, loaded lazily on mount. Shape:
  //   { zebra_mussel: { counts: {...}, stats: {...}, identity_bins: {...} } }
  // File-name convention: ./data/ais_<species.id>.json[.gz].
  let aisPayloads = $state({});
  let aisLoaded = $state(new Set());

  async function loadAisCounts(speciesId) {
    if (aisLoaded.has(speciesId)) return;
    aisLoaded.add(speciesId);
    const base = (import.meta.env.BASE_URL || '/').replace(/\/$/, '');
    const url = `${base}/data/ais_${speciesId}.json`;
    try {
      let res = await fetch(url + '.gz');
      let payload;
      if (res.ok && (res.headers.get('Content-Type') || '').match(/json|gzip|octet/)) {
        const ce = res.headers.get('Content-Encoding') || '';
        if (ce.includes('gzip')) {
          payload = await res.json();
        } else {
          const ds = new DecompressionStream('gzip');
          payload = JSON.parse(await new Response(res.body.pipeThrough(ds)).text());
        }
      } else {
        res = await fetch(url);
        if (!res.ok) throw new Error(`HTTP ${res.status}`);
        payload = await res.json();
      }
      aisPayloads = { ...aisPayloads, [speciesId]: payload };
    } catch (e) {
      console.warn(`AIS counts for ${speciesId} not available:`, e.message);
      aisPayloads = { ...aisPayloads, [speciesId]: { counts: {}, stats: {} } };
    }
  }

  let activePayload = $derived(aisPayloads[species.id] || { counts: {}, stats: {} });
  let activeStats = $derived(activePayload.stats || {});
  let identityBins = $derived(activePayload.identity_bins || { lo_pct: 60, hi_pct: 100, n: 40 });
  let posBins = $derived(activePayload.pos_bins || { n: 200, total_length: 0, n_contigs: 0 });
  // `contigs` is only shipped for chromosome-scale refs (≤64 contigs); for
  // fragmented refs (e.g. CDS catalogs) we render a single linear axis.
  let contigs = $derived(activePayload.contigs || null);

  // Lazy-load each species' payload the first time it's selected.
  $effect(() => { loadAisCounts(species.id); });

  // Cycling button helpers
  function cycle(values, current) {
    const idx = values.indexOf(current);
    return values[(idx + 1) % values.length];
  }
  function getLabel(values, labels, current) {
    const idx = values.indexOf(current);
    return idx >= 0 ? labels[idx] : labels[0];
  }

  // Color-by groups
  const infoGroup =   { values: ['flowcell', 'station'], labels: ['Flowcell', 'Station'] };
  // `cluster` is populated by HeatmapView's k-cut via the sampleClusters store
  // and rendered here as a categorical color. Empty until /heatmap is visited.
  // HQ = high-quality (identity > hqCutoff). The slider drives it live;
  // per-marker hq_hits is recomputed from the per-sample identity_hist in the
  // payload so we never need to redeploy to change the cutoff.
  let hqCutoff = $state(90);
  // Identity histogram bins are 1% wide starting at identity_bins.lo_pct
  // (default 60). `bin i` covers [lo + i, lo + i + 1). For cutoff C, the HQ
  // bins are those with lo+i >= C → i >= C - lo. Snap to bin granularity.
  let hqBinIdx = $derived.by(() => {
    const lo = identityBins.lo_pct ?? 60;
    return Math.max(0, Math.min(identityBins.n ?? 40, Math.round(hqCutoff - lo)));
  });
  const metricGroup = $derived({ values: ['ais_hits', 'ais_hq_hits', 'ais_fraction', 'ais_hq_fraction', 'read_count', 'diversity', 'cluster'], labels: ['AIS hits', `HQ hits (>${hqCutoff}%)`, 'AIS / read', 'HQ / read', 'Reads', 'Diversity', 'Heatmap Cluster'] });
  // Default sizing for /ais is the AIS-hit count — that's the whole point of
  // the page. `ais_fraction` is hits / sample_read_count (relative abundance).
  const sizeGroup =   $derived({ values: ['ais_hits', 'ais_hq_hits', 'ais_fraction', 'ais_hq_fraction', 'fixed', 'read_count', 'total_bases', 'diversity'], labels: ['AIS hits', `HQ hits (>${hqCutoff}%)`, 'AIS / read', 'HQ / read', 'Fixed', 'Reads', 'Bases', 'Diversity'] });

  const BW = { species: '7rem', info: '4.5rem', taxonomy: '4.5rem', metric: '4.5rem', metadata: '5rem', size: '4rem' };

  // AIS species selector — additional species will be wired in as their
  // alignment layers come online. Cycle button pattern matches the rest of
  // the controls.
  // `id` must match the mapping reference name (the .idx basename used by
  // run-nanopore-mapping.sh) AND the JSON filename suffix (./data/ais_<id>.json).
  const AIS_SPECIES = [
    { id: 'zebramussel',  label: 'Zebra mussel',  latin: 'Dreissena polymorpha' },
    { id: 'quaggamussel', label: 'Quagga mussel', latin: 'Dreissena rostriformis bugensis' },
    { id: 'sealamprey',   label: 'Sea lamprey',   latin: 'Petromyzon marinus' },
  ];
  let speciesIdx = $state(0);
  let species = $derived(AIS_SPECIES[speciesIdx]);
  function cycleSpecies() { speciesIdx = (speciesIdx + 1) % AIS_SPECIES.length; }

  // AISView defaults: color by AIS abundance, size by AIS hits.
  let colorMode = $state('metric'); // 'info' | 'taxonomy' | 'metric' | 'metadata'
  let infoField = $state('flowcell');
  let metricField = $state('ais_hits');
  let metaField = $state('');

  let sizeBy = $state('ais_hits');
  let sizeScale = $state(1.0);
  let nudgeIdx = $state(3);
  const NUDGE_STEPS = [0, 1, 10, 100, 1000];
  let nudgeMeters = $derived(NUDGE_STEPS[nudgeIdx]);

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8'];

  // Range filter state
  let log2Min = $state(0);
  let log2Max = $state(21);
  let depthMinPct = $state(0);
  let depthMaxPct = $state(100);

  // Auto-detect metadata columns
  let metaColumns = $derived.by(() => {
    if (!$metadata) return [];
    const cols = new Map();
    for (const m of Object.values($metadata)) {
      for (const [k, v] of Object.entries(m)) {
        if (k === 'lat' || k === 'lon') continue;
        if (!cols.has(k)) cols.set(k, { numeric: 0, total: 0 });
        const c = cols.get(k);
        c.total++;
        if (typeof v === 'number') c.numeric++;
      }
    }
    return [...cols.entries()].map(([k, c]) => ({
      key: k,
      continuous: c.numeric > c.total * 0.8,
    })).sort((a, b) => a.key.localeCompare(b.key));
  });

  // Metadata-filter sidebar: same behavior as /samples — categorical-only
  // columns, value list sorted by descending frequency, click to filter the
  // map markers, click again to clear. Only visible (and only applied) in
  // non-taxonomy modes, where the sidebar slot is otherwise empty.
  let categoricalMetaCols = $derived.by(() => {
    if (!$metadata) return [];
    return metaColumns.map(c => c.key).filter(k => {
      const vals = Object.values($metadata)
        .map(m => m[k])
        .filter(v => v !== undefined && v !== null && v !== '');
      if (vals.length === 0) return false;
      return !vals.every(v => typeof v === 'number');
    });
  });

  let filterCol = $state('');
  let filterVal = $state(null);

  $effect(() => {
    if (!filterCol && categoricalMetaCols.length > 0) filterCol = categoricalMetaCols[0];
    if (filterCol && !categoricalMetaCols.includes(filterCol)) {
      filterCol = categoricalMetaCols[0] || '';
      filterVal = null;
    }
  });

  function cycleFilterCol() {
    if (categoricalMetaCols.length === 0) return;
    filterCol = cycle(categoricalMetaCols, filterCol || categoricalMetaCols[0]);
    filterVal = null;
  }

  let filterValueList = $derived.by(() => {
    if (!$samples || !$metadata || !filterCol) return [];
    const counts = new Map();
    for (const s of $samples) {
      const v = $metadata[s.id]?.[filterCol];
      if (v === undefined || v === null || v === '') continue;
      const key = String(v);
      counts.set(key, (counts.get(key) || 0) + 1);
    }
    return [...counts.entries()]
      .sort((a, b) => b[1] - a[1])
      .map(([val, count]) => ({ val, count }));
  });

  function toggleFilterValue(val) {
    filterVal = filterVal === val ? null : val;
  }

  let metaFilterActive = $derived(colorMode !== 'taxonomy' && filterVal !== null && !!filterCol);

  // Dynamic metadata group (columns not already covered by built-in groups)
  let metaGroup = $derived.by(() => {
    const covered = new Set(['flowcell', 'station', 'read_count', 'diversity', 'date', 'total_bases']);
    const vals = [], labs = [];
    for (const col of metaColumns) {
      if (covered.has(col.key)) continue;
      vals.push(col.key);
      labs.push(col.key);
    }
    return { values: vals, labels: labs };
  });

  // Derived colorBy from cycling state
  let colorBy = $derived(
    colorMode === 'info' ? infoField :
    // Taxonomy mode: marker RINGS are driven by the drill-down store; the
    // base marker fill falls back to flowcell for a legible categorical
    // palette while the rings convey the taxonomy breakdown.
    colorMode === 'taxonomy' ? 'flowcell' :
    colorMode === 'metric' ? metricField :
    colorMode === 'metadata' ? (metaField || metaGroup.values[0] || 'station') :
    'flowcell'
  );

  // Is current colorBy continuous?
  let colorIsContinuous = $derived.by(() => {
    const mc = metaColumns.find(c => c.key === colorBy);
    if (mc) return mc.continuous;
    return ['read_count', 'total_bases', 'diversity', 'ais_hits', 'ais_hq_hits', 'ais_fraction', 'ais_hq_fraction'].includes(colorBy);
  });

  // log2 read count filter ceiling
  let log2Ceil = $derived.by(() => {
    if (!allMarkers.length) return 21;
    const maxReads = Math.max(...allMarkers.map(m => m.read_count || 0));
    return Math.ceil(Math.log2(Math.max(1, maxReads)));
  });

  function readFilterLabel(log2val) {
    const v = Math.round(Math.pow(2, log2val));
    if (v >= 1e6) return `${(v/1e6).toFixed(1)}M`;
    if (v >= 1e3) return `${(v/1e3).toFixed(0)}K`;
    return String(v);
  }

  let depthLimits = $derived.by(() => {
    const vals = allMarkers.filter(m => m.depth_m != null).map(m => m.depth_m).sort((a, b) => a - b);
    if (!vals.length) return { min: 0, max: 1 };
    return { min: vals[0], max: vals[vals.length - 1] };
  });

  function pctToVal(pct, limits) {
    return limits.min + (limits.max - limits.min) * pct / 100;
  }

  // Taxonomy drill-down: shared state in taxNav + activeSubTaxa stores.
  // Only use the sub-taxa index when the map is coloring by taxonomy.
  let taxSubTaxa = $derived(colorMode === 'taxonomy' ? $activeSubTaxa : null);

  // All markers (before filters). In taxonomy drill-down mode, each marker
  // carries a `rings` array describing its per-sub-taxon breakdown at the
  // current (level, filter); LeafletMap renders these as nested concentric
  // circles largest-first so every layer is visible.
  //
  // Ring sizing (cross-sample comparable): each ring's `fraction` is the
  // sub-taxon's count divided by the SAMPLE'S TOTAL reads (s.read_count),
  // not the within-filter sum. This keeps denominators the same across
  // samples so 1% S. enterica in Sample A and 10% S. enterica in Sample B
  // produce visibly different rings (~3.2× radius difference). `radiusScale`
  // is sqrt(fraction / globalMaxFraction) — the single largest (sample,
  // taxon) fraction in the current view fills the marker, everything else
  // scales proportionally by area.
  let allMarkers = $derived.by(() => {
    if (!$samples || !$metadata) return [];

    // First pass: build per-sample fractions, using the sample's own total
    // read count as the denominator for cross-sample comparability.
    const speciesCounts = activePayload.counts || {};
    const speciesStats = activePayload.stats || {};
    const markers = $samples.map(s => {
      const m = $metadata[s.id];
      if (!m || m.lat == null || m.lon == null) return null;
      const ais_hits = speciesCounts[s.id] ?? 0;
      const ais_fraction = s.read_count > 0 ? ais_hits / s.read_count : 0;
      const stat = speciesStats[s.id];
      // Recompute hq_hits live from the identity histogram so the slider
      // takes effect without re-running preprocess.
      let ais_hq_hits = 0;
      if (stat?.identity_hist) {
        for (let i = hqBinIdx; i < stat.identity_hist.length; i++) ais_hq_hits += stat.identity_hist[i];
      } else {
        ais_hq_hits = stat?.hq_hits ?? 0;
      }
      const ais_hq_fraction = s.read_count > 0 ? ais_hq_hits / s.read_count : 0;
      const marker = {
        id: s.id,
        lat: m.lat,
        lon: m.lon,
        read_count: s.read_count,
        total_bases: s.total_bases,
        flowcell: s.flowcell,
        diversity: s.diversity,
        cluster: $sampleClusters?.[s.id] ?? null,
        ais_hits,
        ais_fraction,
        ais_hq_hits,
        ais_hq_fraction,
        ais_mean_identity: stat?.mean_identity ?? null,
        ais_mean_mapq: stat?.mean_mapq ?? null,
        ais_identity_hist: stat?.identity_hist ?? null,
        ...m,
      };
      if (taxSubTaxa) {
        // ALWAYS set `rings` when drill-down is active — even to an empty array.
        // Empty means "sample has zero of the current filter/isolated taxon";
        // LeafletMap renders that as a small hollow placeholder so you can see
        // which samples are truly zero (vs. falling through to the default
        // single-color marker, which would look identical to non-zero samples).
        const agg = taxSubTaxa.perSample[s.id] || {};
        const denom = Math.max(1, s.read_count || 0);
        const entries = Object.entries(agg).sort((a, b) => b[1] - a[1]);
        let rings = entries.map(([name, count]) => ({
          name,
          color: taxSubTaxa.colorMap[name],
          count,
          fraction: count / denom,
        }));
        if ($taxNav.isolated) {
          rings = rings.filter(r => r.name === $taxNav.isolated);
        }
        marker.rings = rings;   // may be empty — that's intentional
      }
      return marker;
    }).filter(Boolean);

    // Second pass: cross-sample size normalization. Find the max fraction of
    // any (sample, sub-taxon) so we can scale ring radii against it.
    let maxFrac = 0;
    for (const m of markers) {
      if (!m.rings) continue;
      for (const r of m.rings) if (r.fraction > maxFrac) maxFrac = r.fraction;
    }
    if (maxFrac > 0) {
      for (const m of markers) {
        if (!m.rings) continue;
        for (const r of m.rings) r.radiusScale = Math.sqrt(r.fraction / maxFrac);
      }
    }
    return markers;
  });

  // ---- Read-selection overlay ----
  // When the cart filter is on AND a selection is active, attach a
  // `selectionRing` to each marker describing its share of that selection's
  // reads. Sized as sqrt(count / globalMaxCount) so the busiest sample's ring
  // hits a fixed outer radius and others shrink proportionally by area.
  let activeSelection = $derived.by(() => {
    if (!$cartActive || !$activeReadSelection) return null;
    return $readSelections.get($activeReadSelection) || null;
  });

  let selectionStats = $derived.by(() => {
    if (!activeSelection) return null;
    return {
      name: $activeReadSelection,
      reads: activeSelection.readIds.size,
      samples: activeSelection.samples.size,
    };
  });

  let selectionMaxCount = $derived.by(() => {
    if (!activeSelection) return 0;
    let mx = 0;
    for (const c of activeSelection.samples.values()) if (c > mx) mx = c;
    return mx;
  });

  function clearActiveSelection() {
    activeReadSelection.set(null);
  }

  // Filtered markers
  let markers = $derived.by(() => {
    let filtered = allMarkers;
    if ($cartActive && $cartItems.size > 0) {
      filtered = filtered.filter(m => $cartItems.has(m.id));
    }
    filtered = filtered.filter(m => {
      const rc = m.read_count || 0;
      if (rc === 0) return log2Min === 0;
      const l = Math.log2(Math.max(1, rc));
      return l >= log2Min && l <= log2Max;
    });
    if (hasDepthData) {
      const dMin = pctToVal(depthMinPct, depthLimits);
      const dMax = pctToVal(depthMaxPct, depthLimits);
      filtered = filtered.filter(m => {
        if (m.depth_m == null) return true;
        return m.depth_m >= dMin && m.depth_m <= dMax;
      });
    }
    if (metaFilterActive) {
      filtered = filtered.filter(m => {
        const v = $metadata?.[m.id]?.[filterCol];
        return v !== undefined && v !== null && String(v) === filterVal;
      });
    }

    // Attach the active read-selection ring (if any) to each marker. Done
    // here (not in allMarkers) so it survives the cartActive sample filter
    // — markers that aren't in the active selection still pass through with
    // selectionRing absent, which LeafletMap interprets as "no ring".
    if (activeSelection && selectionMaxCount > 0) {
      filtered = filtered.map(m => {
        const count = activeSelection.samples.get(m.id) || 0;
        if (count === 0) return m;
        const scale = Math.sqrt(count / selectionMaxCount);
        return {
          ...m,
          selectionRing: { count, scale, color: '#fbbf24' /* amber-400 */ },
        };
      });
    }
    return filtered;
  });

  // Color map
  let markerColorMap = $derived.by(() => {
    if (colorIsContinuous) {
      const vals = markers.map(m => m[colorBy]).filter(v => v != null && typeof v === 'number');
      if (!vals.length) return {};
      const min = Math.min(...vals);
      const max = Math.max(...vals);
      const range = max - min || 1;
      const map = {};
      for (const m of markers) {
        const v = m[colorBy];
        if (v == null || typeof v !== 'number') { map[m.id] = '#475569'; continue; }
        const t = (v - min) / range;
        map[m.id] = viridisColor(t);
      }
      return map;
    }
    const map = {};
    const values = [...new Set(markers.map(m => {
      const v = m[colorBy];
      return v != null ? String(v) : 'Unknown';
    }))].filter(v => v !== 'Unknown').sort();
    // Cluster mode: reuse the heatmap's shared per-label color map so colors
    // are identical across /heatmap, /samples and /map.
    if (colorBy === 'cluster' && Object.keys($sampleClusterColors).length > 0) {
      for (const v of values) map[v] = $sampleClusterColors[v] || '#475569';
      return map;
    }
    values.forEach((v, i) => { map[v] = PALETTE[i % PALETTE.length]; });
    return map;
  });

  function viridisColor(t) {
    const r = Math.round(Math.min(255, Math.max(0, 68 + t * 187)));
    const g = Math.round(Math.min(255, Math.max(0, 1 + t * 230)));
    const b = Math.round(Math.min(255, Math.max(0, 84 + (t < 0.5 ? t * 200 : (1 - t) * 200))));
    return `rgb(${r},${g},${b})`;
  }

  let noGeoData = $derived(!$metadata || allMarkers.length === 0);
  let hasDepthData = $derived(allMarkers.some(m => m.depth_m != null));

  // ---- Search state (taxon/gene/product, plus sample fields) ----
  let searchQuery = $state('');
  let searchDebounced = $state('');
  let searchTimer;
  $effect(() => {
    const q = searchQuery;
    clearTimeout(searchTimer);
    searchTimer = setTimeout(() => { searchDebounced = q; }, 300);
    return () => clearTimeout(searchTimer);
  });

  // Lazy-load read_explorer on first non-empty query so /map opens fast.
  let readsLoadRequested = $state(false);
  $effect(() => {
    if (searchDebounced.trim() && !readsLoadRequested) {
      readsLoadRequested = true;
      loadReadExplorer();
    }
  });

  const READ_SEARCH_FIELDS = [
    'sample',
    'kraken_domain', 'kraken_phylum', 'kraken_class',
    'kraken_order', 'kraken_family', 'kraken_genus', 'kraken_species',
    'gtdb_domain', 'gtdb_phylum', 'gtdb_class',
    'gtdb_order', 'gtdb_family', 'gtdb_genus', 'gtdb_species',
    'genes', 'products',
  ];

  // For each non-empty query, compute (matchingSampleIds, matchingReadCount).
  // Matches against:
  //   - sample id, flowcell, barcode, station, any metadata value
  //   - any read in the sample carrying a matching tax/genes/products field
  let searchResult = $derived.by(() => {
    const q = searchDebounced.trim().toLowerCase();
    if (!q) return null;
    const matchedSamples = new Set();
    let readMatches = 0;

    // (1) Sample-level match (sample object + metadata).
    if ($samples) {
      for (const s of $samples) {
        const m = $metadata?.[s.id] || {};
        const fields = [s.id, s.flowcell, s.barcode, ...Object.values(m)];
        for (const v of fields) {
          if (v != null && String(v).toLowerCase().includes(q)) {
            matchedSamples.add(s.id);
            break;
          }
        }
      }
    }

    // (2) Read-level match — only when read_explorer is loaded.
    const reads = $readExplorer?.reads;
    if (reads) {
      for (const r of reads) {
        let hit = false;
        for (const f of READ_SEARCH_FIELDS) {
          const v = r[f];
          if (v != null && String(v).toLowerCase().includes(q)) { hit = true; break; }
        }
        if (hit) {
          readMatches++;
          if (r.sample) matchedSamples.add(r.sample);
        }
      }
    }

    return { matchedSamples, readMatches, readsLoaded: !!reads };
  });

  // Sample IDs to dim on the map (= visible markers minus matched samples).
  let dimIds = $derived.by(() => {
    if (!searchResult) return null;
    const dim = new Set();
    for (const m of allMarkers) {
      if (!searchResult.matchedSamples.has(m.id)) dim.add(m.id);
    }
    return dim;
  });

  let selectedDetail = $derived.by(() => {
    if (!$selectedSample || !$samples) return null;
    const s = $samples.find(s => s.id === $selectedSample);
    if (!s) return null;
    const m = $metadata?.[s.id];
    return { ...s, ...(m || {}) };
  });

  let selectedAis = $derived(activeStats[$selectedSample ?? ''] || null);

  // Metadata fields explicitly rendered above; everything else falls through
  // into the generic "Metadata" catch-all section.
  const SHOWN_META_KEYS = new Set([
    'lat', 'lon', 'station', 'date', 'depth_m',
    'temperature_c', 'salinity_psu', 'chl_ug_l', 'do_mg_l',
  ]);

  let extraMeta = $derived.by(() => {
    if (!$selectedSample) return [];
    const m = $metadata?.[$selectedSample];
    if (!m) return [];
    return Object.entries(m)
      .filter(([k, v]) => !SHOWN_META_KEYS.has(k) && v != null && v !== '')
      .sort((a, b) => a[0].localeCompare(b[0]));
  });

  function formatMetaValue(v) {
    if (v == null) return '-';
    if (typeof v === 'number') {
      if (Number.isInteger(v)) return v.toLocaleString();
      if (Math.abs(v) >= 1) return v.toLocaleString(undefined, { maximumFractionDigits: 4 });
      return String(v);
    }
    return String(v);
  }

  // ppm = parts per million of total reads. Easier to read than a tiny
  // fraction like 0.000123.
  function fmtFraction(v) {
    if (typeof v !== 'number' || v === 0) return '-';
    const ppm = v * 1e6;
    if (ppm >= 100) return ppm.toFixed(0) + ' ppm';
    if (ppm >= 1)   return ppm.toFixed(1) + ' ppm';
    return ppm.toFixed(2) + ' ppm';
  }

  // Table includes ALL samples — positionless ones (missing lat/lon) too,
  // so contributors without coordinates still show up in the AIS summary.
  // Cart filter is the only filter the table honors; the map-only filters
  // (read-count slider, depth, metadata) intentionally don't apply here.
  let tableRows = $derived.by(() => {
    if (!$samples) return [];
    const speciesCounts = activePayload.counts || {};
    const speciesStats = activePayload.stats || {};
    let rows = $samples.map(s => {
      const m = $metadata?.[s.id] || {};
      const ais_hits = speciesCounts[s.id] ?? 0;
      const ais_fraction = s.read_count > 0 ? ais_hits / s.read_count : 0;
      const stat = speciesStats[s.id];
      let ais_hq_hits = 0;
      if (stat?.identity_hist) {
        for (let i = hqBinIdx; i < stat.identity_hist.length; i++) ais_hq_hits += stat.identity_hist[i];
      } else {
        ais_hq_hits = stat?.hq_hits ?? 0;
      }
      const ais_hq_fraction = s.read_count > 0 ? ais_hq_hits / s.read_count : 0;
      return {
        id: s.id,
        lat: m.lat ?? null,
        lon: m.lon ?? null,
        read_count: s.read_count,
        total_bases: s.total_bases,
        flowcell: s.flowcell,
        diversity: s.diversity,
        ais_hits,
        ais_fraction,
        ais_hq_hits,
        ais_hq_fraction,
        ...m,
      };
    });
    if ($cartActive && $cartItems.size > 0) {
      rows = rows.filter(r => $cartItems.has(r.id));
    }
    return rows;
  });

  // $derived so the "HQ hits (>N%)" column header updates with the slider.
  let tableColumns = $derived([
    { key: 'id', label: 'Sample' },
    { key: 'station', label: 'Station' },
    { key: 'lat', label: 'Lat', render: v => typeof v === 'number' ? v.toFixed(4) : '-' },
    { key: 'lon', label: 'Lon', render: v => typeof v === 'number' ? v.toFixed(4) : '-' },
    { key: 'depth_m', label: 'Depth (m)' },
    { key: 'date', label: 'Date' },
    { key: 'read_count', label: 'Reads', render: v => v?.toLocaleString() ?? '-' },
    { key: 'ais_hits', label: 'AIS hits', render: v => v ? v.toLocaleString() : '-' },
    { key: 'ais_hq_hits', label: `HQ hits (>${hqCutoff}%)`, render: v => v ? v.toLocaleString() : '-' },
    { key: 'ais_fraction', label: 'AIS / read', render: fmtFraction },
    { key: 'ais_hq_fraction', label: 'HQ / read', render: fmtFraction },
  ]);

  // Totals across the table contents (all samples, or cart-filtered).
  let tableTotals = $derived.by(() => {
    let hits = 0, hq = 0, reads = 0;
    for (const r of tableRows) {
      hits += r.ais_hits || 0;
      hq += r.ais_hq_hits || 0;
      reads += r.read_count || 0;
    }
    return {
      hits,
      hq,
      reads,
      n: tableRows.length,
      fraction: reads > 0 ? hits / reads : 0,
      hq_fraction: reads > 0 ? hq / reads : 0,
    };
  });

  let lassoIds = $state(null);

  function handleMarkerClick(id) {
    selectedSample.set(id);
  }

  function handleSelect(ids) {
    lassoIds = ids;
  }

  function addLassoToCart() {
    if (lassoIds) {
      for (const id of lassoIds) addToCart(id);
      lassoIds = null;
    }
  }
</script>

<div class="space-y-6">
  {#if noGeoData}
    <div class="bg-slate-800/60 border border-slate-700 rounded-md px-3 py-2 text-xs text-slate-400">
      Upload metadata with <code class="bg-slate-900 px-1 rounded text-slate-300">lat</code> + <code class="bg-slate-900 px-1 rounded text-slate-300">lon</code> columns to plot samples — or click <span class="text-cyan-400">Template</span> in the NavBar for a starter file.
    </div>
  {/if}
    <!-- Controls -->
    <div class="flex items-center gap-2 flex-wrap text-xs">
      <!-- AIS species selector -->
      <button
        class="px-3 py-1 rounded-md border border-emerald-400 bg-emerald-400/10 text-emerald-300 transition-colors text-center"
        style="min-width: {BW.species}"
        onclick={cycleSpecies}
        title={`AIS species — click to cycle (${AIS_SPECIES.length} available)`}
      >
        <span class="italic">{species.latin}</span> &#x25BE;
      </button>

      <!-- HQ identity cutoff slider. Recomputes hq_hits / hq_fraction live
           from each sample's identity histogram. Range matches the histogram
           bins (60-100%), step 1% to snap to bin boundaries. -->
      <div class="flex items-center gap-1 text-rose-300" title="High-quality identity cutoff. Hits above this threshold are counted as HQ in the totals, table, and histogram cutoff line.">
        <span>HQ &gt;</span>
        <span class="text-slate-200 font-mono w-8 text-right">{hqCutoff}%</span>
        <div class="single-range relative w-28 h-5 flex items-center">
          <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
          <div class="absolute h-1 bg-rose-400/40 rounded" style="left: 0; right: {100 - (hqCutoff - 60) * 2.5}%"></div>
          <input type="range" min="60" max="100" step="1" bind:value={hqCutoff} />
        </div>
      </div>

      <span class="text-slate-600 mx-1">|</span>

      <!-- Color cycling buttons -->
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center
          {colorMode === 'info' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        style="min-width: {BW.info}"
        onclick={() => { if (colorMode === 'info') infoField = cycle(infoGroup.values, infoField); else colorMode = 'info'; }}
        title={`Click to cycle: ${infoGroup.labels.join(' → ')}`}
      >
        {colorMode === 'info' ? getLabel(infoGroup.values, infoGroup.labels, infoField) : 'Run Info'} &#x25BE;
      </button>
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center
          {colorMode === 'taxonomy' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        style="min-width: {BW.taxonomy}"
        onclick={() => colorMode = 'taxonomy'}
        title="Color markers by taxonomy — drill down using the panel at left"
      >
        Taxonomy
      </button>
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center
          {colorMode === 'metric' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        style="min-width: {BW.metric}"
        onclick={() => { if (colorMode === 'metric') metricField = cycle(metricGroup.values, metricField); else colorMode = 'metric'; }}
        title={`Click to cycle: ${metricGroup.labels.join(' → ')}`}
      >
        {colorMode === 'metric' ? getLabel(metricGroup.values, metricGroup.labels, metricField) : 'Metric'} &#x25BE;
      </button>
      {#if metaGroup.values.length > 0}
        <button
          class="px-3 py-1 rounded-md border transition-colors text-center
            {colorMode === 'metadata' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
          style="min-width: {BW.metadata}"
          onclick={() => {
            if (colorMode === 'metadata') {
              metaField = cycle(metaGroup.values, metaField || metaGroup.values[0]);
            } else {
              colorMode = 'metadata';
              if (!metaField) metaField = metaGroup.values[0];
            }
          }}
          title={`Click to cycle: ${metaGroup.labels.join(' → ')}`}
        >
          {colorMode === 'metadata' ? getLabel(metaGroup.values, metaGroup.labels, metaField || metaGroup.values[0]) : 'Metadata'} &#x25BE;
        </button>
      {/if}

      <span class="text-slate-600 mx-1">|</span>

      <!-- Size cycling button -->
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
        style="min-width: {BW.size}"
        onclick={() => sizeBy = cycle(sizeGroup.values, sizeBy)}
        title={`Click to cycle: ${sizeGroup.labels.join(' → ')}`}
      >
        {getLabel(sizeGroup.values, sizeGroup.labels, sizeBy)} &#x25BE;
      </button>
      <div class="text-slate-400 flex items-center gap-1">
        <div class="single-range relative w-16 h-5 flex items-center">
          <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
          <input type="range" min="0.3" max="20" step="0.1" bind:value={sizeScale} />
        </div>
        <span class="text-slate-500 w-8 font-mono">{sizeScale.toFixed(1)}x</span>
      </div>

      <span class="text-slate-600 mx-1">|</span>

      <!-- Nudge -->
      <div class="text-slate-400 flex items-center gap-1">
        Nudge
        <div class="single-range relative w-16 h-5 flex items-center">
          <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
          <input type="range" min="0" max="4" step="1" bind:value={nudgeIdx} />
        </div>
        <span class="text-slate-500 w-12 font-mono">{nudgeMeters === 0 ? 'off' : nudgeMeters >= 1000 ? `${nudgeMeters/1000}km` : `${nudgeMeters}m`}</span>
      </div>

      <!-- Reads dual-range slider -->
      <div class="flex items-center gap-1 text-slate-400">
        <span>Reads</span>
        <span class="text-slate-500 w-10 text-right font-mono">{readFilterLabel(log2Min)}</span>
        <div class="dual-range relative w-28 h-5 flex items-center">
          <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
          <div class="absolute h-1 bg-cyan-400/40 rounded" style="left: {log2Min / log2Ceil * 100}%; right: {100 - log2Max / log2Ceil * 100}%"></div>
          <input type="range" min="0" max={log2Ceil} step="0.5" bind:value={log2Min}
            oninput={() => { if (log2Min > log2Max - 0.5) log2Min = log2Max - 0.5; }} />
          <input type="range" min="0" max={log2Ceil} step="0.5" bind:value={log2Max}
            oninput={() => { if (log2Max < log2Min + 0.5) log2Max = log2Min + 0.5; }} />
        </div>
        <span class="text-slate-500 w-10 font-mono">{readFilterLabel(log2Max)}</span>
      </div>

      <!-- Depth dual-range slider -->
      {#if hasDepthData}
        <div class="flex items-center gap-1 text-slate-400">
          <span>Depth</span>
          <span class="text-slate-500 w-10 text-right font-mono">{Math.round(pctToVal(depthMinPct, depthLimits))}m</span>
          <div class="dual-range relative w-28 h-5 flex items-center">
            <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
            <div class="absolute h-1 bg-cyan-400/40 rounded" style="left: {depthMinPct}%; right: {100 - depthMaxPct}%"></div>
            <input type="range" min="0" max="100" step="1" bind:value={depthMinPct}
              oninput={() => { if (depthMinPct > depthMaxPct - 2) depthMinPct = depthMaxPct - 2; }} />
            <input type="range" min="0" max="100" step="1" bind:value={depthMaxPct}
              oninput={() => { if (depthMaxPct < depthMinPct + 2) depthMaxPct = depthMinPct + 2; }} />
          </div>
          <span class="text-slate-500 w-10 font-mono">{Math.round(pctToVal(depthMaxPct, depthLimits))}m</span>
        </div>
      {/if}

      <span class="text-slate-500">
        {markers.length}/{allMarkers.length} samples
        {#if $cartActive && $cartItems.size > 0}
          <span class="text-cyan-400">(cart filtered)</span>
        {/if}
        <span class="text-slate-600 ml-1">shift-drag to select</span>
      </span>

      <span class="text-slate-600 mx-1">|</span>

      <!-- Search bar (taxon, gene, product, metadata, sample) -->
      <div class="relative flex items-center">
        <input
          type="text"
          bind:value={searchQuery}
          placeholder="taxon, gene, product..."
          class="w-48 px-2 py-1 rounded bg-slate-800 border border-slate-600 text-slate-200 text-xs placeholder-slate-500 focus:border-cyan-400 focus:outline-none"
        />
        {#if searchQuery}
          <button
            class="absolute right-1 text-slate-500 hover:text-slate-300"
            onclick={() => { searchQuery = ''; }}
            title="Clear search"
          >
            <svg class="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
              <path d="M18 6L6 18M6 6l12 12" />
            </svg>
          </button>
        {/if}
      </div>
      {#if searchResult}
        <span class="text-amber-400 font-medium">
          {searchResult.matchedSamples.size} samples / {searchResult.readMatches.toLocaleString()} reads matching
          {#if !searchResult.readsLoaded}
            <span class="text-slate-500 italic ml-1">(loading reads…)</span>
          {/if}
        </span>
      {/if}

      {#if selectionStats}
        <span class="ml-2 inline-flex items-center gap-1 px-2 py-0.5 rounded bg-amber-400/10 border border-amber-400/40 text-amber-300 text-[11px]">
          Showing: {selectionStats.name} ({selectionStats.reads.toLocaleString()} reads across {selectionStats.samples} samples)
          <button
            class="ml-1 text-amber-300/70 hover:text-amber-200"
            onclick={clearActiveSelection}
            title="Clear active selection"
          >&#x2715;</button>
        </span>
      {/if}

      {#if lassoIds}
        <button
          class="px-3 py-1 rounded-md border border-cyan-400 bg-cyan-400/10 text-cyan-400 hover:bg-cyan-400/20 transition-colors"
          onclick={addLassoToCart}
        >
          Add {lassoIds.length} to Cart
        </button>
      {:else if markers.length > 0 && markers.length < allMarkers.length}
        <button
          class="px-3 py-1 rounded-md border border-cyan-400 bg-cyan-400/10 text-cyan-400 hover:bg-cyan-400/20 transition-colors"
          onclick={() => { for (const m of markers) addToCart(m.id); }}
        >
          Add {markers.length} to Cart
        </button>
      {/if}
    </div>

    <!-- Taxonomy drill-down nav (taxonomy mode) OR metadata filter sidebar
         (other modes) + Map + Detail panel. -->
    <div class="flex gap-6">
      {#if colorMode === 'taxonomy'}
        <TaxonomyDrillNav heightClass="h-[500px]" />
      {:else}
        <!-- Mirrors the /samples sidebar: cycle header picks a categorical
             metadata column, list shows unique values sorted by descending
             frequency. Click a value to filter map markers; click again to
             clear. -->
        <div class="w-56 shrink-0 h-[500px] flex flex-col bg-slate-800 rounded-lg border border-slate-700 p-3 text-xs">
          {#if categoricalMetaCols.length === 0}
            <p class="text-slate-500 italic text-center mt-4">
              Upload metadata with categorical columns to filter here.
            </p>
          {:else}
            <div class="flex items-center justify-between gap-2 mb-2">
              <button
                class="flex-1 px-2 py-1 rounded border border-slate-600 text-slate-200 hover:border-slate-500 transition-colors text-left font-semibold truncate"
                onclick={cycleFilterCol}
                title={`Cycle: ${categoricalMetaCols.join(' → ')}`}
              >{(filterCol || '').replace(/_/g, ' ') || '—'} &#x25BE;</button>
              {#if filterVal !== null}
                <button
                  class="px-2 py-1 rounded border border-slate-600 text-slate-400 hover:text-rose-400 hover:border-rose-400/40 transition-colors"
                  onclick={() => (filterVal = null)}
                  title="Clear active filter"
                >✕</button>
              {/if}
            </div>
            {#if filterValueList.length === 0}
              <p class="text-slate-500 italic text-center mt-4">No values for this column.</p>
            {:else}
              <div class="text-[10px] uppercase tracking-wide text-slate-500 mb-1">
                {filterValueList.length} value{filterValueList.length === 1 ? '' : 's'}
              </div>
              <ul class="flex-1 overflow-y-auto divide-y divide-slate-700/50 -mx-1">
                {#each filterValueList as { val, count } (val)}
                  <li>
                    <button
                      class="w-full text-left px-2 py-1 rounded transition-colors flex items-baseline justify-between gap-2
                        {filterVal === val ? 'bg-cyan-400/15 text-cyan-300' : 'text-slate-300 hover:bg-slate-700/50'}"
                      onclick={() => toggleFilterValue(val)}
                    >
                      <span class="truncate" title={val}>{val}</span>
                      <span class="text-[10px] text-slate-500 font-mono shrink-0">{count}</span>
                    </button>
                  </li>
                {/each}
              </ul>
            {/if}
          {/if}
        </div>
      {/if}
      <div class="flex-1 h-[500px] rounded-lg overflow-hidden border border-slate-700">
        <LeafletMap
          {markers}
          colorMap={markerColorMap}
          {colorBy}
          {sizeBy}
          {sizeScale}
          sizeExp={2}
          {nudgeMeters}
          onMarkerClick={handleMarkerClick}
          onSelect={handleSelect}
          {dimIds}
          exportName={`danaseq_ais_${species.id}_color-${colorBy}`}
        />
      </div>

      <!-- Detail panel -->
      <div class="w-72 shrink-0 bg-slate-800 rounded-lg border border-slate-700 p-4 space-y-3 h-fit max-h-[500px] overflow-y-auto">
        {#if selectedDetail}
          <div class="flex items-center justify-between gap-2">
            <h3 class="text-sm font-semibold text-cyan-400 font-mono truncate">{selectedDetail.id}</h3>
            <button
              class="text-xs px-2 py-1 rounded border transition-colors shrink-0
                {$cartItems.has(selectedDetail.id)
                  ? 'bg-cyan-400/20 text-cyan-400 border-cyan-400/40'
                  : 'text-slate-400 border-slate-600 hover:text-cyan-400 hover:border-cyan-400/40'}"
              onclick={() => toggleCart(selectedDetail.id)}
            >
              {$cartItems.has(selectedDetail.id) ? 'In Cart' : '+ Cart'}
            </button>
          </div>

          <!-- AIS hit-quality summary + identity histogram. Bars span the
               identity_bins range from the data file (defaults 60-100%).
               Red vertical line marks the HQ identity cutoff. Bins to the
               right of the line are the HQ hits surfaced in the table. -->
          <div class="border border-emerald-400/30 bg-emerald-400/5 rounded-md p-2 space-y-1">
            <div class="flex items-baseline justify-between text-xs">
              <span class="text-emerald-300 italic font-medium truncate">{species.latin}</span>
              <span class="text-slate-400">
                {selectedAis ? selectedAis.n.toLocaleString() : 0} hits
                {#if selectedAis?.hq_hits}<span class="text-rose-300 ml-1">({selectedAis.hq_hits.toLocaleString()} HQ)</span>{/if}
              </span>
            </div>
            {#if selectedAis && selectedAis.n > 0}
              <div class="text-[10px] text-slate-400 flex justify-between">
                <span>mean id <span class="text-slate-200 font-mono">{selectedAis.mean_identity.toFixed(1)}%</span></span>
                <span>mean mapq <span class="text-slate-200 font-mono">{selectedAis.mean_mapq.toFixed(1)}</span></span>
              </div>
              {@const hist = selectedAis.identity_hist}
              {@const maxBin = Math.max(1, ...hist)}
              {@const lo = identityBins.lo_pct}
              {@const hi = identityBins.hi_pct}
              {@const cutoffX = (hqCutoff - lo) / (hi - lo) * 100}
              <svg viewBox="0 0 100 40" preserveAspectRatio="none" class="w-full h-10 mt-1 block">
                {#each hist as count, i}
                  {#if count > 0}
                    {@const binLo = lo + i * (hi - lo) / hist.length}
                    <rect
                      x={i * 100 / hist.length}
                      y={40 - (count / maxBin) * 36}
                      width={Math.max(0.5, 100 / hist.length - 0.2)}
                      height={(count / maxBin) * 36}
                      fill={binLo >= hqCutoff ? '#fb7185' : '#34d399'}
                    ><title>{binLo.toFixed(0)}–{(lo + (i + 1) * (hi - lo) / hist.length).toFixed(0)}% — {count} hits</title></rect>
                  {/if}
                {/each}
                <!-- HQ cutoff line (e.g. 90% identity) -->
                {#if cutoffX > 0 && cutoffX < 100}
                  <line x1={cutoffX} x2={cutoffX} y1="0" y2="40" stroke="#f43f5e" stroke-width="0.6" stroke-dasharray="1,1" />
                {/if}
                <!-- baseline -->
                <line x1="0" x2="100" y1="40" y2="40" stroke="#475569" stroke-width="0.5" />
              </svg>
              <div class="text-[9px] text-slate-500 font-mono flex justify-between -mt-1">
                <span>{lo}%</span>
                <span class="text-rose-400">{hqCutoff}% (HQ)</span>
                <span>{hi}%</span>
              </div>

              {#if selectedAis.pos_hist}
                {@const ph = selectedAis.pos_hist}
                {@const phMax = Math.max(1, ...ph)}
                {@const totalMb = (posBins.total_length / 1e6).toFixed(0)}
                <div class="text-[10px] text-slate-400 mt-2 flex justify-between">
                  <span>genome position</span>
                  <span class="text-slate-500">spread <span class="text-slate-200 font-mono">{ph.filter(v => v > 0).length}</span>/{ph.length} bins</span>
                </div>
                <svg viewBox="0 0 100 24" preserveAspectRatio="none" class="w-full h-6 block">
                  {#each ph as count, i}
                    {#if count > 0}
                      <rect
                        x={i * 100 / ph.length}
                        y={24 - (count / phMax) * 22}
                        width={Math.max(0.3, 100 / ph.length - 0.05)}
                        height={(count / phMax) * 22}
                        fill="#a78bfa"
                      ><title>bin {i + 1}/{ph.length} — {count} hits</title></rect>
                    {/if}
                  {/each}
                  {#if contigs}
                    {#each contigs as c}
                      {#if c.offset > 0}
                        <line x1={c.offset / posBins.total_length * 100} x2={c.offset / posBins.total_length * 100}
                              y1="0" y2="24" stroke="#475569" stroke-width="0.3" />
                      {/if}
                    {/each}
                  {/if}
                  <line x1="0" x2="100" y1="24" y2="24" stroke="#475569" stroke-width="0.5" />
                </svg>
                <div class="text-[9px] text-slate-500 font-mono flex justify-between -mt-1">
                  <span>0</span>
                  <span class="text-slate-500">
                    {#if posBins.n_contigs > 64}
                      {posBins.n_contigs.toLocaleString()} fragments
                    {:else}
                      {posBins.n_contigs} contigs
                    {/if}
                  </span>
                  <span>{totalMb} Mb</span>
                </div>
              {/if}
            {:else}
              <div class="text-[10px] text-slate-500 italic">no hits passing the filter</div>
            {/if}
          </div>

          <div class="grid grid-cols-2 gap-2 text-xs">
            <div class="text-slate-400">Flowcell</div><div class="text-slate-200 font-mono">{selectedDetail.flowcell ?? '-'}</div>
            <div class="text-slate-400">Barcode</div><div class="text-slate-200 font-mono">{selectedDetail.barcode ?? '-'}</div>
            <div class="text-slate-400">Reads</div><div class="text-slate-200 font-mono">{selectedDetail.read_count?.toLocaleString() ?? '-'}</div>
            <div class="text-slate-400">Bases</div><div class="text-slate-200 font-mono">{selectedDetail.total_bases ? `${(selectedDetail.total_bases/1e6).toFixed(1)}M` : '-'}</div>
            <div class="text-slate-400">Avg Length</div><div class="text-slate-200 font-mono">{selectedDetail.avg_length ? Math.round(selectedDetail.avg_length).toLocaleString() : '-'}</div>
            <div class="text-slate-400">Diversity</div><div class="text-slate-200 font-mono">{selectedDetail.diversity?.toFixed(2) ?? '-'}</div>
          </div>
          {#if selectedDetail.lat != null}
            <div class="border-t border-slate-700 pt-2">
              <div class="text-xs text-slate-400 mb-1">Location</div>
              <div class="grid grid-cols-2 gap-2 text-xs">
                <div class="text-slate-400">Lat</div><div class="text-slate-200 font-mono">{selectedDetail.lat}</div>
                <div class="text-slate-400">Lon</div><div class="text-slate-200 font-mono">{selectedDetail.lon}</div>
                {#if selectedDetail.station}<div class="text-slate-400">Station</div><div class="text-slate-200 font-mono">{selectedDetail.station}</div>{/if}
                {#if selectedDetail.date}<div class="text-slate-400">Date</div><div class="text-slate-200 font-mono">{selectedDetail.date}</div>{/if}
                {#if selectedDetail.depth_m != null}<div class="text-slate-400">Depth</div><div class="text-slate-200 font-mono">{selectedDetail.depth_m}m</div>{/if}
              </div>
            </div>
          {/if}
          {#if selectedDetail.temperature_c != null || selectedDetail.salinity_psu != null}
            <div class="border-t border-slate-700 pt-2">
              <div class="text-xs text-slate-400 mb-1">Environment</div>
              <div class="grid grid-cols-2 gap-2 text-xs">
                {#if typeof selectedDetail.temperature_c === 'number'}<div class="text-slate-400">Temp</div><div class="text-slate-200 font-mono">{selectedDetail.temperature_c.toFixed(1)}&deg;C</div>{/if}
                {#if typeof selectedDetail.salinity_psu === 'number'}<div class="text-slate-400">Salinity</div><div class="text-slate-200 font-mono">{selectedDetail.salinity_psu.toFixed(1)} PSU</div>{/if}
                {#if typeof selectedDetail.chl_ug_l === 'number'}<div class="text-slate-400">Chl-a</div><div class="text-slate-200 font-mono">{selectedDetail.chl_ug_l.toFixed(2)} &mu;g/L</div>{/if}
                {#if typeof selectedDetail.do_mg_l === 'number'}<div class="text-slate-400">DO</div><div class="text-slate-200 font-mono">{selectedDetail.do_mg_l.toFixed(1)} mg/L</div>{/if}
              </div>
            </div>
          {/if}
          {#if extraMeta.length > 0}
            <div class="border-t border-slate-700 pt-2">
              <div class="text-xs text-slate-400 mb-1">Metadata</div>
              <div class="grid grid-cols-2 gap-x-2 gap-y-1 text-xs">
                {#each extraMeta as [k, v] (k)}
                  <div class="text-slate-400 break-words">{k.replace(/_/g, ' ')}</div>
                  <div class="text-slate-200 font-mono break-words">{formatMetaValue(v)}</div>
                {/each}
              </div>
            </div>
          {/if}
        {:else}
          <p class="text-slate-500 text-xs">Click a map marker or table row to view details.</p>
        {/if}
      </div>
    </div>

    <!-- AIS totals across the table contents (all samples, or cart-filtered). -->
    <div class="flex flex-wrap items-baseline gap-x-6 gap-y-1 text-xs px-2">
      <span class="text-slate-400 uppercase tracking-wide text-[10px]">Totals (table):</span>
      <span class="text-slate-300">Samples <span class="text-slate-100 font-mono">{tableTotals.n.toLocaleString()}</span></span>
      <span class="text-slate-300">Reads <span class="text-slate-100 font-mono">{tableTotals.reads.toLocaleString()}</span></span>
      <span class="text-emerald-300">AIS hits <span class="text-slate-100 font-mono">{tableTotals.hits.toLocaleString()}</span></span>
      <span class="text-emerald-300">AIS / read <span class="text-slate-100 font-mono">{fmtFraction(tableTotals.fraction)}</span></span>
      <span class="text-rose-300">HQ hits (&gt;{hqCutoff}%) <span class="text-slate-100 font-mono">{tableTotals.hq.toLocaleString()}</span></span>
      <span class="text-rose-300">HQ / read <span class="text-slate-100 font-mono">{fmtFraction(tableTotals.hq_fraction)}</span></span>
    </div>

    <!-- Location table -->
    <DataTable
      columns={tableColumns}
      rows={tableRows}
      onRowClick={(row) => { selectedSample.set(row.id); }}
      selectedId={$selectedSample}
      idKey="id"
      maxHeight="300px"
      exportFilename={`ais_${species.id}_samples`}
      actionLabel={(row) => $cartItems.has(row.id) ? 'In Cart' : '+ Cart'}
      actionFn={(row) => toggleCart(row.id)}
      actionStyle={(row) => $cartItems.has(row.id)
        ? 'text-[10px] px-2 py-0.5 rounded border bg-cyan-400/20 text-cyan-400 border-cyan-400/40 transition-colors'
        : 'text-[10px] px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-cyan-400 hover:border-cyan-400/40 transition-colors'}
    />
</div>

<style>
  .dual-range input[type="range"] {
    -webkit-appearance: none;
    appearance: none;
    background: transparent;
    pointer-events: none;
    position: absolute;
    width: 100%;
    height: 100%;
    margin: 0;
    padding: 0;
  }
  .dual-range input[type="range"]::-webkit-slider-thumb {
    -webkit-appearance: none;
    pointer-events: all;
    height: 14px;
    width: 14px;
    border-radius: 50%;
    background: #22d3ee;
    cursor: pointer;
    border: 2px solid #0f172a;
    box-shadow: 0 0 3px rgba(0,0,0,0.4);
  }
  .dual-range input[type="range"]::-moz-range-thumb {
    pointer-events: all;
    height: 14px;
    width: 14px;
    border-radius: 50%;
    background: #22d3ee;
    cursor: pointer;
    border: 2px solid #0f172a;
    box-shadow: 0 0 3px rgba(0,0,0,0.4);
  }
  .dual-range input[type="range"]::-webkit-slider-runnable-track { height: 0; }
  .dual-range input[type="range"]::-moz-range-track { height: 0; background: transparent; }

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
