<script>
  import LeafletMap from '../components/charts/LeafletMap.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, metadata } from '../stores/data.js';
  import { cartItems, cartActive, toggleCart, addToCart } from '../stores/cart.js';
  import { selectedSample } from '../stores/selection.js';

  let colorBy = $state('dominant_phylum');
  let sizeBy = $state('fixed');
  let sizeScale = $state(1.0);
  let nudgeIdx = $state(3); // index into NUDGE_STEPS
  const NUDGE_STEPS = [0, 1, 10, 100, 1000];
  let nudgeMeters = $derived(NUDGE_STEPS[nudgeIdx]);

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8'];

  // Range filter state (percentile, 0-100)
  let readMinPct = $state(0);
  let readMaxPct = $state(100);
  let depthMinPct = $state(0);
  let depthMaxPct = $state(100);

  // Auto-detect metadata columns for color/size (exclude lat/lon)
  let metaColumns = $derived.by(() => {
    if (!$metadata) return [];
    const cols = new Map(); // key -> { numeric: count, total: count }
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
      continuous: c.numeric > c.total * 0.8, // >80% numeric = continuous
    })).sort((a, b) => a.key.localeCompare(b.key));
  });

  // Build color and size option lists
  let colorOptions = $derived.by(() => {
    const opts = [
      { value: 'dominant_phylum', label: 'Phylum' },
      { value: 'dominant_class', label: 'Class' },
      { value: 'flowcell', label: 'Flowcell' },
      { value: 'station', label: 'Station' },
      { value: 'read_count', label: 'Read Count' },
      { value: 'diversity', label: 'Diversity' },
    ];
    const existing = new Set(opts.map(o => o.value));
    for (const col of metaColumns) {
      if (existing.has(col.key) || col.key === 'date') continue;
      opts.push({ value: col.key, label: col.key });
    }
    return opts;
  });

  let sizeOptions = $derived.by(() => {
    const opts = [
      { value: 'fixed', label: 'Fixed' },
      { value: 'read_count', label: 'Reads' },
      { value: 'total_bases', label: 'Bases' },
      { value: 'diversity', label: 'Diversity' },
    ];
    const existing = new Set(opts.map(o => o.value));
    for (const col of metaColumns) {
      if (!col.continuous || existing.has(col.key)) continue;
      opts.push({ value: col.key, label: col.key });
    }
    return opts;
  });

  // Is current colorBy continuous?
  let colorIsContinuous = $derived.by(() => {
    const mc = metaColumns.find(c => c.key === colorBy);
    if (mc) return mc.continuous;
    return ['read_count', 'total_bases', 'diversity'].includes(colorBy);
  });

  // Compute absolute limits from data
  let readLimits = $derived.by(() => {
    if (!allMarkers.length) return { min: 0, max: 1 };
    const vals = allMarkers.map(m => m.read_count || 0).sort((a, b) => a - b);
    return { min: vals[0], max: vals[vals.length - 1] };
  });

  let depthLimits = $derived.by(() => {
    const vals = allMarkers.filter(m => m.depth_m != null).map(m => m.depth_m).sort((a, b) => a - b);
    if (!vals.length) return { min: 0, max: 1 };
    return { min: vals[0], max: vals[vals.length - 1] };
  });

  function pctToVal(pct, limits) {
    return limits.min + (limits.max - limits.min) * pct / 100;
  }

  function fmtVal(v) {
    if (v >= 1e6) return `${(v / 1e6).toFixed(1)}M`;
    if (v >= 1e3) return `${(v / 1e3).toFixed(0)}K`;
    return Math.round(v).toString();
  }

  // All markers (before filters) — include all metadata fields
  let allMarkers = $derived.by(() => {
    if (!$samples || !$metadata) return [];
    return $samples.map(s => {
      const m = $metadata[s.id];
      if (!m || m.lat == null || m.lon == null) return null;
      return {
        id: s.id,
        lat: m.lat,
        lon: m.lon,
        read_count: s.read_count,
        dominant_phylum: s.dominant_phylum,
        dominant_class: s.dominant_class,
        diversity: s.diversity,
        ...m, // all metadata fields (station, date, depth_m, temperature_c, etc.)
      };
    }).filter(Boolean);
  });

  // Filtered markers (cart + range filters)
  let markers = $derived.by(() => {
    let filtered = allMarkers;

    // Cart filter
    if ($cartActive && $cartItems.size > 0) {
      filtered = filtered.filter(m => $cartItems.has(m.id));
    }

    // Read count filter
    const rMin = pctToVal(readMinPct, readLimits);
    const rMax = pctToVal(readMaxPct, readLimits);
    filtered = filtered.filter(m => {
      const rc = m.read_count || 0;
      return rc >= rMin && rc <= rMax;
    });

    // Depth filter (only if we have depth data)
    if (hasDepthData) {
      const dMin = pctToVal(depthMinPct, depthLimits);
      const dMax = pctToVal(depthMaxPct, depthLimits);
      filtered = filtered.filter(m => {
        if (m.depth_m == null) return true;
        return m.depth_m >= dMin && m.depth_m <= dMax;
      });
    }

    return filtered;
  });

  // Color map: categorical = discrete palette, continuous = not used (LeafletMap handles by value)
  let markerColorMap = $derived.by(() => {
    if (colorIsContinuous) {
      // For continuous, generate a color scale from marker values
      const vals = markers.map(m => m[colorBy]).filter(v => v != null && typeof v === 'number');
      if (!vals.length) return {};
      const min = Math.min(...vals);
      const max = Math.max(...vals);
      const range = max - min || 1;
      const map = {};
      // Use Viridis-like interpolation via simple gradient
      for (const m of markers) {
        const v = m[colorBy];
        if (v == null || typeof v !== 'number') { map[m.id] = '#475569'; continue; }
        const t = (v - min) / range;
        map[m.id] = viridisColor(t);
      }
      return map;
    }
    // Categorical
    const map = {};
    const values = [...new Set(markers.map(m => {
      const v = m[colorBy];
      return v != null ? String(v) : 'Unknown';
    }))].filter(v => v !== 'Unknown').sort();
    values.forEach((v, i) => { map[v] = PALETTE[i % PALETTE.length]; });
    return map;
  });

  function viridisColor(t) {
    // Simplified viridis: dark purple -> blue -> teal -> green -> yellow
    const r = Math.round(Math.min(255, Math.max(0, 68 + t * 187)));
    const g = Math.round(Math.min(255, Math.max(0, 1 + t * 230)));
    const b = Math.round(Math.min(255, Math.max(0, 84 + (t < 0.5 ? t * 200 : (1 - t) * 200))));
    return `rgb(${r},${g},${b})`;
  }

  let highlightIds = $derived.by(() => {
    if ($cartItems.size === 0) return null;
    return $cartItems;
  });

  let noGeoData = $derived(!$metadata || allMarkers.length === 0);
  let hasDepthData = $derived(allMarkers.some(m => m.depth_m != null));

  // Selected sample detail (merges sample data + metadata)
  let selectedDetail = $derived.by(() => {
    if (!$selectedSample || !$samples) return null;
    const s = $samples.find(s => s.id === $selectedSample);
    if (!s) return null;
    const m = $metadata?.[s.id];
    return { ...s, ...(m || {}) };
  });

  const tableColumns = [
    { key: 'id', label: 'Sample' },
    { key: 'station', label: 'Station' },
    { key: 'lat', label: 'Lat', render: v => v?.toFixed(4) ?? '-' },
    { key: 'lon', label: 'Lon', render: v => v?.toFixed(4) ?? '-' },
    { key: 'depth_m', label: 'Depth (m)' },
    { key: 'date', label: 'Date' },
    { key: 'read_count', label: 'Reads', render: v => v?.toLocaleString() ?? '-' },
    { key: 'dominant_phylum', label: 'Dominant' },
  ];

  function handleMarkerClick(id) {
    selectedSample.set(id);
  }
</script>

<div class="space-y-6">
  {#if noGeoData}
    <div class="bg-amber-900/30 border border-amber-700 rounded-lg p-4 text-amber-300">
      <h3 class="font-semibold mb-1">No geographic metadata</h3>
      <p class="text-sm">Upload a metadata TSV with <code class="bg-slate-800 px-1 rounded">lat</code> and <code class="bg-slate-800 px-1 rounded">lon</code> columns using the Metadata button above.</p>
    </div>
  {:else}
    <!-- Controls -->
    <div class="flex items-center gap-3 flex-wrap text-xs">
      <select bind:value={colorBy}
        class="px-2 py-1 rounded-md border border-slate-600 bg-slate-800 text-slate-300 text-xs focus:border-cyan-400 focus:outline-none cursor-pointer">
        {#each colorOptions as opt}
          <option value={opt.value}>Color: {opt.label}</option>
        {/each}
      </select>
      <select bind:value={sizeBy}
        class="px-2 py-1 rounded-md border border-slate-600 bg-slate-800 text-slate-300 text-xs focus:border-cyan-400 focus:outline-none cursor-pointer">
        {#each sizeOptions as opt}
          <option value={opt.value}>Size: {opt.label}</option>
        {/each}
      </select>
      <div class="text-slate-400 flex items-center gap-1">
        Size
        <div class="single-range relative w-16 h-5 flex items-center">
          <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
          <input type="range" min="0.3" max="3" step="0.1" bind:value={sizeScale} />
        </div>
        <span class="text-slate-500 w-8 font-mono">{sizeScale.toFixed(1)}x</span>
      </div>
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
        <span class="text-slate-500 w-10 text-right font-mono">{fmtVal(pctToVal(readMinPct, readLimits))}</span>
        <div class="dual-range relative w-28 h-5 flex items-center">
          <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
          <div class="absolute h-1 bg-cyan-400/40 rounded" style="left: {readMinPct}%; right: {100 - readMaxPct}%"></div>
          <input type="range" min="0" max="100" step="1" bind:value={readMinPct}
            oninput={() => { if (readMinPct > readMaxPct - 2) readMinPct = readMaxPct - 2; }} />
          <input type="range" min="0" max="100" step="1" bind:value={readMaxPct}
            oninput={() => { if (readMaxPct < readMinPct + 2) readMaxPct = readMinPct + 2; }} />
        </div>
        <span class="text-slate-500 w-10 font-mono">{fmtVal(pctToVal(readMaxPct, readLimits))}</span>
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
      </span>
      {#if markers.length > 0 && markers.length < allMarkers.length}
        <button
          class="px-3 py-1 rounded-md border border-cyan-400 bg-cyan-400/10 text-cyan-400 hover:bg-cyan-400/20 transition-colors"
          onclick={() => { for (const m of markers) addToCart(m.id); }}
        >
          Add {markers.length} to Cart
        </button>
      {/if}
    </div>

    <!-- Map + Detail panel -->
    <div class="flex gap-6">
      <div class="flex-1 h-[500px] rounded-lg overflow-hidden border border-slate-700">
        <LeafletMap
          {markers}
          colorMap={markerColorMap}
          {colorBy}
          {sizeBy}
          {sizeScale}
          {nudgeMeters}
          onMarkerClick={handleMarkerClick}
          {highlightIds}
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
          <div class="grid grid-cols-2 gap-2 text-xs">
            <div class="text-slate-400">Flowcell</div><div class="text-slate-200 font-mono">{selectedDetail.flowcell ?? '-'}</div>
            <div class="text-slate-400">Barcode</div><div class="text-slate-200 font-mono">{selectedDetail.barcode ?? '-'}</div>
            <div class="text-slate-400">Reads</div><div class="text-slate-200 font-mono">{selectedDetail.read_count?.toLocaleString() ?? '-'}</div>
            <div class="text-slate-400">Bases</div><div class="text-slate-200 font-mono">{selectedDetail.total_bases ? `${(selectedDetail.total_bases/1e6).toFixed(1)}M` : '-'}</div>
            <div class="text-slate-400">Avg Length</div><div class="text-slate-200 font-mono">{selectedDetail.avg_length ? Math.round(selectedDetail.avg_length).toLocaleString() : '-'}</div>
            <div class="text-slate-400">Diversity</div><div class="text-slate-200 font-mono">{selectedDetail.diversity?.toFixed(2) ?? '-'}</div>
            <div class="text-slate-400">Phylum</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.dominant_phylum ?? '-'}</div>
            <div class="text-slate-400">Class</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.dominant_class ?? '-'}</div>
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
        {:else}
          <p class="text-slate-500 text-xs">Click a map marker or table row to view details.</p>
        {/if}
      </div>
    </div>

    <!-- Location table -->
    <DataTable
      columns={tableColumns}
      rows={markers}
      onRowClick={(row) => { selectedSample.set(row.id); }}
      selectedId={$selectedSample}
      idKey="id"
      maxHeight="300px"
      exportFilename="sample_locations"
      actionLabel={(row) => $cartItems.has(row.id) ? 'In Cart' : '+ Cart'}
      actionFn={(row) => toggleCart(row.id)}
      actionStyle={(row) => $cartItems.has(row.id)
        ? 'text-[10px] px-2 py-0.5 rounded border bg-cyan-400/20 text-cyan-400 border-cyan-400/40 transition-colors'
        : 'text-[10px] px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-cyan-400 hover:border-cyan-400/40 transition-colors'}
    />
  {/if}
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
