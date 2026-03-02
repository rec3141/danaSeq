<script>
  import LeafletMap from '../components/charts/LeafletMap.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, metadata } from '../stores/data.js';
  import { cartItems, cartActive, toggleCart, addToCart } from '../stores/cart.js';
  import { selectedSample } from '../stores/selection.js';

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
  const taxGroup =    { values: ['dominant_phylum', 'dominant_class'], labels: ['Phylum', 'Class'] };
  const metricGroup = { values: ['read_count', 'diversity'], labels: ['Reads', 'Diversity'] };
  const sizeGroup =   { values: ['fixed', 'read_count', 'total_bases', 'diversity'], labels: ['Fixed', 'Reads', 'Bases', 'Diversity'] };

  const BW = { info: '4.5rem', taxonomy: '4.5rem', metric: '4.5rem', metadata: '5rem', size: '4rem' };

  let colorMode = $state('taxonomy'); // 'info' | 'taxonomy' | 'metric' | 'metadata'
  let infoField = $state('flowcell');
  let taxField = $state('dominant_phylum');
  let metricField = $state('read_count');
  let metaField = $state('');

  let sizeBy = $state('fixed');
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

  // Dynamic metadata group (columns not already covered by built-in groups)
  let metaGroup = $derived.by(() => {
    const covered = new Set(['dominant_phylum', 'dominant_class', 'flowcell', 'station', 'read_count', 'diversity', 'date', 'total_bases']);
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
    colorMode === 'taxonomy' ? taxField :
    colorMode === 'metric' ? metricField :
    colorMode === 'metadata' ? (metaField || metaGroup.values[0] || 'station') :
    'dominant_phylum'
  );

  // Is current colorBy continuous?
  let colorIsContinuous = $derived.by(() => {
    const mc = metaColumns.find(c => c.key === colorBy);
    if (mc) return mc.continuous;
    return ['read_count', 'total_bases', 'diversity'].includes(colorBy);
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

  // All markers (before filters)
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
        total_bases: s.total_bases,
        flowcell: s.flowcell,
        dominant_phylum: s.dominant_phylum,
        dominant_class: s.dominant_class,
        diversity: s.diversity,
        ...m,
      };
    }).filter(Boolean);
  });

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
    values.forEach((v, i) => { map[v] = PALETTE[i % PALETTE.length]; });
    return map;
  });

  function viridisColor(t) {
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
    { key: 'lat', label: 'Lat', render: v => typeof v === 'number' ? v.toFixed(4) : '-' },
    { key: 'lon', label: 'Lon', render: v => typeof v === 'number' ? v.toFixed(4) : '-' },
    { key: 'depth_m', label: 'Depth (m)' },
    { key: 'date', label: 'Date' },
    { key: 'read_count', label: 'Reads', render: v => v?.toLocaleString() ?? '-' },
    { key: 'dominant_phylum', label: 'Dominant' },
  ];

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
    <div class="bg-amber-900/30 border border-amber-700 rounded-lg p-4 text-amber-300">
      <h3 class="font-semibold mb-1">No geographic metadata</h3>
      <p class="text-sm">Upload a metadata TSV with <code class="bg-slate-800 px-1 rounded">lat</code> and <code class="bg-slate-800 px-1 rounded">lon</code> columns using the Metadata button above.</p>
    </div>
  {:else}
    <!-- Controls -->
    <div class="flex items-center gap-2 flex-wrap text-xs">
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
        onclick={() => { if (colorMode === 'taxonomy') taxField = cycle(taxGroup.values, taxField); else colorMode = 'taxonomy'; }}
        title={`Click to cycle: ${taxGroup.labels.join(' → ')}`}
      >
        {colorMode === 'taxonomy' ? getLabel(taxGroup.values, taxGroup.labels, taxField) : 'Taxonomy'} &#x25BE;
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
          <input type="range" min="0.3" max="3" step="0.1" bind:value={sizeScale} />
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
          onSelect={handleSelect}
          {highlightIds}
          exportName={`danaseq_map_color-${colorBy}`}
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
