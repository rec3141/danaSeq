<script>
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, metadata } from '../stores/data.js';
  import { cartItems, cartActive, toggleCart, addToCart } from '../stores/cart.js';
  import { selectedSample } from '../stores/selection.js';

  let xVar = $state('temperature_c');
  let yVar = $state('diversity');
  let colorBy = $state('station');
  let sizeBy = $state('fixed');
  let sizeScale = $state(1.0);
  let lassoIds = $state(null);
  let boxXVar = $derived(categoricalCols.includes(xVar) ? xVar : categoricalCols[0] || 'station');

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8'];

  // Discover all metadata columns, classify as numeric or categorical
  let allColumns = $derived.by(() => {
    if (!$metadata) return { numeric: [], categorical: [], all: [] };
    const stats = new Map();
    for (const m of Object.values($metadata)) {
      for (const [k, v] of Object.entries(m)) {
        if (k === 'lat' || k === 'lon') continue;
        if (!stats.has(k)) stats.set(k, { num: 0, total: 0 });
        const s = stats.get(k);
        s.total++;
        if (typeof v === 'number') s.num++;
      }
    }
    const numeric = [], categorical = [], all = [];
    for (const [k, s] of stats) {
      const isNum = s.num > s.total * 0.8;
      all.push(k);
      if (isNum) numeric.push(k);
      else categorical.push(k);
    }
    return { numeric: numeric.sort(), categorical: categorical.sort(), all: all.sort() };
  });

  // Continuous columns for axes (metadata numeric + sample numeric)
  let continuousCols = $derived([...allColumns.numeric, 'read_count', 'total_bases', 'gc', 'diversity'].filter((v, i, a) => a.indexOf(v) === i));

  // Categorical columns for color-by and boxplot x-axis
  let categoricalCols = $derived([...allColumns.categorical, 'station', 'dominant_phylum', 'dominant_class', 'flowcell'].filter((v, i, a) => a.indexOf(v) === i));

  // Merge samples with metadata for plotting
  let mergedData = $derived.by(() => {
    if (!$samples) return [];
    return $samples.map(s => {
      const m = $metadata?.[s.id] || {};
      return { ...s, ...m };
    }).filter(d => {
      if ($cartActive && $cartItems.size > 0) return $cartItems.has(d.id);
      return true;
    });
  });

  // Selected sample detail
  let selectedDetail = $derived.by(() => {
    if (!$selectedSample || !$samples) return null;
    const s = $samples.find(s => s.id === $selectedSample);
    if (!s) return null;
    const m = $metadata?.[s.id];
    return { ...s, ...(m || {}) };
  });

  // Color map for categorical coloring
  let colorMap = $derived.by(() => {
    const map = {};
    const vals = [...new Set(mergedData.map(d => d[colorBy]).filter(v => v != null).map(String))].sort();
    vals.forEach((v, i) => { map[v] = PALETTE[i % PALETTE.length]; });
    return map;
  });

  // Scatter traces with color-by and size-by
  let scatterTraces = $derived.by(() => {
    if (!mergedData.length) return [];
    const valid = mergedData.filter(d => d[xVar] != null && d[yVar] != null);
    if (!valid.length) return [];

    // Group by colorBy value for categorical coloring
    const groups = {};
    for (const d of valid) {
      const key = d[colorBy] != null ? String(d[colorBy]) : 'Unknown';
      if (!groups[key]) groups[key] = { x: [], y: [], text: [], ids: [], sizes: [] };
      const g = groups[key];
      g.x.push(d[xVar]);
      g.y.push(d[yVar]);
      g.text.push(`${d.id}<br>${xVar}: ${d[xVar]}<br>${yVar}: ${d[yVar]}<br>${colorBy}: ${key}`);
      g.ids.push(d.id);
      const sVal = typeof d[sizeBy] === 'number' ? d[sizeBy] : 0;
      const sz = (sizeBy === 'fixed' ? 8 : (sVal === 0 ? 4 : Math.max(4, Math.min(24, 4 + Math.log1p(sVal) * 1.2)))) * sizeScale;
      g.sizes.push(sz);
    }

    return Object.entries(groups).sort((a, b) => a[0].localeCompare(b[0])).map(([name, g]) => ({
      type: 'scattergl', mode: 'markers',
      name: name.length > 20 ? name.slice(0, 18) + '..' : name,
      x: g.x, y: g.y, text: g.text, customdata: g.ids,
      marker: {
        size: g.sizes,
        color: colorMap[name] || '#94a3b8',
        opacity: 0.7,
        line: { width: 1, color: '#0f172a' },
      },
      hoverinfo: 'text',
    }));
  });

  let scatterLayout = $derived({
    xaxis: { title: { text: xVar, font: { color: '#94a3b8' } } },
    yaxis: { title: { text: yVar, font: { color: '#94a3b8' } } },
    legend: { bgcolor: 'rgba(15,23,42,0.7)', font: { color: '#94a3b8', size: 10 },
              x: 1, y: 1, xanchor: 'right', yanchor: 'top' },
    height: 420,
  });

  // Boxplot traces: categorical boxXVar vs continuous yVar, with jittered points
  let boxTraces = $derived.by(() => {
    if (!mergedData.length) return [];
    const valid = mergedData.filter(d => d[boxXVar] != null && typeof d[yVar] === 'number');
    if (!valid.length) return [];

    // Group by boxXVar (categorical only)
    const groups = {};
    for (const d of valid) {
      const key = String(d[boxXVar]);
      if (!groups[key]) groups[key] = [];
      groups[key].push(d);
    }

    const keys = Object.keys(groups).sort();

    // Build box traces (whiskers only) + scatter traces (jittered points)
    // Use numeric x positions so jitter math works, then map to category labels via layout
    const traces = [];
    keys.forEach((key, ki) => {
      const pts = groups[key];
      const boxColor = '#22d3ee';
      const shortName = key.length > 15 ? key.slice(0, 13) + '..' : key;

      // Box trace — no points, just box/whisker
      traces.push({
        type: 'box',
        name: shortName,
        y: pts.map(d => d[yVar]),
        boxpoints: false,
        line: { color: boxColor, width: 1.5 },
        fillcolor: 'rgba(34,211,238,0.08)',
        showlegend: false,
        x0: ki, x: pts.map(() => ki),
      });

      // Scatter overlay — jittered points with per-point colorBy colors
      traces.push({
        type: 'scattergl', mode: 'markers',
        name: shortName,
        showlegend: false,
        x: pts.map(() => ki + (Math.random() - 0.5) * 0.3),
        y: pts.map(d => d[yVar]),
        text: pts.map(d => `${d.id}<br>${colorBy}: ${d[colorBy] ?? '-'}<br>${yVar}: ${d[yVar]}`),
        customdata: pts.map(d => d.id),
        marker: {
          size: pts.map(d => {
            const sVal = typeof d[sizeBy] === 'number' ? d[sizeBy] : 0;
            const sz = sizeBy === 'fixed' ? 5 : (sVal === 0 ? 3 : Math.max(3, Math.min(18, 3 + Math.log1p(sVal) * 0.9)));
            return sz * sizeScale;
          }),
          color: pts.map(d => {
            const cv = d[colorBy] != null ? String(d[colorBy]) : 'Unknown';
            return colorMap[cv] || '#94a3b8';
          }),
          opacity: 0.7,
          line: { width: 0 },
        },
        hoverinfo: 'text',
      });
    });
    // Stash keys for layout tick mapping
    traces._keys = keys;
    return traces;
  });

  let boxLayout = $derived.by(() => {
    const keys = boxTraces._keys || [];
    return {
      yaxis: { title: { text: yVar, font: { color: '#94a3b8' } } },
      xaxis: {
        title: { text: boxXVar, font: { color: '#94a3b8' } },
        tickangle: -45,
        tickvals: keys.map((_, i) => i),
        ticktext: keys.map(k => k.length > 15 ? k.slice(0, 13) + '..' : k),
      },
      margin: { b: 100 },
      showlegend: false,
      height: 420,
    };
  });

  // Correlation heatmap
  let corrTraces = $derived.by(() => {
    const numCols = continuousCols.filter(c => mergedData.some(d => d[c] != null));
    if (numCols.length < 2) return [];

    const n = numCols.length;
    const z = [];
    for (let i = 0; i < n; i++) {
      const row = [];
      for (let j = 0; j < n; j++) {
        const pairs = mergedData.filter(d => d[numCols[i]] != null && d[numCols[j]] != null);
        if (pairs.length < 3) { row.push(0); continue; }
        const xi = pairs.map(d => d[numCols[i]]);
        const yj = pairs.map(d => d[numCols[j]]);
        row.push(pearson(xi, yj));
      }
      z.push(row);
    }
    return [{
      type: 'heatmap', z, x: numCols, y: numCols,
      colorscale: 'RdBu', zmin: -1, zmax: 1, reversescale: true,
      hovertemplate: '%{x} vs %{y}: %{z:.2f}<extra></extra>',
    }];
  });

  let corrLayout = $derived({ height: 350, margin: { l: 100, b: 80 } });

  function pearson(x, y) {
    const n = x.length;
    const mx = x.reduce((s, v) => s + v, 0) / n;
    const my = y.reduce((s, v) => s + v, 0) / n;
    let num = 0, dx2 = 0, dy2 = 0;
    for (let i = 0; i < n; i++) {
      const dx = x[i] - mx, dy = y[i] - my;
      num += dx * dy; dx2 += dx * dx; dy2 += dy * dy;
    }
    const denom = Math.sqrt(dx2 * dy2);
    return denom === 0 ? 0 : num / denom;
  }

  function handlePointClick(id) {
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

  // Metadata table
  let tableColumns = $derived.by(() => {
    const cols = [{ key: 'id', label: 'Sample' }];
    for (const c of allColumns.numeric) {
      cols.push({ key: c, label: c, render: v => v != null ? (typeof v === 'number' ? v.toFixed(2) : v) : '-' });
    }
    cols.push({ key: 'diversity', label: 'Shannon H', render: v => v != null ? v.toFixed(2) : '-' });
    return cols;
  });
</script>

<div class="space-y-6">
  {#if allColumns.numeric.length === 0}
    <div class="bg-amber-900/30 border border-amber-700 rounded-lg p-4 text-amber-300">
      <h3 class="font-semibold mb-1">No environmental metadata</h3>
      <p class="text-sm">Provide a metadata TSV with numeric columns (temperature_c, salinity_psu, depth_m, etc.) to enable environmental plots.</p>
    </div>
  {:else}
    <!-- Variable selectors -->
    <div class="flex items-center gap-3 flex-wrap text-xs">
      <select bind:value={xVar}
        class="px-2 py-1 rounded-md border border-slate-600 bg-slate-800 text-slate-300 text-xs focus:border-cyan-400 focus:outline-none cursor-pointer">
        {#each [...categoricalCols, ...continuousCols] as col}
          <option value={col}>X: {col}</option>
        {/each}
      </select>
      <select bind:value={yVar}
        class="px-2 py-1 rounded-md border border-slate-600 bg-slate-800 text-slate-300 text-xs focus:border-cyan-400 focus:outline-none cursor-pointer">
        {#each continuousCols as col}
          <option value={col}>Y: {col}</option>
        {/each}
      </select>
      <select bind:value={colorBy}
        class="px-2 py-1 rounded-md border border-slate-600 bg-slate-800 text-slate-300 text-xs focus:border-cyan-400 focus:outline-none cursor-pointer">
        {#each categoricalCols as col}
          <option value={col}>Color: {col}</option>
        {/each}
      </select>
      <select bind:value={sizeBy}
        class="px-2 py-1 rounded-md border border-slate-600 bg-slate-800 text-slate-300 text-xs focus:border-cyan-400 focus:outline-none cursor-pointer">
        <option value="fixed">Size: Fixed</option>
        {#each continuousCols as col}
          <option value={col}>Size: {col}</option>
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
      <span class="text-slate-500">
        {mergedData.length} samples
        {#if $cartActive && $cartItems.size > 0}
          <span class="text-cyan-400">(cart filtered)</span>
        {/if}
      </span>
      {#if lassoIds}
        <button
          class="px-3 py-1 rounded-md border border-cyan-400 bg-cyan-400/10 text-cyan-400 hover:bg-cyan-400/20 transition-colors"
          onclick={addLassoToCart}
        >
          Add {lassoIds.length} to Cart
        </button>
      {/if}
    </div>

    <!-- Top row: Scatter | Info panel | Boxplot -->
    <div class="grid grid-cols-1 lg:grid-cols-[1fr,auto,1fr] gap-6">
      <!-- Scatter plot -->
      <div class="bg-slate-800 rounded-lg border border-slate-700 p-4 overflow-hidden">
        <h3 class="text-sm font-semibold text-slate-300 mb-2">{xVar} vs {yVar}</h3>
        {#if scatterTraces.length}
          <PlotlyChart traces={scatterTraces} layout={scatterLayout} onclick={handlePointClick} onselect={handleSelect} exportName={`danaseq_env_scatter_${xVar}_vs_${yVar}`} />
        {:else}
          <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">No data for selected variables</div>
        {/if}
      </div>

      <!-- Detail panel (center) -->
      <div class="w-64 shrink-0 bg-slate-800 rounded-lg border border-slate-700 p-4 space-y-3 h-fit max-h-[500px] overflow-y-auto">
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
            <div class="text-slate-400">Reads</div><div class="text-slate-200 font-mono">{selectedDetail.read_count?.toLocaleString() ?? '-'}</div>
            <div class="text-slate-400">Bases</div><div class="text-slate-200 font-mono">{selectedDetail.total_bases ? `${(selectedDetail.total_bases/1e6).toFixed(1)}M` : '-'}</div>
            <div class="text-slate-400">Diversity</div><div class="text-slate-200 font-mono">{selectedDetail.diversity?.toFixed(2) ?? '-'}</div>
            <div class="text-slate-400">Phylum</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.dominant_phylum ?? '-'}</div>
            <div class="text-slate-400">Class</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.dominant_class ?? '-'}</div>
          </div>
          {#if allColumns.numeric.length > 0}
            <div class="border-t border-slate-700 pt-2">
              <div class="text-xs text-slate-400 mb-1">Environment</div>
              <div class="grid grid-cols-2 gap-2 text-xs">
                {#each allColumns.numeric as col}
                  {#if selectedDetail[col] != null}
                    <div class="text-slate-400">{col}</div>
                    <div class="text-slate-200 font-mono">{typeof selectedDetail[col] === 'number' ? selectedDetail[col].toFixed(2) : selectedDetail[col]}</div>
                  {/if}
                {/each}
              </div>
            </div>
          {/if}
          {#if selectedDetail.lat != null}
            <div class="border-t border-slate-700 pt-2">
              <div class="grid grid-cols-2 gap-2 text-xs">
                {#if selectedDetail.station}<div class="text-slate-400">Station</div><div class="text-slate-200 font-mono">{selectedDetail.station}</div>{/if}
                <div class="text-slate-400">Lat</div><div class="text-slate-200 font-mono">{selectedDetail.lat}</div>
                <div class="text-slate-400">Lon</div><div class="text-slate-200 font-mono">{selectedDetail.lon}</div>
                {#if selectedDetail.depth_m != null}<div class="text-slate-400">Depth</div><div class="text-slate-200 font-mono">{selectedDetail.depth_m}m</div>{/if}
              </div>
            </div>
          {/if}
        {:else}
          <p class="text-slate-500 text-xs">Click a scatter point or table row to view details.</p>
        {/if}
      </div>

      <!-- Boxplot -->
      <div class="bg-slate-800 rounded-lg border border-slate-700 p-4 overflow-hidden">
        <h3 class="text-sm font-semibold text-slate-300 mb-2">{boxXVar} vs {yVar} (box)</h3>
        {#if boxTraces.length}
          <PlotlyChart traces={boxTraces} layout={boxLayout} onclick={handlePointClick} onselect={handleSelect} exportName={`danaseq_env_boxplot_${boxXVar}_vs_${yVar}`} />
        {:else}
          <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">
            {#if mergedData.length === 0}No data
            {:else}X variable has too many unique values for boxplot (max 30){/if}
          </div>
        {/if}
      </div>
    </div>

    <!-- Full metadata table -->
    <DataTable
      columns={tableColumns}
      rows={mergedData}
      onRowClick={(row) => { selectedSample.set(row.id); }}
      selectedId={$selectedSample}
      idKey="id"
      maxHeight="300px"
      exportFilename="environmental_metadata"
      actionLabel={(row) => $cartItems.has(row.id) ? 'In Cart' : '+ Cart'}
      actionFn={(row) => toggleCart(row.id)}
      actionStyle={(row) => $cartItems.has(row.id)
        ? 'text-[10px] px-2 py-0.5 rounded border bg-cyan-400/20 text-cyan-400 border-cyan-400/40 transition-colors'
        : 'text-[10px] px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-cyan-400 hover:border-cyan-400/40 transition-colors'}
    />

    <!-- Correlation heatmap (full width, bottom) -->
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4 overflow-hidden">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">Correlation Matrix</h3>
      {#if corrTraces.length}
        <PlotlyChart traces={corrTraces} layout={corrLayout} exportName="danaseq_env_correlation_matrix" />
      {:else}
        <div class="h-[300px] flex items-center justify-center text-slate-500 text-sm">Not enough numeric columns for correlation</div>
      {/if}
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
