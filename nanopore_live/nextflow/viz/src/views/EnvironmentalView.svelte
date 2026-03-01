<script>
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, metadata } from '../stores/data.js';
  import { cartItems, cartActive } from '../stores/cart.js';

  let xVar = $state('temperature_c');
  let yVar = $state('diversity');

  // Discover numeric metadata columns
  let envColumns = $derived.by(() => {
    if (!$metadata) return [];
    const cols = new Set();
    for (const m of Object.values($metadata)) {
      for (const [k, v] of Object.entries(m)) {
        if (typeof v === 'number' && k !== 'lat' && k !== 'lon') cols.add(k);
      }
    }
    return [...cols].sort();
  });

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

  // Scatter trace for x vs y
  let scatterTraces = $derived.by(() => {
    if (!mergedData.length) return [];
    const valid = mergedData.filter(d => d[xVar] != null && d[yVar] != null);
    if (!valid.length) return [];
    return [{
      type: 'scattergl', mode: 'markers',
      x: valid.map(d => d[xVar]),
      y: valid.map(d => d[yVar]),
      text: valid.map(d => `${d.id}<br>${xVar}: ${d[xVar]}<br>${yVar}: ${d[yVar]}`),
      customdata: valid.map(d => d.id),
      marker: { size: 10, color: '#22d3ee', opacity: 0.7, line: { width: 1, color: '#0f172a' } },
      hoverinfo: 'text',
    }];
  });

  let scatterLayout = $derived({
    xaxis: { title: { text: xVar, font: { color: '#94a3b8' } } },
    yaxis: { title: { text: yVar, font: { color: '#94a3b8' } } },
    height: 450,
  });

  // Correlation heatmap
  let corrTraces = $derived.by(() => {
    const numCols = [...envColumns, 'read_count', 'total_bases', 'gc', 'diversity'].filter(c =>
      mergedData.some(d => d[c] != null)
    );
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

  let corrLayout = $derived({ height: 400, margin: { l: 100, b: 80 } });

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

  // Metadata table
  let tableColumns = $derived.by(() => {
    const cols = [{ key: 'id', label: 'Sample' }];
    for (const c of envColumns) {
      cols.push({ key: c, label: c, render: v => v != null ? (typeof v === 'number' ? v.toFixed(2) : v) : '-' });
    }
    cols.push({ key: 'diversity', label: 'Shannon H', render: v => v != null ? v.toFixed(2) : '-' });
    return cols;
  });
</script>

<div class="space-y-6">
  {#if envColumns.length === 0}
    <div class="bg-amber-900/30 border border-amber-700 rounded-lg p-4 text-amber-300">
      <h3 class="font-semibold mb-1">No environmental metadata</h3>
      <p class="text-sm">Provide a metadata TSV with numeric columns (temperature_c, salinity_psu, depth_m, etc.) to enable environmental plots.</p>
    </div>
  {:else}
    <!-- Variable selectors -->
    <div class="flex items-center gap-4 flex-wrap">
      <label class="text-xs text-slate-400">
        X axis:
        <select bind:value={xVar} class="ml-1 bg-slate-800 border border-slate-600 rounded px-2 py-1 text-xs text-slate-200">
          {#each [...envColumns, 'read_count', 'total_bases', 'gc', 'diversity'] as col}
            <option value={col}>{col}</option>
          {/each}
        </select>
      </label>
      <label class="text-xs text-slate-400">
        Y axis:
        <select bind:value={yVar} class="ml-1 bg-slate-800 border border-slate-600 rounded px-2 py-1 text-xs text-slate-200">
          {#each [...envColumns, 'read_count', 'total_bases', 'gc', 'diversity'] as col}
            <option value={col}>{col}</option>
          {/each}
        </select>
      </label>
    </div>

    <div class="grid grid-cols-1 lg:grid-cols-2 gap-6">
      <!-- Scatter plot -->
      <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
        <h3 class="text-sm font-semibold text-slate-300 mb-2">{xVar} vs {yVar}</h3>
        {#if scatterTraces.length}
          <PlotlyChart traces={scatterTraces} layout={scatterLayout} />
        {:else}
          <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">No data for selected variables</div>
        {/if}
      </div>

      <!-- Correlation heatmap -->
      <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
        <h3 class="text-sm font-semibold text-slate-300 mb-2">Correlation Matrix</h3>
        {#if corrTraces.length}
          <PlotlyChart traces={corrTraces} layout={corrLayout} />
        {:else}
          <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">Not enough numeric columns for correlation</div>
        {/if}
      </div>
    </div>

    <!-- Full metadata table -->
    <DataTable columns={tableColumns} rows={mergedData} maxHeight="300px" />
  {/if}
</div>
