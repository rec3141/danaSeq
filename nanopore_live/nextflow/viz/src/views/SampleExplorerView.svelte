<script>
  import PlotlyScatter from '../components/charts/PlotlyScatter.svelte';
  import StatCard from '../components/layout/StatCard.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, sampleTsne, sampleTaxonomy, metadata, overview } from '../stores/data.js';
  import { cartItems, cartActive, addToCart, toggleCart } from '../stores/cart.js';
  import { selectedSample } from '../stores/selection.js';

  let colorBy = $state('dominant_phylum');
  let sizeBy = $state('total_bases');
  let sizeScale = $state(1.0);
  let lassoIds = $state(null);

  // log2 read count filter
  let log2Min = $state(0);
  let log2Max = $state(21);
  let log2Ceil = $derived.by(() => {
    if (!$samples) return 21;
    const maxReads = Math.max(...$samples.map(s => s.read_count || 0));
    return Math.ceil(Math.log2(Math.max(1, maxReads)));
  });

  function readFilterLabel(log2val) {
    const v = Math.round(Math.pow(2, log2val));
    if (v >= 1e6) return `${(v/1e6).toFixed(1)}M`;
    if (v >= 1e3) return `${(v/1e3).toFixed(0)}K`;
    return String(v);
  }

  // Date filter (initialized once from data range)
  let dateMin = $state(null);
  let dateMax = $state(null);
  let dateRange = $derived.by(() => {
    if (!$samples) return { min: '', max: '' };
    const dates = $samples.map(s => s.start_time?.slice(0, 10)).filter(Boolean).sort();
    if (!dates.length) return { min: '', max: '' };
    return { min: dates[0], max: dates[dates.length - 1] };
  });

  // Default date inputs to full range once data loads
  $effect(() => {
    if (dateRange.min && dateMin === null) dateMin = dateRange.min;
    if (dateRange.max && dateMax === null) dateMax = dateRange.max;
  });

  function passesFilters(s) {
    const rc = s.read_count || 0;
    if (rc === 0) return false;
    const l = Math.log2(Math.max(1, rc));
    if (l < log2Min || l > log2Max) return false;
    if (s.start_time) {
      const d = s.start_time.slice(0, 10);
      if (dateMin && d < dateMin) return false;
      if (dateMax && d > dateMax) return false;
    }
    return true;
  }

  // Embedding mode: cycle between available embeddings (derived from data)
  const EMBED_LABELS = { grid: 'Grid', tsne: 't-SNE' };
  let embedIdx = $state(0);
  let availableModes = $derived.by(() => {
    if (!$sampleTsne) return ['grid'];
    return Object.keys($sampleTsne).filter(k => $sampleTsne[k] && Object.keys($sampleTsne[k]).length > 0);
  });
  let embedMode = $derived(availableModes[embedIdx % availableModes.length] || 'grid');

  function cycleEmbed() {
    embedIdx++;
  }

  const colorOptions = [
    { value: 'dominant_phylum', label: 'Dominant Phylum' },
    { value: 'dominant_class', label: 'Dominant Class' },
    { value: 'flowcell', label: 'Flowcell' },
    { value: 'read_count', label: 'Read Count' },
    { value: 'total_bases', label: 'Total Bases' },
    { value: 'diversity', label: 'Shannon Diversity' },
  ];

  // Get coordinates for current embedding mode
  function getCoords(sampleId) {
    if (!$sampleTsne) return [0, 0];
    const embedData = $sampleTsne[embedMode];
    if (!embedData) return [0, 0];
    return embedData[sampleId] ?? [0, 0];
  }

  // Fixed extents from ALL samples (not filtered) so axes don't jump
  let coordExtents = $derived.by(() => {
    if (!$samples || !$sampleTsne) return null;
    let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
    for (const s of $samples) {
      const [x, y] = getCoords(s.id);
      if (x < xMin) xMin = x;
      if (x > xMax) xMax = x;
      if (y < yMin) yMin = y;
      if (y > yMax) yMax = y;
    }
    if (!isFinite(xMin)) return null;
    return { xMin, xMax, yMin, yMax };
  });

  // Build scatter data from samples + embedding coordinates
  let scatterData = $derived.by(() => {
    if (!$samples || !$sampleTsne) return null;
    const points = $samples.filter(passesFilters).map(s => ({
      ...s,
      tsne_x: getCoords(s.id)[0],
      tsne_y: getCoords(s.id)[1],
    }));

    // Merge metadata if available
    if ($metadata) {
      for (const p of points) {
        const m = $metadata[p.id];
        if (m) Object.assign(p, m);
      }
    }

    // Filter by cart if active
    const filtered = $cartActive && $cartItems.size > 0
      ? points.filter(p => $cartItems.has(p.id))
      : points;

    return { points: filtered };
  });

  // Stable color map
  let colorMap = $derived.by(() => {
    if (!$samples) return {};
    const palette = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                     '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8',
                     '#94a3b8','#d4d4d8','#78716c','#fb7185','#facc15','#4ade80',
                     '#60a5fa','#c084fc','#f97316','#14b8a6','#e879f9','#a3e635'];
    const map = {};
    const values = [...new Set($samples.map(s => s[colorBy] || 'Unknown'))].filter(v => v !== 'Unknown').sort();
    values.forEach((v, i) => { map[v] = palette[i % palette.length]; });
    return map;
  });

  // Selected sample detail
  let selectedDetail = $derived.by(() => {
    if (!$selectedSample || !$samples) return null;
    const s = $samples.find(s => s.id === $selectedSample);
    if (!s) return null;
    const m = $metadata?.[s.id];
    return { ...s, ...(m || {}) };
  });

  // Table columns — separate flowcell + barcode
  const tableColumns = [
    { key: 'flowcell', label: 'Flowcell' },
    { key: 'barcode', label: 'Barcode' },
    { key: 'read_count', label: 'Reads', render: v => v?.toLocaleString() ?? '-' },
    { key: 'total_bases', label: 'Bases', render: v => v ? `${(v/1e6).toFixed(1)}M` : '-' },
    { key: 'avg_length', label: 'Avg Len', render: v => v ? Math.round(v).toLocaleString() : '-' },
    { key: 'n_classified', label: 'Classified', render: v => v?.toLocaleString() ?? '-' },
    { key: 'dominant_class', label: 'Dominant Class' },
    { key: 'diversity', label: 'Shannon H', render: v => v != null ? v.toFixed(2) : '-' },
  ];

  function handleSelect(ids) {
    lassoIds = ids;
  }

  function handleClick(id) {
    selectedSample.set(id);
  }

  function addLassoToCart() {
    if (lassoIds) {
      for (const id of lassoIds) addToCart(id);
      lassoIds = null;
    }
  }
</script>

<div class="space-y-6">
  <!-- Header stats -->
  {#if $overview}
    <div class="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-4">
      <StatCard label="Flowcells" value={$overview.n_flowcells ?? 0} color="purple" />
      <StatCard label="Barcodes" value={$overview.n_samples ?? 0} color="cyan" />
      <StatCard label="Total Reads" value={($overview.total_reads ?? 0).toLocaleString()} color="emerald" />
      <StatCard label="Total Bases" value={$overview.total_bases ? `${($overview.total_bases / 1e9).toFixed(1)} Gbp` : '-'} color="amber" />
      <StatCard label="In Cart" value={$cartItems.size} color="slate" />
    </div>
  {/if}

  <!-- Controls -->
  <div class="flex items-center gap-4 flex-wrap">
    <button
      class="text-xs px-3 py-1.5 rounded border border-slate-600 text-slate-300 transition-colors font-medium
        {availableModes.length > 1 ? 'hover:text-cyan-400 hover:border-cyan-400/40 cursor-pointer' : 'opacity-60 cursor-default'}"
      onclick={availableModes.length > 1 ? cycleEmbed : undefined}
      title={availableModes.length > 1 ? `Cycle layout (${availableModes.map(m => EMBED_LABELS[m] || m).join(' / ')})` : 'Only grid layout available (run sketching for t-SNE)'}
    >
      {EMBED_LABELS[embedMode] || embedMode}
    </button>
    <label class="text-xs text-slate-400">
      Color by:
      <select bind:value={colorBy} class="ml-1 bg-slate-800 border border-slate-600 rounded px-2 py-1 text-xs text-slate-200">
        {#each colorOptions as opt}
          <option value={opt.value}>{opt.label}</option>
        {/each}
      </select>
    </label>
    <label class="text-xs text-slate-400">
      Size by:
      <select bind:value={sizeBy} class="ml-1 bg-slate-800 border border-slate-600 rounded px-2 py-1 text-xs text-slate-200">
        <option value="fixed">Fixed</option>
        <option value="read_count">Read Count</option>
        <option value="total_bases">Total Bases</option>
      </select>
    </label>
    <label class="text-xs text-slate-400 flex items-center gap-2">
      Size:
      <input type="range" min="0.2" max="3" step="0.1" bind:value={sizeScale}
        class="w-20 h-1 accent-cyan-400" />
      <span class="text-slate-500 w-8">{sizeScale.toFixed(1)}x</span>
    </label>
    <!-- log2 read count filter -->
    <label class="text-xs text-slate-400 flex items-center gap-2">
      Reads:
      <span class="text-slate-500 w-10 text-right font-mono">{readFilterLabel(log2Min)}</span>
      <div class="relative w-28 h-4 flex items-center">
        <input type="range" min="0" max={log2Ceil} step="0.5" bind:value={log2Min}
          oninput={(e) => { if (log2Min > log2Max - 0.5) log2Min = log2Max - 0.5; }}
          class="absolute w-full h-1 accent-cyan-400 pointer-events-none [&::-webkit-slider-thumb]:pointer-events-auto [&::-moz-range-thumb]:pointer-events-auto appearance-none bg-slate-700 rounded" />
        <input type="range" min="0" max={log2Ceil} step="0.5" bind:value={log2Max}
          oninput={(e) => { if (log2Max < log2Min + 0.5) log2Max = log2Min + 0.5; }}
          class="absolute w-full h-1 accent-cyan-400 pointer-events-none [&::-webkit-slider-thumb]:pointer-events-auto [&::-moz-range-thumb]:pointer-events-auto appearance-none bg-transparent rounded" />
      </div>
      <span class="text-slate-500 w-10 font-mono">{readFilterLabel(log2Max)}</span>
    </label>
    <!-- Date filter -->
    {#if dateRange.min}
      <label class="text-xs text-slate-400 flex items-center gap-2">
        Run Date:
        <input type="date" bind:value={dateMin}
          min={dateRange.min} max={dateRange.max}
          class="bg-slate-800 border border-slate-600 rounded px-1.5 py-1 text-xs text-slate-200 w-[130px]" />
        <span class="text-slate-500">–</span>
        <input type="date" bind:value={dateMax}
          min={dateRange.min} max={dateRange.max}
          class="bg-slate-800 border border-slate-600 rounded px-1.5 py-1 text-xs text-slate-200 w-[130px]" />
      </label>
    {/if}
    {#if lassoIds}
      <button
        class="text-xs px-3 py-1.5 rounded bg-cyan-400/20 text-cyan-400 border border-cyan-400/40 hover:bg-cyan-400/30 transition-colors"
        onclick={addLassoToCart}
      >
        Add {lassoIds.length} to Cart
      </button>
    {/if}
  </div>

  <!-- Main layout: scatter + detail -->
  <div class="flex gap-6">
    <div class="flex-1 h-[600px] flex flex-col">
      {#if scatterData}
        <PlotlyScatter
          data={scatterData}
          {colorBy}
          {sizeBy}
          {sizeScale}
          mode="tsne"
          {colorMap}
          {coordExtents}
          onselect={handleSelect}
          onclick={handleClick}
        />
      {:else}
        <div class="flex items-center justify-center h-full text-slate-500">No sample data available</div>
      {/if}
    </div>

    <!-- Detail panel -->
    <div class="w-72 shrink-0 bg-slate-800 rounded-lg border border-slate-700 p-4 space-y-3 h-fit max-h-[600px] overflow-y-auto">
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
          {#if selectedDetail.start_time}
            <div class="text-slate-400">Run Date</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.start_time.slice(0, 10)}</div>
          {/if}
        </div>
        {#if selectedDetail.lat != null}
          <div class="border-t border-slate-700 pt-2">
            <div class="text-xs text-slate-400 mb-1">Location</div>
            <div class="grid grid-cols-2 gap-2 text-xs">
              <div class="text-slate-400">Lat</div><div class="text-slate-200 font-mono">{selectedDetail.lat}</div>
              <div class="text-slate-400">Lon</div><div class="text-slate-200 font-mono">{selectedDetail.lon}</div>
              {#if selectedDetail.station}<div class="text-slate-400">Station</div><div class="text-slate-200 font-mono">{selectedDetail.station}</div>{/if}
              {#if selectedDetail.depth_m != null}<div class="text-slate-400">Depth</div><div class="text-slate-200 font-mono">{selectedDetail.depth_m}m</div>{/if}
            </div>
          </div>
        {/if}
      {:else}
        <p class="text-slate-500 text-xs">Click a sample point or table row to view details.</p>
      {/if}
    </div>
  </div>

  <!-- Sample table -->
  {#if $samples}
    {@const filteredSamples = $samples.filter(passesFilters)}
    {@const tableRows = lassoIds
      ? filteredSamples.filter(s => lassoIds.includes(s.id))
      : $cartActive && $cartItems.size > 0
        ? filteredSamples.filter(s => $cartItems.has(s.id))
        : filteredSamples}
    <div>
      <h3 class="text-sm font-semibold text-slate-300 mb-2">
        {lassoIds ? `Selected (${lassoIds.length})` : $cartActive && $cartItems.size > 0 ? `Cart (${$cartItems.size})` : `Barcodes (${filteredSamples.length}/${$samples.length})`}
      </h3>
      <DataTable
        columns={tableColumns}
        rows={tableRows}
        onRowClick={(row) => { selectedSample.set(row.id); }}
        selectedId={$selectedSample}
        idKey="id"
        maxHeight="300px"
        actionLabel={(row) => $cartItems.has(row.id) ? 'In Cart' : '+ Cart'}
        actionFn={(row) => toggleCart(row.id)}
        actionStyle={(row) => $cartItems.has(row.id)
          ? 'text-[10px] px-2 py-0.5 rounded border bg-cyan-400/20 text-cyan-400 border-cyan-400/40 transition-colors'
          : 'text-[10px] px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-cyan-400 hover:border-cyan-400/40 transition-colors'}
      />
    </div>
  {/if}
</div>
