<script>
  import { onMount } from 'svelte';
  import ReglScatter from '../components/charts/ReglScatter.svelte';
  import StatCard from '../components/layout/StatCard.svelte';
  import { readExplorer, loadReadExplorer } from '../stores/data.js';
  import { cartItems, cartActive } from '../stores/cart.js';
  import { selectedRead } from '../stores/selection.js';

  let loading = $state(true);

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
  const taxGroup =    { values: ['kraken_phylum', 'kraken_class', 'kraken_order', 'kraken_family', 'kraken_domain'],
                        labels: ['Phylum', 'Class', 'Order', 'Family', 'Domain'] };
  const metricGroup = { values: ['gc', 'length'], labels: ['GC%', 'Length'] };

  const BW = { sample: '4.5rem', taxonomy: '4.5rem', metric: '4.5rem' };

  let colorMode = $state('sample'); // 'sample' | 'taxonomy' | 'metric'
  let taxRank = $state('kraken_phylum');
  let metric = $state('gc');

  let sizeScale = $state(0.6);

  let colorBy = $derived(
    colorMode === 'sample' ? 'sample' :
    colorMode === 'taxonomy' ? taxRank :
    metric
  );

  onMount(async () => {
    await loadReadExplorer();
    loading = false;
  });

  // Filter reads by cart
  let scatterData = $derived.by(() => {
    if (!$readExplorer?.reads) return null;
    let reads = $readExplorer.reads;

    if ($cartActive && $cartItems.size > 0) {
      reads = reads.filter(r => $cartItems.has(r.sample));
    }

    return { points: reads };
  });

  // Cart-based highlighting (dim non-carted)
  let searchMatchIds = $derived.by(() => {
    if (!$cartActive || $cartItems.size === 0 || !$readExplorer?.reads) return null;
    const ids = new Set();
    for (const r of $readExplorer.reads) {
      if ($cartItems.has(r.sample)) ids.add(r.id);
    }
    return ids;
  });

  let selectedDetail = $derived.by(() => {
    if (!$selectedRead || !$readExplorer?.reads) return null;
    return $readExplorer.reads.find(r => r.id === $selectedRead);
  });

  let stats = $derived.by(() => {
    const reads = scatterData?.points;
    if (!reads?.length) return null;
    const n = reads.length;
    const samples = new Set(reads.map(r => r.sample));
    const avgLen = reads.reduce((s, r) => s + (r.length || 0), 0) / n;
    const avgGc = reads.reduce((s, r) => s + (r.gc || 0), 0) / n;
    return { n, nSamples: samples.size, avgLen, avgGc };
  });
</script>

<div class="space-y-6">
  {#if loading}
    <div class="flex flex-col items-center justify-center py-20">
      <div class="w-10 h-10 border-4 border-slate-700 border-t-cyan-400 rounded-full animate-spin"></div>
      <p class="text-slate-400 text-sm mt-4">Loading read data (this may take a moment)...</p>
    </div>
  {:else}
    <!-- Stats -->
    {#if stats}
      <div class="grid grid-cols-2 md:grid-cols-4 gap-4">
        <StatCard label="Reads" value={stats.n.toLocaleString()} color="cyan" />
        <StatCard label="Samples" value={stats.nSamples} color="emerald" />
        <StatCard label="Avg Length" value={`${Math.round(stats.avgLen).toLocaleString()} bp`} color="amber" />
        <StatCard label="Avg GC" value={`${stats.avgGc.toFixed(1)}%`} color="purple" />
      </div>
    {/if}

    <!-- Controls -->
    <div class="flex items-center gap-2 text-xs">
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center
          {colorMode === 'sample' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        style="min-width: {BW.sample}"
        onclick={() => colorMode = 'sample'}
      >
        Sample
      </button>
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center
          {colorMode === 'taxonomy' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        style="min-width: {BW.taxonomy}"
        onclick={() => { if (colorMode === 'taxonomy') taxRank = cycle(taxGroup.values, taxRank); else colorMode = 'taxonomy'; }}
        title={`Click to cycle: ${taxGroup.labels.join(' → ')}`}
      >
        {colorMode === 'taxonomy' ? getLabel(taxGroup.values, taxGroup.labels, taxRank) : 'Taxonomy'} &#x25BE;
      </button>
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center
          {colorMode === 'metric' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        style="min-width: {BW.metric}"
        onclick={() => { if (colorMode === 'metric') metric = cycle(metricGroup.values, metric); else colorMode = 'metric'; }}
        title={`Click to cycle: ${metricGroup.labels.join(' → ')}`}
      >
        {colorMode === 'metric' ? getLabel(metricGroup.values, metricGroup.labels, metric) : 'Metric'} &#x25BE;
      </button>
      <span class="text-slate-600 mx-1">|</span>

      <div class="text-slate-400 flex items-center gap-1">
        Size
        <div class="single-range relative w-16 h-5 flex items-center">
          <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
          <input type="range" min="0.2" max="3" step="0.1" bind:value={sizeScale} />
        </div>
        <span class="text-slate-500 w-8 font-mono">{sizeScale.toFixed(1)}x</span>
      </div>

      {#if scatterData?.points}
        <span class="text-slate-500 ml-1">
          {scatterData.points.length.toLocaleString()} reads
        </span>
      {/if}
    </div>

    <!-- Scatter -->
    <div class="flex gap-6">
      <div class="flex-1 h-[700px] flex flex-col">
        {#if scatterData}
          <ReglScatter
            data={scatterData}
            {colorBy}
            sizeBy="fixed"
            {sizeScale}
            mode="tsne"
            {searchMatchIds}
            exportName={`danaseq_read_tsne_color-${colorBy}`}
          />
        {:else}
          <div class="flex items-center justify-center h-full text-slate-500">
            No read data available. Run preprocess to generate read_explorer.json.
          </div>
        {/if}
      </div>

      <!-- Read detail -->
      {#if selectedDetail}
        <div class="w-64 bg-slate-800 rounded-lg border border-slate-700 p-4 space-y-3 h-fit">
          <h3 class="text-sm font-semibold text-cyan-400 font-mono truncate" title={selectedDetail.id}>
            {selectedDetail.id.slice(0, 20)}...
          </h3>
          <div class="grid grid-cols-2 gap-2 text-xs">
            <div class="text-slate-400">Sample</div><div class="text-slate-200 font-mono">{selectedDetail.sample}</div>
            <div class="text-slate-400">Length</div><div class="text-slate-200 font-mono">{selectedDetail.length?.toLocaleString() ?? '-'} bp</div>
            <div class="text-slate-400">GC%</div><div class="text-slate-200 font-mono">{selectedDetail.gc?.toFixed(1) ?? '-'}</div>
            {#if selectedDetail.kraken_phylum}
              <div class="text-slate-400">Phylum</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.kraken_phylum}</div>
            {/if}
            {#if selectedDetail.kraken_class}
              <div class="text-slate-400">Class</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.kraken_class}</div>
            {/if}
            {#if selectedDetail.sketch_hit}
              <div class="text-slate-400">Sketch Hit</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.sketch_hit}</div>
            {/if}
          </div>
        </div>
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
