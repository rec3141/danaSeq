<script>
  import { onMount } from 'svelte';
  import ReglScatter from '../components/charts/ReglScatter.svelte';
  import StatCard from '../components/layout/StatCard.svelte';
  import { readExplorer, loadReadExplorer } from '../stores/data.js';
  import { cartItems, cartActive } from '../stores/cart.js';
  import { selectedRead } from '../stores/selection.js';

  let colorBy = $state('sample');
  let loading = $state(true);

  const colorOptions = [
    { value: 'sample', label: 'Sample' },
    { value: 'kraken_domain', label: 'Domain' },
    { value: 'kraken_phylum', label: 'Phylum' },
    { value: 'kraken_class', label: 'Class' },
    { value: 'kraken_order', label: 'Order' },
    { value: 'kraken_family', label: 'Family' },
    { value: 'gc', label: 'GC%' },
    { value: 'length', label: 'Read Length' },
  ];

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
    <div class="flex items-center gap-3 text-xs">
      <select bind:value={colorBy}
        class="px-2 py-1 rounded-md border border-slate-600 bg-slate-800 text-slate-300 text-xs focus:border-cyan-400 focus:outline-none cursor-pointer">
        {#each colorOptions as opt}
          <option value={opt.value}>{opt.label}</option>
        {/each}
      </select>
      {#if scatterData?.points}
        <span class="text-xs text-slate-500">
          {scatterData.points.length.toLocaleString()} reads displayed
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
            sizeScale={0.6}
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
