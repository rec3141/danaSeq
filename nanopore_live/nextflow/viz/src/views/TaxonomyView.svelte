<script>
  import D3Sunburst from '../components/charts/D3Sunburst.svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, sampleTaxonomy, taxonomySunburst } from '../stores/data.js';
  import { cartItems, cartActive } from '../stores/cart.js';

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8',
                   '#94a3b8','#d4d4d8','#78716c'];

  let showPct = $state(false);

  // Filtered sample list
  let activeSamples = $derived.by(() => {
    if (!$samples) return [];
    if ($cartActive && $cartItems.size > 0) return $samples.filter(s => $cartItems.has(s.id));
    return $samples;
  });

  // Helper to get phylum-level taxonomy for a sample
  function getPhylumTax(sampleId) {
    const tax = $sampleTaxonomy?.[sampleId];
    if (!tax) return {};
    return tax.phylum || tax;
  }

  // Per-sample totals for percentage mode
  function sampleTotal(id) {
    const tax = getPhylumTax(id);
    return Object.values(tax).reduce((s, v) => s + v, 0) || 1;
  }

  // Stacked bar chart: per-sample phylum composition
  let stackedTraces = $derived.by(() => {
    if (!$sampleTaxonomy || !activeSamples.length) return [];
    const allPhyla = new Set();
    for (const s of activeSamples) {
      const tax = getPhylumTax(s.id);
      Object.keys(tax).forEach(k => allPhyla.add(k));
    }
    const phylaList = [...allPhyla].sort();
    const sampleIds = activeSamples.map(s => s.id);

    return phylaList.map((phylum, i) => ({
      type: 'bar',
      name: phylum.length > 20 ? phylum.slice(0, 18) + '..' : phylum,
      x: sampleIds,
      y: sampleIds.map(id => {
        const count = getPhylumTax(id)[phylum] ?? 0;
        return showPct ? (count / sampleTotal(id)) * 100 : count;
      }),
      marker: { color: PALETTE[i % PALETTE.length] },
    }));
  });

  let stackedLayout = $derived({
    barmode: 'stack',
    height: 400,
    xaxis: { title: { text: 'Sample', font: { color: '#94a3b8' } }, tickangle: -45, tickfont: { size: 9 } },
    yaxis: { title: { text: showPct ? '% of Classified' : 'Read Count', font: { color: '#94a3b8' } },
             range: showPct ? [0, 100] : undefined },
    legend: { font: { size: 9 }, bgcolor: 'rgba(15,23,42,0.7)' },
    margin: { b: 80 },
  });

  // Top taxa table (aggregated across active samples)
  let topTaxa = $derived.by(() => {
    if (!$sampleTaxonomy) return [];
    const counts = {};
    for (const s of activeSamples) {
      const tax = getPhylumTax(s.id);
      for (const [phylum, count] of Object.entries(tax)) {
        counts[phylum] = (counts[phylum] || 0) + count;
      }
    }
    const total = Object.values(counts).reduce((s, v) => s + v, 0);
    return Object.entries(counts)
      .map(([phylum, count]) => ({
        phylum,
        count,
        pct: total > 0 ? Math.round((count / total) * 1000) / 10 : 0,
      }))
      .sort((a, b) => b.count - a.count);
  });

  const taxColumns = [
    { key: 'phylum', label: 'Phylum' },
    { key: 'count', label: 'Reads', render: v => v?.toLocaleString() ?? '-' },
    { key: 'pct', label: '%', render: v => `${v.toFixed(1)}%` },
  ];
</script>

<div class="space-y-6">
  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6">
    <!-- Sunburst -->
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">
        Taxonomy Overview
        <span class="text-xs text-slate-500 font-normal ml-2">click to zoom, click center to zoom out</span>
      </h3>
      {#if $taxonomySunburst}
        <D3Sunburst data={$taxonomySunburst} colorDepth={3} />
      {:else}
        <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">No taxonomy data</div>
      {/if}
    </div>

    <!-- Top taxa table -->
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">
        Top Phyla
        {#if $cartActive && $cartItems.size > 0}
          <span class="text-xs text-cyan-400 ml-2">(cart filtered)</span>
        {/if}
      </h3>
      <DataTable columns={taxColumns} rows={topTaxa} maxHeight="450px" idKey="phylum" />
    </div>
  </div>

  <!-- Stacked bar chart -->
  <div class="bg-slate-800 rounded-lg border border-slate-700 p-4 overflow-hidden">
    <div class="flex items-center justify-between mb-2">
      <h3 class="text-sm font-semibold text-slate-300">
        Per-Sample Phylum Composition
        {#if $cartActive && $cartItems.size > 0}
          <span class="text-xs text-cyan-400 ml-2">(cart filtered)</span>
        {/if}
      </h3>
      <button
        class="text-xs px-2 py-1 rounded border transition-colors
          {showPct
            ? 'bg-cyan-400/20 text-cyan-400 border-cyan-400/40'
            : 'text-slate-400 border-slate-600 hover:text-slate-200 hover:border-slate-500'}"
        onclick={() => showPct = !showPct}
      >
        {showPct ? '%' : '#'}
      </button>
    </div>
    {#if stackedTraces.length}
      <PlotlyChart traces={stackedTraces} layout={stackedLayout} />
    {:else}
      <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">No per-sample taxonomy data</div>
    {/if}
  </div>
</div>
