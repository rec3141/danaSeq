<script>
  import D3Sunburst from '../components/charts/D3Sunburst.svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, sampleTaxonomy, taxonomySunburst } from '../stores/data.js';
  import { cartItems, cartActive } from '../stores/cart.js';

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8',
                   '#94a3b8','#d4d4d8','#78716c'];

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
    // New format: {phylum: {...}, class: {...}}; fallback to flat dict
    return tax.phylum || tax;
  }

  // Stacked bar chart: per-sample phylum composition
  let stackedTraces = $derived.by(() => {
    if (!$sampleTaxonomy || !activeSamples.length) return [];
    // Collect all phyla
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
      y: sampleIds.map(id => getPhylumTax(id)[phylum] ?? 0),
      marker: { color: PALETTE[i % PALETTE.length] },
    }));
  });

  let stackedLayout = $derived({
    barmode: 'stack',
    height: 400,
    xaxis: { title: { text: 'Sample', font: { color: '#94a3b8' } }, tickangle: -45, tickfont: { size: 9 } },
    yaxis: { title: { text: 'Read Count', font: { color: '#94a3b8' } } },
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
    return Object.entries(counts)
      .map(([phylum, count]) => ({ phylum, count, pct: 0 }))
      .sort((a, b) => b.count - a.count)
      .map((row, _, arr) => {
        const total = arr.reduce((s, r) => s + r.count, 0);
        row.pct = total > 0 ? ((row.count / total) * 100).toFixed(1) : '0.0';
        return row;
      });
  });

  const taxColumns = [
    { key: 'phylum', label: 'Phylum' },
    { key: 'count', label: 'Reads', render: v => v?.toLocaleString() ?? '-' },
    { key: 'pct', label: '%', render: v => `${v}%` },
  ];

  // Build filtered sunburst from sample taxonomy when cart is active
  let sunburstData = $derived.by(() => {
    if ($cartActive && $cartItems.size > 0 && $sampleTaxonomy) {
      // Build from per-sample taxonomy (phylum level only for now)
      const counts = {};
      for (const s of activeSamples) {
        const tax = getPhylumTax(s.id);
        for (const [phylum, count] of Object.entries(tax)) {
          counts[phylum] = (counts[phylum] || 0) + count;
        }
      }
      return {
        name: 'Community',
        children: Object.entries(counts).map(([name, value]) => ({ name, value })),
      };
    }
    return $taxonomySunburst;
  });
</script>

<div class="space-y-6">
  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6">
    <!-- Sunburst -->
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">
        Taxonomy Overview
        {#if $cartActive && $cartItems.size > 0}
          <span class="text-xs text-cyan-400 ml-2">(cart filtered)</span>
        {/if}
      </h3>
      {#if sunburstData}
        <D3Sunburst data={sunburstData} />
      {:else}
        <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">No taxonomy data</div>
      {/if}
    </div>

    <!-- Top taxa table -->
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">Top Taxa</h3>
      <DataTable columns={taxColumns} rows={topTaxa} maxHeight="450px" idKey="phylum" />
    </div>
  </div>

  <!-- Stacked bar chart -->
  <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
    <h3 class="text-sm font-semibold text-slate-300 mb-2">Per-Sample Taxonomy</h3>
    {#if stackedTraces.length}
      <PlotlyChart traces={stackedTraces} layout={stackedLayout} />
    {:else}
      <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">No per-sample taxonomy data</div>
    {/if}
  </div>
</div>
