<script>
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, sampleFunction } from '../stores/data.js';
  import { cartItems, cartActive } from '../stores/cart.js';

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c'];

  let activeSamples = $derived.by(() => {
    if (!$samples) return [];
    if ($cartActive && $cartItems.size > 0) return $samples.filter(s => $cartItems.has(s.id));
    return $samples;
  });

  // Gene type bar chart (CDS, rRNA, tRNA per sample)
  let geneBarTraces = $derived.by(() => {
    if (!$sampleFunction || !activeSamples.length) return [];
    const ids = activeSamples.map(s => s.id);
    const types = ['cds_count', 'rrna_count', 'trna_count'];
    const labels = ['CDS', 'rRNA', 'tRNA'];
    return types.map((t, i) => ({
      type: 'bar', name: labels[i],
      x: ids,
      y: ids.map(id => $sampleFunction[id]?.[t] ?? 0),
      marker: { color: PALETTE[i] },
    }));
  });

  let geneBarLayout = $derived({
    barmode: 'stack', height: 400,
    xaxis: { title: { text: 'Sample', font: { color: '#94a3b8' } }, tickangle: -45, tickfont: { size: 9 } },
    yaxis: { title: { text: 'Feature Count', font: { color: '#94a3b8' } } },
    legend: { font: { size: 10 } },
    margin: { b: 80 },
  });

  // Hypothetical protein ratio per sample
  let hypotheticalTraces = $derived.by(() => {
    if (!$sampleFunction || !activeSamples.length) return [];
    const ids = activeSamples.map(s => s.id);
    return [{
      type: 'bar',
      x: ids,
      y: ids.map(id => {
        const f = $sampleFunction[id];
        return f?.hypothetical_pct != null ? f.hypothetical_pct : 0;
      }),
      marker: { color: '#fb923c' },
    }];
  });

  let hypotheticalLayout = $derived({
    height: 300,
    xaxis: { title: { text: 'Sample', font: { color: '#94a3b8' } }, tickangle: -45, tickfont: { size: 9 } },
    yaxis: { title: { text: '% Hypothetical', font: { color: '#94a3b8' } }, range: [0, 100] },
    margin: { b: 80 },
    showlegend: false,
  });

  // Top EC numbers across samples
  let topEc = $derived.by(() => {
    if (!$sampleFunction) return [];
    const counts = {};
    for (const s of activeSamples) {
      const ec = $sampleFunction[s.id]?.ec_counts;
      if (!ec) continue;
      for (const [ec_num, count] of Object.entries(ec)) {
        counts[ec_num] = (counts[ec_num] || 0) + count;
      }
    }
    return Object.entries(counts)
      .map(([ec, count]) => ({ ec, count }))
      .sort((a, b) => b.count - a.count)
      .slice(0, 30);
  });

  // Function summary table
  let funcTableRows = $derived.by(() => {
    if (!$sampleFunction) return [];
    return activeSamples.map(s => {
      const f = $sampleFunction[s.id] || {};
      return {
        id: s.id,
        cds_count: f.cds_count ?? 0,
        rrna_count: f.rrna_count ?? 0,
        trna_count: f.trna_count ?? 0,
        hypothetical_pct: f.hypothetical_pct ?? 0,
        n_ec: f.ec_counts ? Object.keys(f.ec_counts).length : 0,
        n_genes: f.n_genes ?? 0,
      };
    });
  });

  const funcColumns = [
    { key: 'id', label: 'Sample' },
    { key: 'n_genes', label: 'Total Genes', render: v => v.toLocaleString() },
    { key: 'cds_count', label: 'CDS', render: v => v.toLocaleString() },
    { key: 'rrna_count', label: 'rRNA' },
    { key: 'trna_count', label: 'tRNA' },
    { key: 'hypothetical_pct', label: '% Hypothetical', render: v => v.toFixed(1) },
    { key: 'n_ec', label: 'EC Numbers' },
  ];

  const ecColumns = [
    { key: 'ec', label: 'EC Number' },
    { key: 'count', label: 'Occurrences', render: v => v.toLocaleString() },
  ];
</script>

<div class="space-y-6">
  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6">
    <!-- Gene type breakdown -->
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">Gene Features per Sample</h3>
      {#if geneBarTraces.length}
        <PlotlyChart traces={geneBarTraces} layout={geneBarLayout} />
      {:else}
        <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">No functional data</div>
      {/if}
    </div>

    <!-- Hypothetical protein ratio -->
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">Hypothetical Protein %</h3>
      {#if hypotheticalTraces.length}
        <PlotlyChart traces={hypotheticalTraces} layout={hypotheticalLayout} />
      {:else}
        <div class="h-[300px] flex items-center justify-center text-slate-500 text-sm">No data</div>
      {/if}
    </div>
  </div>

  <!-- Function summary table -->
  <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
    <h3 class="text-sm font-semibold text-slate-300 mb-2">Functional Summary</h3>
    <DataTable columns={funcColumns} rows={funcTableRows} maxHeight="300px" exportFilename="functional_summary" />
  </div>

  <!-- Top EC numbers -->
  {#if topEc.length}
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">Top EC Numbers</h3>
      <DataTable columns={ecColumns} rows={topEc} maxHeight="300px" idKey="ec" exportFilename="top_ec_numbers" />
    </div>
  {/if}
</div>
