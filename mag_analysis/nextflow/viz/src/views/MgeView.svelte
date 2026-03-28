<script>
  import StatCard from '../components/layout/StatCard.svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { mgeSummary, mgePerBin, mags } from '../stores/data.js';
  import { selectedMag } from '../stores/selection.js';

  let summary = $derived($mgeSummary);
  let perBin = $derived($mgePerBin);
  let magsData = $derived($mags);
  let selected = $derived($selectedMag);

  let defenseBarData = $derived.by(() => {
    if (!summary?.defense_types) return [];
    const entries = Object.entries(summary.defense_types).sort((a, b) => b[1] - a[1]);
    return [{
      type: 'bar',
      x: entries.map(([k]) => k),
      y: entries.map(([, v]) => v),
      marker: { color: '#22d3ee', opacity: 0.8 },
      hovertemplate: '%{x}: %{y}<extra></extra>',
    }];
  });

  let viralDonutData = $derived.by(() => {
    if (!summary?.viral_taxonomy) return [];
    const entries = Object.entries(summary.viral_taxonomy).sort((a, b) => b[1] - a[1]);
    const palette = ['#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa', '#fb923c', '#818cf8', '#f472b6'];
    return [{
      type: 'pie',
      labels: entries.map(([k]) => k),
      values: entries.map(([, v]) => v),
      marker: { colors: palette },
      hole: 0.5,
      textinfo: 'label+percent',
      textfont: { color: '#e2e8f0', size: 10 },
      hovertemplate: '%{label}: %{value}<extra></extra>',
    }];
  });

  let tableColumns = [
    { key: 'name', label: 'MAG' },
    { key: 'n_virus', label: 'Virus' },
    { key: 'n_plasmid', label: 'Plasmid' },
    { key: 'n_defense', label: 'Defense' },
    { key: 'n_integron', label: 'Integron' },
    { key: 'total', label: 'Total' },
  ];

  let tableRows = $derived.by(() => {
    if (!magsData) return [];
    return magsData.map(m => ({
      name: m.name,
      n_virus: m.n_virus,
      n_plasmid: m.n_plasmid,
      n_defense: m.n_defense,
      n_integron: m.n_integron,
      total: m.n_virus + m.n_plasmid + m.n_defense + m.n_integron,
    })).sort((a, b) => b.total - a.total);
  });

  let expandedMag = $state(null);

  let expandedDetail = $derived.by(() => {
    if (!expandedMag || !perBin) return null;
    return perBin[expandedMag] || null;
  });
</script>

{#if summary}
  <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
    <StatCard label="Viral Contigs" value={summary.n_virus.toLocaleString()} color="rose" />
    <StatCard label="Plasmid Contigs" value={summary.n_plasmid.toLocaleString()} color="purple" />
    <StatCard label="Defense Systems" value={summary.n_defense} color="cyan" />
    <StatCard label="Integron Elements" value={summary.n_integron_elements.toLocaleString()} color="amber" />
  </div>

  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-6">
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart
        data={defenseBarData}
        layout={{
          title: { text: 'Defense System Types', font: { size: 14 } },
          xaxis: { tickangle: -45 },
          yaxis: { title: 'Count' },
        }}
      />
    </div>
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart
        data={viralDonutData}
        layout={{
          title: { text: 'Viral Taxonomy (geNomad)', font: { size: 14 } },
          showlegend: true,
          legend: { font: { size: 9 } },
        }}
      />
    </div>
  </div>

  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
    <h3 class="text-sm font-medium text-slate-400 mb-2">Per-MAG MGE Summary</h3>
    <DataTable
      columns={tableColumns}
      rows={tableRows}
      onRowClick={(row) => { expandedMag = expandedMag === row.name ? null : row.name; }}
      selectedId={expandedMag}
      maxHeight="400px"
    />
  </div>

  {#if expandedDetail}
    <div class="mt-4 bg-slate-800 rounded-lg p-4 border border-cyan-400/30">
      <h4 class="text-cyan-400 text-sm font-semibold mb-2">{expandedMag} — MGE Details</h4>
      <div class="grid grid-cols-1 md:grid-cols-2 gap-4 text-xs font-mono">
        {#if expandedDetail.viruses.length}
          <div>
            <h5 class="text-rose-400 mb-1">Viruses ({expandedDetail.viruses.length})</h5>
            {#each expandedDetail.viruses.slice(0, 10) as v}
              <div class="text-slate-300">{v.contig} — {v.length} bp, score {v.score}</div>
            {/each}
            {#if expandedDetail.viruses.length > 10}
              <div class="text-slate-500">...and {expandedDetail.viruses.length - 10} more</div>
            {/if}
          </div>
        {/if}
        {#if expandedDetail.plasmids.length}
          <div>
            <h5 class="text-purple-400 mb-1">Plasmids ({expandedDetail.plasmids.length})</h5>
            {#each expandedDetail.plasmids.slice(0, 10) as p}
              <div class="text-slate-300">{p.contig} — {p.length} bp, score {p.score}</div>
            {/each}
          </div>
        {/if}
        {#if expandedDetail.defense.length}
          <div>
            <h5 class="text-cyan-400 mb-1">Defense ({expandedDetail.defense.length})</h5>
            {#each expandedDetail.defense as d}
              <div class="text-slate-300">{d.type} ({d.subtype}) — {d.contig}</div>
            {/each}
          </div>
        {/if}
      </div>
    </div>
  {/if}
{/if}
