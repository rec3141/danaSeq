<script>
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { eukaryotic } from '../stores/data.js';

  let eukData = $derived($eukaryotic);

  // Tiara domain donut (size-weighted)
  let tiaraDonutData = $derived.by(() => {
    if (!eukData?.tiara_size_weighted) return [];
    const entries = Object.entries(eukData.tiara_size_weighted).sort((a, b) => b[1] - a[1]);
    const colors = {
      bacteria: '#22d3ee', archaea: '#fbbf24', eukaryota: '#34d399',
      organelle: '#a78bfa', unknown: '#64748b',
    };
    return [{
      type: 'pie',
      labels: entries.map(([k]) => k),
      values: entries.map(([, v]) => v),
      marker: { colors: entries.map(([k]) => colors[k] || '#64748b') },
      hole: 0.5,
      textinfo: 'label+percent',
      textfont: { color: '#e2e8f0', size: 11 },
      hovertemplate: '%{label}: %{value:.0f} bp (%{percent})<extra></extra>',
    }];
  });

  // MarFERReT phylum bar
  let marferretBarData = $derived.by(() => {
    if (!eukData?.marferret_taxonomy) return [];
    const entries = Object.entries(eukData.marferret_taxonomy)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 15);
    const palette = ['#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa', '#fb923c',
                     '#2dd4bf', '#818cf8', '#f472b6', '#4ade80', '#e879f9', '#38bdf8'];
    return [{
      type: 'bar',
      x: entries.map(([k]) => k),
      y: entries.map(([, v]) => v),
      marker: { color: entries.map((_, i) => palette[i % palette.length]), opacity: 0.85 },
      hovertemplate: '%{x}: %{y} contigs<extra></extra>',
    }];
  });

  // Whokaryote counts as simple bar
  let whokaryoteBarData = $derived.by(() => {
    if (!eukData?.whokaryote_counts) return [];
    const entries = Object.entries(eukData.whokaryote_counts).sort((a, b) => b[1] - a[1]);
    return [{
      type: 'bar',
      x: entries.map(([k]) => k),
      y: entries.map(([, v]) => v),
      marker: { color: '#a78bfa', opacity: 0.8 },
      hovertemplate: '%{x}: %{y}<extra></extra>',
    }];
  });

  // MarFERReT contig table
  let tableColumns = [
    { key: 'contig', label: 'Contig' },
    { key: 'n_proteins', label: 'Proteins' },
    { key: 'n_classified', label: 'Classified' },
    { key: 'taxonomy', label: 'Top Taxonomy' },
    { key: 'pfam', label: 'Pfam Domains' },
  ];
</script>

{#if eukData}
  <div class="grid grid-cols-1 lg:grid-cols-3 gap-6 mb-6">
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart
        data={tiaraDonutData}
        layout={{
          title: { text: 'Tiara Classification (size-weighted)', font: { size: 13 } },
          showlegend: true,
          legend: { font: { size: 10 } },
        }}
      />
    </div>
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart
        data={whokaryoteBarData}
        layout={{
          title: { text: 'Whokaryote Classification', font: { size: 13 } },
          yaxis: { title: 'Count' },
        }}
      />
    </div>
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart
        data={marferretBarData}
        layout={{
          title: { text: 'MarFERReT Eukaryotic Phyla', font: { size: 13 } },
          xaxis: { tickangle: -45 },
          yaxis: { title: 'Contigs' },
        }}
      />
    </div>
  </div>

  {#if eukData.marferret_contigs?.length}
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <h3 class="text-sm font-medium text-slate-400 mb-2">
        Eukaryotic Contigs (MarFERReT, top 500)
      </h3>
      <DataTable
        columns={tableColumns}
        rows={eukData.marferret_contigs}
        idKey="contig"
        maxHeight="400px"
      />
    </div>
  {/if}
{/if}
