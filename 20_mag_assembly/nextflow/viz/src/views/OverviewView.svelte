<script>
  import StatCard from '../components/layout/StatCard.svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import { overview, contigLengths, mags } from '../stores/data.js';

  let overviewData = $derived($overview);
  let lenData = $derived($contigLengths);
  let magsData = $derived($mags);

  function formatBp(n) {
    if (n >= 1e9) return (n / 1e9).toFixed(1) + ' Gbp';
    if (n >= 1e6) return (n / 1e6).toFixed(1) + ' Mbp';
    if (n >= 1e3) return (n / 1e3).toFixed(1) + ' Kbp';
    return n + ' bp';
  }

  let histogramData = $derived.by(() => {
    if (!lenData) return [];
    const edges = lenData.bin_edges;
    const counts = lenData.counts;
    const x = [];
    for (let i = 0; i < counts.length; i++) {
      x.push((edges[i] + edges[i + 1]) / 2);
    }
    return [{
      type: 'bar',
      x, y: counts,
      marker: { color: '#22d3ee', opacity: 0.8 },
      hovertemplate: '%{x:.0f} bp: %{y} contigs<extra></extra>',
    }];
  });

  let histogramLayout = $derived.by(() => {
    if (!lenData) return {};
    return {
      title: { text: 'Contig Length Distribution', font: { size: 14 } },
      xaxis: { title: 'Contig length (bp)', type: 'log' },
      yaxis: { title: 'Count' },
      shapes: [{
        type: 'line',
        x0: lenData.n50, x1: lenData.n50,
        y0: 0, y1: 1, yref: 'paper',
        line: { color: '#fbbf24', width: 2, dash: 'dash' },
      }],
      annotations: [{
        x: Math.log10(lenData.n50), y: 1, yref: 'paper',
        text: `N50: ${formatBp(lenData.n50)}`,
        showarrow: false,
        font: { color: '#fbbf24', size: 11 },
        yanchor: 'bottom',
      }],
    };
  });

  let donutData = $derived.by(() => {
    if (!overviewData?.binner_counts) return [];
    const entries = Object.entries(overviewData.binner_counts);
    const colors = {
      semibin: '#22d3ee', metabat: '#34d399', maxbin: '#fbbf24',
      lorbin: '#a78bfa', comebin: '#fb923c', other: '#64748b',
    };
    return [{
      type: 'pie',
      labels: entries.map(([k]) => k),
      values: entries.map(([, v]) => v),
      marker: { colors: entries.map(([k]) => colors[k] || '#64748b') },
      hole: 0.5,
      textinfo: 'label+value',
      textfont: { color: '#e2e8f0', size: 12 },
      hovertemplate: '%{label}: %{value} MAGs (%{percent})<extra></extra>',
    }];
  });
</script>

{#if overviewData}
  <div class="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-6 gap-4 mb-6">
    <StatCard label="Assembly Size" value={formatBp(overviewData.assembly_size)} color="cyan" />
    <StatCard label="N50" value={formatBp(overviewData.n50)} color="cyan" />
    <StatCard label="Contigs" value={overviewData.n_contigs.toLocaleString()} color="slate" />
    <StatCard label="Total MAGs" value={overviewData.n_mags} sub={`${overviewData.hq} HQ / ${overviewData.mq} MQ / ${overviewData.lq} LQ`} color="emerald" />
    <StatCard label="Viral Contigs" value={overviewData.n_virus.toLocaleString()} color="rose" />
    <StatCard label="Plasmid Contigs" value={overviewData.n_plasmid.toLocaleString()} color="purple" />
  </div>

  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6">
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart data={histogramData} layout={histogramLayout} />
    </div>
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart
        data={donutData}
        layout={{ title: { text: 'Binner Origin (DAS Tool)', font: { size: 14 } }, showlegend: true }}
      />
    </div>
  </div>
{/if}
