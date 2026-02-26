<script>
  import StatCard from '../components/layout/StatCard.svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import PipelineDAG from '../components/charts/PipelineDAG.svelte';
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

  // Format a log10 value as human-readable
  function logLabel(logVal) {
    const v = Math.pow(10, logVal);
    if (v >= 1e6) return (v / 1e6).toFixed(0) + 'M';
    if (v >= 1e3) return (v / 1e3).toFixed(0) + 'K';
    if (v >= 1) return v.toFixed(0);
    return v.toFixed(2);
  }

  function logTicksForRange(lo, hi, suffix) {
    const vals = [];
    const text = [];
    for (let e = Math.ceil(lo); e <= Math.floor(hi); e++) {
      vals.push(e);
      const v = Math.pow(10, e);
      let label;
      if (v >= 1e6) label = (v / 1e6).toFixed(0) + 'M';
      else if (v >= 1e3) label = (v / 1e3).toFixed(0) + 'K';
      else label = v.toFixed(v < 1 ? 2 : 0);
      text.push(label + (suffix || ''));
    }
    return { vals, text };
  }

  // Length histogram: use log10(center) as x, label axis with real bp values
  let histogramData = $derived.by(() => {
    if (!lenData) return [];
    const edges = lenData.bin_edges;
    const counts = lenData.counts;
    const x = [];
    const hoverText = [];
    for (let i = 0; i < counts.length; i++) {
      const center = (edges[i] + edges[i + 1]) / 2;
      x.push(Math.log10(center));
      hoverText.push(`${logLabel(Math.log10(center))} bp: ${counts[i]} contigs`);
    }
    return [{
      type: 'bar',
      x, y: counts,
      text: hoverText,
      hoverinfo: 'text',
      textposition: 'none',
      marker: { color: '#22d3ee', opacity: 0.8 },
    }];
  });

  let histogramLayout = $derived.by(() => {
    if (!lenData) return {};
    const edges = lenData.bin_edges;
    const lo = Math.log10(edges[0]);
    const hi = Math.log10(edges[edges.length - 1]);
    const ticks = logTicksForRange(lo, hi, ' bp');
    return {
      title: { text: 'Contig Length Distribution', font: { size: 14 } },
      xaxis: { title: 'Contig length', tickvals: ticks.vals, ticktext: ticks.text, fixedrange: true },
      yaxis: { title: 'Count', fixedrange: true },
      shapes: [{
        type: 'line',
        x0: Math.log10(lenData.n50), x1: Math.log10(lenData.n50),
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

  // Coverage histogram: log10-space bins
  let covHistData = $derived.by(() => {
    if (!lenData?.cov_log_edges) return [];
    const edges = lenData.cov_log_edges;
    const counts = lenData.cov_counts;
    const x = [];
    const hoverText = [];
    for (let i = 0; i < counts.length; i++) {
      const center = (edges[i] + edges[i + 1]) / 2;
      x.push(center);
      hoverText.push(`${logLabel(center)}x: ${counts[i]} contigs`);
    }
    return [{
      type: 'bar',
      x, y: counts,
      text: hoverText,
      hoverinfo: 'text',
      textposition: 'none',
      marker: { color: '#34d399', opacity: 0.8 },
    }];
  });

  let covHistLayout = $derived.by(() => {
    if (!lenData?.cov_log_edges) return {};
    const edges = lenData.cov_log_edges;
    const lo = edges[0];
    const hi = edges[edges.length - 1];
    const ticks = logTicksForRange(lo, hi, 'x');
    return {
      title: { text: 'Coverage Distribution', font: { size: 14 } },
      xaxis: { title: 'Coverage (depth)', tickvals: ticks.vals, ticktext: ticks.text, fixedrange: true },
      yaxis: { title: 'Count', fixedrange: true },
    };
  });

  // Length-coverage scatter, colored by GC%
  let scatterData = $derived.by(() => {
    if (!lenData?.scatter_log_length) return [];
    const hasGc = lenData.scatter_gc && lenData.scatter_gc.length === lenData.scatter_log_length.length;
    return [{
      type: 'scattergl',
      mode: 'markers',
      x: lenData.scatter_log_length,
      y: lenData.scatter_log_depth,
      marker: hasGc ? {
        color: lenData.scatter_gc,
        colorscale: 'Viridis',
        cmin: 20,
        cmax: 80,
        opacity: 0.5,
        size: 3,
        colorbar: { title: 'GC%', thickness: 12, len: 0.5, tickfont: { size: 9 } },
      } : {
        color: '#a78bfa',
        opacity: 0.3,
        size: 3,
      },
      hovertemplate: hasGc
        ? 'Length: %{customdata[0]}<br>Depth: %{customdata[1]}x<br>GC: %{customdata[2]}%<extra></extra>'
        : 'Length: %{customdata[0]}<br>Depth: %{customdata[1]}x<extra></extra>',
      customdata: lenData.scatter_log_length.map((lx, i) => {
        const row = [logLabel(lx) + ' bp', logLabel(lenData.scatter_log_depth[i])];
        if (hasGc) row.push(lenData.scatter_gc[i]);
        return row;
      }),
    }];
  });

  let scatterLayout = $derived.by(() => {
    if (!lenData?.scatter_log_length) return {};
    const lx = lenData.scatter_log_length;
    const ly = lenData.scatter_log_depth;
    let xlo = Infinity, xhi = -Infinity, ylo = Infinity, yhi = -Infinity;
    for (const v of lx) { if (v < xlo) xlo = v; if (v > xhi) xhi = v; }
    for (const v of ly) { if (v < ylo) ylo = v; if (v > yhi) yhi = v; }
    const xpad = (xhi - xlo) * 0.05 || 0.5;
    const ypad = (yhi - ylo) * 0.05 || 0.5;
    const xticks = logTicksForRange(xlo, xhi, ' bp');
    const yticks = logTicksForRange(ylo, yhi, 'x');
    return {
      title: { text: 'Length vs Coverage', font: { size: 14 } },
      xaxis: { title: 'Contig length', tickvals: xticks.vals, ticktext: xticks.text,
               range: [xlo - xpad, xhi + xpad], autorange: false },
      yaxis: { title: 'Coverage (depth)', tickvals: yticks.vals, ticktext: yticks.text,
               range: [ylo - ypad, yhi + ypad], autorange: false },
    };
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

  <div class="grid grid-cols-1 lg:grid-cols-3 gap-6 mb-6">
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart data={histogramData} layout={histogramLayout} />
    </div>
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      {#if scatterData.length}
        <PlotlyChart data={scatterData} layout={scatterLayout} />
      {:else}
        <p class="text-slate-500 text-sm py-8 text-center">No coverage data available</p>
      {/if}
    </div>
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      {#if covHistData.length}
        <PlotlyChart data={covHistData} layout={covHistLayout} />
      {:else}
        <p class="text-slate-500 text-sm py-8 text-center">No coverage data available</p>
      {/if}
    </div>
  </div>

  {#if overviewData?.processes}
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <div class="flex items-center justify-between mb-3">
        <h3 class="text-sm font-medium text-slate-400">
          Pipeline Status ({overviewData.pipeline_completed}/{overviewData.pipeline_total} processes)
        </h3>
        <div class="flex gap-3 text-xs">
          <a href="pipeline_info/report.html" target="_blank" rel="noopener noreferrer"
             class="text-slate-500 hover:text-cyan-400 transition-colors">Nextflow Report</a>
          <a href="pipeline_info/timeline.html" target="_blank" rel="noopener noreferrer"
             class="text-slate-500 hover:text-cyan-400 transition-colors">Timeline</a>
        </div>
      </div>
      <PipelineDAG processes={overviewData.processes} />
    </div>
  {/if}
{/if}
