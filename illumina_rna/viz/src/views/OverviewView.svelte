<script>
  import PlotlyChart from '../components/PlotlyChart.svelte';
  import { overview, readFlow } from '../stores/data.js';

  // ── Mapping-rate matrix
  let matrixTraces = $derived.by(() => {
    const o = $overview;
    if (!o || !o.samples?.length || !o.references?.length) return [];
    return [{
      type: 'heatmap',
      z: o.mapping_rate_matrix,
      x: o.references, y: o.samples,
      colorscale: 'Viridis',
      zmin: 0, zmax: 1,
      colorbar: { title: 'mapping rate' },
      hovertemplate: '%{y} → %{x}<br>%{z:.1%}<extra></extra>',
    }];
  });
  let matrixLayout = {
    xaxis: { tickangle: -45, automargin: true },
    yaxis: { automargin: true },
    margin: { l: 140, r: 60, t: 20, b: 90 },
  };

  // ── Read-flow funnel
  let funnelTraces = $derived.by(() => {
    const rf = $readFlow;
    if (!rf || !rf.samples?.length) return [];
    const order = ['clumped', 'trimmed', 'filtered', 'nohuman', 'norrna'];
    const present = order.filter(s => rf.samples.some(sid => rf.stages[sid]?.[s] != null));
    if (!present.length) return [];
    const medians = present.map(stage => {
      const vals = rf.samples.map(sid => rf.stages[sid]?.[stage]).filter(v => v != null).sort((a, b) => a - b);
      return vals.length ? vals[Math.floor(vals.length / 2)] : 0;
    });
    return [{
      type: 'funnel', y: present, x: medians,
      textinfo: 'value+percent initial',
      marker: { color: '#22d3ee' },
    }];
  });
  let funnelLayout = { margin: { l: 100, r: 40, t: 20, b: 40 } };
</script>

<div class="space-y-6">
  <div class="grid grid-cols-3 gap-4">
    <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
      <div class="text-xs uppercase text-slate-500">Samples</div>
      <div class="text-3xl font-semibold text-cyan-400 mt-1">{$overview?.n_samples ?? '—'}</div>
    </div>
    <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
      <div class="text-xs uppercase text-slate-500">References</div>
      <div class="text-3xl font-semibold text-cyan-400 mt-1">{$overview?.n_references ?? '—'}</div>
    </div>
    <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
      <div class="text-xs uppercase text-slate-500">Pipeline</div>
      <div class="text-sm font-medium text-slate-200 mt-2">illumina_rna · BBmap + featureCounts</div>
    </div>
  </div>

  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
    <div class="text-sm font-semibold text-slate-200 mb-2">Mapping rate — samples × references</div>
    <PlotlyChart data={matrixTraces} layout={matrixLayout} height="520px" />
  </div>

  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
    <div class="text-sm font-semibold text-slate-200 mb-2">Read flow (median across samples)</div>
    <PlotlyChart data={funnelTraces} layout={funnelLayout} height="320px" />
  </div>
</div>
