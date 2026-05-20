<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import { overview, readFlow } from '../stores/data.js';

  let matrixEl, funnelEl;

  // ── Mapping-rate matrix (samples × references)
  function drawMatrix(o) {
    if (!o || !o.samples?.length || !o.references?.length) return;
    Plotly.purge(matrixEl);
    Plotly.newPlot(matrixEl, [{
      type: 'heatmap',
      z: o.mapping_rate_matrix,
      x: o.references,
      y: o.samples,
      colorscale: 'Viridis',
      zmin: 0, zmax: 1,
      colorbar: { title: 'mapping rate' },
      hovertemplate: '%{y} → %{x}<br>%{z:.1%}<extra></extra>',
    }], {
      paper_bgcolor: '#0f172a', plot_bgcolor: '#0f172a',
      font: { color: '#cbd5e1' },
      margin: { l: 140, r: 60, t: 30, b: 90 },
      xaxis: { tickangle: -45, automargin: true },
      yaxis: { automargin: true },
    }, { responsive: true, displaylogo: false });
  }

  // ── Read-flow funnel (median per stage across samples)
  function drawFunnel(rf) {
    if (!rf || !rf.samples?.length) return;
    const stageOrder = ['clumped', 'trimmed', 'filtered', 'nohuman', 'norrna'];
    const present = stageOrder.filter(s =>
      rf.samples.some(sid => rf.stages[sid]?.[s] != null)
    );
    if (!present.length) return;
    const medians = present.map(stage => {
      const vals = rf.samples
        .map(sid => rf.stages[sid]?.[stage])
        .filter(v => v != null)
        .sort((a, b) => a - b);
      return vals.length ? vals[Math.floor(vals.length / 2)] : 0;
    });
    Plotly.purge(funnelEl);
    Plotly.newPlot(funnelEl, [{
      type: 'funnel',
      y: present,
      x: medians,
      textinfo: 'value+percent initial',
      marker: { color: '#22d3ee' },
    }], {
      paper_bgcolor: '#0f172a', plot_bgcolor: '#0f172a',
      font: { color: '#cbd5e1' },
      margin: { l: 100, r: 40, t: 30, b: 40 },
    }, { responsive: true, displaylogo: false });
  }

  $effect(() => { drawMatrix($overview); });
  $effect(() => { drawFunnel($readFlow); });
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
    <div bind:this={matrixEl} style="width:100%;height:520px"></div>
  </div>

  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
    <div class="text-sm font-semibold text-slate-200 mb-2">Read flow (median across samples)</div>
    <div bind:this={funnelEl} style="width:100%;height:320px"></div>
  </div>
</div>
