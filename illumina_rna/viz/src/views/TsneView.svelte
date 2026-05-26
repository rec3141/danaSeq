<script>
  import PlotlyChart from '../components/PlotlyChart.svelte';
  import { tsne } from '../stores/data.js';

  let selectedRef = $state('all');
  let colorBy = $state('mag');
  let sizeByExpr = $state(true);

  let refKeys = $derived(Object.keys($tsne || {}).sort());

  const COLORS = [
    '#22d3ee', '#34d399', '#a78bfa', '#f472b6', '#fbbf24',
    '#fb7185', '#60a5fa', '#facc15', '#84cc16', '#06b6d4',
    '#c084fc', '#f87171', '#2dd4bf', '#a3e635', '#fb923c',
    '#818cf8', '#ec4899',
  ];
  function palette(keys) {
    const m = {};
    keys.forEach((k, i) => { m[k] = COLORS[i % COLORS.length]; });
    return m;
  }

  let traces = $derived.by(() => {
    if (!$tsne) return [];
    const pal = palette(refKeys);
    const targets = selectedRef === 'all' ? refKeys : [selectedRef];
    const out = [];
    for (const ref of targets) {
      const d = $tsne[ref];
      if (!d) continue;
      const totals = d.totals || [];
      const maxT = Math.max(1, ...totals);
      let marker;
      if (colorBy === 'mag') {
        marker = {
          color: pal[ref], opacity: 0.7, line: { width: 0 },
          size: sizeByExpr
            ? totals.map(x => 4 + 14 * Math.log2(x + 1) / Math.log2(maxT + 1))
            : 7,
        };
      } else {
        const logT = totals.map(x => Math.log2(x + 1));
        marker = {
          color: logT, colorscale: 'Viridis',
          showscale: targets.length === 1,
          colorbar: targets.length === 1 ? { title: 'log₂(Σ+1)' } : undefined,
          opacity: 0.8, line: { width: 0 },
          size: sizeByExpr
            ? totals.map(x => 4 + 14 * Math.log2(x + 1) / Math.log2(maxT + 1))
            : 7,
        };
      }
      out.push({
        x: d.coords.map(c => c[0]), y: d.coords.map(c => c[1]),
        mode: 'markers', type: 'scattergl', name: ref, marker,
        text: d.gene_ids.map((g, i) =>
          `<b>${g}</b><br>${d.products[i] || '—'}<br>Σ=${totals[i]?.toLocaleString() ?? 0}`),
        hovertemplate: '%{text}<extra></extra>',
      });
    }
    return out;
  });

  let layout = $derived({
    xaxis: { title: 't-SNE 1' }, yaxis: { title: 't-SNE 2' },
    showlegend: (selectedRef === 'all'),
    legend: { font: { size: 10 } },
  });
</script>

<div class="space-y-4">
  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4 flex flex-wrap items-center gap-6">
    <label class="text-sm text-slate-400">
      MAG:
      <select bind:value={selectedRef} class="ml-2 bg-slate-800 border border-slate-700 rounded px-2 py-1 text-slate-200 text-sm">
        <option value="all">all MAGs combined</option>
        {#each refKeys as k}<option value={k}>{k}</option>{/each}
      </select>
    </label>
    <label class="text-sm text-slate-400">
      Color by:
      <select bind:value={colorBy} class="ml-2 bg-slate-800 border border-slate-700 rounded px-2 py-1 text-slate-200 text-sm">
        <option value="mag">MAG</option>
        <option value="expression">Total expression</option>
      </select>
    </label>
    <label class="text-sm text-slate-400 flex items-center gap-2">
      <input type="checkbox" bind:checked={sizeByExpr}/> Size by expression
    </label>
    <div class="text-xs text-slate-500">{refKeys.length} MAGs · top-500-variance genes each · scroll to zoom, drag to pan</div>
  </div>

  <div class="bg-slate-900 border border-slate-800 rounded-lg p-2">
    {#if refKeys.length === 0}
      <div class="text-slate-500 text-sm py-12 text-center">No gene t-SNE data. Run augment_viz.py with sklearn installed.</div>
    {:else}
      <PlotlyChart data={traces} {layout} height="680px" />
    {/if}
  </div>
</div>
