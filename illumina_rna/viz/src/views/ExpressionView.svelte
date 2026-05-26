<script>
  import PlotlyChart from '../components/PlotlyChart.svelte';
  import { expression } from '../stores/data.js';

  let selectedRef = $state(null);
  let topN = $state(50);

  let refKeys = $derived(Object.keys($expression || {}));

  $effect(() => { if (!selectedRef && refKeys.length) selectedRef = refKeys[0]; });

  function logTransform(matrix) {
    return matrix.map(row => row.map(x => Math.log2((x || 0) + 1)));
  }
  function topByVariance(genes, matrix, n) {
    if (!matrix.length) return { genes: [], matrix: [] };
    const stats = matrix.map((row, i) => {
      const m = row.reduce((a, b) => a + b, 0) / row.length;
      const v = row.reduce((a, b) => a + (b - m) ** 2, 0) / row.length;
      return { idx: i, v };
    }).sort((a, b) => b.v - a.v).slice(0, n);
    return {
      genes:  stats.map(s => genes[s.idx]),
      matrix: stats.map(s => matrix[s.idx]),
    };
  }

  let traces = $derived.by(() => {
    if (!selectedRef || !$expression?.[selectedRef]) return [];
    const { genes, samples, matrix } = $expression[selectedRef];
    if (!matrix.length) return [];
    const lg = logTransform(matrix);
    const { genes: top_genes, matrix: top_mat } = topByVariance(genes, lg, topN);
    return [{
      type: 'heatmap', z: top_mat, x: samples, y: top_genes,
      colorscale: 'Cividis',
      colorbar: { title: 'log₂(count+1)' },
      hovertemplate: '%{y} · %{x}<br>%{z:.2f}<extra></extra>',
    }];
  });

  let layout = {
    margin: { l: 200, r: 40, t: 20, b: 100 },
    xaxis: { tickangle: -45, automargin: true },
    yaxis: { automargin: true },
  };
</script>

<div class="space-y-4">
  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
    <div class="flex items-center gap-4">
      <label class="text-sm text-slate-400">
        Reference:
        <select bind:value={selectedRef} class="ml-2 bg-slate-800 border border-slate-700 rounded px-2 py-1 text-slate-200 text-sm">
          {#each refKeys as ref}<option value={ref}>{ref}</option>{/each}
        </select>
      </label>
      <label class="text-sm text-slate-400">
        Top genes by variance:
        <input type="number" min="10" max="500" step="10" bind:value={topN}
               class="ml-2 bg-slate-800 border border-slate-700 rounded px-2 py-1 text-slate-200 text-sm w-20"/>
      </label>
      <div class="text-xs text-slate-500 ml-auto">scroll to zoom, drag to pan</div>
    </div>
  </div>

  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
    {#if refKeys.length === 0}
      <div class="text-slate-500 text-sm py-8 text-center">
        No expression data available. featureCounts requires a <code>&lt;ref&gt;.gff</code> alongside each reference FASTA.
      </div>
    {:else}
      <PlotlyChart data={traces} {layout} height="640px" />
    {/if}
  </div>
</div>
