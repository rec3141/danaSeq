<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import { expression } from '../stores/data.js';

  let heatmapEl;
  let selectedRef = $state(null);
  let topN = $state(50);

  let refKeys = $derived(Object.keys($expression || {}));

  $effect(() => {
    if (!selectedRef && refKeys.length) selectedRef = refKeys[0];
  });

  // log2(count + 1) — keeps low-counts visible without taking actual log of zeros.
  function logTransform(matrix) {
    return matrix.map(row => row.map(x => Math.log2((x || 0) + 1)));
  }

  // Rank genes by row variance, take top N
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

  function drawHeatmap(ref) {
    if (!heatmapEl || !ref || !$expression?.[ref]) return;
    const { genes, samples, matrix } = $expression[ref];
    if (!matrix.length) {
      Plotly.purge(heatmapEl);
      return;
    }
    const lg = logTransform(matrix);
    const { genes: top_genes, matrix: top_mat } = topByVariance(genes, lg, topN);
    Plotly.purge(heatmapEl);
    Plotly.newPlot(heatmapEl, [{
      type: 'heatmap',
      z: top_mat,
      x: samples,
      y: top_genes,
      colorscale: 'Cividis',
      colorbar: { title: 'log₂(count+1)' },
      hovertemplate: '%{y} · %{x}<br>%{z:.2f}<extra></extra>',
    }], {
      paper_bgcolor: '#0f172a', plot_bgcolor: '#0f172a',
      font: { color: '#cbd5e1', size: 10 },
      margin: { l: 200, r: 40, t: 20, b: 100 },
      xaxis: { tickangle: -45, automargin: true },
      yaxis: { automargin: true },
    }, { responsive: true, displaylogo: false });
  }

  $effect(() => { drawHeatmap(selectedRef); });
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
    </div>
  </div>

  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
    {#if refKeys.length === 0}
      <div class="text-slate-500 text-sm py-8 text-center">
        No expression data available. featureCounts requires a <code>&lt;ref&gt;.gff</code> alongside each reference FASTA.
      </div>
    {:else}
      <div bind:this={heatmapEl} style="width:100%;height:640px"></div>
    {/if}
  </div>
</div>
