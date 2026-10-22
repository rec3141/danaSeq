<script>
  import PlotlyChart from '../components/PlotlyChart.svelte';
  import { contigMap, refName, loadContigMap, loadContigLayer, SEQ_BLUE, fmtInt } from '../stores/data.js';

  // The heatmap ramp minus its surface-coloured floor (a point that matched
  // the background would vanish), rescaled to start at 0.
  const POINT_RAMP = SEQ_BLUE.slice(1).map(([p, c], _, a) => [(p - a[0][0]) / (1 - a[0][0]), c]);

  let status = $state('loading');
  let sampleIdx = $state(0);
  let layer = $state(null);
  let lockScale = $state(false);
  let lockedMax = $state(null);

  $effect(() => {
    $refName;   // reload when the reference changes
    status = 'loading';
    loadContigMap().then(m => { status = m ? 'ready' : 'none'; sampleIdx = 0; })
                   .catch(e => { status = 'error: ' + e.message; });
  });

  $effect(() => {
    const m = $contigMap;
    const j = sampleIdx;
    if (!m || status !== 'ready') return;
    loadContigLayer(m, j).then(l => { if (j === sampleIdx) layer = l; });
  });

  // Reads per kb per million reads on the mapped contigs, log10(+1)
  let values = $derived.by(() => {
    const m = $contigMap;
    if (!m || !layer) return null;
    const total = m.sample_totals[sampleIdx] || 1;
    const v = layer.idx.map((i, k) => Math.log10(1 + layer.reads[k] / (m.length[i] / 1e3) / (total / 1e6)));
    // draw high values last so they sit on top
    const order = v.map((_, k) => k).sort((a, b) => v[a] - v[b]);
    return { order, v };
  });

  let p99 = $derived.by(() => {
    if (!values) return 1;
    const s = [...values.v].sort((a, b) => a - b);
    return s[Math.floor(s.length * 0.99)] || 1;
  });
  let cmax = $derived(lockScale && lockedMax != null ? lockedMax : p99);
  $effect(() => { if (lockScale && lockedMax == null) lockedMax = p99; if (!lockScale) lockedMax = null; });

  let traces = $derived.by(() => {
    const m = $contigMap;
    if (!m) return [];
    const bg = {
      type: 'scattergl', mode: 'markers', name: 'RNA-active contigs',
      x: m.x, y: m.y, hoverinfo: 'skip',
      marker: { size: 2, color: '#475569', opacity: 0.7 },
    };
    if (!values) return [bg];
    const { order, v } = values;
    const pick = order.map(k => layer.idx[k]);
    return [bg, {
      type: 'scattergl', mode: 'markers', name: m.samples[sampleIdx],
      x: pick.map(i => m.x[i]), y: pick.map(i => m.y[i]),
      customdata: order.map(k => [layer.reads[k], m.length[layer.idx[k]], m.bin[layer.idx[k]] >= 0 ? m.bins[m.bin[layer.idx[k]]] : 'unbinned']),
      marker: {
        size: 4, color: order.map(k => v[k]),
        colorscale: POINT_RAMP, cmin: 0, cmax,
        colorbar: { title: { text: 'log₁₀(RPKM+1)', side: 'right' }, thickness: 10, outlinewidth: 0, len: 0.8, y: 0.4 },
      },
      hovertemplate: '%{customdata[2]}<br>%{customdata[0]:,} reads · %{customdata[1]:,} bp<extra></extra>',
    }];
  });

  const layout = {
    xaxis: { visible: false }, yaxis: { visible: false, scaleanchor: 'x' },
    showlegend: false, margin: { l: 10, r: 10, t: 10, b: 10 },
  };

  const step = d => {
    const n = $contigMap?.samples?.length ?? 0;
    if (n) sampleIdx = (sampleIdx + d + n) % n;
  };
</script>

{#if status === 'none'}
  <div class="bg-slate-900 border border-slate-800 rounded-lg p-6 text-sm text-slate-400">
    No contig embedding for this reference. The map reuses the contig t-SNE from a mag_analysis run's viz; run with <code class="text-slate-300">--from_mag</code> on an output that has <code class="text-slate-300">viz/data/contig_explorer_*</code>.
  </div>
{:else if status.startsWith('error')}
  <div class="text-rose-300 text-sm">{status}</div>
{:else}
  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
    <div class="flex items-center gap-3 flex-wrap mb-2">
      <div class="mr-auto">
        <div class="text-sm font-semibold text-slate-200">Contig t-SNE · expression by sample</div>
        <div class="text-xs text-slate-400">
          Contig layout from the MAG analysis. Grey: the {fmtInt($contigMap?.n)} contigs with RNA reads in any sample; colour: reads per kb per million for the selected sample.
        </div>
      </div>
      <button class="px-2 py-1 rounded border border-slate-700 text-slate-300 text-xs" onclick={() => step(-1)} aria-label="Previous sample">◀</button>
      <select class="bg-slate-950 border border-slate-700 rounded px-2 py-1 text-xs text-slate-200 font-mono" bind:value={sampleIdx}>
        {#each $contigMap?.samples ?? [] as s, j}<option value={j}>{s}</option>{/each}
      </select>
      <button class="px-2 py-1 rounded border border-slate-700 text-slate-300 text-xs" onclick={() => step(1)} aria-label="Next sample">▶</button>
      <label class="text-xs text-slate-400 flex items-center gap-1">
        <input type="checkbox" bind:checked={lockScale} class="rounded bg-slate-950 border-slate-700" /> lock scale
      </label>
    </div>
    {#if status === 'loading'}
      <div class="text-slate-400 text-center py-24 text-sm">Loading contig map…</div>
    {:else}
      <PlotlyChart data={traces} {layout} height="760px" />
      <div class="text-xs text-slate-500 mt-1">
        {#if layer}{fmtInt(layer.idx.length)} contigs with reads in {layer.sample}{/if}
      </div>
    {/if}
  </div>
{/if}
