<script>
  /**
   * Sample t-SNE beside a heatmap: Bray-Curtis on fourth-root relative values
   * of the entities that view shows (MAGs, functions, genes, taxa), as in the
   * MAG dashboard's Abundance tab. Colour by a QC metric or sequencing run.
   *
   * samples: names; matrix: entities × samples (non-negative); what: label
   */
  import PlotlyChart from './charts/PlotlyChart.svelte';
  import { samples as sampleStore, fmtPct } from '../stores/data.js';
  import { tsne, brayCurtis } from '../lib/tsne.js';

  let { samples = [], matrix = [], what = 'shown rows' } = $props();

  const colorDefs = [
    { key: 'run', label: 'Sequencing run' },
    { key: 'pct_cds', label: '% on CDS' },
    { key: 'pct_rrna_total', label: '% rRNA' },
    { key: 'pct_mapped', label: '% mapped' },
  ];
  let colorBy = $state('run');

  let meta = $derived(Object.fromEntries(($sampleStore?.rows ?? []).map(r => [r.sample, r])));

  // Only rerun the embedding when the data changes, not the colouring
  let embedding = $derived.by(() => {
    if (samples.length < 4 || !matrix.length) return null;
    const rows = [];
    samples.forEach((name, j) => {
      const col = matrix.map(r => r[j] ?? 0);
      const tot = col.reduce((a, b) => a + b, 0);
      if (tot > 0) rows.push({ name, v: col.map(x => Math.pow(x / tot, 0.25)) });
    });
    if (rows.length < 4) return null;
    return { names: rows.map(r => r.name), coords: tsne(brayCurtis(rows.map(r => r.v))) };
  });

  const runOf = n => n.split('_').slice(0, -1).join('_') || n;
  const PAL = ['#22d3ee', '#f472b6', '#fbbf24', '#34d399', '#a78bfa', '#fb923c',
               '#60a5fa', '#f87171', '#a3e635', '#e879f9', '#2dd4bf', '#facc15'];

  let traces = $derived.by(() => {
    if (!embedding) return [];
    const marker = { size: 10, line: { color: '#0f172a', width: 1 } };
    const pt = idx => ({
      type: 'scatter', mode: 'markers',
      x: idx.map(i => embedding.coords[i][0]), y: idx.map(i => embedding.coords[i][1]),
      hovertext: idx.map(i => {
        const m = meta[embedding.names[i]] ?? {};
        return `<b>${embedding.names[i]}</b><br>rRNA ${fmtPct(m.pct_rrna_total ?? m.pct_rrna)} · mapped ${fmtPct(m.pct_mapped)} · CDS ${fmtPct(m.pct_cds)}`;
      }),
      hovertemplate: '%{hovertext}<extra></extra>',
    });
    const all = embedding.names.map((_, i) => i);
    if (colorBy === 'run') {
      const groups = new Map();
      for (const i of all) {
        const k = runOf(embedding.names[i]);
        if (!groups.has(k)) groups.set(k, []);
        groups.get(k).push(i);
      }
      return [...groups.keys()].sort().map((k, g) => ({ ...pt(groups.get(k)), name: k,
        marker: { ...marker, color: g < PAL.length ? PAL[g] : '#475569' } }));
    }
    return [{ ...pt(all), marker: { ...marker,
      color: all.map(i => meta[embedding.names[i]]?.[colorBy] ?? null),
      colorscale: 'Viridis',
      colorbar: { title: { text: colorDefs.find(d => d.key === colorBy)?.label, side: 'right' }, thickness: 10 },
    } }];
  });
</script>

<div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
  <div class="flex items-center gap-3 text-xs text-slate-400 mb-1">
    <label class="flex items-center gap-2">
      Color by:
      <select class="bg-slate-800 border border-slate-600 text-slate-200 text-xs rounded px-2 py-1" bind:value={colorBy}>
        {#each colorDefs as d}<option value={d.key}>{d.label}</option>{/each}
      </select>
    </label>
  </div>
  {#if embedding}
    <PlotlyChart data={traces} layout={{
      title: { text: `Sample t-SNE (${embedding.names.length} samples, Bray-Curtis on ⁴√ ${what})`, font: { size: 14 } },
      xaxis: { showticklabels: false, zeroline: false, title: 't-SNE 1' },
      yaxis: { showticklabels: false, zeroline: false, title: 't-SNE 2' },
      showlegend: colorBy === 'run',
      legend: { font: { size: 10 } },
    }} />
  {:else}
    <p class="text-slate-500 text-sm py-12 text-center">Needs at least 4 samples with values in the {what}</p>
  {/if}
</div>
