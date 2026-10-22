<script>
  /**
   * Sample t-SNE beside a heatmap: Bray-Curtis on fourth-root relative values
   * of the entities that view shows (MAGs, functions, genes, taxa), as in the
   * MAG dashboard's Abundance tab. Coloured by the shared sample annotation (run, QC or metadata).
   *
   * samples: names; matrix: entities × samples (non-negative); what: label
   */
  import PlotlyChart from './charts/PlotlyChart.svelte';
  import SampleAnnotSelect from './SampleAnnotSelect.svelte';
  import AnnotLegend from './AnnotLegend.svelte';
  import { samples as sampleStore, fmtPct } from '../stores/data.js';
  import { sampleAnnot } from '../stores/metadata.js';
  import { tsne, brayCurtis } from '../lib/tsne.js';

  let { samples = [], matrix = [], what = 'shown rows' } = $props();

  let meta = $derived(Object.fromEntries(($sampleStore?.rows ?? []).map(r => [r.sample, r])));

  // The embedding (exact t-SNE, ~0.2 s for 130 samples) reruns only once the
  // shown rows stop changing, and only if they changed; colouring never reruns it.
  let embedding = $state(null);
  let lastKey = '';
  let tsneTimer;
  $effect(() => {
    const m = matrix, names = samples;
    let sum = 0;
    for (const r of m) for (const v of r) sum += v;
    const key = `${names.length}:${m.length}:${sum.toPrecision(12)}`;
    if (key === lastKey) return;
    clearTimeout(tsneTimer);
    tsneTimer = setTimeout(() => { lastKey = key; embedding = embed(names, m); }, embedding ? 350 : 0);
    return () => clearTimeout(tsneTimer);
  });

  function embed(samples, matrix) {
    if (samples.length < 4 || !matrix.length) return null;
    const rows = [];
    samples.forEach((name, j) => {
      const col = matrix.map(r => r[j] ?? 0);
      const tot = col.reduce((a, b) => a + b, 0);
      if (tot > 0) rows.push({ name, v: col.map(x => Math.pow(x / tot, 0.25)) });
    });
    if (rows.length < 4) return null;
    return { names: rows.map(r => r.name), coords: tsne(brayCurtis(rows.map(r => r.v))) };
  }

  // Coloured by the sample annotation every view shares (run, QC, metadata)
  let traces = $derived.by(() => {
    const an = $sampleAnnot;
    if (!embedding || !an) return [];
    const marker = { size: 10, line: { color: '#0f172a', width: 1 } };
    const pt = idx => ({
      type: 'scatter', mode: 'markers',
      x: idx.map(i => embedding.coords[i][0]), y: idx.map(i => embedding.coords[i][1]),
      hovertext: idx.map(i => {
        const n = embedding.names[i], m = meta[n] ?? {};
        return `<b>${n}</b><br>${an.label}: ${an.display(n)}<br>rRNA ${fmtPct(m.pct_rrna_total ?? m.pct_rrna)} · mapped ${fmtPct(m.pct_mapped)} · CDS ${fmtPct(m.pct_cds)}`;
      }),
      hovertemplate: '%{hovertext}<extra></extra>',
    });
    const all = embedding.names.map((_, i) => i);
    const missing = all.filter(i => an.value(embedding.names[i]) == null);
    const has = all.filter(i => an.value(embedding.names[i]) != null);
    const missingTrace = missing.length ? [{ ...pt(missing), name: 'No value', marker: { ...marker, color: '#334155' } }] : [];
    if (an.type === 'category') {
      const groups = new Map();
      for (const i of has) {
        const k = String(an.value(embedding.names[i]));
        if (!groups.has(k)) groups.set(k, []);
        groups.get(k).push(i);
      }
      return [...[...groups.keys()].sort((x, y) => x.localeCompare(y, undefined, { numeric: true })).map(k => ({
        ...pt(groups.get(k)), name: k, marker: { ...marker, color: an.legend.colorFor(k) } })), ...missingTrace];
    }
    return [...missingTrace, { ...pt(has), name: an.label, marker: { ...marker,
      color: has.map(i => an.colorOf(embedding.names[i])) } }];
  });
</script>

<div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
  <div class="flex items-center gap-3 text-xs text-slate-400 mb-1 flex-wrap">
    <SampleAnnotSelect label="Color by:" />
    {#if $sampleAnnot?.type !== 'category'}<AnnotLegend />{/if}
  </div>
  {#if embedding}
    <PlotlyChart data={traces} layout={{
      title: { text: `Sample t-SNE (${embedding.names.length} samples, Bray-Curtis on ⁴√ ${what})`, font: { size: 14 } },
      xaxis: { showticklabels: false, zeroline: false, title: 't-SNE 1' },
      yaxis: { showticklabels: false, zeroline: false, title: 't-SNE 2' },
      showlegend: $sampleAnnot?.type === 'category',
      legend: { font: { size: 10 } },
    }} />
  {:else}
    <p class="text-slate-500 text-sm py-12 text-center">Needs at least 4 samples with values in the {what}</p>
  {/if}
</div>
