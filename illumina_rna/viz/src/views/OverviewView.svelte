<script>
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import StatCard from '../components/layout/StatCard.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, bins, median, fmtInt, fmtPct, qualityTier } from '../stores/data.js';
  import { tsne, brayCurtis } from '../lib/tsne.js';

  // Where each sample's raw reads went; every read lands in exactly one fate.
  // House hues (MAG dashboard palette) in an order validated for adjacent
  // stacks: CVD dE 10.5, normal-vision dE 21.2, all >= 3:1 on slate-800. They
  // sit above the dark-mode lightness band on purpose, to match the siblings.
  const FATES = [
    { key: 'assigned',  label: 'Assigned to a gene', color: '#22d3ee' },
    { key: 'nofeature', label: 'Mapped, no gene',    color: '#fb923c' },
    { key: 'unmapped',  label: 'Not in assembly',    color: '#34d399' },
    { key: 'human',     label: 'Human',              color: '#a78bfa' },
    { key: 'rrna',      label: 'rRNA',               color: '#fbbf24' },
    { key: 'qc',        label: 'Removed by QC',      color: '#f87171' },
  ];

  let rows = $derived(($samples?.rows ?? []).map((r, i) => {
    const raw = r.reads_raw ?? 0;
    const filtered = r.reads_filtered ?? raw;
    const norrna = r.reads_norrna ?? filtered;
    const final = r.reads_nohuman ?? norrna;
    const mapped = r.reads_mapped ?? 0;
    const assigned = 2 * (r.pairs_assigned ?? 0);
    const f = {
      qc: raw - filtered, rrna: filtered - norrna, human: norrna - final,
      unmapped: Math.max(0, final - mapped), nofeature: Math.max(0, mapped - assigned), assigned,
    };
    return { ...r, raw, f, genes_detected: $samples.genes_detected?.[i] };
  }));

  // ── Read fate chart ──
  const sortDefs = [
    { key: 'name', label: 'Sample' },
    { key: 'assigned', label: 'Assigned' },
    { key: 'rrna', label: 'rRNA' },
    { key: 'mapped', label: 'Mapped' },
  ];
  let sortBy = $state('name');
  let normalize = $state(true);
  const cycleSort = () => { const k = sortDefs.map(d => d.key); sortBy = k[(k.indexOf(sortBy) + 1) % k.length]; };

  let ordered = $derived.by(() => {
    const by = {
      name: (a, b) => a.sample.localeCompare(b.sample),
      assigned: (a, b) => (b.f.assigned / b.raw) - (a.f.assigned / a.raw),
      rrna: (a, b) => (b.pct_rrna ?? 0) - (a.pct_rrna ?? 0),
      mapped: (a, b) => (b.pct_mapped ?? 0) - (a.pct_mapped ?? 0),
    }[sortBy];
    return [...rows].sort(by);
  });

  let fateTraces = $derived(FATES.map(fate => ({
    type: 'bar',
    name: fate.label,
    x: ordered.map(r => r.sample),
    y: ordered.map(r => normalize ? (r.raw ? 100 * r.f[fate.key] / r.raw : 0) : r.f[fate.key]),
    customdata: ordered.map(r => [r.f[fate.key], r.raw ? 100 * r.f[fate.key] / r.raw : 0]),
    marker: { color: fate.color, line: { color: '#1e293b', width: ordered.length > 80 ? 0.5 : 1 } },
    hovertemplate: `<b>%{x}</b><br>${fate.label}: %{customdata[0]:,} reads (%{customdata[1]:.1f}%)<extra></extra>`,
  })));

  // ── Stat cards ──
  let magCounts = $derived.by(() => {
    const list = ($bins?.bins ?? []).filter(b => b.bin !== 'unbinned');
    const c = { HQ: 0, MQ: 0, LQ: 0, none: 0 };
    for (const b of list) c[qualityTier(b) ?? 'none']++;
    return { n: list.length, ...c };
  });
  let binnedShare = $derived.by(() => {
    const u = $bins?.bins?.find(b => b.bin === 'unbinned');
    return u ? 100 - median(u.tpm) / 1e4 : null;
  });

  // ── Sample t-SNE: Bray-Curtis on fourth-root MAG expression shares ──
  const colorDefs = [
    { key: 'pct_assigned', label: '% assigned' },
    { key: 'pct_rrna', label: '% rRNA' },
    { key: 'pct_mapped', label: '% mapped' },
    { key: 'run', label: 'Sequencing run' },
  ];
  let tsneColor = $state('pct_assigned');

  let embedding = $derived.by(() => {
    const B = ($bins?.bins ?? []).filter(b => b.bin !== 'unbinned');
    const names = $bins?.samples ?? [];
    if (names.length < 4 || !B.length) return null;
    const vecs = names.map((_, j) => {
      const col = B.map(b => b.tpm[j]);
      const tot = col.reduce((a, b) => a + b, 0) || 1;
      return col.map(v => Math.pow(v / tot, 0.25));
    });
    return { names, coords: tsne(brayCurtis(vecs)) };
  });

  let tsneTraces = $derived.by(() => {
    if (!embedding) return [];
    const byName = Object.fromEntries(rows.map(r => [r.sample, r]));
    const marker = { size: 10, line: { color: '#0f172a', width: 1 } };
    const pt = idx => ({
      type: 'scatter', mode: 'markers',
      x: idx.map(i => embedding.coords[i][0]), y: idx.map(i => embedding.coords[i][1]),
      hovertext: idx.map(i => {
        const r = byName[embedding.names[i]] ?? {};
        return `<b>${embedding.names[i]}</b><br>rRNA ${fmtPct(r.pct_rrna)} · mapped ${fmtPct(r.pct_mapped)} · assigned ${fmtPct(r.pct_assigned)}`;
      }),
      hovertemplate: '%{hovertext}<extra></extra>',
    });
    const all = embedding.names.map((_, i) => i);
    if (tsneColor === 'run') {
      // Run = flowcell (+ lane) prefix of the sample name
      const runOf = n => n.split('_').slice(0, -1).join('_') || n;
      const groups = new Map();
      for (const i of all) { const k = runOf(embedding.names[i]); if (!groups.has(k)) groups.set(k, []); groups.get(k).push(i); }
      const pal = ['#22d3ee', '#f472b6', '#fbbf24', '#34d399', '#a78bfa', '#fb923c', '#60a5fa', '#f87171', '#a3e635', '#e879f9', '#2dd4bf', '#facc15'];
      return [...groups.keys()].sort().map((k, g) => ({ ...pt(groups.get(k)), name: k,
        marker: { ...marker, color: g < pal.length ? pal[g] : '#475569' } }));
    }
    return [{ ...pt(all), marker: { ...marker,
      color: all.map(i => byName[embedding.names[i]]?.[tsneColor] ?? null),
      colorscale: 'Viridis',
      colorbar: { title: { text: colorDefs.find(d => d.key === tsneColor)?.label, side: 'right' }, thickness: 10 },
    } }];
  });

  const sampleColumns = [
    { key: 'sample', label: 'Sample' },
    { key: 'reads_raw', label: 'Raw reads', render: v => fmtInt(v) },
    { key: 'pct_rrna', label: 'rRNA', render: v => fmtPct(v) },
    { key: 'pct_mapped', label: 'Mapped', render: v => fmtPct(v) },
    { key: 'pct_assigned', label: 'Assigned', render: v => fmtPct(v) },
    { key: 'genes_detected', label: 'Genes detected', render: v => fmtInt(v) },
  ];
</script>

{#if $samples}
  <div class="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-6 gap-4 mb-6">
    <StatCard label="Samples" value={rows.length.toLocaleString()} color="cyan" />
    <StatCard label="Genes" value={fmtInt($samples.n_genes)} sub="{fmtInt($samples.n_contigs)} contigs" color="slate" />
    <StatCard label="MAGs" value={magCounts.n.toLocaleString()}
      sub={magCounts.none === magCounts.n ? 'no CheckM2 yet' : `${magCounts.HQ} HQ / ${magCounts.MQ} MQ / ${magCounts.LQ} LQ`} color="emerald" />
    <StatCard label="Median rRNA" value={fmtPct(median(rows.map(r => r.pct_rrna)))} sub="of reads after QC" color="amber" />
    <StatCard label="Median Mapped" value={fmtPct(median(rows.map(r => r.pct_mapped)))} sub="of non-rRNA, non-human" color="cyan" />
    <StatCard label="Expression in MAGs" value={fmtPct(binnedShare)} sub="median share of TPM" color="purple" />
  </div>

  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700 mb-6">
    <div class="flex items-center justify-between mb-2 gap-3 flex-wrap">
      <h3 class="text-sm font-medium text-slate-400">Read Fate ({rows.length} samples)</h3>
      <div class="flex gap-2 text-xs">
        <button
          class="px-2 py-0.5 rounded border border-cyan-400 bg-cyan-400/10 text-cyan-400"
          onclick={cycleSort}
          title={`Sort: ${sortDefs.map(d => d.label).join(' → ')}`}
        >Sort: {sortDefs.find(d => d.key === sortBy)?.label} &#x25BE;</button>
        <button
          class="px-2 py-0.5 rounded border transition-colors border-cyan-400 bg-cyan-400/10 text-cyan-400"
          onclick={() => normalize = !normalize}
        >{normalize ? 'Percent' : 'Reads'} &#x25BE;</button>
      </div>
    </div>
    <PlotlyChart data={fateTraces} layout={{
      barmode: 'stack', bargap: rows.length > 80 ? 0.05 : 0.15,
      height: 460,
      xaxis: { tickangle: -60, tickfont: { size: rows.length > 80 ? 7 : 9 }, automargin: true },
      yaxis: { title: normalize ? '% of raw reads' : 'Reads', range: normalize ? [0, 100] : undefined },
      legend: { orientation: 'h', y: 1.12, x: 0, traceorder: 'normal' },
      margin: { t: 40, r: 20, b: 100, l: 60 },
      hovermode: 'closest',
    }} />
  </div>

  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6">
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <div class="flex items-center gap-3 text-xs text-slate-400 mb-1">
        <label class="flex items-center gap-2">
          Color by:
          <select class="bg-slate-800 border border-slate-600 text-slate-200 text-xs rounded px-2 py-1" bind:value={tsneColor}>
            {#each colorDefs as d}<option value={d.key}>{d.label}</option>{/each}
          </select>
        </label>
      </div>
      {#if embedding}
        <PlotlyChart data={tsneTraces} layout={{
          title: { text: `Sample t-SNE (${embedding.names.length} samples, Bray-Curtis on ⁴√ MAG expression)`, font: { size: 14 } },
          xaxis: { showticklabels: false, zeroline: false, title: 't-SNE 1' },
          yaxis: { showticklabels: false, zeroline: false, title: 't-SNE 2' },
          showlegend: tsneColor === 'run',
          legend: { font: { size: 10 } },
        }} />
      {:else}
        <p class="text-slate-500 text-sm py-12 text-center">Needs at least 4 samples and MAG bins</p>
      {/if}
    </div>
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <h3 class="text-sm font-medium text-slate-400 mb-2">Samples</h3>
      <DataTable columns={sampleColumns} rows={rows} idKey="sample" maxHeight="420px" />
    </div>
  </div>
{/if}
