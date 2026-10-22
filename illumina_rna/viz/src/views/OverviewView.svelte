<script>
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import StatCard from '../components/layout/StatCard.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, bins, median, fmtInt, fmtPct, qualityTier } from '../stores/data.js';

  // Where each sample's raw reads went; every read lands in exactly one fate.
  // House hues (MAG dashboard palette) in an order validated for adjacent
  // stacks: CVD dE 10.5, normal-vision dE 21.2, all >= 3:1 on slate-800. They
  // sit above the dark-mode lightness band on purpose, to match the siblings.
  const FATES = [
    { key: 'assigned',  label: 'Protein-coding genes', color: '#22d3ee' },
    { key: 'nofeature', label: 'Mapped, other',      color: '#fb923c' },
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
    // CDS pairs when the run reports them; older runs only have all features
    const assigned = 2 * (r.pairs_cds ?? r.pairs_assigned ?? 0);
    // rRNA = what the filter removed + mapped reads on rRNA loci it missed
    const residual = 2 * (r.pairs_rrna_features ?? 0);
    const f = {
      qc: raw - filtered, rrna: filtered - norrna + residual, human: norrna - final,
      unmapped: Math.max(0, final - mapped), nofeature: Math.max(0, mapped - assigned - residual), assigned,
    };
    return { ...r, raw, f, genes_detected: $samples.genes_detected?.[i] };
  }));

  // ── Read fate chart ──
  const sortDefs = [
    { key: 'name', label: 'Sample' },
    { key: 'assigned', label: 'Protein-coding' },
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
      rrna: (a, b) => (b.f.rrna / b.raw) - (a.f.rrna / a.raw),
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

  const sampleColumns = [
    { key: 'sample', label: 'Sample' },
    { key: 'reads_raw', label: 'Raw reads', render: v => fmtInt(v) },
    { key: 'pct_rrna_total', label: 'rRNA', render: (v, r) => fmtPct(v ?? r.pct_rrna) },
    { key: 'pct_mapped', label: 'Mapped', render: v => fmtPct(v) },
    { key: 'pct_cds', label: 'On CDS', render: (v, r) => fmtPct(v ?? r.pct_assigned) },
    { key: 'genes_detected', label: 'Genes detected', render: v => fmtInt(v) },
  ];
</script>

{#if $samples}
  <div class="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-6 gap-4 mb-6">
    <StatCard label="Samples" value={rows.length.toLocaleString()} color="cyan" />
    <StatCard label="Genes" value={fmtInt($samples.n_genes)} sub="{fmtInt($samples.n_contigs)} contigs" color="slate" />
    <StatCard label="MAGs" value={magCounts.n.toLocaleString()}
      sub={magCounts.none === magCounts.n ? 'no CheckM2 yet' : `${magCounts.HQ} HQ / ${magCounts.MQ} MQ / ${magCounts.LQ} LQ`} color="emerald" />
    <StatCard label="Median rRNA" value={fmtPct(median(rows.map(r => r.pct_rrna_total ?? r.pct_rrna)))} sub="filtered + residual on rRNA loci" color="amber" />
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

  <div class="grid grid-cols-1 gap-6">
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <h3 class="text-sm font-medium text-slate-400 mb-2">Samples</h3>
      <DataTable columns={sampleColumns} rows={rows} idKey="sample" maxHeight="420px" />
    </div>
  </div>
{/if}
