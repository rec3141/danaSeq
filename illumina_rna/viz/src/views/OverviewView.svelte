<script>
  import PlotlyChart from '../components/PlotlyChart.svelte';
  import { samples, bins, median, fmtInt, fmtPct } from '../stores/data.js';

  // Where each sample's raw reads went. Every read lands in exactly one
  // fate, so the six segments sum to the raw count.
  // Colors: categorical slots 1-6, dark steps (validated adjacent on #0f172a).
  const FATES = [
    { key: 'assigned',   label: 'Assigned to a gene',   color: '#3987e5' },
    { key: 'nofeature',  label: 'Mapped, no gene',      color: '#d95926' },
    { key: 'unmapped',   label: 'Not in assembly',      color: '#199e70' },
    { key: 'human',      label: 'Human',                color: '#c98500' },
    { key: 'rrna',       label: 'rRNA',                 color: '#d55181' },
    { key: 'qc',         label: 'Removed by QC',        color: '#008300' },
  ];

  let rows = $derived(($samples?.rows ?? []).map(r => {
    const raw = r.reads_raw ?? 0;
    const filtered = r.reads_filtered ?? raw;
    const norrna = r.reads_norrna ?? filtered;
    const final = r.reads_nohuman ?? norrna;
    const mapped = r.reads_mapped ?? 0;
    const assigned = 2 * (r.pairs_assigned ?? 0);
    const f = {
      qc: raw - filtered,
      rrna: filtered - norrna,
      human: norrna - final,
      unmapped: Math.max(0, final - mapped),
      nofeature: Math.max(0, mapped - assigned),
      assigned,
    };
    return { sample: r.sample, raw, f, pct: r };
  }));

  let showTable = $state(false);
  let normalize = $state(true);

  let traces = $derived(FATES.map(fate => ({
    type: 'bar',
    name: fate.label,
    x: rows.map(r => r.sample),
    y: rows.map(r => normalize ? (r.raw ? 100 * r.f[fate.key] / r.raw : 0) : r.f[fate.key]),
    customdata: rows.map(r => [r.f[fate.key], r.raw ? 100 * r.f[fate.key] / r.raw : 0]),
    marker: { color: fate.color, line: { color: '#0f172a', width: 1 } },
    hovertemplate: `%{x}<br>${fate.label}: %{customdata[0]:,} reads (%{customdata[1]:.1f}%)<extra></extra>`,
  })));

  let layout = $derived({
    barmode: 'stack',
    bargap: 0.15,
    hovermode: 'closest',
    xaxis: { tickangle: -60, automargin: true, tickfont: { size: 9 } },
    yaxis: { title: normalize ? '% of raw reads' : 'reads', range: normalize ? [0, 100] : undefined, gridcolor: '#1e293b' },
    legend: { orientation: 'h', y: 1.12, x: 0, traceorder: 'normal' },
    margin: { l: 60, r: 20, t: 40, b: 110 },
  });

  let unbinnedShare = $derived.by(() => {
    const u = $bins?.bins?.find(b => b.bin === 'unbinned');
    return u ? median(u.tpm) / 1e4 : null;
  });

  let tiles = $derived([
    { label: 'Samples', value: fmtInt($samples?.samples?.length) },
    { label: 'Genes in reference', value: fmtInt($samples?.n_genes) },
    { label: 'MAGs', value: fmtInt(($bins?.bins?.length ?? 1) - 1) },
    { label: 'Median rRNA', value: fmtPct(median(rows.map(r => r.pct.pct_rrna))) },
    { label: 'Median mapped', value: fmtPct(median(rows.map(r => r.pct.pct_mapped))), hint: 'of non-rRNA, non-human reads' },
    { label: 'Expression in MAGs', value: unbinnedShare == null ? '—' : fmtPct(100 - unbinnedShare), hint: 'median share of TPM on binned contigs' },
  ]);
</script>

<div class="space-y-6">
  <div class="grid grid-cols-2 md:grid-cols-3 xl:grid-cols-6 gap-4">
    {#each tiles as t}
      <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
        <div class="text-xs uppercase tracking-wide text-slate-400">{t.label}</div>
        <div class="text-2xl font-semibold text-slate-100 mt-1 tabular-nums">{t.value}</div>
        {#if t.hint}<div class="text-[11px] text-slate-500 mt-1">{t.hint}</div>{/if}
      </div>
    {/each}
  </div>

  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
    <div class="flex items-center justify-between mb-2 gap-4 flex-wrap">
      <div>
        <div class="text-sm font-semibold text-slate-200">Where the reads went</div>
        <div class="text-xs text-slate-400">Each sample's raw reads, by fate through QC, rRNA and human removal, mapping and gene assignment</div>
      </div>
      <div class="flex gap-2 text-xs">
        <button class="px-2 py-1 rounded border {normalize ? 'border-cyan-500 text-cyan-300' : 'border-slate-700 text-slate-400'}" onclick={() => (normalize = !normalize)}>
          {normalize ? 'Percent' : 'Reads'}
        </button>
        <button class="px-2 py-1 rounded border border-slate-700 text-slate-400 hover:text-slate-200" onclick={() => (showTable = !showTable)}>
          {showTable ? 'Chart' : 'Table'}
        </button>
      </div>
    </div>
    {#if showTable}
      <div class="overflow-auto max-h-[560px]">
        <table class="w-full text-xs tabular-nums">
          <thead class="text-slate-400 sticky top-0 bg-slate-900">
            <tr><th class="text-left p-1">Sample</th><th class="text-right p-1">Raw reads</th>
              {#each FATES as f}<th class="text-right p-1">{f.label}</th>{/each}</tr>
          </thead>
          <tbody class="text-slate-200">
            {#each rows as r}
              <tr class="border-t border-slate-800">
                <td class="p-1 font-mono">{r.sample}</td><td class="text-right p-1">{fmtInt(r.raw)}</td>
                {#each FATES as f}<td class="text-right p-1">{r.raw ? (100 * r.f[f.key] / r.raw).toFixed(1) : '—'}%</td>{/each}
              </tr>
            {/each}
          </tbody>
        </table>
      </div>
    {:else}
      <PlotlyChart data={traces} {layout} height="520px" />
    {/if}
  </div>
</div>
