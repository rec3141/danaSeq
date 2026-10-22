<script>
  import { onMount } from 'svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import ExpressionMatrix from '../components/ExpressionMatrix.svelte';
  import { index, refName, fmtPct } from '../stores/data.js';
  import { fetchJSON } from '../stores/fetchJson.js';
  import { wardOrder } from '../lib/cluster.js';

  // Loaded on first visit (not every reference has it)
  let tax = $state(null);
  let state = $state('loading');
  async function load() {
    const ref = $index?.references?.find(r => r.name === $refName);
    const f = ref?.files?.taxonomy;
    if (!f) { state = 'none'; return; }
    try { tax = await fetchJSON(`./data/${f.replace(/\.gz$/, '')}`); state = 'ready'; }
    catch (e) { state = 'error: ' + e.message; }
  }
  onMount(load);

  const SOURCES = [
    { key: 'rrna', label: 'rRNA reads', title: 'Who is there: rRNA reads classified against SILVA 138.2' },
    { key: 'mrna_kaiju', label: 'mRNA · contigs (Kaiju)', title: 'Who is active: CDS read pairs by the contig\'s Kaiju lineage' },
    { key: 'mrna_gtdb', label: 'mRNA · MAGs (GTDB)', title: 'Who is active: CDS read pairs by the MAG\'s GTDB-Tk lineage' },
  ];
  const RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus'];

  let sources = $derived(SOURCES.filter(s => tax?.sources?.[s.key]));
  let source = $state(null);
  $effect(() => { if (!source || !tax?.sources?.[source]) source = sources[0]?.key ?? null; });
  let rank = $state('phylum');
  const cycleRank = () => { rank = RANKS[(RANKS.indexOf(rank) + 1) % RANKS.length]; };

  let src = $derived(source ? tax.sources[source] : null);
  let samples = $derived(tax?.samples ?? []);
  // Relative abundance (%) of each taxon in each sample
  let rows = $derived((src?.ranks?.[rank] ?? []).map(t => ({
    id: t.taxon,
    name: t.taxon.split(';').filter(x => x !== 'unclassified').pop() ?? 'unclassified',
    unclassified: t.taxon.endsWith('unclassified') || t.taxon === 'unbinned',
    tpm: t.counts.map((c, j) => 100 * c / (src.totals[j] || 1)),
    reads: t.counts.reduce((a, b) => a + b, 0),
  })));

  // ── Composition bars: top taxa + Other, unclassified in grey ──
  const PALETTE = ['#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa', '#fb923c',
                   '#2dd4bf', '#818cf8', '#f472b6', '#4ade80', '#e879f9', '#38bdf8'];
  let barSort = $state('clustered');
  const barSorts = [{ key: 'clustered', label: 'Clustered' }, { key: 'name', label: 'Sample' }, { key: 'top', label: 'Top taxon' }];
  const cycleBarSort = () => { const k = barSorts.map(d => d.key); barSort = k[(k.indexOf(barSort) + 1) % k.length]; };

  let bars = $derived.by(() => {
    if (!rows.length) return { traces: [], order: [] };
    const classified = rows.filter(r => !r.unclassified).sort((a, b) =>
      b.tpm.reduce((x, y) => x + y, 0) - a.tpm.reduce((x, y) => x + y, 0));
    const top = classified.slice(0, PALETTE.length);
    const unc = rows.filter(r => r.unclassified);
    const other = samples.map((_, j) => Math.max(0, 100 - top.reduce((s, r) => s + r.tpm[j], 0) - unc.reduce((s, r) => s + r.tpm[j], 0)));
    const uncl = samples.map((_, j) => unc.reduce((s, r) => s + r.tpm[j], 0));
    let order = samples.map((_, j) => j);
    if (barSort === 'name') order.sort((a, b) => samples[a].localeCompare(samples[b]));
    else if (barSort === 'top' && top.length) order.sort((a, b) => top[0].tpm[b] - top[0].tpm[a]);
    else if (barSort === 'clustered') order = wardOrder(samples.map((_, j) => [...top.map(r => Math.sqrt(r.tpm[j])), Math.sqrt(other[j]), Math.sqrt(uncl[j])]));
    const x = order.map(j => samples[j]);
    const tr = (name, y, color) => ({
      type: 'bar', name, x, y: order.map(j => y[j]),
      marker: { color, line: { color: '#1e293b', width: samples.length > 80 ? 0.5 : 1 } },
      hovertemplate: `<b>%{x}</b><br>${name}: %{y:.1f}%<extra></extra>`,
    });
    return { traces: [...top.map((r, i) => tr(r.name, r.tpm, PALETTE[i])), tr('Other', other, '#94a3b8'), tr('Unclassified', uncl, '#475569')] };
  });
</script>

{#if state === 'none'}
  <div class="flex flex-col items-center justify-center py-24 text-center">
    <p class="text-slate-400 text-lg mb-2">No taxonomy for this reference</p>
    <p class="text-slate-500 text-sm max-w-xl">rRNA-read taxonomy needs <code class="bg-slate-800 px-1 rounded">--rrna_method bbduk</code> with SILVA-style rRNA references; mRNA taxonomy needs a contig embedding (Kaiju) or GTDB-classified MAGs from mag_analysis.</p>
  </div>
{:else if state !== 'ready'}
  <div class="flex flex-col items-center justify-center py-20">
    <div class="w-10 h-10 border-4 border-slate-700 border-t-cyan-400 rounded-full animate-spin"></div>
    <p class="text-slate-400 text-sm mt-4">{state === 'loading' ? 'Loading taxonomy...' : state}</p>
  </div>
{:else}
  <div class="flex items-center gap-3 text-xs mb-2 flex-wrap">
    <span class="text-slate-400">Source:</span>
    {#each sources as s}
      <button
        class="px-3 py-1 rounded-md border transition-colors
          {source === s.key ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        onclick={() => source = s.key} title={s.title}
      >{s.label}</button>
    {/each}
    <span class="text-slate-600">|</span>
    <span class="text-slate-400">Rank:</span>
    <button class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
      style="min-width: 5.5rem" onclick={cycleRank} title={`Click to cycle: ${RANKS.join(' → ')}`}
    >{rank[0].toUpperCase() + rank.slice(1)} &#x25BE;</button>
    {#if src}<span class="text-slate-500">{SOURCES.find(s => s.key === source)?.title}; % of {src.unit} per sample</span>{/if}
  </div>

  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700 mb-6">
    <div class="flex items-center justify-between mb-2">
      <h3 class="text-sm font-medium text-slate-400">Composition by {rank} ({samples.length} samples)</h3>
      <button class="text-xs px-2 py-0.5 rounded border border-cyan-400 bg-cyan-400/10 text-cyan-400" onclick={cycleBarSort}
        title={`Sort samples: ${barSorts.map(d => d.label).join(' → ')}`}>Sort: {barSorts.find(d => d.key === barSort)?.label} &#x25BE;</button>
    </div>
    <PlotlyChart data={bars.traces} layout={{
      barmode: 'stack', bargap: samples.length > 80 ? 0.05 : 0.15, height: 460,
      xaxis: { tickangle: -60, tickfont: { size: samples.length > 80 ? 7 : 9 }, automargin: true },
      yaxis: { title: `% of ${src?.unit}`, range: [0, 100] },
      legend: { font: { size: 10 }, traceorder: 'normal' },
      margin: { t: 20, r: 20, b: 100, l: 60 },
    }} />
  </div>

  {#key `${source}:${rank}`}
    <ExpressionMatrix {rows} {samples} noun="taxa" title="Taxa ({rank})" unit="%"
      idLabel={rank[0].toUpperCase() + rank.slice(1)} idWidth={170}
      searchFields={['id']} searchPlaceholder="taxon..."
      labelOf={r => r.name} rowLabel={r => r.name}
      columns={[{ key: 'name', label: rank[0].toUpperCase() + rank.slice(1) }, { key: 'id', label: 'Lineage', render: v => v }]}
      annotations={shown => [{ label: 'Lineage', values: shown.map(r => r.id), align: 'start', width: 600 }]}
      subtitleOf={r => r.id}
      facts={r => [{ label: src.unit === 'reads' ? 'rRNA reads' : 'CDS read pairs', value: r.reads.toLocaleString() }]}
      note="Relative abundance per sample, log-scaled colour." />
  {/key}
{/if}
