<script>
  /**
   * Toolbar + heatmap + detail panel + table for any entity × samples TPM
   * matrix (functions, genes). Same controls and layout as the MAG view.
   *
   * rows:         [{ id, tpm[], ...fields }]
   * groups:       optional { label, field, values: [{key, label}] } toggle row
   * searchFields: fields matched by the search box
   * annotations:  (shown) => D3Heatmap rowAnnotations
   * facts:        (row) => ExpressionDetail facts
   * links:        (row) => [{label, href}]
   * labelOf:      (row) => heatmap row label
   */
  import * as d3 from 'd3';
  import D3Heatmap from './charts/D3Heatmap.svelte';
  import DataTable from './ui/DataTable.svelte';
  import ExpressionDetail from './ExpressionDetail.svelte';
  import SampleTsne from './SampleTsne.svelte';
  import PlotlyChart from './charts/PlotlyChart.svelte';
  import { clusterMatrix } from '../lib/cluster.js';
  import { mean, fmtTpm, downloadTsv } from '../stores/data.js';

  let { rows = [], samples = [], noun = 'rows', title = '', note = '',
        groups = null, searchFields = ['id'], searchPlaceholder = 'search...',
        annotations = null, columns = [], facts = () => [], links = () => [],
        labelOf = r => r.id, subtitleOf = () => '', idLabel = 'ID', idWidth = 120, unit = 'TPM', rowLabel = r => r.id } = $props();

  let all = $derived(rows.map(r => ({
    ...r,
    meanTpm: mean(r.tpm),
    maxTpm: Math.max(...r.tpm),
    detected: r.tpm.filter(v => v > 0).length,
  })));

  // Group toggles (e.g. feature type), all on by default
  let activeGroups = $state(null);
  let groupCounts = $derived.by(() => {
    if (!groups) return {};
    const c = {};
    for (const r of all) c[r[groups.field] ?? ''] = (c[r[groups.field] ?? ''] || 0) + 1;
    return c;
  });
  let active = $derived(activeGroups ?? new Set(Object.keys(groupCounts)));
  function toggleGroup(k) {
    const next = new Set(active);
    if (next.has(k)) { if (next.size > 1) next.delete(k); } else next.add(k);
    activeGroups = next;
  }

  let minDetected = $state(1);
  // applied once the slider pauses (the heatmap, clustering and t-SNE follow it)
  let minDetectedApplied = $state(1);
  let detTimer;
  $effect(() => { const v = minDetected; clearTimeout(detTimer); detTimer = setTimeout(() => { minDetectedApplied = v; }, 150); return () => clearTimeout(detTimer); });

  let searchQuery = $state('');
  let searchDebounced = $state('');
  let timer;
  $effect(() => {
    const q = searchQuery;
    clearTimeout(timer);
    timer = setTimeout(() => { searchDebounced = q.trim().toLowerCase(); }, 300);
    return () => clearTimeout(timer);
  });

  const topDefs = [50, 100, 200];
  let topN = $state(50);
  const cycleTop = () => { topN = topDefs[(topDefs.indexOf(topN) + 1) % topDefs.length]; };

  const sortDefs = [
    { key: 'clustered', label: 'Clustered' },
    { key: 'mean', label: `Mean ${unit}` },
    { key: 'max', label: `Max ${unit}` },
    { key: 'detected', label: 'Prevalence' },
    { key: 'cv', label: 'Variability' },
    { key: 'name', label: 'A-Z' },
  ];
  let sortBy = $state('clustered');
  const cycleSort = () => { const k = sortDefs.map(d => d.key); sortBy = k[(k.indexOf(sortBy) + 1) % k.length]; };

  let filtered = $derived(all.filter(r => {
    if (groups && !active.has(r[groups.field] ?? '')) return false;
    if (r.detected < minDetectedApplied) return false;
    if (searchDebounced && !searchFields.some(f => r[f] != null && String(r[f]).toLowerCase().includes(searchDebounced))) return false;
    return true;
  }));
  const cv = r => { const m = r.meanTpm; if (!m) return 0; return Math.sqrt(mean(r.tpm.map(x => (x - m) ** 2))) / m; };
  let sorted = $derived.by(() => {
    const by = {
      clustered: (a, b) => b.meanTpm - a.meanTpm,   // top N by mean, then Ward order
      mean: (a, b) => b.meanTpm - a.meanTpm,
      max: (a, b) => b.maxTpm - a.maxTpm,
      detected: (a, b) => b.detected - a.detected || b.meanTpm - a.meanTpm,
      cv: (a, b) => cv(b) - cv(a),
      name: (a, b) => String(a.id).localeCompare(String(b.id)),
    }[sortBy];
    return [...filtered].sort(by);
  });
  // Clustered: Ward on log1p TPM orders both the shown rows and the samples
  let clust = $derived(sortBy === 'clustered' && sorted.length > 2 ? clusterMatrix(sorted.slice(0, topN).map(r => r.tpm)) : null);
  let shown = $derived(clust ? clust.rowOrder.map(i => sorted[i]) : sorted.slice(0, topN));

  let maxTpm = $derived(Math.max(1, ...shown.flatMap(r => r.tpm)));
  let colorScale = $derived.by(() => {
    const viridis = d3.scaleSequential(d3.interpolateViridis).domain([0, Math.log10(maxTpm + 1)]);
    return v => (v <= 0 ? '#0f172a' : viridis(Math.log10(v + 1)));
  });
  let heatmapData = $derived(shown.length ? {
    mag_ids: shown.map(r => r.id),
    module_ids: samples, module_names: samples,
    matrix: shown.map(r => r.tpm),
    row_order: shown.map((_, i) => i),
    col_order: clust ? clust.colOrder : samples.map((_, i) => i),
  } : null);
  // Column 0 is the row label (clickable, bold when selected)
  let rowAnnotations = $derived(heatmapData ? [
    { label: idLabel, values: shown.map(r => rowLabel(r)), width: idWidth },
    ...(annotations ? annotations(shown) ?? [] : []),
  ] : null);
  const tooltip = (val, id, sample) => {
    const r = shown.find(x => x.id === id);
    return `<strong>${sample}</strong><br>${r ? labelOf(r) : id}<br>${unit}: ${fmtTpm(val)}`;
  };

  // Many samples: full-width heatmap, detail panel below it
  let wide = $derived(samples.length > 60);

  let selected = $state(null);
  let selectedRow = $derived(all.find(r => r.id === selected) ?? null);
  let autoOpened = false;
  $effect(() => { if (!autoOpened && sorted.length) { autoOpened = true; selected = sorted[0].id; } });

  let tableColumns = $derived([
    ...columns,
    { key: 'meanTpm', label: `Mean ${unit}`, render: v => fmtTpm(v) },
    { key: 'maxTpm', label: `Max ${unit}`, render: v => fmtTpm(v) },
    { key: 'detected', label: 'Samples' },
  ]);
  let tableRows = $derived(sorted.slice(0, 500));

  // Every filtered row with its TPM in each sample
  function exportMatrix() {
    const cols = columns.map(c => c.key).filter(k => k !== 'id');
    downloadTsv(`${noun.replace(/\W+/g, '_')}_tpm_${sorted.length}.tsv`,
      ['id', ...cols, ...samples],
      sorted.map(r => [r.id, ...cols.map(k => r[k]), ...r.tpm]));
  }
</script>

<div class="sticky top-14 z-10 bg-slate-950 pb-2 -mx-4 px-4 pt-1">
<div class="flex items-center gap-3 flex-wrap text-xs">
  {#if groups}
    <span class="text-slate-400">{groups.label}:</span>
    {#each groups.values as { key, label }}
      {#if groupCounts[key]}
        <button
          class="px-3 py-1 rounded-md border transition-colors text-center
            {active.has(key) ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
          onclick={() => toggleGroup(key)}
        >{label} <span class="text-slate-500">({groupCounts[key].toLocaleString()})</span></button>
      {/if}
    {/each}
    <span class="text-slate-600">|</span>
  {/if}
  <label class="text-slate-400 flex items-center gap-1">
    In &ge;
    <input type="range" min="1" max={Math.max(1, samples.length)} step="1" bind:value={minDetected} class="w-28 h-1 accent-cyan-400" />
    <span class="text-cyan-400 font-mono w-16">{minDetected} sample{minDetected === 1 ? '' : 's'}</span>
  </label>
  <span class="text-slate-600">|</span>
  <div class="relative flex items-center">
    <input type="text" bind:value={searchQuery} placeholder={searchPlaceholder}
      class="w-56 px-2 py-1 rounded bg-slate-800 border border-slate-600 text-slate-200 text-xs placeholder-slate-500 focus:border-cyan-400 focus:outline-none" />
    {#if searchQuery}
      <button class="absolute right-1 text-slate-500 hover:text-slate-300" onclick={() => { searchQuery = ''; }} title="Clear search">
        <svg class="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M18 6L6 18M6 6l12 12" /></svg>
      </button>
    {/if}
  </div>
  {#if searchDebounced}<span class="text-amber-400 font-medium">{filtered.length.toLocaleString()} matches</span>{/if}
  <span class="text-slate-500 ml-auto">{filtered.length.toLocaleString()} of {all.length.toLocaleString()} {noun}</span>
</div>
</div>

{#if shown.length}
<div class="grid grid-cols-1 {selectedRow && !wide ? 'lg:grid-cols-3' : ''} gap-6 mb-6">
  <div class="{selectedRow && !wide ? 'lg:col-span-2' : ''} bg-slate-800 rounded-lg p-4 border border-slate-700 min-w-0">
    <div class="flex items-center justify-between mb-2 gap-2 flex-wrap">
      <h3 class="text-sm font-medium text-slate-400">{title} ({shown.length} of {filtered.length.toLocaleString()} x {samples.length} samples)</h3>
      <div class="flex gap-2 text-xs">
        <button class="px-2 py-0.5 rounded border border-cyan-400 bg-cyan-400/10 text-cyan-400" onclick={cycleTop}
          title="Rows shown: Top 50 → Top 100 → Top 200">Top {topN} &#x25BE;</button>
        <button class="px-2 py-0.5 rounded border border-cyan-400 bg-cyan-400/10 text-cyan-400" onclick={cycleSort}
          title={`Sort rows: ${sortDefs.map(d => d.label).join(' → ')}`}>Sort: {sortDefs.find(d => d.key === sortBy)?.label} &#x25BE;</button>
      </div>
    </div>
    <div class="overflow-auto" style="max-height: 700px;">
      <D3Heatmap data={heatmapData} onRowClick={id => selected = id} tooltipFormat={tooltip}
        {colorScale} legendLabel={['0', `${fmtTpm(maxTpm)} ${unit}`]} selectedRow={selected} {rowAnnotations} />
    </div>
    <div class="mt-2 text-xs text-slate-500">{note}{sortBy === 'clustered' ? ' Ward clustering of rows and samples on log(value+1).' : ''}</div>
  </div>

  {#if selectedRow}
    <ExpressionDetail title={selectedRow.id} subtitle={subtitleOf(selectedRow)}
      facts={facts(selectedRow)} tpm={selectedRow.tpm} {samples} links={links(selectedRow)} {unit}
      horizontal={wide} onclose={() => selected = null} />
  {/if}
</div>

<div class="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-6">
  <SampleTsne {samples} matrix={shown.map(r => r.tpm)} what="shown {noun}" />
  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
    <PlotlyChart data={[{
      type: 'bar', x: shown.map(r => labelOf(r)), y: shown.map(r => r.meanTpm),
      marker: { color: '#22d3ee', opacity: 0.8 },
      hovertemplate: `%{x}: %{y:,.2f} ${unit}<extra></extra>`,
    }]} layout={{
      title: { text: `Mean ${unit} (shown ${noun})`, font: { size: 14 } },
      xaxis: { tickangle: -45, tickfont: { size: 8 }, showticklabels: shown.length <= 60, automargin: true },
      yaxis: { title: `Mean ${unit}` },
    }} />
  </div>
</div>
{:else}
  <div class="flex flex-col items-center justify-center py-24 text-center">
    <p class="text-slate-400 text-lg mb-2">Nothing matches the current filters</p>
    <p class="text-slate-500 text-sm">Try clearing the search or lowering the prevalence filter.</p>
  </div>
{/if}

<div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
  <div class="flex items-center justify-between mb-2">
  <h3 class="text-sm font-medium text-slate-400">
    {searchDebounced ? `${filtered.length.toLocaleString()} matching ${noun}` : `Top ${Math.min(500, sorted.length)} ${noun} by ${sortBy === 'clustered' ? `Mean ${unit}` : sortDefs.find(d => d.key === sortBy)?.label}`}
    {#if sorted.length > 500}<span class="text-slate-500 text-xs">({sorted.length.toLocaleString()} total)</span>{/if}
  </h3>
  <button
    class="text-xs px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-slate-200 hover:border-slate-500 transition-colors"
    onclick={exportMatrix}
    title="Export all {sorted.length.toLocaleString()} rows with {unit} per sample"
  >TSV + samples</button>
  </div>
  <DataTable columns={tableColumns} rows={tableRows} idKey="id" maxHeight="400px" selectedId={selected} onRowClick={row => selected = row.id} />
</div>
