<script>
  import * as d3 from 'd3';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import D3Heatmap from '../components/charts/D3Heatmap.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import ExpressionDetail from '../components/ExpressionDetail.svelte';
  import SampleTsne from '../components/SampleTsne.svelte';
  import { clusterMatrix } from '../lib/cluster.js';
  import { bins, index, mean, fmtTpm, fmtInt, fmtBp, qualityTier, gtdbRanks, lastRank, downloadTsv } from '../stores/data.js';

  let names = $derived($bins?.samples ?? []);
  let all = $derived(($bins?.bins ?? []).map(b => ({
    ...b,
    id: b.bin,
    tier: b.bin === 'unbinned' ? null : qualityTier(b),
    taxon: lastRank(b.classification),
    meanTpm: mean(b.tpm),
    maxTpm: Math.max(...b.tpm),
    totalReads: b.reads.reduce((a, v) => a + v, 0),
    detected: b.tpm.filter(v => v > 0).length,
  })));

  // ── Toolbar state (same controls as the MAG dashboard's Abundance tab) ──
  const tierDefs = [
    { key: 'HQ', label: 'HQ' }, { key: 'MQ', label: 'MQ' },
    { key: 'LQ', label: 'LQ' }, { key: 'none', label: 'No CheckM2' },
  ];
  let activeTiers = $state(new Set(['HQ', 'MQ', 'LQ', 'none']));
  let tierCounts = $derived.by(() => {
    const c = {};
    for (const b of all) if (b.id !== 'unbinned') c[b.tier ?? 'none'] = (c[b.tier ?? 'none'] || 0) + 1;
    return c;
  });
  function toggleTier(k) {
    const next = new Set(activeTiers);
    if (next.has(k)) { if (next.size > 1) next.delete(k); } else next.add(k);
    activeTiers = next;
  }
  let minCompleteness = $state(0);
  let maxContamination = $state(100);
  let showUnbinned = $state(false);

  let searchQuery = $state('');
  let searchDebounced = $state('');
  let timer;
  $effect(() => {
    const q = searchQuery;
    clearTimeout(timer);
    timer = setTimeout(() => { searchDebounced = q.trim().toLowerCase(); }, 300);
    return () => clearTimeout(timer);
  });

  const topDefs = [50, 100, 200, 0];
  let topN = $state(50);
  const cycleTop = () => { topN = topDefs[(topDefs.indexOf(topN) + 1) % topDefs.length]; };

  const sortDefs = [
    { key: 'clustered', label: 'Clustered' },
    { key: 'mean', label: 'Mean TPM' },
    { key: 'max', label: 'Max TPM' },
    { key: 'detected', label: 'Prevalence' },
    { key: 'cv', label: 'Variability' },
    { key: 'name', label: 'A-Z' },
    { key: 'completeness', label: 'Completeness' },
    { key: 'taxonomy', label: 'Taxonomy' },
    { key: 'size', label: 'Size' },
  ];
  let sortBy = $state('clustered');
  const cycleSort = () => { const k = sortDefs.map(d => d.key); sortBy = k[(k.indexOf(sortBy) + 1) % k.length]; };

  let filtered = $derived(all.filter(b => {
    if (b.id === 'unbinned') return showUnbinned;
    if (!activeTiers.has(b.tier ?? 'none')) return false;
    if (b.completeness != null && b.completeness < minCompleteness) return false;
    if (b.contamination != null && b.contamination > maxContamination) return false;
    if ((minCompleteness > 0 || maxContamination < 100) && b.completeness == null) return false;
    if (searchDebounced && ![b.id, b.classification].some(v => v && v.toLowerCase().includes(searchDebounced))) return false;
    return true;
  }));

  const cv = b => { const m = b.meanTpm; if (!m) return 0; const v = mean(b.tpm.map(x => (x - m) ** 2)); return Math.sqrt(v) / m; };
  let sorted = $derived.by(() => {
    const by = {
      clustered: (a, b) => b.meanTpm - a.meanTpm,   // top N by mean, then Ward order
      mean: (a, b) => b.meanTpm - a.meanTpm,
      max: (a, b) => b.maxTpm - a.maxTpm,
      detected: (a, b) => b.detected - a.detected || b.meanTpm - a.meanTpm,
      cv: (a, b) => cv(b) - cv(a),
      name: (a, b) => a.id.localeCompare(b.id),
      completeness: (a, b) => (b.completeness ?? -1) - (a.completeness ?? -1),
      taxonomy: (a, b) => (a.classification ?? '~').localeCompare(b.classification ?? '~'),
      size: (a, b) => (b.length ?? 0) - (a.length ?? 0),
    }[sortBy];
    return [...filtered].sort(by);
  });
  let topRows = $derived(topN ? sorted.slice(0, topN) : sorted);
  // Clustered: Ward on log1p TPM orders both the shown MAGs and the samples
  let clust = $derived(sortBy === 'clustered' && topRows.length > 2 ? clusterMatrix(topRows.map(b => b.tpm)) : null);
  let shown = $derived(clust ? clust.rowOrder.map(i => topRows[i]) : topRows);

  // ── Heatmap ──
  let maxTpm = $derived(Math.max(1, ...shown.flatMap(b => b.tpm)));
  let colorScale = $derived.by(() => {
    const viridis = d3.scaleSequential(d3.interpolateViridis).domain([0, Math.log10(maxTpm + 1)]);
    return v => (v <= 0 ? '#0f172a' : viridis(Math.log10(v + 1)));
  });
  let heatmapData = $derived(shown.length ? {
    mag_ids: shown.map(b => b.id),
    module_ids: names, module_names: names,
    matrix: shown.map(b => b.tpm),
    row_order: shown.map((_, i) => i),
    col_order: clust ? clust.colOrder : names.map((_, i) => i),
  } : null);

  const tierColor = v => v === 'HQ' ? '#34d399' : v === 'MQ' ? '#fbbf24' : '#64748b';
  // Column 0 is the row label (clickable, bold when selected)
  let rowAnnotations = $derived(heatmapData ? [
    { label: 'MAG', values: shown.map(b => b.id), width: 170 },
    { label: 'Quality', values: shown.map(b => b.tier ?? '-'), colorFn: tierColor },
    { label: 'Comp%', values: shown.map(b => b.completeness?.toFixed(0) ?? '-'),
      colorFn: v => { const n = parseFloat(v); return n >= 90 ? '#34d399' : n >= 50 ? '#fbbf24' : '#64748b'; } },
    { label: 'Cont%', values: shown.map(b => b.contamination?.toFixed(1) ?? '-'),
      colorFn: v => { const n = parseFloat(v); return n < 5 ? '#34d399' : n < 10 ? '#fbbf24' : '#ef4444'; } },
    { label: 'Size', values: shown.map(b => b.length ? (b.length / 1e6).toFixed(1) + 'M' : '-') },
    { label: 'Mean TPM', values: shown.map(b => fmtTpm(b.meanTpm)) },
    { label: 'Tax (GTDB)', values: shown.map(b => b.classification ? b.classification.replace(/;s__$/, '') : '-'),
      align: 'start', width: 700 },
  ] : null);

  const tooltip = (val, magId, sample) => `<strong>${sample}</strong><br>MAG: ${magId}<br>TPM: ${fmtTpm(val)}`;

  // ── Selection / detail ──
  // Many samples: full-width heatmap, detail panel below it
  let wide = $derived(names.length > 60);

  let selected = $state(null);
  let selectedBin = $derived(all.find(b => b.id === selected) ?? null);
  // Open the top MAG once, as the MAG dashboard opens its first contig
  let autoOpened = false;
  $effect(() => { if (!autoOpened && sorted.length) { autoOpened = true; selected = sorted[0].id; } });

  let detailFacts = $derived.by(() => {
    const b = selectedBin;
    if (!b) return [];
    const r = gtdbRanks(b.classification);
    return [
      { label: 'Quality', value: b.tier },
      { label: 'Completeness', value: b.completeness != null ? `${b.completeness}%` : null, cls: 'text-emerald-400' },
      { label: 'Contamination', value: b.contamination != null ? `${b.contamination}%` : null, cls: 'text-amber-400' },
      { label: 'Size', value: b.length ? fmtBp(b.length) : null },
      { label: 'Genes', value: b.n_genes ? fmtInt(b.n_genes) : null },
      { label: 'Read pairs', value: fmtInt(b.totalReads) },
      ...['phylum', 'class', 'order', 'family', 'genus', 'species'].map(k => ({ label: k[0].toUpperCase() + k.slice(1), value: r[k] })),
    ];
  });

  // ── Summary charts ──
  let boxData = $derived(names.length && shown.length ? names.map((s, j) => ({
    type: 'box', name: s,
    y: shown.filter(b => b.id !== 'unbinned').map(b => b.tpm[j]).filter(v => v > 0),
    marker: { color: '#22d3ee' }, line: { width: 1 }, boxpoints: false,
  })) : []);

  const columns = [
    { key: 'id', label: 'MAG' },
    { key: 'tier', label: 'Quality' },
    { key: 'completeness', label: 'Comp%' },
    { key: 'contamination', label: 'Cont%' },
    { key: 'length', label: 'Size', render: v => (v ? fmtBp(v) : '-') },
    { key: 'taxon', label: 'Taxon' },
    { key: 'meanTpm', label: 'Mean TPM', render: v => fmtTpm(v) },
    { key: 'maxTpm', label: 'Max TPM', render: v => fmtTpm(v) },
    { key: 'detected', label: 'Samples' },
    { key: 'totalReads', label: 'Read pairs', render: v => fmtInt(v) },
  ];
</script>

{#if $bins?.bins?.length}
  <div class="sticky top-14 z-10 bg-slate-950 pb-2 -mx-4 px-4 pt-1">
  <div class="flex items-center gap-3 flex-wrap text-xs">
    <span class="text-slate-400">Quality:</span>
    {#each tierDefs as { key, label }}
      {#if tierCounts[key]}
        <button
          class="px-3 py-1 rounded-md border transition-colors text-center
            {activeTiers.has(key) ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
          onclick={() => toggleTier(key)}
        >{label} <span class="text-slate-500">({tierCounts[key]})</span></button>
      {/if}
    {/each}
    <button
      class="px-3 py-1 rounded-md border transition-colors text-center
        {showUnbinned ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
      onclick={() => showUnbinned = !showUnbinned}
      title="Expression on contigs that are in no MAG"
    >Unbinned</button>
    <span class="text-slate-500">{filtered.length} bins</span>
    <span class="text-slate-600">|</span>
    <label class="text-slate-400 flex items-center gap-1">
      Comp &ge;
      <input type="range" min="0" max="100" step="5" bind:value={minCompleteness} class="w-28 h-1 accent-cyan-400" />
      <span class="text-cyan-400 font-mono w-8">{minCompleteness}%</span>
    </label>
    <label class="text-slate-400 flex items-center gap-1">
      Cont &le;
      <input type="range" min="0" max="100" step="1" bind:value={maxContamination} class="w-28 h-1 accent-cyan-400" />
      <span class="text-cyan-400 font-mono w-8">{maxContamination}%</span>
    </label>
    <span class="text-slate-600">|</span>
    <div class="relative flex items-center">
      <input type="text" bind:value={searchQuery} placeholder="MAG or taxon..."
        class="w-48 px-2 py-1 rounded bg-slate-800 border border-slate-600 text-slate-200 text-xs placeholder-slate-500 focus:border-cyan-400 focus:outline-none" />
      {#if searchQuery}
        <button class="absolute right-1 text-slate-500 hover:text-slate-300" onclick={() => { searchQuery = ''; }} title="Clear search">
          <svg class="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M18 6L6 18M6 6l12 12" /></svg>
        </button>
      {/if}
    </div>
    {#if searchDebounced}<span class="text-amber-400 font-medium">{filtered.length.toLocaleString()} matches</span>{/if}
  </div>
  </div>

  {#if shown.length}
  <div class="grid grid-cols-1 {selectedBin && !wide ? 'lg:grid-cols-3' : ''} gap-6 mb-6">
    <div class="{selectedBin && !wide ? 'lg:col-span-2' : ''} bg-slate-800 rounded-lg p-4 border border-slate-700 min-w-0">
      <div class="flex items-center justify-between mb-2 gap-2 flex-wrap">
        <h3 class="text-sm font-medium text-slate-400">
          MAG Expression ({shown.length}{shown.length < filtered.length ? ` of ${filtered.length}` : ''} bins x {names.length} samples)
        </h3>
        <div class="flex gap-2 text-xs">
          <button class="px-2 py-0.5 rounded border border-cyan-400 bg-cyan-400/10 text-cyan-400" onclick={cycleTop}
            title="Rows shown: Top 50 → Top 100 → Top 200 → All">{topN ? `Top ${topN}` : 'All'} &#x25BE;</button>
          <button class="px-2 py-0.5 rounded border border-cyan-400 bg-cyan-400/10 text-cyan-400" onclick={cycleSort}
            title={`Sort rows: ${sortDefs.map(d => d.label).join(' → ')}`}>Sort: {sortDefs.find(d => d.key === sortBy)?.label} &#x25BE;</button>
        </div>
      </div>
      <div class="overflow-auto" style="max-height: 700px;">
        <D3Heatmap data={heatmapData} onRowClick={id => selected = id} tooltipFormat={tooltip}
          {colorScale} legendLabel={['0', `${fmtTpm(maxTpm)} TPM`]} selectedRow={selected} {rowAnnotations} />
      </div>
      <div class="mt-2 text-xs text-slate-500">Sum of protein-coding gene TPM per MAG, log-scaled colour.{sortBy === 'clustered' ? ' Ward clustering of MAGs and samples on log(TPM+1).' : ''} Click a row for its detail.</div>
    </div>

    {#if selectedBin}
      <ExpressionDetail
        title={selectedBin.id}
        subtitle={selectedBin.id === 'unbinned' ? 'Contigs in no MAG' : selectedBin.classification ?? 'No GTDB-Tk classification'}
        facts={detailFacts} tpm={selectedBin.tpm} samples={names}
        links={$index?.mag_url ? [{ label: 'MAG Dashboard', href: `${$index.mag_url}#quality` }] : []}
        horizontal={wide} onclose={() => selected = null} />
    {/if}
  </div>

  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-6">
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart data={boxData} layout={{
        title: { text: 'MAG TPM Distribution per Sample (shown bins)', font: { size: 14 } },
        yaxis: { title: 'TPM', type: 'log' },
        xaxis: { showticklabels: names.length <= 60, tickangle: -60, tickfont: { size: 8 } },
        showlegend: false,
      }} />
    </div>
    <SampleTsne samples={names} matrix={shown.filter(b => b.id !== 'unbinned').map(b => b.tpm)} what="shown MAGs" />
  </div>
  {:else}
    <div class="flex flex-col items-center justify-center py-24 text-center">
      <p class="text-slate-400 text-lg mb-2">No MAGs match the current filters</p>
      <p class="text-slate-500 text-sm">Try adjusting the quality toggles, thresholds or search above.</p>
    </div>
  {/if}

  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
    <div class="flex items-center justify-between mb-2">
      <h3 class="text-sm font-medium text-slate-400">All MAGs ({sorted.length.toLocaleString()})</h3>
      <button
        class="text-xs px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-slate-200 hover:border-slate-500 transition-colors"
        onclick={() => downloadTsv(`mag_tpm_${sorted.length}.tsv`,
          ['mag', 'quality', 'completeness', 'contamination', 'length', 'classification', ...names],
          sorted.map(b => [b.id, b.tier, b.completeness, b.contamination, b.length, b.classification, ...b.tpm]))}
        title="Export all {sorted.length} MAGs with TPM per sample"
      >TSV + samples</button>
    </div>
    <DataTable {columns} rows={sorted} idKey="id" maxHeight="400px" selectedId={selected} onRowClick={row => selected = row.id} />
  </div>
{:else}
  <div class="flex flex-col items-center justify-center py-24 text-center">
    <p class="text-slate-400 text-lg mb-2">No MAG bins for this reference</p>
    <p class="text-slate-500 text-sm">Run with <code class="bg-slate-800 px-1 rounded">--from_mag</code> or give a <code class="bg-slate-800 px-1 rounded">&lt;name&gt;.bins.tsv</code> next to the reference.</p>
  </div>
{/if}
