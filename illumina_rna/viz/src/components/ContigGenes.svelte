<script>
  /**
   * Contig panel for the Explorer, after the MAG dashboard's ContigDetail:
   * header, taxonomy and MAG, gene arrow map, and a gene table with
   * expression bars. Clicking a gene (row or arrow) shows its TPM in every
   * sample; otherwise the contig's own profile is shown.
   */
  import GeneArrowMap from './charts/GeneArrowMap.svelte';
  import PlotlyChart from './charts/PlotlyChart.svelte';
  import { fmtInt, downloadTsv } from '../stores/data.js';

  let { contig = null, genes = null, loading = false, onclose = null, samples = [],
        sampleIdx = null, contigProfile = null, query = '' } = $props();

  let colorMode = $state('type');
  let colorModes = $derived([
    { id: 'type', label: 'Type' },
    { id: 'pooled', label: 'Pooled TPM' },
    ...(sampleIdx != null ? [{ id: 'sample', label: 'Sample TPM' }] : []),
  ]);
  const cycleColor = () => { const ids = colorModes.map(m => m.id); colorMode = ids[(ids.indexOf(colorMode) + 1) % ids.length]; };
  $effect(() => { if (!colorModes.some(m => m.id === colorMode)) colorMode = 'type'; });

  let selectedGene = $state(null);
  $effect(() => { contig; selectedGene = null; });

  // Arrow-map features: the 'dp' channel carries the chosen expression
  let features = $derived((genes ?? []).map(f => ({
    ...f,
    dp: colorMode === 'pooled' ? f.mean : colorMode === 'sample' ? f.tpm[sampleIdx] : null,
  })));

  let stats = $derived.by(() => {
    if (!genes?.length) return null;
    const byType = {};
    let named = 0, expressed = 0;
    for (const f of genes) {
      byType[f.t] = (byType[f.t] || 0) + 1;
      if (f.g || (f.p && f.p !== 'hypothetical protein')) named++;
      if (f.reads > 0) expressed++;
    }
    return { byType, named, expressed, total: genes.length };
  });

  let q = $derived(query.trim().toLowerCase());
  const matches = f => q && [f.id, f.g, f.p, f.ko].some(v => v && v.toLowerCase().includes(q));

  // Bars: log-scaled within this contig (TPM spans orders of magnitude)
  let maxPooled = $derived(Math.max(1e-9, ...(genes ?? []).map(f => f.mean)));
  let maxSample = $derived(sampleIdx == null ? 1e-9 : Math.max(1e-9, ...(genes ?? []).map(f => f.tpm[sampleIdx])));
  const barW = (v, mx) => v > 0 ? Math.max(2, Math.round(40 * Math.log1p(v) / Math.log1p(mx))) : 0;
  const fmt = v => v >= 100 ? Math.round(v).toLocaleString() : v >= 1 ? v.toFixed(1) : v > 0 ? v.toFixed(2) : '0';

  const TYPE_CLS = { cds: 'text-cyan-400', rrna: 'text-amber-400', trna: 'text-emerald-400', tmrna: 'text-violet-400', ncrna: 'text-pink-400' };

  function exportGenes() {
    downloadTsv(`${contig.id}_genes.tsv`,
      ['contig_id', 'locus_tag', 'type', 'gene', 'product', 'ko', 'start', 'end', 'strand', 'read_pairs', 'mean_tpm', ...samples],
      genes.map(f => [contig.id, f.id, f.t, f.g, f.p, f.ko, f.s, f.e, f.d === 1 ? '+' : '-', f.reads, f.mean.toFixed(3), ...f.tpm.map(v => v.toFixed(3))]));
  }

  let profile = $derived.by(() => {
    const y = selectedGene ? selectedGene.tpm : contigProfile;
    if (!y) return null;
    return [{
      type: 'bar', x: samples, y,
      marker: { color: samples.map((_, j) => j === sampleIdx ? '#fbbf24' : '#22d3ee'), opacity: 0.85 },
      hovertemplate: `%{x}<br>${selectedGene ? 'TPM' : 'RPKM'} %{y:.2f}<extra></extra>`,
    }];
  });
</script>

<div class="flex flex-col bg-slate-800 rounded-lg border border-slate-700 overflow-hidden h-[700px]">
  <div class="flex items-center justify-between px-4 py-2 border-b border-slate-700 flex-shrink-0">
    {#if contig}
      <div class="min-w-0">
        <h3 class="text-sm font-semibold text-slate-200 truncate">{contig.id}</h3>
        <div class="flex gap-3 text-xs text-slate-400 mt-0.5">
          <span>{contig.length.toLocaleString()} bp</span>
          <span>{fmtInt(contig.reads)} CDS read pairs</span>
          <span>Pooled RPKM: {fmt(contig.depth)}</span>
          <span>In {contig.n_detected}/{samples.length} samples</span>
        </div>
      </div>
    {:else}
      <div class="text-sm text-slate-500">No contig selected</div>
    {/if}
    {#if onclose}
      <button class="text-slate-500 hover:text-slate-300 p-1 -mr-1" onclick={onclose} title="Close panel">
        <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24"><path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M6 18L18 6M6 6l12 12" /></svg>
      </button>
    {/if}
  </div>

  {#if !contig}
    <div class="flex-1 flex items-center justify-center text-slate-500 text-sm p-8 text-center">Click a point in the scatter plot to view its genes</div>
  {:else if loading}
    <div class="flex-1 flex flex-col items-center justify-center">
      <div class="w-6 h-6 border-2 border-slate-700 border-t-cyan-400 rounded-full animate-spin"></div>
      <p class="text-slate-500 text-xs mt-2">Loading genes...</p>
    </div>
  {:else}
    <div class="px-4 py-2 border-b border-slate-700/50 flex-shrink-0 space-y-1 text-xs">
      <div>
        <span class="text-slate-500">Kaiju:</span>
        <span class="text-slate-200">{[contig.kaiju_phylum, contig.kaiju_genus].filter(Boolean).join('; ') || 'unclassified'}</span>
      </div>
      <div class="flex gap-x-4 flex-wrap">
        <span><span class="text-slate-500">MAG:</span> <span class="text-slate-300">{contig.bin || 'unbinned'}</span></span>
        {#if contig.rrna}<span class="text-amber-400" title="Reads here are mostly rRNA the filter missed; their CDS calls may be spurious">carries rRNA loci</span>{/if}
      </div>
    </div>

    {#if genes?.length}
      <div class="px-3 pt-3 flex-shrink-0">
        <GeneArrowMap {features} contigLength={contig.length} height={140} colorMode={colorMode === 'type' ? 'type' : 'depth'}
          dpLabel={colorMode === 'sample' ? `TPM in ${samples[sampleIdx]}` : 'pooled TPM'} dpUnit="TPM"
          onclick={f => selectedGene = genes.find(g => g.id === f.id) ?? null} />
      </div>

      {#if stats}
        <div class="flex items-center gap-3 px-4 py-1 text-xs text-slate-400 flex-wrap flex-shrink-0">
          {#each Object.entries(stats.byType) as [type, count]}<span>{count} {type.toUpperCase()}</span>{/each}
          <span class="text-slate-600">|</span>
          <span>{stats.named} / {stats.total} annotated</span>
          <span>{stats.expressed} expressed</span>
          <button
            class="ml-auto px-2 py-0.5 rounded border transition-colors {colorMode !== 'type' ? 'border-cyan-400 text-cyan-400' : 'border-slate-600 text-slate-400 hover:text-slate-200 hover:border-slate-500'}"
            onclick={cycleColor} title="Cycle arrow map coloring"
          >Color: {colorModes.find(m => m.id === colorMode)?.label}</button>
          <button class="px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-slate-200 hover:border-slate-500 transition-colors"
            onclick={exportGenes} title="Export genes with TPM per sample">TSV</button>
        </div>
      {/if}

      <div class="flex-1 overflow-y-auto min-h-0 px-3">
        <table class="w-full text-xs">
          <thead class="sticky top-0 bg-slate-800">
            <tr class="text-slate-500 border-b border-slate-700">
              <th class="text-left py-1 px-1 font-medium">ID</th>
              <th class="text-left py-1 px-1 font-medium">Type</th>
              <th class="text-left py-1 px-1 font-medium">Gene</th>
              <th class="text-left py-1 px-1 font-medium">Product</th>
              <th class="text-right py-1 px-1 font-medium">bp</th>
              <th class="text-center py-1 px-1 font-medium">Str</th>
              <th class="text-left py-1 px-1 font-medium" style="min-width: 80px" title="Mean TPM over all samples">Pooled</th>
              {#if sampleIdx != null}<th class="text-left py-1 px-1 font-medium" style="min-width: 80px" title="TPM in {samples[sampleIdx]}">Sample</th>{/if}
              <th class="text-right py-1 px-1 font-medium" title="Samples with reads">n</th>
            </tr>
          </thead>
          <tbody>
            {#each genes as f}
              <tr class="border-b border-slate-700/50 hover:bg-slate-700/30 cursor-pointer
                  {selectedGene?.id === f.id ? 'bg-cyan-400/10' : matches(f) ? 'bg-amber-400/10' : ''}"
                  onclick={() => selectedGene = selectedGene?.id === f.id ? null : f}>
                <td class="py-0.5 px-1 text-slate-400 font-mono text-[10px]">{f.id}</td>
                <td class="py-0.5 px-1 {TYPE_CLS[f.t] ?? 'text-slate-400'} font-mono text-[10px]">{f.t.toUpperCase()}</td>
                <td class="py-0.5 px-1 text-cyan-400 font-mono">{f.g || '-'}</td>
                <td class="py-0.5 px-1 text-slate-300 truncate max-w-[180px]" title={f.ko ? `${f.p} (${f.ko})` : f.p}>{f.p || 'hypothetical protein'}</td>
                <td class="py-0.5 px-1 text-right text-slate-400 font-mono">{(f.e - f.s + 1).toLocaleString()}</td>
                <td class="py-0.5 px-1 text-center text-slate-400">{f.d === 1 ? '+' : '-'}</td>
                <td class="py-0.5 px-1">
                  <div class="inline-flex items-center gap-1">
                    <span class="text-slate-400 font-mono text-[10px] w-9 text-right">{fmt(f.mean)}</span>
                    <div class="h-2.5 rounded-sm bg-cyan-400" style="width: {barW(f.mean, maxPooled)}px"></div>
                  </div>
                </td>
                {#if sampleIdx != null}
                  <td class="py-0.5 px-1">
                    <div class="inline-flex items-center gap-1">
                      <span class="text-slate-400 font-mono text-[10px] w-9 text-right">{fmt(f.tpm[sampleIdx])}</span>
                      <div class="h-2.5 rounded-sm bg-amber-400" style="width: {barW(f.tpm[sampleIdx], maxSample)}px"></div>
                    </div>
                  </td>
                {/if}
                <td class="py-0.5 px-1 text-right text-slate-400 font-mono">{f.detected}</td>
              </tr>
            {/each}
          </tbody>
        </table>
      </div>
    {:else}
      <div class="flex-1 flex items-center justify-center text-slate-600 text-xs italic">No predicted genes on this contig</div>
    {/if}

    <div class="border-t border-slate-700 px-3 pt-1 flex-shrink-0">
      {#if profile}
        <PlotlyChart minHeight="170px" data={profile} layout={{
          title: { text: selectedGene ? `${selectedGene.id}${selectedGene.g ? ' ' + selectedGene.g : ''}: TPM per sample` : 'Contig RPKM per sample (click a gene for its own)', font: { size: 11 } },
          height: 170, margin: { t: 26, r: 10, b: 20, l: 45 },
          xaxis: { showticklabels: false }, bargap: 0.15,
        }} />
      {:else}
        <p class="text-slate-500 text-xs py-6 text-center">Loading sample profile...</p>
      {/if}
    </div>
  {/if}
</div>
