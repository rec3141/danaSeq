<script>
  import GeneArrowMap from './charts/GeneArrowMap.svelte';

  let { contig = null, genes = null, loading = false, onclose = null, activeBin = null, taxSource = 'kaiju' } = $props();

  const RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus'];
  const SOURCE_LABELS = { kaiju: 'Kaiju', kraken2: 'Kraken2', rrna: 'rRNA', sendsketch: 'BBSketch' };
  const BIN_FIELDS = ['bin', 'semibin_bin', 'metabat_bin', 'maxbin_bin', 'lorbin_bin', 'comebin_bin'];

  let fullTaxString = $derived.by(() => {
    if (!contig) return '';
    return RANKS.map(r => contig[`${taxSource}_${r}`]).filter(Boolean).join('; ') || 'unclassified';
  });

  // All bins this contig belongs to (deduplicated names)
  let binList = $derived.by(() => {
    if (!contig) return [];
    const bins = [];
    for (const field of BIN_FIELDS) {
      const val = contig[field];
      if (val && !bins.includes(val)) bins.push(val);
    }
    return bins;
  });

  // Feature summary stats
  let stats = $derived.by(() => {
    if (!genes?.length) return null;
    const byType = {};
    let named = 0;
    for (const f of genes) {
      byType[f.t] = (byType[f.t] || 0) + 1;
      if (f.g || (f.p && f.p !== 'hypothetical protein')) named++;
    }
    return { byType, named, total: genes.length };
  });

  const TYPE_ORDER = { cds: 0, rrna: 1, trna: 2, tmrna: 3, ncrna: 4, crispr: 5 };

  let geneTable = $derived.by(() => {
    if (!genes?.length) return [];
    return [...genes].sort((a, b) => a.s - b.s);
  });

  function exportGenesTSV() {
    if (!contig || !geneTable.length) return;
    const header = ['contig_id', 'locus_tag', 'type', 'gene', 'product', 'start', 'end', 'length_bp', 'strand'].join('\t');
    const lines = geneTable.map(f =>
      [contig.id, f.id || '', f.t, f.g || '', f.p || (f.t === 'cds' ? 'hypothetical protein' : ''), f.s, f.e, f.e - f.s + 1, f.d === 1 ? '+' : '-'].join('\t')
    );
    const blob = new Blob([header + '\n' + lines.join('\n') + '\n'], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${contig.id}_genes.tsv`;
    a.click();
    URL.revokeObjectURL(url);
  }
</script>

<div class="flex flex-col bg-slate-800 rounded-lg border border-slate-700 overflow-hidden h-[700px]">
  <!-- Header -->
  <div class="flex items-center justify-between px-4 py-2 border-b border-slate-700 flex-shrink-0">
    {#if contig}
      <div class="min-w-0">
        <h3 class="text-sm font-semibold text-slate-200 truncate">{contig.id}</h3>
        <div class="flex gap-3 text-xs text-slate-400 mt-0.5">
          <span>{contig.length?.toLocaleString()} bp</span>
          <span>Depth: {contig.depth}</span>
          <span>GC: {contig.gc ?? '?'}%</span>
        </div>
      </div>
    {:else}
      <div class="text-sm text-slate-500">No contig selected</div>
    {/if}
    {#if onclose}
      <button
        class="text-slate-500 hover:text-slate-300 p-1 -mr-1"
        onclick={onclose}
        title="Close panel"
      >
        <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
          <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M6 18L18 6M6 6l12 12" />
        </svg>
      </button>
    {/if}
  </div>

  {#if !contig}
    <div class="flex-1 flex items-center justify-center text-slate-500 text-sm p-8 text-center">
      Click a point in the scatter plot to view contig details
    </div>
  {:else if loading}
    <div class="flex-1 flex flex-col items-center justify-center">
      <div class="w-6 h-6 border-2 border-slate-700 border-t-cyan-400 rounded-full animate-spin"></div>
      <p class="text-slate-500 text-xs mt-2">Loading gene features...</p>
    </div>
  {:else}
    <!-- Taxonomy + bins + classification (always shown) -->
    <div class="px-4 py-2 border-b border-slate-700/50 flex-shrink-0 space-y-1">
      <!-- Full taxonomy for selected classifier -->
      <div class="text-xs">
        <span class="text-slate-500">{SOURCE_LABELS[taxSource] || taxSource}:</span>
        <span class="text-slate-200">{fullTaxString}</span>
      </div>
      <!-- Other classifiers (phylum only) + classification tags -->
      <div class="flex gap-x-4 gap-y-0.5 text-xs flex-wrap">
        {#each ['kaiju', 'kraken2', 'rrna', 'sendsketch'] as src}
          {#if src !== taxSource}
            {@const phylum = contig[`${src}_phylum`]}
            {#if phylum}
              <span><span class="text-slate-500">{SOURCE_LABELS[src] || src}:</span> <span class="text-slate-400">{phylum}</span></span>
            {/if}
          {/if}
        {/each}
        {#if contig.replicon}
          <span><span class="text-slate-500">Replicon:</span> <span class="text-slate-400">{contig.replicon}</span></span>
        {/if}
        {#if contig.tiara}
          <span><span class="text-slate-500">Tiara:</span> <span class="text-slate-400">{contig.tiara}</span></span>
        {/if}
      </div>
      <!-- Bins (simple list) -->
      {#if binList.length}
        <div class="text-xs">
          <span class="text-slate-500">Bins:</span>
          <span class="text-slate-300">{binList.join(', ')}</span>
        </div>
      {/if}
    </div>

    {#if genes?.length}
      <!-- Gene arrow map (fixed, not scrolled) -->
      <div class="px-3 pt-3 flex-shrink-0">
        <GeneArrowMap features={genes} contigLength={contig.length} height={160} />
      </div>

      <!-- Stats bar -->
      {#if stats}
        <div class="flex items-center gap-3 px-4 py-1 text-xs text-slate-400 flex-wrap flex-shrink-0">
          {#each Object.entries(stats.byType) as [type, count]}
            <span>{count} {type.toUpperCase()}</span>
          {/each}
          <span class="text-slate-600">|</span>
          <span>{stats.named} / {stats.total} annotated</span>
          <button
            class="ml-auto px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-slate-200 hover:border-slate-500 transition-colors"
            onclick={exportGenesTSV}
            title="Export gene annotations as TSV"
          >TSV</button>
        </div>
      {/if}

      <!-- Gene table (scrollable) -->
      <div class="flex-1 overflow-y-auto min-h-0 px-3 pb-3">
        <table class="w-full text-xs">
          <thead class="sticky top-0 bg-slate-800">
            <tr class="text-slate-500 border-b border-slate-700">
              <th class="text-left py-1 px-1 font-medium">ID</th>
              <th class="text-left py-1 px-1 font-medium">Type</th>
              <th class="text-left py-1 px-1 font-medium">Gene</th>
              <th class="text-left py-1 px-1 font-medium">Product</th>
              <th class="text-right py-1 px-1 font-medium">bp</th>
              <th class="text-center py-1 px-1 font-medium">Str</th>
            </tr>
          </thead>
          <tbody>
            {#each geneTable as f}
              {@const typeColor = f.t === 'cds' ? 'text-cyan-400' : f.t === 'rrna' ? 'text-amber-400' : f.t === 'trna' ? 'text-emerald-400' : f.t === 'tmrna' ? 'text-violet-400' : f.t === 'ncrna' ? 'text-pink-400' : f.t === 'crispr' ? 'text-red-400' : 'text-slate-400'}
              <tr class="border-b border-slate-700/50 hover:bg-slate-700/30">
                <td class="py-0.5 px-1 text-slate-400 font-mono text-[10px]">{f.id || '-'}</td>
                <td class="py-0.5 px-1 {typeColor} font-mono text-[10px]">{f.t.toUpperCase()}</td>
                <td class="py-0.5 px-1 text-cyan-400 font-mono">{f.g || '-'}</td>
                <td class="py-0.5 px-1 text-slate-300 truncate max-w-[200px]" title={f.p}>{f.p || (f.t === 'cds' ? 'hypothetical protein' : '')}</td>
                <td class="py-0.5 px-1 text-right text-slate-400 font-mono">{(f.e - f.s + 1).toLocaleString()}</td>
                <td class="py-0.5 px-1 text-center text-slate-400">{f.d === 1 ? '+' : '-'}</td>
              </tr>
            {/each}
          </tbody>
        </table>
      </div>
    {:else}
      <!-- No genes for this contig -->
      <div class="flex-1 flex items-center justify-center text-slate-600 text-xs italic">
        No predicted genes on this contig
      </div>
    {/if}
  {/if}
</div>
