<script>
  import ContigScatter from '../components/charts/ContigScatter.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { contigExplorer, loadContigExplorer } from '../stores/data.js';
  import { onMount } from 'svelte';

  let explorerData = $derived($contigExplorer);
  let sizeBy = $state('length');
  let sizeScale = $state(1.0);
  let mode = $state('tsne');

  // Color state: separate mode, binner source, taxonomy source + rank, metric, replicon
  let colorMode = $state('bins');   // 'bins' | 'taxa' | 'metric' | 'replicon'
  let binSource = $state('bin');
  let taxSource = $state('kaiju');
  let taxRank = $state('phylum');
  let metric = $state('depth');
  let replicon = $state('replicon');

  // Derived colorBy from the active mode + sub-selections
  let colorBy = $derived(
    colorMode === 'bins' ? binSource :
    colorMode === 'taxa' ? `${taxSource}_${taxRank}` :
    colorMode === 'replicon' ? replicon :
    metric
  );

  onMount(() => {
    loadContigExplorer();
  });

  // Cycle group definitions
  const binGroup =     { values: ['bin', 'semibin_bin', 'metabat_bin', 'maxbin_bin', 'lorbin_bin', 'comebin_bin'],
                         labels: ['DAS Tool', 'SemiBin2', 'MetaBAT2', 'MaxBin2', 'LorBin', 'COMEBin'] };
  const sourceGroup =  { values: ['kaiju', 'kraken2', 'rrna', 'genomad'],
                         labels: ['Kaiju', 'Kraken2', 'rRNA', 'geNomad'] };
  const rankGroup =    { values: ['domain', 'phylum', 'class', 'order', 'family', 'genus'],
                         labels: ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus'] };
  const metricGroup =  { values: ['depth', 'length', 'gc'],
                         labels: ['Depth', 'Length', 'GC%'] };
  const repliconGroup ={ values: ['replicon', 'tiara', 'whokaryote'],
                         labels: ['Replicon', 'Tiara', 'Whokaryote'] };
  const sizeGroup =    { values: ['length', 'depth', 'fixed'],
                         labels: ['Length', 'Depth', 'Fixed'] };
  const modeGroup =    { values: ['tsne', 'umap', 'pca'],
                         labels: ['t-SNE', 'UMAP', 'PCA'] };

  function cycle(values, current) {
    const idx = values.indexOf(current);
    return values[(idx + 1) % values.length];
  }

  function getLabel(values, labels, current) {
    const idx = values.indexOf(current);
    return idx >= 0 ? labels[idx] : labels[0];
  }

  // Fixed pixel widths per button to prevent layout shift when cycling
  const BW = { bin: '5.5rem', source: '5rem', rank: '4.5rem', metric: '4rem', replicon: '6rem', size: '4rem', mode: '4rem' };

  // Selection state
  let selectedIds = $state(null);  // null = no selection, [] would be empty lasso

  function handleSelection(ids) {
    selectedIds = ids;
  }

  const RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus'];

  function taxString(row, source) {
    return RANKS.map(r => row[`${source}_${r}`]).filter(Boolean).join('; ') || '-';
  }

  let tableColumns = $derived([
    { key: 'id', label: 'Contig' },
    { key: 'length', label: 'Length', render: (v) => v.toLocaleString() },
    { key: 'bin', label: 'Bin', render: (v) => v || '-' },
    { key: 'depth', label: 'Depth' },
    { key: 'gc', label: 'GC%' },
    { key: `${taxSource}_domain`, label: getLabel(sourceGroup.values, sourceGroup.labels, taxSource),
      render: (v, row) => taxString(row, taxSource) },
  ]);

  // Selected contigs lookup + table rows
  let selectedContigs = $derived.by(() => {
    if (!explorerData?.contigs || !selectedIds?.length) return null;
    const idSet = new Set(selectedIds);
    return explorerData.contigs.filter(c => idSet.has(c.id));
  });

  let selectionStats = $derived.by(() => {
    if (!selectedContigs) return null;
    const contigs = selectedContigs;
    const totalLen = contigs.reduce((s, c) => s + c.length, 0);
    const meanDepth = contigs.reduce((s, c) => s + c.depth, 0) / contigs.length;
    const meanGc = contigs.reduce((s, c) => s + (c.gc || 0), 0) / contigs.length;
    const bins = new Set(contigs.map(c => c.bin).filter(Boolean));
    const phyla = new Set(contigs.map(c => c.kaiju_phylum).filter(Boolean));
    return {
      count: contigs.length,
      totalLen,
      meanDepth: meanDepth.toFixed(1),
      meanGc: meanGc.toFixed(1),
      nBins: bins.size,
      nPhyla: phyla.size,
    };
  });

  let tableRows = $derived.by(() => {
    if (!explorerData?.contigs) return [];
    const source = selectedContigs || explorerData.contigs;
    return [...source]
      .sort((a, b) => b.length - a.length)
      .slice(0, 500);
  });
</script>

<div class="mb-4 flex items-center gap-3 flex-wrap text-xs">
  <span class="text-slate-400">Color:</span>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center
      {colorMode === 'bins' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
    style="min-width: {BW.bin}"
    onclick={() => { if (colorMode === 'bins') binSource = cycle(binGroup.values, binSource); else colorMode = 'bins'; }}
    title={`Click to cycle: ${binGroup.labels.join(' → ')}`}
  >
    {colorMode === 'bins' ? getLabel(binGroup.values, binGroup.labels, binSource) : 'Bins'} &#x25BE;
  </button>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center
      {colorMode === 'taxa' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
    style="min-width: {BW.source}"
    onclick={() => { if (colorMode === 'taxa') taxSource = cycle(sourceGroup.values, taxSource); else colorMode = 'taxa'; }}
    title={`Click to cycle: ${sourceGroup.labels.join(' → ')}`}
  >
    {colorMode === 'taxa' ? getLabel(sourceGroup.values, sourceGroup.labels, taxSource) : 'Source'} &#x25BE;
  </button>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center
      {colorMode === 'taxa' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
    style="min-width: {BW.rank}"
    onclick={() => { if (colorMode === 'taxa') taxRank = cycle(rankGroup.values, taxRank); else { colorMode = 'taxa'; } }}
    title={`Click to cycle: ${rankGroup.labels.join(' → ')}`}
  >
    {colorMode === 'taxa' ? getLabel(rankGroup.values, rankGroup.labels, taxRank) : 'Rank'} &#x25BE;
  </button>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center
      {colorMode === 'replicon' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
    style="min-width: {BW.replicon}"
    onclick={() => { if (colorMode === 'replicon') replicon = cycle(repliconGroup.values, replicon); else colorMode = 'replicon'; }}
    title={`Click to cycle: ${repliconGroup.labels.join(' → ')}`}
  >
    {colorMode === 'replicon' ? getLabel(repliconGroup.values, repliconGroup.labels, replicon) : 'Replicon'} &#x25BE;
  </button>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center
      {colorMode === 'metric' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
    style="min-width: {BW.metric}"
    onclick={() => { if (colorMode === 'metric') metric = cycle(metricGroup.values, metric); else colorMode = 'metric'; }}
    title={`Click to cycle: ${metricGroup.labels.join(' → ')}`}
  >
    {colorMode === 'metric' ? getLabel(metricGroup.values, metricGroup.labels, metric) : 'Metric'} &#x25BE;
  </button>
  <span class="text-slate-600">|</span>
  <span class="text-slate-400">Size:</span>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
    style="min-width: {BW.size}"
    onclick={() => sizeBy = cycle(sizeGroup.values, sizeBy)}
    title={`Click to cycle: ${sizeGroup.labels.join(' → ')}`}
  >
    {getLabel(sizeGroup.values, sizeGroup.labels, sizeBy)} &#x25BE;
  </button>
  <input
    type="range" min="0.2" max="3" step="0.1"
    bind:value={sizeScale}
    class="w-20 h-1 accent-cyan-400 cursor-pointer"
  />
  <span class="text-slate-500 w-8">{sizeScale.toFixed(1)}x</span>
  <span class="text-slate-600">|</span>
  <span class="text-slate-400">Projection:</span>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
    style="min-width: {BW.mode}"
    onclick={() => {
      let next = cycle(modeGroup.values, mode);
      if (next === 'tsne' && !explorerData?.has_tsne) next = cycle(modeGroup.values, next);
      if (next === 'umap' && !explorerData?.has_umap) next = cycle(modeGroup.values, next);
      mode = next;
    }}
    title={`Click to cycle: ${modeGroup.labels.join(' → ')}`}
  >
    {getLabel(modeGroup.values, modeGroup.labels, mode)} &#x25BE;
  </button>
</div>

{#if !explorerData}
  <div class="flex flex-col items-center justify-center py-20">
    <div class="w-10 h-10 border-4 border-slate-700 border-t-cyan-400 rounded-full animate-spin"></div>
    <p class="text-slate-400 text-sm mt-4">Loading contig data (2.4 MB)...</p>
  </div>
{:else}
  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700 mb-6">
    <ContigScatter data={explorerData} {colorBy} {sizeBy} {sizeScale} {mode} onselect={handleSelection} />
    <div class="text-xs text-slate-500 mt-2">
      {explorerData.contigs.length.toLocaleString()} contigs | WebGL rendering |
      Fourth-root transformed TNF |
      PCA variance explained: {explorerData.pca_variance_explained?.slice(0, 3).map(v => (v * 100).toFixed(1) + '%').join(', ')}
    </div>
  </div>

  {#if selectionStats}
    <div class="flex items-center gap-4 mb-4 text-xs flex-wrap">
      <span class="text-slate-400 font-medium">Selection:</span>
      <span class="text-cyan-400">{selectionStats.count.toLocaleString()} contigs</span>
      <span class="text-slate-400">{(selectionStats.totalLen / 1e6).toFixed(2)} Mb</span>
      <span class="text-slate-400">Mean depth: {selectionStats.meanDepth}x</span>
      <span class="text-slate-400">Mean GC: {selectionStats.meanGc}%</span>
      <span class="text-slate-400">{selectionStats.nBins} bins</span>
      <span class="text-slate-400">{selectionStats.nPhyla} phyla</span>
      <button
        class="text-slate-500 hover:text-slate-300 underline"
        onclick={() => selectedIds = null}
      >clear</button>
    </div>
  {/if}

  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
    <h3 class="text-sm font-medium text-slate-400 mb-2">
      {selectedContigs ? `${selectedContigs.length.toLocaleString()} selected contigs` : `Top 500 contigs by length`}
    </h3>
    <DataTable columns={tableColumns} rows={tableRows} idKey="id" maxHeight="350px" />
  </div>
{/if}
