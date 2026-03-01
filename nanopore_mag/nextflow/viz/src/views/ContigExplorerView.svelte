<script>
  import ContigScatter from '../components/charts/ContigScatter.svelte';
  import ReglScatter from '../components/charts/ReglScatter.svelte';
  import ContigDetail from '../components/ContigDetail.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { contigExplorer, loadContigExplorer, loadContigGenes, loadAllGenes, buildGeneSearchIndex, sampleDepths, loadSampleDepths } from '../stores/data.js';
  import { onMount } from 'svelte';

  let explorerData = $derived($contigExplorer);
  let sizeBy = $state('length');
  let sizeScale = $state(1.0);
  let mode = $state('tsne');
  let renderer = $state('plotly');  // 'plotly' | 'regl'

  // Color state: separate mode, binner source, taxonomy source + rank, metric, replicon
  let colorMode = $state('metric');   // 'bins' | 'taxa' | 'metric' | 'replicon'
  let binSource = $state('bin');
  let taxSource = $state('kaiju');
  let taxRank = $state('phylum');
  let metric = $state('gc');
  let replicon = $state('replicon');

  // Per-sample depth state
  let selectedSample = $state('');
  let sampleDepthData = $derived($sampleDepths);
  let sampleNames = $derived(sampleDepthData?.samples?.length ? sampleDepthData.samples : []);
  let hasSamples = $derived(sampleNames.length > 0);

  // Per-sample nonzero contig counts for dropdown labels
  let sampleCounts = $derived.by(() => {
    if (!sampleDepthData?.depths || !sampleNames.length) return {};
    const depths = sampleDepthData.depths;
    const keys = Object.keys(depths);
    const counts = {};
    for (let si = 0; si < sampleNames.length; si++) {
      let nz = 0;
      for (const k of keys) {
        if (depths[k][si] > 0) nz++;
      }
      counts[sampleNames[si]] = nz;
    }
    return counts;
  });

  // Trigger lazy load of per-sample depths when user enters depth metric mode
  $effect(() => {
    if (colorMode === 'metric' && metric === 'depth') {
      loadSampleDepths();
    }
  });

  // Derived colorBy from the active mode + sub-selections
  let colorBy = $derived(
    colorMode === 'bins' ? binSource :
    colorMode === 'taxa' ? `${taxSource}_${taxRank}` :
    colorMode === 'replicon' ? replicon :
    (metric === 'depth' && selectedSample) ? 'sample_depth' : metric
  );

  onMount(() => {
    loadContigExplorer();
  });

  // Cycle group definitions
  const binGroup =     { values: ['bin', 'semibin_bin', 'metabat_bin', 'maxbin_bin', 'lorbin_bin', 'comebin_bin'],
                         labels: ['DAS Tool', 'SemiBin2', 'MetaBAT2', 'MaxBin2', 'LorBin', 'COMEBin'] };
  const sourceGroup =  { values: ['kaiju', 'kraken2', 'sendsketch', 'rrna', 'genomad'],
                         labels: ['Kaiju', 'Kraken2', 'BBSketch', 'rRNA', 'geNomad'] };
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

  // Auto-correct mode if the initial choice isn't available
  $effect(() => {
    if (!explorerData) return;
    if (mode === 'tsne' && !explorerData.has_tsne && explorerData.has_umap) mode = 'umap';
    else if (mode === 'umap' && !explorerData.has_umap && explorerData.has_tsne) mode = 'tsne';
    else if (mode !== 'pca' && !explorerData.has_tsne && !explorerData.has_umap) mode = 'pca';
  });

  const rendererGroup = { values: ['regl', 'plotly'], labels: ['WebGL', 'Plotly'] };

  // Fixed pixel widths per button to prevent layout shift when cycling
  const BW = { bin: '5.5rem', source: '5rem', rank: '4.5rem', metric: '4rem', replicon: '6rem', size: '4rem', mode: '4rem', renderer: '4.5rem' };

  // Stable color map: always derived from FULL dataset so filtering doesn't shift colors
  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8',
                   '#94a3b8','#d4d4d8','#78716c'];

  let stableColorMap = $derived.by(() => {
    if (!explorerData?.contigs) return {};
    const isCont = ['depth', 'length', 'gc', 'sample_depth'].includes(colorBy);
    if (isCont) return {};
    const isBin = colorBy === 'bin' || colorBy.endsWith('_bin');
    const bgLabel = isBin ? 'unbinned' : 'Unknown';
    const names = new Set();
    for (const c of explorerData.contigs) names.add(c[colorBy] || bgLabel);
    const sorted = [...names].filter(n => n !== bgLabel).sort();
    const map = { [bgLabel]: '#475569' };
    for (let i = 0; i < sorted.length; i++) map[sorted[i]] = PALETTE[i % PALETTE.length];
    return map;
  });

  // Stable coordinate extents: full-dataset bounds so filtering doesn't re-zoom
  let coordExtents = $derived.by(() => {
    if (!explorerData?.contigs?.length) return null;
    const contigs = explorerData.contigs;
    const modes = [['pca', 'pca_x', 'pca_y'], ['tsne', 'tsne_x', 'tsne_y'], ['umap', 'umap_x', 'umap_y']];
    const result = {};
    for (const [name, xKey, yKey] of modes) {
      let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
      for (const c of contigs) {
        const x = c[xKey] || 0, y = c[yKey] || 0;
        if (x < xMin) xMin = x; if (x > xMax) xMax = x;
        if (y < yMin) yMin = y; if (y > yMax) yMax = y;
      }
      result[name] = { xMin, xMax, yMin, yMax };
    }
    return result;
  });

  // Stable size extents: pre-transformed min/max from full dataset so sizes don't jump on filter
  let fullSizeRange = $derived.by(() => {
    if (!explorerData?.contigs) return null;
    let sqrtDepMin = Infinity, sqrtDepMax = -Infinity;
    let powLenMin = Infinity, powLenMax = -Infinity;
    for (const c of explorerData.contigs) {
      const sv = Math.sqrt(c.depth + 0.01);
      if (sv < sqrtDepMin) sqrtDepMin = sv;
      if (sv > sqrtDepMax) sqrtDepMax = sv;
      const pv = Math.pow(c.length, 0.3);
      if (pv < powLenMin) powLenMin = pv;
      if (pv > powLenMax) powLenMax = pv;
    }
    return { depth: [sqrtDepMin, sqrtDepMax], length: [powLenMin, powLenMax] };
  });

  // ---- Filter state (log-scale dual-range sliders) ----
  let lenMinPct = $state(0);
  let lenMaxPct = $state(100);
  let depMinPct = $state(0);
  let depMaxPct = $state(100);

  let dataExtents = $derived.by(() => {
    if (!explorerData?.contigs?.length) return null;
    let lMin = Infinity, lMax = -Infinity, dMin = Infinity, dMax = -Infinity;
    for (const c of explorerData.contigs) {
      if (c.length < lMin) lMin = c.length;
      if (c.length > lMax) lMax = c.length;
      if (c.depth < dMin) dMin = c.depth;
      if (c.depth > dMax) dMax = c.depth;
    }
    return {
      lMin: Math.max(1, lMin), lMax: Math.max(2, lMax),
      dMin: Math.max(0.001, dMin), dMax: Math.max(0.01, dMax),
    };
  });

  function pctToLog(pct, logMin, logMax) {
    return Math.pow(10, logMin + (logMax - logMin) * pct / 100);
  }

  let lenRange = $derived.by(() => {
    if (!dataExtents) return [0, Infinity];
    const lo = Math.log10(dataExtents.lMin), hi = Math.log10(dataExtents.lMax);
    // Use exact extent values at boundaries to avoid log↔pow floating-point drift
    return [
      lenMinPct <= 0 ? dataExtents.lMin : pctToLog(lenMinPct, lo, hi),
      lenMaxPct >= 100 ? dataExtents.lMax : pctToLog(lenMaxPct, lo, hi),
    ];
  });
  let depRange = $derived.by(() => {
    if (!dataExtents) return [0, Infinity];
    const lo = Math.log10(dataExtents.dMin), hi = Math.log10(dataExtents.dMax);
    return [
      depMinPct <= 0 ? dataExtents.dMin : pctToLog(depMinPct, lo, hi),
      depMaxPct >= 100 ? dataExtents.dMax : pctToLog(depMaxPct, lo, hi),
    ];
  });

  let isFiltered = $derived(lenMinPct > 0 || lenMaxPct < 100 || depMinPct > 0 || depMaxPct < 100);

  let filteredData = $derived.by(() => {
    if (!explorerData?.contigs) return null;
    if (!isFiltered) return explorerData;
    const [lMin, lMax] = lenRange;
    const [dMin, dMax] = depRange;
    const filtered = explorerData.contigs.filter(c =>
      c.length >= lMin && c.length <= lMax &&
      c.depth >= dMin && c.depth <= dMax
    );
    return { ...explorerData, contigs: filtered };
  });

  function fmtLen(v) {
    if (v >= 1e6) return (v / 1e6).toFixed(1) + ' Mb';
    if (v >= 1e3) return (v / 1e3).toFixed(1) + ' Kb';
    return Math.round(v) + ' bp';
  }
  function fmtDep(v) {
    if (v >= 100) return Math.round(v) + 'x';
    if (v >= 10) return v.toFixed(1) + 'x';
    if (v >= 1) return v.toFixed(2) + 'x';
    return v.toFixed(3) + 'x';
  }

  function resetFilters() {
    lenMinPct = 0; lenMaxPct = 100;
    depMinPct = 0; depMaxPct = 100;
    searchQuery = '';
  }

  // ---- Search state ----
  let searchQuery = $state('');
  let searchDebounced = $state('');
  let geneIndex = $state(null);
  let geneIndexLoading = $state(false);
  let debounceTimer;

  // Debounce search input (300ms)
  $effect(() => {
    const q = searchQuery;
    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(() => { searchDebounced = q; }, 300);
    return () => clearTimeout(debounceTimer);
  });

  // Trigger gene loading on first keystroke
  $effect(() => {
    if (searchQuery.length > 0 && !geneIndex && !geneIndexLoading) {
      geneIndexLoading = true;
      loadAllGenes().then(allGenes => {
        geneIndex = buildGeneSearchIndex(allGenes);
        geneIndexLoading = false;
      });
    }
  });

  // Search fields on each contig to check (contig ID + 24 taxonomy + 6 bin + replicon/tiara/whokaryote)
  const SEARCH_FIELDS = [
    'id',
    ...['kaiju', 'kraken2', 'sendsketch', 'rrna', 'genomad'].flatMap(s =>
      ['domain', 'phylum', 'class', 'order', 'family', 'genus'].map(r => `${s}_${r}`)
    ),
    'bin', 'semibin_bin', 'metabat_bin', 'maxbin_bin', 'lorbin_bin', 'comebin_bin',
    'replicon', 'tiara', 'whokaryote',
  ];

  // Compute matching contig IDs
  let searchMatchIds = $derived.by(() => {
    const q = searchDebounced.trim().toLowerCase();
    if (!q || !filteredData?.contigs) return null;

    const matched = new Set();
    const contigs = filteredData.contigs;

    // Phase 1: search contig fields
    for (const c of contigs) {
      for (const field of SEARCH_FIELDS) {
        const val = c[field];
        if (val && String(val).toLowerCase().includes(q)) {
          matched.add(c.id);
          break; // short-circuit on first match
        }
      }
    }

    // Phase 2: search gene index (for contigs not already matched)
    if (geneIndex) {
      for (const c of contigs) {
        if (matched.has(c.id)) continue;
        const geneStr = geneIndex.get(c.id);
        if (geneStr && geneStr.includes(q)) {
          matched.add(c.id);
        }
      }
    }

    return matched;
  });

  // Selection state
  let selectedIds = $state(null);  // null = no selection, [] would be empty lasso

  function handleSelection(ids) {
    selectedIds = ids;
  }

  // ---- Detail panel state ----
  let detailContigId = $state(null);
  let detailGenes = $state(null);
  let detailLoading = $state(false);

  let detailContig = $derived.by(() => {
    if (!detailContigId || !explorerData?.contigs) return null;
    return explorerData.contigs.find(c => c.id === detailContigId) || null;
  });

  function handleContigClick(contigId) {
    if (detailContigId === contigId) return; // already showing
    detailContigId = contigId;
    detailGenes = null;
    detailLoading = true;
    loadContigGenes(contigId).then(genes => {
      // Only update if still showing this contig
      if (detailContigId === contigId) {
        detailGenes = genes;
        detailLoading = false;
      }
    });
  }

  // Active bin info for tooltips and detail panel (null when not coloring by binner)
  let activeBin = $derived.by(() => {
    if (colorMode !== 'bins') return null;
    const idx = binGroup.values.indexOf(binSource);
    return { field: binSource, label: binGroup.labels[idx] || binSource };
  });

  // Auto-open detail panel with the longest contig on first load
  let autoOpened = false;
  $effect(() => {
    if (explorerData?.contigs?.length && !autoOpened) {
      autoOpened = true;
      const longest = explorerData.contigs.reduce((a, b) => a.length > b.length ? a : b);
      handleContigClick(longest.id);
    }
  });

  const RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus'];

  function taxString(row, source) {
    return RANKS.map(r => row[`${source}_${r}`]).filter(Boolean).join('; ') || '-';
  }

  let tableColumns = $derived([
    { key: 'id', label: 'Contig' },
    { key: 'length', label: 'Length', render: (v) => v.toLocaleString() },
    ...(activeBin ? [{ key: activeBin.field, label: activeBin.label, render: (v) => v || '-' }] : []),
    { key: 'depth', label: 'Depth' },
    { key: 'gc', label: 'GC%' },
    { key: `${taxSource}_domain`, label: getLabel(sourceGroup.values, sourceGroup.labels, taxSource),
      render: (v, row) => taxString(row, taxSource) },
  ]);

  // Selected contigs lookup + table rows (use filteredData)
  let selectedContigs = $derived.by(() => {
    if (!filteredData?.contigs || !selectedIds?.length) return null;
    const idSet = new Set(selectedIds);
    return filteredData.contigs.filter(c => idSet.has(c.id));
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

  // Full sorted source (for export — no 500 limit)
  // Priority: lasso selection > search matches > all filtered contigs
  let tableSrc = $derived.by(() => {
    if (!filteredData?.contigs) return [];
    let source;
    if (selectedContigs) {
      source = selectedContigs;
    } else if (searchMatchIds) {
      source = filteredData.contigs.filter(c => searchMatchIds.has(c.id));
    } else {
      source = filteredData.contigs;
    }
    return [...source].sort((a, b) => b.length - a.length);
  });

  let tableRows = $derived(tableSrc.slice(0, 500));

  function exportTSV() {
    if (!tableSrc.length) return;
    const binField = activeBin?.field || 'bin';
    const binLabel = activeBin?.label || 'DAS Tool';
    const header = ['contig_id', 'length', binLabel, 'depth', 'gc', `${taxSource}_taxonomy`].join('\t');
    const lines = tableSrc.map(row =>
      [row.id, row.length, row[binField] || '', row.depth, row.gc ?? '', taxString(row, taxSource)].join('\t')
    );
    const blob = new Blob([header + '\n' + lines.join('\n') + '\n'], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `contigs_${selectedContigs ? 'selection' : 'all'}_${tableSrc.length}.tsv`;
    a.click();
    URL.revokeObjectURL(url);
  }
</script>

<div class="sticky top-14 z-10 bg-slate-950 pb-2 -mx-4 px-4 pt-1">
<div class="mb-2 flex items-center gap-3 flex-wrap text-xs">
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
  {#if colorMode === 'metric' && metric === 'depth' && hasSamples}
    <select
      class="px-2 py-1 rounded-md border border-slate-600 bg-slate-800 text-slate-300 text-xs focus:border-cyan-400 focus:outline-none cursor-pointer"
      bind:value={selectedSample}
    >
      <option value="">Total depth</option>
      {#each sampleNames as name}
        <option value={name}>{name} ({(sampleCounts[name] || 0).toLocaleString()})</option>
      {/each}
    </select>
  {/if}
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
  <span class="text-slate-600">|</span>
  <span class="text-slate-400">Renderer:</span>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
    style="min-width: {BW.renderer}"
    onclick={() => renderer = cycle(rendererGroup.values, renderer)}
    title={`Click to cycle: ${rendererGroup.labels.join(' → ')}`}
  >
    {getLabel(rendererGroup.values, rendererGroup.labels, renderer)} &#x25BE;
  </button>
</div>

{#if explorerData}
  <div class="flex items-center gap-3 flex-wrap text-xs">
    <span class="text-slate-400">Filter:</span>

    <span class="text-slate-500">Length</span>
    <span class="text-slate-500 w-14 text-right font-mono">{fmtLen(lenRange[0])}</span>
    <div class="dual-range relative w-36 h-5 flex items-center">
      <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
      <div class="absolute h-1 bg-cyan-400/40 rounded" style="left: {lenMinPct}%; right: {100 - lenMaxPct}%"></div>
      <input type="range" min="0" max="100" step="1"
        bind:value={lenMinPct}
        oninput={() => { if (lenMinPct > lenMaxPct - 2) lenMinPct = lenMaxPct - 2; }}
      />
      <input type="range" min="0" max="100" step="1"
        bind:value={lenMaxPct}
        oninput={() => { if (lenMaxPct < lenMinPct + 2) lenMaxPct = lenMinPct + 2; }}
      />
    </div>
    <span class="text-slate-500 w-14 font-mono">{fmtLen(lenRange[1])}</span>

    <span class="text-slate-600">|</span>

    <span class="text-slate-500">Depth</span>
    <span class="text-slate-500 w-14 text-right font-mono">{fmtDep(depRange[0])}</span>
    <div class="dual-range relative w-36 h-5 flex items-center">
      <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
      <div class="absolute h-1 bg-cyan-400/40 rounded" style="left: {depMinPct}%; right: {100 - depMaxPct}%"></div>
      <input type="range" min="0" max="100" step="1"
        bind:value={depMinPct}
        oninput={() => { if (depMinPct > depMaxPct - 2) depMinPct = depMaxPct - 2; }}
      />
      <input type="range" min="0" max="100" step="1"
        bind:value={depMaxPct}
        oninput={() => { if (depMaxPct < depMinPct + 2) depMaxPct = depMinPct + 2; }}
      />
    </div>
    <span class="text-slate-500 w-14 font-mono">{fmtDep(depRange[1])}</span>

    <span class="text-slate-600">|</span>

    <div class="relative flex items-center">
      <input
        type="text"
        bind:value={searchQuery}
        placeholder="contig, taxon, bin, gene..."
        class="w-48 px-2 py-1 rounded bg-slate-800 border border-slate-600 text-slate-200 text-xs placeholder-slate-500 focus:border-cyan-400 focus:outline-none"
      />
      {#if searchQuery}
        <button
          class="absolute right-1 text-slate-500 hover:text-slate-300"
          onclick={() => { searchQuery = ''; }}
          title="Clear search"
        >
          <svg class="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
            <path d="M18 6L6 18M6 6l12 12" />
          </svg>
        </button>
      {/if}
    </div>
    {#if geneIndexLoading}
      <span class="text-slate-500 text-xs italic">loading genes...</span>
    {/if}
    {#if searchMatchIds}
      <span class="text-amber-400 font-medium">{searchMatchIds.size.toLocaleString()} matches</span>
    {/if}

    {#if isFiltered || searchMatchIds}
      {#if isFiltered}
        <span class="text-cyan-400 font-medium">{filteredData.contigs.length.toLocaleString()} / {explorerData.contigs.length.toLocaleString()}</span>
      {/if}
      <button class="text-slate-500 hover:text-slate-300 underline" onclick={resetFilters}>reset</button>
    {/if}
  </div>
{/if}
</div>

{#if !explorerData}
  <div class="flex flex-col items-center justify-center py-20">
    <div class="w-10 h-10 border-4 border-slate-700 border-t-cyan-400 rounded-full animate-spin"></div>
    <p class="text-slate-400 text-sm mt-4">Loading contig data...</p>
  </div>
{:else}
  <div class="flex gap-4 mb-6">
    <!-- Scatter plot -->
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700 h-[700px] flex flex-col flex-1 min-w-0">
      {#if renderer === 'regl'}
        <ReglScatter data={filteredData} {colorBy} {sizeBy} {sizeScale} {mode} colorMap={stableColorMap} sizeRange={fullSizeRange} {coordExtents} onselect={handleSelection} onclick={handleContigClick} {searchMatchIds} {sampleDepthData} {selectedSample} />
      {:else}
        <ContigScatter data={filteredData} {colorBy} {sizeBy} {sizeScale} {mode} colorMap={stableColorMap} {coordExtents} onselect={handleSelection} onclick={handleContigClick} {searchMatchIds} {sampleDepthData} {selectedSample} />
      {/if}
      <div class="text-xs text-slate-500 mt-2 flex-shrink-0">
        {filteredData.contigs.length.toLocaleString()} contigs |
        {renderer === 'regl' ? 'regl-scatterplot' : 'Plotly scattergl'} |
        Fourth-root transformed TNF |
        {renderer === 'regl' ? 'Shift+drag to lasso | Click point for details' : 'Lasso via toolbar | Click point for details'}
      </div>
    </div>

    <!-- Detail panel (right side, always visible) -->
    <div class="w-1/2 min-w-[400px] max-w-[600px]">
      <ContigDetail
        contig={detailContig}
        genes={detailGenes}
        loading={detailLoading}
        {activeBin}
        {taxSource}
      />
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
    <div class="flex items-center justify-between mb-2">
      <h3 class="text-sm font-medium text-slate-400">
        {selectedContigs ? `${selectedContigs.length.toLocaleString()} selected contigs` : searchMatchIds ? `${searchMatchIds.size.toLocaleString()} matching contigs` : `Top 500 contigs by length`}
        {#if tableSrc.length > 500}
          <span class="text-slate-500 text-xs">({tableSrc.length.toLocaleString()} total)</span>
        {/if}
      </h3>
      <button
        class="text-xs px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-slate-200 hover:border-slate-500 transition-colors"
        onclick={exportTSV}
        title="Export all {tableSrc.length.toLocaleString()} contigs as TSV"
      >TSV</button>
    </div>
    <DataTable columns={tableColumns} rows={tableRows} idKey="id" maxHeight="350px" onRowClick={(row) => handleContigClick(row.id)} selectedId={detailContigId} hideExport={true} />
  </div>
{/if}

<style>
  .dual-range input[type="range"] {
    -webkit-appearance: none;
    appearance: none;
    background: transparent;
    pointer-events: none;
    position: absolute;
    width: 100%;
    height: 100%;
    margin: 0;
    padding: 0;
  }
  .dual-range input[type="range"]::-webkit-slider-thumb {
    -webkit-appearance: none;
    pointer-events: all;
    height: 14px;
    width: 14px;
    border-radius: 50%;
    background: #22d3ee;
    cursor: pointer;
    border: 2px solid #0f172a;
    box-shadow: 0 0 3px rgba(0,0,0,0.4);
  }
  .dual-range input[type="range"]::-moz-range-thumb {
    pointer-events: all;
    height: 14px;
    width: 14px;
    border-radius: 50%;
    background: #22d3ee;
    cursor: pointer;
    border: 2px solid #0f172a;
    box-shadow: 0 0 3px rgba(0,0,0,0.4);
  }
  .dual-range input[type="range"]::-webkit-slider-runnable-track {
    height: 0;
  }
  .dual-range input[type="range"]::-moz-range-track {
    height: 0;
    background: transparent;
  }
</style>
