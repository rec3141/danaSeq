<script>
  import ReglScatter from '../components/charts/ReglScatter.svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { onMount } from 'svelte';
  import ContigGenes from '../components/ContigGenes.svelte';
  import { contigMap, contigMapState, loadContigMap, loadContigLayer, refName,
           loadContigGenes, ensureGeneSearch, searchGenes, geneSearchState,
           fmtBp, fmtInt, downloadTsv } from '../stores/data.js';

  // Contig map from the MAG analysis' t-SNE, with RNA expression laid over it.
  // Same controls as the MAG dashboard's Explorer; ReglScatter's "depth"
  // channels carry expression (total RPKM, or the selected sample's RPKM).
  onMount(() => { loadContigMap(); });
  $effect(() => { $refName; loadContigMap(); });

  let map = $derived($contigMap);
  let mapState = $derived($contigMapState);

  // ── Contig objects in the shape ReglScatter expects ──
  let explorerData = $derived.by(() => {
    if (!map) return null;
    const grand = map.sample_totals.reduce((a, b) => a + b, 0) || 1;
    const ph = map.kaiju_phylum, ge = map.kaiju_genus;
    const contigs = map.id.map((id, i) => ({
      id, idx: i,
      length: map.length[i],
      tsne_x: map.x[i], tsne_y: map.y[i],
      bin: map.bin[i] >= 0 ? map.bins[map.bin[i]] : '',
      kaiju_phylum: ph ? map.kaiju_phylum_vocab[ph[i]] : '',
      kaiju_genus: ge ? map.kaiju_genus_vocab[ge[i]] : '',
      reads: map.total_reads[i],
      n_detected: map.n_detected[i],
      rrna: map.rrna?.[i] === 1,
      // pooled RPKM across all samples
      depth: map.total_reads[i] / (map.length[i] / 1e3) / (grand / 1e6),
    }));
    return { contigs, has_tsne: true, has_umap: false, n_total: contigs.length };
  });

  // ── Colour / size controls ──
  let colorMode = $state('sample');   // 'sample' | 'total' | 'bins' | 'taxa' | 'length'
  let taxRank = $state('phylum');
  let sizeBy = $state('length');
  let sizeScale = $state(1.0);
  let sizeExponent = $state(1);
  let sampleIdx = $state(0);
  let layer = $state(null);
  let layerFor = $state(-1);

  $effect(() => {
    const m = map, j = sampleIdx;
    if (!m || colorMode !== 'sample') return;
    loadContigLayer(m, j).then(l => { if (j === sampleIdx) { layer = l; layerFor = j; } });
  });

  // Fixed colour ceiling across samples (the largest per-sample 99th percentile)
  let rpkmCeil = $derived(map ? Math.max(1, ...map.rpkm_p99) : 1);
  let sampleDepthData = $derived.by(() => {
    if (!map || !layer || layerFor !== sampleIdx) return null;
    const total = map.sample_totals[sampleIdx] || 1;
    const l = layer;
    return {
      samples: map.samples,
      maxDepth: rpkmCeil,
      depthOf: (c) => l[c.idx] / (c.length / 1e3) / (total / 1e6),
    };
  });
  let selectedSampleName = $derived(map?.samples?.[sampleIdx] ?? '');

  let colorBy = $derived(
    colorMode === 'sample' ? (sampleDepthData ? 'sample_depth' : 'depth') :
    colorMode === 'total' ? 'depth' :
    colorMode === 'bins' ? 'bin' :
    colorMode === 'taxa' ? `kaiju_${taxRank}` : 'length'
  );
  let scatterSize = $derived(sizeBy === 'expression' ? 'depth' : sizeBy);

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8',
                   '#94a3b8','#d4d4d8','#78716c'];
  let stableColorMap = $derived.by(() => {
    if (!explorerData || !['bin', 'kaiju_phylum', 'kaiju_genus'].includes(colorBy)) return {};
    const bg = colorBy === 'bin' ? 'unbinned' : 'Unknown';
    const names = colorBy === 'bin' ? map.bins
      : (map[`${colorBy}_vocab`] ?? []).filter(Boolean);
    const out = { [bg]: '#475569' };
    [...names].sort().forEach((n, i) => { out[n] = PALETTE[i % PALETTE.length]; });
    return out;
  });
  let coordExtents = $derived.by(() => {
    if (!map) return null;
    let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
    for (let i = 0; i < map.n; i++) {
      const x = map.x[i], y = map.y[i];
      if (x < xMin) xMin = x; if (x > xMax) xMax = x;
      if (y < yMin) yMin = y; if (y > yMax) yMax = y;
    }
    return { tsne: { xMin, xMax, yMin, yMax } };
  });

  const modeDefs = [
    { key: 'sample', label: 'Sample' }, { key: 'total', label: 'Total' },
    { key: 'bins', label: 'MAGs' }, { key: 'taxa', label: 'Taxon' }, { key: 'length', label: 'Length' },
  ];
  const rankDefs = ['phylum', 'genus'];
  const sizeDefs = [{ key: 'length', label: 'Length' }, { key: 'expression', label: 'Expression' }, { key: 'fixed', label: 'Fixed' }];
  const cycleSize = () => { const k = sizeDefs.map(d => d.key); sizeBy = k[(k.indexOf(sizeBy) + 1) % k.length]; };
  const step = d => { const n = map?.samples?.length ?? 0; if (n) sampleIdx = (sampleIdx + d + n) % n; };
  let hasTaxa = $derived(!!map?.kaiju_phylum_vocab?.some(Boolean));

  // ── Filters: log-scale dual-range sliders on length and pooled expression ──
  let lenMinPct = $state(0), lenMaxPct = $state(100);
  let expMinPct = $state(0), expMaxPct = $state(100);
  let ext = $derived.by(() => {
    if (!explorerData) return null;
    let lMin = Infinity, lMax = 0, eMin = Infinity, eMax = 0;
    for (const c of explorerData.contigs) {
      if (c.length < lMin) lMin = c.length; if (c.length > lMax) lMax = c.length;
      if (c.depth > 0 && c.depth < eMin) eMin = c.depth; if (c.depth > eMax) eMax = c.depth;
    }
    return { lMin: Math.max(1, lMin), lMax: Math.max(2, lMax), eMin: Math.max(1e-4, eMin), eMax: Math.max(1e-3, eMax) };
  });
  const pctToLog = (pct, lo, hi) => Math.pow(10, lo + (hi - lo) * pct / 100);
  let lenRange = $derived(ext ? [
    lenMinPct <= 0 ? ext.lMin : pctToLog(lenMinPct, Math.log10(ext.lMin), Math.log10(ext.lMax)),
    lenMaxPct >= 100 ? ext.lMax : pctToLog(lenMaxPct, Math.log10(ext.lMin), Math.log10(ext.lMax)),
  ] : [0, Infinity]);
  let expRange = $derived(ext ? [
    expMinPct <= 0 ? ext.eMin : pctToLog(expMinPct, Math.log10(ext.eMin), Math.log10(ext.eMax)),
    expMaxPct >= 100 ? ext.eMax : pctToLog(expMaxPct, Math.log10(ext.eMin), Math.log10(ext.eMax)),
  ] : [0, Infinity]);
  // Slider labels follow the drag; the filter itself applies once it pauses
  let committed = $state({ lenMinPct: 0, lenMaxPct: 100, expMinPct: 0, expMaxPct: 100, lenRange: [0, Infinity], expRange: [0, Infinity] });
  let commitTimer;
  $effect(() => {
    const next = { lenMinPct, lenMaxPct, expMinPct, expMaxPct, lenRange, expRange };
    clearTimeout(commitTimer);
    commitTimer = setTimeout(() => { committed = next; }, 150);
    return () => clearTimeout(commitTimer);
  });
  // Fixed size scale from the whole map: filtering never resizes points
  let sizeMax = $derived(ext ? { length: ext.lMax, depth: ext.eMax } : null);
  let isFiltered = $derived(committed.lenMinPct > 0 || committed.lenMaxPct < 100 || committed.expMinPct > 0 || committed.expMaxPct < 100);
  let passes = $derived.by(() => {
    const c0 = committed;
    const lMin = c0.lenMinPct > 0 ? c0.lenRange[0] : -Infinity, lMax = c0.lenMaxPct < 100 ? c0.lenRange[1] : Infinity;
    const eMin = c0.expMinPct > 0 ? c0.expRange[0] : -Infinity, eMax = c0.expMaxPct < 100 ? c0.expRange[1] : Infinity;
    return c => c.length >= lMin && c.length <= lMax && c.depth >= eMin && c.depth <= eMax;
  });

  let filteredData = $derived.by(() => {
    if (!explorerData) return null;
    if (!isFiltered) return explorerData;
    return { ...explorerData, contigs: explorerData.contigs.filter(passes) };
  });
  // Sorted by pooled expression once; filters keep this order for the table
  let byExpression = $derived(explorerData ? [...explorerData.contigs].sort((a, b) => b.depth - a.depth) : []);

  const fmtExp = v => v >= 100 ? Math.round(v).toLocaleString() : v >= 1 ? v.toFixed(1) : v.toFixed(3);

  // ── Search ──
  let searchQuery = $state('');
  let searchDebounced = $state('');
  let timer;
  $effect(() => {
    const q = searchQuery;
    clearTimeout(timer);
    timer = setTimeout(() => { searchDebounced = q.trim().toLowerCase(); }, 300);
    return () => clearTimeout(timer);
  });
  // Gene names, products and KOs: the index loads on the first keystroke
  $effect(() => { if (searchQuery.length > 0 && map) ensureGeneSearch(map); });
  let geneState = $derived($geneSearchState);
  let searchMatchIds = $derived.by(() => {
    const q = searchDebounced;
    if (!q || !filteredData) return null;
    const out = new Set();
    const genePos = geneState.ready ? searchGenes(q) : null;
    if (genePos) for (const c of filteredData.contigs) if (genePos.has(c.idx)) out.add(c.id);
    for (const c of filteredData.contigs) {
      if (c.id.toLowerCase().includes(q) || c.bin.toLowerCase().includes(q)
          || c.kaiju_phylum.toLowerCase().includes(q) || c.kaiju_genus.toLowerCase().includes(q)) out.add(c.id);
    }
    return out;
  });
  function resetFilters() { lenMinPct = 0; lenMaxPct = 100; expMinPct = 0; expMaxPct = 100; searchQuery = ''; }

  // ── Lasso selection ──
  let selectedIds = $state(null);
  let selectedContigs = $derived.by(() => {
    if (!filteredData || !selectedIds?.length) return null;
    const s = new Set(selectedIds);
    return filteredData.contigs.filter(c => s.has(c.id));
  });
  let selectionStats = $derived.by(() => {
    const cs = selectedContigs;
    if (!cs) return null;
    return {
      count: cs.length,
      totalLen: cs.reduce((s, c) => s + c.length, 0),
      reads: cs.reduce((s, c) => s + c.reads, 0),
      nBins: new Set(cs.map(c => c.bin).filter(Boolean)).size,
      nPhyla: new Set(cs.map(c => c.kaiju_phylum).filter(Boolean)).size,
    };
  });

  // ── Contig detail: expression across every sample (loads all layers once) ──
  let detailId = $state(null);
  let detail = $derived(explorerData?.contigs?.find(c => c.id === detailId) ?? null);
  let detailProfile = $state(null);
  let profileProgress = $state(0);
  let detailGenes = $state(null);
  let genesLoading = $state(false);
  async function openDetail(id) {
    detailId = id;
    detailProfile = null;
    const c = explorerData.contigs.find(x => x.id === id);
    if (!c) return;
    genesLoading = true;
    detailGenes = null;
    loadContigGenes(map, c.idx).then(g => { if (detailId === id) { detailGenes = g; genesLoading = false; } })
                               .catch(() => { if (detailId === id) { detailGenes = []; genesLoading = false; } });
    const n = map.samples.length;
    const out = new Array(n).fill(null);
    profileProgress = 0;
    for (let j = 0; j < n; j += 8) {
      await Promise.all(Array.from({ length: Math.min(8, n - j) }, (_, k) => j + k).map(async jj => {
        const l = await loadContigLayer(map, jj);
        out[jj] = l[c.idx] / (c.length / 1e3) / ((map.sample_totals[jj] || 1) / 1e6);
      }));
      profileProgress = Math.min(n, j + 8);
      if (detailId !== id) return;
    }
    detailProfile = out;
  }
  let autoOpened = false;
  $effect(() => {
    if (!autoOpened && explorerData?.contigs?.length) {
      autoOpened = true;
      let best = explorerData.contigs[0];
      for (const c of explorerData.contigs) if (c.reads > best.reads) best = c;
      openDetail(best.id);
    }
  });

  // ── Table ──
  let tableSrc = $derived.by(() => {
    if (!filteredData) return [];
    if (selectedContigs) return [...selectedContigs].sort((a, b) => b.depth - a.depth);
    const p = isFiltered ? passes : null;
    return byExpression.filter(c => (!p || p(c)) && (!searchMatchIds || searchMatchIds.has(c.id)));
  });
  let tableRows = $derived(tableSrc.slice(0, 500).map(c => ({
    ...c,
    sample_rpkm: sampleDepthData && colorMode === 'sample' ? sampleDepthData.depthOf(c) : null,
  })));
  let tableColumns = $derived([
    { key: 'id', label: 'Contig' },
    { key: 'length', label: 'Length', render: v => v.toLocaleString() },
    { key: 'bin', label: 'MAG', render: v => v || '-' },
    ...(hasTaxa ? [{ key: 'kaiju_phylum', label: 'Phylum', render: v => v || '-' }, { key: 'kaiju_genus', label: 'Genus', render: v => v || '-' }] : []),
    { key: 'reads', label: 'Read pairs', render: v => fmtInt(v) },
    { key: 'depth', label: 'Pooled RPKM', render: v => fmtExp(v) },
    { key: 'n_detected', label: 'Samples' },
    ...(colorMode === 'sample' && sampleDepthData ? [{ key: 'sample_rpkm', label: `RPKM ${selectedSampleName}`, render: v => v == null ? '-' : fmtExp(v) }] : []),
  ]);
  function exportTSV() {
    downloadTsv(`contigs_${selectedContigs ? 'selection' : 'all'}_${tableSrc.length}.tsv`,
      ['contig_id', 'length', 'mag', 'phylum', 'genus', 'read_pairs', 'pooled_rpkm', 'n_samples'],
      tableSrc.map(c => [c.id, c.length, c.bin, c.kaiju_phylum, c.kaiju_genus, c.reads, c.depth.toFixed(4), c.n_detected]));
  }

  const tip = c => {
    const s = sampleDepthData && colorMode === 'sample' ? ` | ${selectedSampleName}: ${fmtExp(sampleDepthData.depthOf(c))} RPKM` : '';
    const tax = c.kaiju_genus || c.kaiju_phylum;
    return `${c.id} | ${c.length.toLocaleString()} bp | ${fmtInt(c.reads)} read pairs | pooled ${fmtExp(c.depth)} RPKM${s}${c.bin ? ` | ${c.bin}` : ''}${tax ? ` | ${tax}` : ''}`;
  };
</script>

{#if !mapState.available}
  <div class="flex flex-col items-center justify-center py-24 text-center">
    <p class="text-slate-400 text-lg mb-2">No contig embedding for this reference</p>
    <p class="text-slate-500 text-sm max-w-xl">The Explorer reuses the contig t-SNE from a mag_analysis run's viz. Run with <code class="bg-slate-800 px-1 rounded">--from_mag</code> on an output that has <code class="bg-slate-800 px-1 rounded">viz/data/contig_explorer_*</code>, or pass <code class="bg-slate-800 px-1 rounded">--contig_explorer DIR</code>.</p>
  </div>
{:else if !explorerData}
  <div class="flex flex-col items-center justify-center py-20">
    <div class="w-10 h-10 border-4 border-slate-700 border-t-cyan-400 rounded-full animate-spin"></div>
    <p class="text-slate-400 text-sm mt-4">{mapState.error ? `Failed to load contig map: ${mapState.error}` : 'Loading contig map...'}</p>
  </div>
{:else}
<div class="sticky top-14 z-10 bg-slate-950 pb-2 -mx-4 px-4 pt-1">
<div class="mb-2 flex items-center gap-3 flex-wrap text-xs">
  <span class="text-slate-400">Color:</span>
  {#each modeDefs as d}
    {#if d.key !== 'taxa' || hasTaxa}
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center
          {colorMode === d.key ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        onclick={() => { if (d.key === 'taxa' && colorMode === 'taxa') taxRank = rankDefs[(rankDefs.indexOf(taxRank) + 1) % rankDefs.length]; else colorMode = d.key; }}
        title={d.key === 'taxa' ? 'Kaiju taxonomy from the MAG analysis; click again to cycle Phylum → Genus' : ''}
      >{d.key === 'taxa' && colorMode === 'taxa' ? `Kaiju ${taxRank}` : d.label}{d.key === 'taxa' ? ' ▾' : ''}</button>
    {/if}
  {/each}
  {#if colorMode === 'sample'}
    <button class="px-2 py-1 rounded-md border border-slate-600 text-slate-400 hover:border-slate-500" onclick={() => step(-1)} title="Previous sample">&#x25C0;</button>
    <select
      class="px-2 py-1 rounded-md border border-slate-600 bg-slate-800 text-slate-300 text-xs focus:border-cyan-400 focus:outline-none cursor-pointer font-mono"
      bind:value={sampleIdx}
    >
      {#each map.samples as s, j}<option value={j}>{s}</option>{/each}
    </select>
    <button class="px-2 py-1 rounded-md border border-slate-600 text-slate-400 hover:border-slate-500" onclick={() => step(1)} title="Next sample">&#x25B6;</button>
    {#if layerFor !== sampleIdx}<span class="text-slate-500 italic">loading sample...</span>{/if}
  {/if}
  <span class="text-slate-600">|</span>
  <span class="text-slate-400">Size:</span>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
    style="min-width: 5.5rem"
    onclick={cycleSize}
    title={`Click to cycle: ${sizeDefs.map(d => d.label).join(' → ')}`}
  >{sizeDefs.find(d => d.key === sizeBy)?.label} &#x25BE;</button>
  <input type="range" min="0.2" max="20" step="0.1" bind:value={sizeScale} class="w-20 h-1 accent-cyan-400 cursor-pointer" />
  <span class="text-slate-500 w-8">{sizeScale.toFixed(1)}x</span>
  <span class="text-slate-500" title="Size exponent: higher makes the largest contigs stand out">exp</span>
  <input type="range" min="1" max="4" step="0.25" bind:value={sizeExponent} class="w-16 h-1 accent-cyan-400 cursor-pointer"
    title="Size exponent: higher makes the largest contigs stand out" />
  <span class="text-slate-500 w-8">{sizeExponent.toFixed(2).replace(/0$/, '')}</span>
</div>

<div class="flex items-center gap-3 flex-wrap text-xs">
  <span class="text-slate-400">Filter:</span>

  <span class="text-slate-500">Length</span>
  <span class="text-slate-500 w-16 text-right font-mono whitespace-nowrap">{fmtBp(lenRange[0])}</span>
  <div class="dual-range relative w-36 h-5 flex items-center">
    <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
    <div class="absolute h-1 bg-cyan-400/40 rounded" style="left: {lenMinPct}%; right: {100 - lenMaxPct}%"></div>
    <input type="range" min="0" max="100" step="1" bind:value={lenMinPct}
      oninput={() => { if (lenMinPct > lenMaxPct - 2) lenMinPct = lenMaxPct - 2; }} />
    <input type="range" min="0" max="100" step="1" bind:value={lenMaxPct}
      oninput={() => { if (lenMaxPct < lenMinPct + 2) lenMaxPct = lenMinPct + 2; }} />
  </div>
  <span class="text-slate-500 w-16 font-mono whitespace-nowrap">{fmtBp(lenRange[1])}</span>

  <span class="text-slate-600">|</span>

  <span class="text-slate-500" title="Reads per kb per million, pooled over all samples">Expression</span>
  <span class="text-slate-500 w-16 text-right font-mono whitespace-nowrap">{fmtExp(expRange[0])}</span>
  <div class="dual-range relative w-36 h-5 flex items-center">
    <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
    <div class="absolute h-1 bg-cyan-400/40 rounded" style="left: {expMinPct}%; right: {100 - expMaxPct}%"></div>
    <input type="range" min="0" max="100" step="1" bind:value={expMinPct}
      oninput={() => { if (expMinPct > expMaxPct - 2) expMinPct = expMaxPct - 2; }} />
    <input type="range" min="0" max="100" step="1" bind:value={expMaxPct}
      oninput={() => { if (expMaxPct < expMinPct + 2) expMaxPct = expMinPct + 2; }} />
  </div>
  <span class="text-slate-500 w-16 font-mono whitespace-nowrap">{fmtExp(expRange[1])}</span>

  <span class="text-slate-600">|</span>

  <div class="relative flex items-center">
    <input type="text" bind:value={searchQuery} placeholder={hasTaxa ? 'contig, MAG, taxon, gene, product, KO...' : 'contig, MAG, gene, product, KO...'}
      class="w-48 px-2 py-1 rounded bg-slate-800 border border-slate-600 text-slate-200 text-xs placeholder-slate-500 focus:border-cyan-400 focus:outline-none" />
    {#if searchQuery}
      <button class="absolute right-1 text-slate-500 hover:text-slate-300" onclick={() => { searchQuery = ''; }} title="Clear search">
        <svg class="w-3.5 h-3.5" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M18 6L6 18M6 6l12 12" /></svg>
      </button>
    {/if}
  </div>
  {#if geneState.loading}<span class="text-slate-500 italic">loading gene index...</span>{/if}
  {#if searchMatchIds}<span class="text-amber-400 font-medium">{searchMatchIds.size.toLocaleString()} matching contigs</span>{/if}
  {#if isFiltered || searchMatchIds}
    {#if isFiltered}<span class="text-cyan-400 font-medium">{filteredData.contigs.length.toLocaleString()} / {explorerData.n_total.toLocaleString()}</span>{/if}
    <button class="text-slate-500 hover:text-slate-300 underline" onclick={resetFilters}>reset</button>
  {/if}
</div>
</div>

<div class="flex gap-4 mb-6">
  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700 h-[820px] flex flex-col flex-1 min-w-0">
    <ReglScatter data={filteredData} {colorBy} sizeBy={scatterSize} {sizeScale} {sizeExponent} sizeRange={sizeMax} mode="tsne"
      colorMap={stableColorMap} {coordExtents} onselect={ids => selectedIds = ids} onclick={openDetail}
      {searchMatchIds} {sampleDepthData} selectedSample={colorMode === 'sample' ? selectedSampleName : ''} tipFormat={tip} />
    <div class="text-xs text-slate-500 mt-2 flex-shrink-0">
      {filteredData.contigs.length.toLocaleString()} of {map.n_contigs_total.toLocaleString()} contigs (those with RNA reads) |
      regl-scatterplot | t-SNE from the MAG analysis |
      {colorMode === 'sample' ? `colour: RPKM in ${selectedSampleName}, fixed scale to ${fmtExp(rpkmCeil)}` : colorMode === 'total' ? 'colour: pooled RPKM (rank)' : ''} |
      Shift+drag to lasso | Click point for details
    </div>
  </div>

  <div class="w-1/2 min-w-[420px] max-w-[640px]">
    <ContigGenes contig={detail} genes={detailGenes} loading={genesLoading} samples={map.samples}
      sampleIdx={colorMode === 'sample' ? sampleIdx : null} contigProfile={detailProfile}
      query={searchDebounced} onclose={() => detailId = null} />
  </div>
</div>

{#if selectionStats}
  <div class="flex items-center gap-4 mb-4 text-xs flex-wrap">
    <span class="text-slate-400 font-medium">Selection:</span>
    <span class="text-cyan-400">{selectionStats.count.toLocaleString()} contigs</span>
    <span class="text-slate-400">{(selectionStats.totalLen / 1e6).toFixed(2)} Mb</span>
    <span class="text-slate-400">{fmtInt(selectionStats.reads)} read pairs</span>
    <span class="text-slate-400">{selectionStats.nBins} MAGs</span>
    {#if hasTaxa}<span class="text-slate-400">{selectionStats.nPhyla} phyla</span>{/if}
    <button class="text-slate-500 hover:text-slate-300 underline" onclick={() => selectedIds = null}>clear</button>
  </div>
{/if}

<div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
  <div class="flex items-center justify-between mb-2">
    <h3 class="text-sm font-medium text-slate-400">
      {selectedContigs ? `${selectedContigs.length.toLocaleString()} selected contigs` : searchMatchIds ? `${searchMatchIds.size.toLocaleString()} matching contigs` : 'Top 500 contigs by pooled expression'}
      {#if tableSrc.length > 500}<span class="text-slate-500 text-xs">({tableSrc.length.toLocaleString()} total)</span>{/if}
    </h3>
    <button
      class="text-xs px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-slate-200 hover:border-slate-500 transition-colors"
      onclick={exportTSV}
      title="Export all {tableSrc.length.toLocaleString()} contigs as TSV"
    >TSV</button>
  </div>
  <DataTable columns={tableColumns} rows={tableRows} idKey="id" maxHeight="350px" onRowClick={row => openDetail(row.id)} selectedId={detailId} hideExport={true} />
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
  .dual-range input[type="range"]::-webkit-slider-runnable-track { height: 0; }
  .dual-range input[type="range"]::-moz-range-track { height: 0; background: transparent; }
</style>
