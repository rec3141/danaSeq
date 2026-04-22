<script>
  import { onMount } from 'svelte';
  import ReglScatter from '../components/charts/ReglScatter.svelte';
  import StatCard from '../components/layout/StatCard.svelte';
  import { readExplorer, loadReadExplorer, sampleTaxonomy } from '../stores/data.js';
  import { cartItems, cartActive } from '../stores/cart.js';
  import { selectedRead } from '../stores/selection.js';
  import { taxonomySource } from '../stores/taxonomySource.js';
  import { taxNav, activeSubTaxa, rankOrder, RANK_LABELS } from '../stores/taxonomy.js';
  import TaxonomyDrillNav from '../components/layout/TaxonomyDrillNav.svelte';

  // Map rank code → the per-read field suffix emitted by compute_read_tsne.py.
  // R2 (Kraken superkingdom) and D (GTDB domain) both surface as *_domain
  // on reads; S surfaces as *_species. Reads now carry the full chain so
  // drilling down to species colors/filters without a fallback.
  const RANK_FIELD = { D: 'domain', R2: 'domain', P: 'phylum', C: 'class',
                       O: 'order', F: 'family', G: 'genus', S: 'species' };

  let loading = $state(true);

  // Cycling button helpers
  function cycle(values, current) {
    const idx = values.indexOf(current);
    return values[(idx + 1) % values.length];
  }
  function getLabel(values, labels, current) {
    const idx = values.indexOf(current);
    return idx >= 0 ? labels[idx] : labels[0];
  }

  // Taxonomy label/cycle state. The cycle button walks the active rank order
  // via taxNav.setLevel() — this shares the single source of truth with the
  // drill-down sidebar, so both cycling and drilling move the same cursor
  // and stay in sync.
  const metricGroup = { values: ['gc', 'length'], labels: ['GC%', 'Length'] };

  const BW = { sample: '4.5rem', taxonomy: '4.5rem', metric: '4.5rem' };

  let colorMode = $state('metric');  // 'sample' | 'taxonomy' | 'metric'
  let metric = $state('gc');         // default metric → GC%

  // Scatter color field = `${source}_${rank_suffix}` from the active drill
  // level. Missing entries (shouldn't happen now that RANK_FIELD covers all
  // drillable ranks) fall back to phylum so the scatter still has something
  // to color by.
  let taxField = $derived(
    `${$taxonomySource}_${RANK_FIELD[$taxNav.level] ?? 'phylum'}`
  );

  let sizeScale = $state(0.6);

  let colorBy = $derived(
    colorMode === 'sample' ? 'sample' :
    colorMode === 'taxonomy' ? taxField :
    metric
  );

  // Color map at the current drill rank — shared with the sidebar.
  let taxColorMap = $derived(
    colorMode === 'taxonomy' && $activeSubTaxa ? $activeSubTaxa.colorMap : {}
  );

  // Cycle forward through the active rank order (R2→P→C→O→F→G→S→R2 for Kraken,
  // D→…→S→D for GTDB). Overrides the drill state without wiping the filter.
  function cycleTaxRank() {
    const order = $rankOrder;
    const i = order.indexOf($taxNav.level);
    const next = order[(i + 1) % order.length];
    taxNav.setLevel(next);
  }

  onMount(async () => {
    await loadReadExplorer();
    loading = false;
  });

  // Filter reads by cart + (in taxonomy mode) by the drill-down filter.
  // The drill filter scopes to the filter taxon's rank: when drilled to a
  // Phylum, we keep reads whose `${source}_phylum` matches; at Class, we
  // keep reads whose class matches; etc. Isolation (leaf-rank) also filters
  // to that taxon at its own rank.
  let scatterData = $derived.by(() => {
    if (!$readExplorer?.reads) return null;
    let reads = $readExplorer.reads;

    if ($cartActive && $cartItems.size > 0) {
      reads = reads.filter(r => $cartItems.has(r.sample));
    }

    if (colorMode === 'taxonomy' && $sampleTaxonomy?.ranks) {
      const tax = $sampleTaxonomy;
      // Isolation takes precedence when set (species click at S → isolate
      // that species); otherwise filter by the drill breadcrumb. Both use
      // the taxon's rank to pick which per-read field to match against.
      const scope = $taxNav.isolated || $taxNav.filter;
      if (scope) {
        const scopeRank = tax.ranks[scope];
        const field = scopeRank ? RANK_FIELD[scopeRank] : null;
        if (field) {
          const readField = `${$taxonomySource}_${field}`;
          reads = reads.filter(r => r[readField] === scope);
        }
      }
    }

    return { points: reads };
  });

  // ---- Search state ----
  let searchQuery = $state('');
  let searchDebounced = $state('');
  let debounceTimer;

  $effect(() => {
    const q = searchQuery;
    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(() => { searchDebounced = q; }, 300);
    return () => clearTimeout(debounceTimer);
  });

  // Search across both taxonomies so a user can find reads by either
  // Kraken or GTDB names regardless of the current toggle.
  const SEARCH_FIELDS = [
    'sample',
    'kraken_domain', 'kraken_phylum', 'kraken_class',
    'kraken_order', 'kraken_family', 'kraken_genus', 'kraken_species',
    'gtdb_domain', 'gtdb_phylum', 'gtdb_class',
    'gtdb_order', 'gtdb_family', 'gtdb_genus', 'gtdb_species',
    'genes', 'products',
  ];

  // Combined cart + search highlighting
  let searchMatchIds = $derived.by(() => {
    const hasCart = $cartActive && $cartItems.size > 0;
    const q = searchDebounced.trim().toLowerCase();
    if (!hasCart && !q) return null;
    if (!$readExplorer?.reads) return null;

    const ids = new Set();
    for (const r of $readExplorer.reads) {
      if (hasCart && !$cartItems.has(r.sample)) continue;
      if (q) {
        let found = false;
        for (const f of SEARCH_FIELDS) {
          const val = r[f];
          if (val && String(val).toLowerCase().includes(q)) { found = true; break; }
        }
        if (!found) continue;
      }
      ids.add(r.id);
    }
    return ids;
  });

  let searchMatchCount = $derived(searchMatchIds?.size ?? null);

  let selectedIds = $state(null);  // Set of selected read IDs from lasso

  function handleSelect(ids) {
    selectedIds = ids ? new Set(ids) : null;
  }

  function exportSelection() {
    if (!selectedIds || !$readExplorer?.reads) return;
    const selected = $readExplorer.reads.filter(r => selectedIds.has(r.id));
    const cols = [
      'id', 'sample', 'length', 'gc',
      'kraken_domain', 'kraken_phylum', 'kraken_class', 'kraken_order', 'kraken_family', 'kraken_genus', 'kraken_species',
      'gtdb_domain', 'gtdb_phylum', 'gtdb_class', 'gtdb_order', 'gtdb_family', 'gtdb_genus', 'gtdb_species',
      'genes', 'products',
    ];
    const header = cols.join('\t');
    const rows = selected.map(r => cols.map(c => r[c] ?? '').join('\t'));
    const tsv = header + '\n' + rows.join('\n') + '\n';
    const blob = new Blob([tsv], { type: 'text/tab-separated-values' });
    const a = document.createElement('a');
    a.href = URL.createObjectURL(blob);
    a.download = `selected_reads_${selected.length}.tsv`;
    a.click();
    URL.revokeObjectURL(a.href);
  }

  let selectedDetail = $derived.by(() => {
    if (!$selectedRead || !$readExplorer?.reads) return null;
    return $readExplorer.reads.find(r => r.id === $selectedRead);
  });

  // Selection-panel rank picker. Tracks the active classifier's rank order
  // so the cycle walks the same levels as the main Taxonomy button and
  // the drill-down sidebar.
  let selTaxRanks = $derived(
    $rankOrder
      .filter(r => RANK_FIELD[r])
      .map(r => ({ key: `${$taxonomySource}_${RANK_FIELD[r]}`, label: RANK_LABELS[r] ?? r }))
  );
  let selTaxIdx = $state(1); // default phylum (index 1 in both orderings)

  let selectionStats = $derived.by(() => {
    if (!selectedIds || !$readExplorer?.reads) return null;
    const selected = $readExplorer.reads.filter(r => selectedIds.has(r.id));
    if (!selected.length) return null;
    const n = selected.length;
    const totalBp = selected.reduce((s, r) => s + (r.length || 0), 0);
    const avgLen = totalBp / n;
    const avgGc = selected.reduce((s, r) => s + (r.gc || 0), 0) / n;
    const samples = new Set(selected.map(r => r.sample));

    // Taxonomy at selected rank
    const rankKey = selTaxRanks[selTaxIdx].key;
    const taxCounts = {};
    for (const r of selected) {
      const t = r[rankKey] || 'Unclassified';
      taxCounts[t] = (taxCounts[t] || 0) + 1;
    }
    const topTaxa = Object.entries(taxCounts)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 12)
      .map(([name, count]) => ({ name, count, pct: (100 * count / n).toFixed(1) }));

    // Top gene products
    const prodCounts = {};
    for (const r of selected) {
      if (!r.products) continue;
      for (const p of r.products.split('; ')) {
        if (p) prodCounts[p] = (prodCounts[p] || 0) + 1;
      }
    }
    const topProducts = Object.entries(prodCounts)
      .sort((a, b) => b[1] - a[1])
      .slice(0, 12)
      .map(([name, count]) => ({ name, count }));

    return { n, totalBp, avgLen, avgGc, nSamples: samples.size, topTaxa, topProducts };
  });

  let stats = $derived.by(() => {
    const reads = scatterData?.points;
    if (!reads?.length) return null;
    const n = reads.length;
    const samples = new Set(reads.map(r => r.sample));
    const totalBp = reads.reduce((s, r) => s + (r.length || 0), 0);
    const avgLen = totalBp / n;
    const avgGc = reads.reduce((s, r) => s + (r.gc || 0), 0) / n;
    return { n, nSamples: samples.size, totalBp, avgLen, avgGc };
  });
</script>

<div class="space-y-6">
  {#if loading}
    <div class="flex flex-col items-center justify-center py-20">
      <div class="w-10 h-10 border-4 border-slate-700 border-t-cyan-400 rounded-full animate-spin"></div>
      <p class="text-slate-400 text-sm mt-4">Loading read data (this may take a moment)...</p>
    </div>
  {:else}
    <!-- Stats (update to selection when active) -->
    {#if stats}
      {@const s = selectionStats || stats}
      {@const label = selectionStats ? 'Selected' : 'Reads Shown'}
      <div class="grid grid-cols-2 md:grid-cols-5 gap-4">
        <StatCard label={label} value={s.n.toLocaleString()} color="cyan" />
        <StatCard label="Total bp" value={s.totalBp >= 1e9 ? `${(s.totalBp/1e9).toFixed(1)} Gb` : s.totalBp >= 1e6 ? `${(s.totalBp/1e6).toFixed(1)} Mb` : s.totalBp.toLocaleString()} color="teal" />
        <StatCard label="Samples" value={s.nSamples} color="emerald" />
        <StatCard label="Avg Length" value={`${Math.round(s.avgLen).toLocaleString()} bp`} color="amber" />
        <StatCard label="Avg GC" value={`${s.avgGc.toFixed(1)}%`} color="purple" />
      </div>
    {/if}

    <!-- Controls -->
    <div class="flex items-center gap-2 text-xs">
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center
          {colorMode === 'sample' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        style="min-width: {BW.sample}"
        onclick={() => colorMode = 'sample'}
      >
        Sample
      </button>
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center
          {colorMode === 'taxonomy' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        style="min-width: {BW.taxonomy}"
        onclick={() => { if (colorMode === 'taxonomy') cycleTaxRank(); else colorMode = 'taxonomy'; }}
        title={`Click to cycle rank (${$taxonomySource.toUpperCase()}): ${$rankOrder.map(r => RANK_LABELS[r] ?? r).join(' → ')}`}
      >
        {colorMode === 'taxonomy' ? (RANK_LABELS[$taxNav.level] ?? $taxNav.level) : 'Taxonomy'} &#x25BE;
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
      <span class="text-slate-600 mx-1">|</span>

      <div class="text-slate-400 flex items-center gap-1">
        Size
        <div class="single-range relative w-16 h-5 flex items-center">
          <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
          <input type="range" min="0.2" max="20" step="0.1" bind:value={sizeScale} />
        </div>
        <span class="text-slate-500 w-8 font-mono">{sizeScale.toFixed(1)}x</span>
      </div>

      {#if scatterData?.points}
        <span class="text-slate-500 ml-1">
          {scatterData.points.length.toLocaleString()} reads
        </span>
      {/if}

      <span class="text-slate-600 mx-1">|</span>

      <div class="relative flex items-center">
        <input
          type="text"
          bind:value={searchQuery}
          placeholder="taxon, gene, product..."
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
      {#if searchMatchCount !== null}
        <span class="text-amber-400 font-medium">{searchMatchCount.toLocaleString()} matches</span>
      {/if}
    </div>

    <!-- Scatter -->
    <div class="flex gap-6">
      {#if colorMode === 'taxonomy'}
        <TaxonomyDrillNav heightClass="h-[700px]" />
      {/if}
      <div class="flex-1 h-[700px] flex flex-col">
        {#if scatterData}
          <ReglScatter
            data={scatterData}
            {colorBy}
            colorMap={taxColorMap}
            sizeBy="fixed"
            {sizeScale}
            mode="tsne"
            {searchMatchIds}
            onselect={handleSelect}
            onclick={(id) => { $selectedRead = id; }}
            exportName={`danaseq_read_tsne_color-${colorBy}`}
          />
        {:else}
          <div class="flex items-center justify-center h-full text-slate-500">
            No read data available. Run preprocess to generate read_explorer.json.
          </div>
        {/if}
      </div>

      <!-- Right panel: lasso-selection summary AND single-read detail
           stack together so clicking a read inside an active selection
           doesn't hide either view. Both panels share the same w-96
           column; the wrapper handles width/flex-shrink. -->
      {#if selectionStats || selectedDetail}
        <div class="w-96 flex-shrink-0 space-y-3 h-fit">
          {#if selectionStats}
            <div class="bg-slate-800 rounded-lg border border-cyan-900 p-4 space-y-3">
          <div class="flex items-center justify-between">
            <h3 class="text-sm font-semibold text-cyan-400">{selectionStats.n.toLocaleString()} reads selected</h3>
            <button
              class="text-slate-500 hover:text-slate-300 text-xs"
              onclick={() => { selectedIds = null; }}
              title="Clear selection"
            >&#x2715;</button>
          </div>
          <div class="space-y-1">
            <div class="flex items-center justify-between">
              <button
                class="text-xs text-cyan-400 hover:text-cyan-300 transition-colors"
                onclick={() => { selTaxIdx = (selTaxIdx + 1) % selTaxRanks.length; }}
                title="Click to cycle rank"
              >{selTaxRanks[selTaxIdx].label} &#x25BE;</button>
              <span class="text-[10px] text-slate-600">{selectionStats.topTaxa.length} shown</span>
            </div>
            {#each selectionStats.topTaxa as t}
              <div class="flex justify-between text-xs gap-1">
                <span class="text-slate-300 truncate" title={t.name}>{t.name}</span>
                <span class="text-slate-500 font-mono flex-shrink-0">{t.count.toLocaleString()} <span class="text-slate-600">({t.pct}%)</span></span>
              </div>
            {/each}
          </div>
          {#if selectionStats.topProducts.length > 0}
            <div class="space-y-1">
              <div class="text-xs text-slate-400">Functions</div>
              {#each selectionStats.topProducts as p}
                <div class="flex justify-between text-xs gap-1">
                  <span class="text-slate-300 truncate" title={p.name}>{p.name}</span>
                  <span class="text-slate-500 font-mono flex-shrink-0">{p.count.toLocaleString()}</span>
                </div>
              {/each}
            </div>
          {/if}
          <button
            class="w-full px-3 py-1.5 text-xs bg-cyan-600 hover:bg-cyan-500 text-white rounded transition-colors"
            onclick={exportSelection}
          >Export TSV</button>
        </div>
      {/if}
      {#if selectedDetail}
        <div class="bg-slate-800 rounded-lg border border-slate-700 p-4 space-y-3">
          <h3 class="text-sm font-semibold text-cyan-400 font-mono truncate" title={selectedDetail.id}>
            {selectedDetail.id.slice(0, 20)}...
          </h3>
          <div class="grid grid-cols-2 gap-2 text-xs">
            <div class="text-slate-400">Sample</div><div class="text-slate-200 font-mono">{selectedDetail.sample}</div>
            <div class="text-slate-400">Length</div><div class="text-slate-200 font-mono">{selectedDetail.length?.toLocaleString() ?? '-'} bp</div>
            <div class="text-slate-400">GC%</div><div class="text-slate-200 font-mono">{selectedDetail.gc?.toFixed(1) ?? '-'}</div>
            {#if selectedDetail.kraken_domain}
              <div class="text-slate-400">Domain</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.kraken_domain}</div>
            {/if}
            {#if selectedDetail.kraken_phylum}
              <div class="text-slate-400">Phylum</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.kraken_phylum}</div>
            {/if}
            {#if selectedDetail.kraken_class}
              <div class="text-slate-400">Class</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.kraken_class}</div>
            {/if}
            {#if selectedDetail.kraken_order}
              <div class="text-slate-400">Order</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.kraken_order}</div>
            {/if}
            {#if selectedDetail.kraken_family}
              <div class="text-slate-400">Family</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.kraken_family}</div>
            {/if}
            {#if selectedDetail.kraken_genus}
              <div class="text-slate-400">Genus</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.kraken_genus}</div>
            {/if}
            {#if selectedDetail.sketch_hit}
              <div class="text-slate-400">Sketch Hit</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.sketch_hit}</div>
            {/if}
            {#if selectedDetail.genes}
              <div class="text-slate-400">Genes</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.genes}</div>
            {/if}
            {#if selectedDetail.products}
              <div class="text-slate-400">Products</div><div class="text-slate-200 font-mono text-[11px] whitespace-normal">{selectedDetail.products}</div>
            {/if}
          </div>
        </div>
      {/if}
        </div>
      {/if}
    </div>
  {/if}
</div>

<style>
  .single-range input[type="range"] {
    -webkit-appearance: none;
    appearance: none;
    background: transparent;
    position: absolute;
    width: 100%;
    height: 100%;
    margin: 0;
    padding: 0;
    cursor: pointer;
  }
  .single-range input[type="range"]::-webkit-slider-thumb {
    -webkit-appearance: none;
    height: 14px;
    width: 14px;
    border-radius: 50%;
    background: #22d3ee;
    cursor: pointer;
    border: 2px solid #0f172a;
    box-shadow: 0 0 3px rgba(0,0,0,0.4);
  }
  .single-range input[type="range"]::-moz-range-thumb {
    height: 14px;
    width: 14px;
    border-radius: 50%;
    background: #22d3ee;
    cursor: pointer;
    border: 2px solid #0f172a;
    box-shadow: 0 0 3px rgba(0,0,0,0.4);
  }
  .single-range input[type="range"]::-webkit-slider-runnable-track { height: 0; }
  .single-range input[type="range"]::-moz-range-track { height: 0; background: transparent; }
</style>
