<script>
  import PlotlyScatter from '../components/charts/PlotlyScatter.svelte';
  import StatCard from '../components/layout/StatCard.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, sampleTsne, sampleTaxonomy, metadata, overview } from '../stores/data.js';
  import { cartItems, cartActive, addToCart, toggleCart } from '../stores/cart.js';
  import { selectedSample } from '../stores/selection.js';
  import { taxNav, activeSubTaxa, RANK_LABELS } from '../stores/taxonomy.js';
  import TaxonomyDrillNav from '../components/layout/TaxonomyDrillNav.svelte';

  // Taxonomy-overlay view: one disc per (sample, sub-taxon) at the current
  // drill level. Size = sqrt(fraction_of_sample_reads) × scale, so abundances
  // are cross-sample comparable (a 10% class in Sample A looks √10 bigger
  // than a 1% class in Sample B). Replaces the prior `colorBy = dominant_*`
  // hack which could only show a single lineage's dominant taxon per sample.
  let sizeScale = $state(1.0);
  let lassoIds = $state(null);

  // log2 read count filter
  let log2Min = $state(0);
  let log2Max = $state(21);
  let log2Ceil = $derived.by(() => {
    if (!$samples) return 21;
    const maxReads = Math.max(...$samples.map(s => s.read_count || 0));
    return Math.ceil(Math.log2(Math.max(1, maxReads)));
  });

  function readFilterLabel(log2val) {
    const v = Math.round(Math.pow(2, log2val));
    if (v >= 1e6) return `${(v/1e6).toFixed(1)}M`;
    if (v >= 1e3) return `${(v/1e3).toFixed(0)}K`;
    return String(v);
  }

  // Date filter (initialized once from data range)
  let dateMin = $state(null);
  let dateMax = $state(null);
  let dateRange = $derived.by(() => {
    if (!$samples) return { min: '', max: '' };
    const dates = $samples.map(s => s.start_time?.slice(0, 10)).filter(Boolean).sort();
    if (!dates.length) return { min: '', max: '' };
    return { min: dates[0], max: dates[dates.length - 1] };
  });

  // Default date inputs to full range once data loads
  $effect(() => {
    if (dateRange.min && dateMin === null) dateMin = dateRange.min;
    if (dateRange.max && dateMax === null) dateMax = dateRange.max;
  });

  function passesFilters(s) {
    const rc = s.read_count || 0;
    if (rc === 0) return false;
    const l = Math.log2(Math.max(1, rc));
    if (l < log2Min || l > log2Max) return false;
    if (s.start_time) {
      const d = s.start_time.slice(0, 10);
      if (dateMin && d < dateMin) return false;
      if (dateMax && d > dateMax) return false;
    }
    // Metadata-filter-sidebar selection (only applies in non-taxonomy
    // modes where the sidebar is visible).
    if (colorMode !== 'taxonomy' && filterVal !== null && filterCol) {
      const v = $metadata?.[s.id]?.[filterCol];
      if (v === undefined || v === null || String(v) !== filterVal) return false;
    }
    return true;
  }

  // Stats for the currently filtered (and optionally cart-filtered) dataset
  let filteredStats = $derived.by(() => {
    if (!$samples) return null;
    const passedFilter = $samples.filter(passesFilters);
    const visible = $cartActive && $cartItems.size > 0
      ? passedFilter.filter(s => $cartItems.has(s.id))
      : passedFilter;
    const flowcells = new Set(visible.map(s => s.flowcell).filter(Boolean));
    const totalReads = visible.reduce((sum, s) => sum + (s.read_count || 0), 0);
    const totalBases = visible.reduce((sum, s) => sum + (s.total_bases || 0), 0);
    return {
      nFlowcells: flowcells.size,
      nSamples: visible.length,
      totalReads,
      totalBases,
      isFiltered: visible.length < $samples.length,
      totalAll: $samples.length,
    };
  });

  // Embedding mode: cycle between available embeddings (derived from data)
  const EMBED_LABELS = { grid: 'Grid', tsne: 't-SNE' };
  let embedIdx = $state(-1); // -1 = uninitialized, will auto-pick on data load
  let availableModes = $derived.by(() => {
    if (!$sampleTsne) return ['grid'];
    return Object.keys($sampleTsne).filter(k => $sampleTsne[k] && Object.keys($sampleTsne[k]).length > 0);
  });

  // Default to t-SNE when available
  $effect(() => {
    if (embedIdx === -1 && availableModes.length > 0) {
      const tsneIdx = availableModes.indexOf('tsne');
      embedIdx = tsneIdx >= 0 ? tsneIdx : 0;
    }
  });

  let embedMode = $derived(availableModes[Math.max(0, embedIdx) % availableModes.length] || 'grid');

  function cycleEmbed() {
    embedIdx++;
  }

  // ------------------------------------------------------------------
  // Search + color-mode cycling.
  //
  // Search: free-text substring match across id / flowcell / barcode /
  // dominant_phylum / dominant_class / any uploaded metadata value.
  // Result is a Set<sampleId> passed to PlotlyScatter via searchMatchIds
  // (dims non-matches there) and also filters the sample table below.
  //
  // Color mode: cycles between the default taxonomy overlay (colorBy
  // 'taxon', using the shared drill-down palette) and coloring every disc
  // of a sample by a sample-level field — a cycled metadata column, or a
  // cycled numeric metric (reads / bases / avg length / Shannon H).
  // ------------------------------------------------------------------

  let searchQuery = $state('');
  let searchMatchIds = $derived.by(() => {
    const q = searchQuery.trim().toLowerCase();
    if (!q || !$samples) return null;
    const out = new Set();
    for (const s of $samples) {
      const meta = $metadata?.[s.id];
      const parts = [
        s.id, s.flowcell, s.barcode, s.dominant_phylum, s.dominant_class,
        ...(meta ? Object.values(meta) : []),
      ];
      const blob = parts.filter(v => v != null).map(String).join(' ').toLowerCase();
      if (blob.includes(q)) out.add(s.id);
    }
    return out;
  });

  let colorMode = $state('taxonomy');  // 'taxonomy' | 'metric' | 'metadata'
  let metaField = $state('');
  let metricField = $state('read_count');

  const METRIC_VALUES = ['read_count', 'total_bases', 'avg_length', 'diversity'];
  const METRIC_LABELS = {
    read_count: 'Reads', total_bases: 'Bases',
    avg_length: 'Avg Length', diversity: 'Shannon H',
  };

  // Discover uploaded metadata columns. lat/lon are excluded because
  // coloring by a lat/lon gradient on a t-SNE is meaningless; users who
  // want geographic coloring use the /map view.
  let metaGroup = $derived.by(() => {
    if (!$metadata) return { values: [], labels: [] };
    const cols = new Set();
    for (const m of Object.values($metadata)) {
      for (const k of Object.keys(m)) {
        if (k === 'lat' || k === 'lon') continue;
        cols.add(k);
      }
    }
    const values = [...cols].sort();
    return { values, labels: values.map(k => k.replace(/_/g, ' ')) };
  });

  // Categorical subset of the metadata columns — any column whose values
  // aren't uniformly numeric. Used as the left-sidebar filter's column
  // cycle in non-taxonomy modes; a numeric column (e.g. depth_m) doesn't
  // make sense as a "pick a value" filter.
  let categoricalMetaCols = $derived.by(() => {
    if (!$metadata) return [];
    return metaGroup.values.filter(k => {
      const vals = Object.values($metadata)
        .map(m => m[k])
        .filter(v => v !== undefined && v !== null && v !== '');
      if (vals.length === 0) return false;
      return !vals.every(v => typeof v === 'number');
    });
  });

  // Filter sidebar state: which categorical column are we listing values
  // from, and (optionally) which value is selected as an active filter.
  let filterCol = $state('');
  let filterVal = $state(null);

  $effect(() => {
    // Initialize / re-initialize when the available categorical set changes.
    if (!filterCol && categoricalMetaCols.length > 0) filterCol = categoricalMetaCols[0];
    if (filterCol && !categoricalMetaCols.includes(filterCol)) {
      filterCol = categoricalMetaCols[0] || '';
      filterVal = null;
    }
  });

  function cycleFilterCol() {
    if (categoricalMetaCols.length === 0) return;
    filterCol = cycle(categoricalMetaCols, filterCol || categoricalMetaCols[0]);
    filterVal = null;  // clear the selected value when switching columns
  }

  // Unique values of the active filter column, sorted by descending
  // frequency across the currently-loaded samples (ignoring all other
  // filters — this panel shows the full universe so an operator can see
  // "what's in my data" even when the current filter/drill hides most).
  let filterValueList = $derived.by(() => {
    if (!$samples || !$metadata || !filterCol) return [];
    const counts = new Map();
    for (const s of $samples) {
      const v = $metadata[s.id]?.[filterCol];
      if (v === undefined || v === null || v === '') continue;
      const key = String(v);
      counts.set(key, (counts.get(key) || 0) + 1);
    }
    return [...counts.entries()]
      .sort((a, b) => b[1] - a[1])
      .map(([val, count]) => ({ val, count }));
  });

  function toggleFilterValue(val) {
    filterVal = filterVal === val ? null : val;
  }

  // Initialize metaField to first discovered column once metadata arrives.
  $effect(() => {
    if (!metaField && metaGroup.values.length > 0) metaField = metaGroup.values[0];
  });

  // Next-in-list helper used by the cycling buttons. Wraps to the start.
  function cycle(list, current) {
    if (!list || list.length === 0) return current;
    const i = list.indexOf(current);
    return list[(i + 1) % list.length];
  }

  // Active colorBy passed to PlotlyScatter. 'taxon' defers to the
  // shared drill-down palette via colorMap; metric/metadata modes let
  // PlotlyScatter auto-assign (continuous for numeric, categorical otherwise).
  let activeColorBy = $derived(
    colorMode === 'metric' ? metricField :
    colorMode === 'metadata' ? (metaField || metaGroup.values[0] || 'taxon') :
    'taxon'
  );

  // Get coordinates for current embedding mode
  function getCoords(sampleId) {
    if (!$sampleTsne) return [0, 0];
    const embedData = $sampleTsne[embedMode];
    if (!embedData) return [0, 0];
    return embedData[sampleId] ?? [0, 0];
  }

  // Fixed extents from ALL samples (not filtered) so axes don't jump
  let coordExtents = $derived.by(() => {
    if (!$samples || !$sampleTsne) return null;
    let xMin = Infinity, xMax = -Infinity, yMin = Infinity, yMax = -Infinity;
    for (const s of $samples) {
      const [x, y] = getCoords(s.id);
      if (x < xMin) xMin = x;
      if (x > xMax) xMax = x;
      if (y < yMin) yMin = y;
      if (y > yMax) yMax = y;
    }
    if (!isFinite(xMin)) return null;
    return { xMin, xMax, yMin, yMax };
  });

  // Taxonomy-overlay scatter: one disc per (sample × sub-taxon) at the
  // current drill level. perSample[sid][subTaxon] holds the read count at
  // that sub-taxon; we normalize against sample.read_count so disc areas
  // are proportional to "% of this sample's reads" — same semantics as the
  // /map rings — making sizes directly comparable across samples.
  //
  // Samples with zero discs (none of their reads fell under the drilled
  // filter) disappear entirely, per spec.
  let scatterData = $derived.by(() => {
    if (!$samples || !$sampleTsne) return null;
    const basePass = $samples.filter(passesFilters);
    const cartOn = $cartActive && $cartItems.size > 0;

    // Per-sample fields available as color-by targets — every disc gets
    // the sample-level fields copied on so colorBy='read_count' /
    // 'station' / 'depth_m' etc. all work without the chart component
    // needing to know about the sample→taxon fan-out.
    const fieldsFor = (s) => {
      const m = $metadata?.[s.id];
      return {
        flowcell: s.flowcell,
        barcode: s.barcode,
        read_count: s.read_count,
        total_bases: s.total_bases,
        avg_length: s.avg_length,
        diversity: s.diversity,
        dominant_phylum: s.dominant_phylum,
        dominant_class: s.dominant_class,
        start_time: s.start_time,
        ...(m || {}),
      };
    };

    const out = [];

    if (colorMode === 'taxonomy') {
      // Taxonomy mode: one disc per (sample × sub-taxon) at the current
      // drill level. Samples with no taxa at this level are hidden.
      if (!$activeSubTaxa) return null;
      const perSample = $activeSubTaxa.perSample;
      const isolatedTaxon = $taxNav.isolated;
      for (const s of basePass) {
        if (cartOn && !$cartItems.has(s.id)) continue;
        const [x, y] = getCoords(s.id);
        const subCounts = perSample[s.id];
        if (!subCounts) continue;
        const denom = s.read_count > 0 ? s.read_count : null;
        if (!denom) continue;
        const sampleFields = fieldsFor(s);
        for (const taxon in subCounts) {
          if (isolatedTaxon && taxon !== isolatedTaxon) continue;
          const cnt = subCounts[taxon];
          if (!cnt) continue;
          const fraction = cnt / denom;
          out.push({
            id: s.id,
            taxon,
            ...sampleFields,
            tsne_x: x,
            tsne_y: y,
            size: Math.sqrt(fraction) * 28,
            fraction,
            hover_text:
              `<b>${s.id}</b><br>` +
              `${taxon}: ${(fraction * 100).toFixed(1)}% ` +
              `(${cnt.toLocaleString()} / ${s.read_count.toLocaleString()} reads)`,
          });
        }
      }
    } else {
      // Metric / Metadata mode: one disc per sample. Every sample shows
      // regardless of whether taxonomy data is present for it — that
      // fixes the "missing dots" for samples where GTDB hasn't landed.
      // PlotlyScatter's default "Unknown" bucket handles samples with no
      // value for the current colorBy field.
      for (const s of basePass) {
        if (cartOn && !$cartItems.has(s.id)) continue;
        const [x, y] = getCoords(s.id);
        out.push({
          id: s.id,
          taxon: null,
          ...fieldsFor(s),
          tsne_x: x,
          tsne_y: y,
          size: 14,           // fixed; sizeScale still applies downstream
          fraction: 1,
          hover_text:
            `<b>${s.id}</b><br>` +
            `${(s.read_count ?? 0).toLocaleString()} reads`,
        });
      }
    }

    // Largest first — big discs render underneath so smaller overlays stay
    // visible. (In non-taxonomy mode all discs are equal; the sort is a
    // no-op but harmless.)
    out.sort((a, b) => b.fraction - a.fraction);
    if (typeof window !== 'undefined') {
      const taxa = new Set(out.map(p => p.taxon).filter(t => t));
      const palette = colorMode === 'taxonomy' ? ($activeSubTaxa?.colorMap ?? {}) : {};
      window.__sampleScatterDiag = {
        mode: colorMode,
        level: $taxNav.level, filter: $taxNav.filter, isolated: $taxNav.isolated,
        discs: out.length, uniqueTaxa: taxa.size,
        colorMapSize: Object.keys(palette).length,
      };
    }
    return { points: out };
  });

  // Shared color map with the drill-down sidebar and the /map view.
  let taxColorMap = $derived($activeSubTaxa?.colorMap ?? {});

  // Selected sample detail
  let selectedDetail = $derived.by(() => {
    if (!$selectedSample || !$samples) return null;
    const s = $samples.find(s => s.id === $selectedSample);
    if (!s) return null;
    const m = $metadata?.[s.id];
    return { ...s, ...(m || {}) };
  });

  // Table columns — separate flowcell + barcode
  const tableColumns = [
    { key: 'flowcell', label: 'Flowcell' },
    { key: 'barcode', label: 'Barcode' },
    { key: 'read_count', label: 'Reads', render: v => typeof v === 'number' ? v.toLocaleString() : '-' },
    { key: 'total_bases', label: 'Bases', render: v => typeof v === 'number' ? `${(v/1e6).toFixed(1)}M` : '-' },
    { key: 'avg_length', label: 'Avg Len', render: v => typeof v === 'number' ? Math.round(v).toLocaleString() : '-' },
    { key: 'n_classified', label: 'Classified', render: v => typeof v === 'number' ? v.toLocaleString() : '-' },
    { key: 'dominant_class', label: 'Dominant Class' },
    { key: 'diversity', label: 'Shannon H', render: v => typeof v === 'number' ? v.toFixed(2) : '-' },
  ];

  function handleSelect(ids) {
    lassoIds = ids;
  }

  function handleClick(id) {
    selectedSample.set(id);
  }

  function addLassoToCart() {
    if (lassoIds) {
      for (const id of lassoIds) addToCart(id);
      lassoIds = null;
    }
  }
</script>

<div class="space-y-6">
  <!-- Header stats -->
  {#if filteredStats}
    <div class="grid grid-cols-2 md:grid-cols-4 lg:grid-cols-5 gap-4">
      <StatCard label="Flowcells" value={filteredStats.nFlowcells} color="purple" />
      <StatCard label={filteredStats.isFiltered ? `Barcodes (of ${filteredStats.totalAll})` : 'Barcodes'} value={filteredStats.nSamples} color="cyan" />
      <StatCard label="Total Reads" value={filteredStats.totalReads.toLocaleString()} color="emerald" />
      <StatCard label="Total Bases" value={filteredStats.totalBases ? `${(filteredStats.totalBases / 1e9).toFixed(1)} Gbp` : '-'} color="amber" />
      <StatCard label="In Cart" value={$cartItems.size} color="slate" />
    </div>
  {/if}

  <!-- Controls -->
  <div class="flex items-center gap-3 flex-wrap text-xs">
    <button
      class="px-3 py-1 rounded-md border transition-colors text-center
        {availableModes.length > 1 ? 'border-slate-600 text-slate-400 hover:border-slate-500 cursor-pointer' : 'border-slate-700 text-slate-500 opacity-60 cursor-default'}"
      onclick={availableModes.length > 1 ? cycleEmbed : undefined}
      title={availableModes.length > 1 ? `Cycle layout (${availableModes.map(m => EMBED_LABELS[m] || m).join(' / ')})` : 'Only grid layout available (run sketching for t-SNE)'}
    >
      {EMBED_LABELS[embedMode] || embedMode} &#x25BE;
    </button>

    <!-- Color-mode cycle: Taxonomy (default) / Metric / Metadata.
         Each button is a toggle+cycle: first click switches the mode;
         subsequent clicks while in that mode cycle the inner field. -->
    <div class="flex items-center gap-1">
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center min-w-[5rem]
          {colorMode === 'taxonomy' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        onclick={() => { colorMode = 'taxonomy'; }}
        title="Color by active taxonomy level (shared palette with /map)"
      >Taxonomy</button>
      <button
        class="px-3 py-1 rounded-md border transition-colors text-center min-w-[5rem]
          {colorMode === 'metric' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        onclick={() => {
          if (colorMode === 'metric') metricField = cycle(METRIC_VALUES, metricField);
          else colorMode = 'metric';
        }}
        title={`Click to cycle: ${METRIC_VALUES.map(v => METRIC_LABELS[v]).join(' → ')}`}
      >{colorMode === 'metric' ? METRIC_LABELS[metricField] : 'Metric'} &#x25BE;</button>
      {#if metaGroup.values.length > 0}
        <button
          class="px-3 py-1 rounded-md border transition-colors text-center min-w-[6rem]
            {colorMode === 'metadata' ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
          onclick={() => {
            if (colorMode === 'metadata') metaField = cycle(metaGroup.values, metaField || metaGroup.values[0]);
            else colorMode = 'metadata';
          }}
          title={`Click to cycle: ${metaGroup.labels.join(' → ')}`}
        >{colorMode === 'metadata' ? (metaField || metaGroup.values[0]).replace(/_/g, ' ') : 'Metadata'} &#x25BE;</button>
      {/if}
    </div>

    <!-- Search: free-text match across id/flowcell/barcode/taxonomy/metadata.
         Matches stay colored; non-matches are dimmed by PlotlyScatter. -->
    <div class="flex items-center gap-1">
      <input
        type="search"
        bind:value={searchQuery}
        placeholder="Search id / flowcell / metadata…"
        class="bg-slate-800 border border-slate-600 rounded-md px-2 py-1 text-xs text-slate-200 w-[220px] focus:border-cyan-400 focus:outline-none"
      />
      {#if searchMatchIds}
        <span class="text-[10px] text-slate-500 font-mono">
          {searchMatchIds.size}/{$samples?.length ?? 0}
        </span>
        <button class="text-slate-500 hover:text-slate-200 text-xs" onclick={() => (searchQuery = '')} title="Clear search">✕</button>
      {/if}
    </div>

    <span class="text-slate-500 text-[11px] uppercase tracking-wide">
      {RANK_LABELS[$taxNav.level] ?? $taxNav.level}
      {#if $taxNav.filter} <span class="text-cyan-400 normal-case tracking-normal">· {$taxNav.filter}</span>{/if}
    </span>
    <div class="text-slate-400 flex items-center gap-1">
      Size
      <div class="single-range relative w-20 h-5 flex items-center">
        <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
        <!-- Max is clamped at 6 because raw disc size reaches 28 (sqrt(1)*28
             for a 100%-dominant taxon) and Plotly's scattergl marker
             rasterization wraps around above ~256 px — sizeScale*28 above
             that point visually collapses the biggest discs back to tiny.
             6*28 = 168 stays well clear of the wrap threshold. -->
        <input type="range" min="0.2" max="6" step="0.1" bind:value={sizeScale} />
      </div>
      <span class="text-slate-500 w-8 font-mono">{sizeScale.toFixed(1)}x</span>
    </div>
    <!-- log2 read count filter (MAG-style dual-range) -->
    <div class="text-slate-400 flex items-center gap-1">
      Reads
      <span class="text-slate-500 w-10 text-right font-mono">{readFilterLabel(log2Min)}</span>
      <div class="dual-range relative w-28 h-5 flex items-center">
        <div class="absolute h-1 w-full bg-slate-700 rounded"></div>
        <div class="absolute h-1 bg-cyan-400/40 rounded" style="left: {log2Min / log2Ceil * 100}%; right: {100 - log2Max / log2Ceil * 100}%"></div>
        <input type="range" min="0" max={log2Ceil} step="0.5" bind:value={log2Min}
          oninput={() => { if (log2Min > log2Max - 0.5) log2Min = log2Max - 0.5; }} />
        <input type="range" min="0" max={log2Ceil} step="0.5" bind:value={log2Max}
          oninput={() => { if (log2Max < log2Min + 0.5) log2Max = log2Min + 0.5; }} />
      </div>
      <span class="text-slate-500 w-10 font-mono">{readFilterLabel(log2Max)}</span>
    </div>
    <!-- Date filter -->
    {#if dateRange.min}
      <div class="text-slate-400 flex items-center gap-1">
        Date
        <input type="date" bind:value={dateMin}
          min={dateRange.min} max={dateRange.max}
          class="bg-slate-800 border border-slate-600 rounded-md px-1.5 py-1 text-xs text-slate-200 w-[130px] focus:border-cyan-400 focus:outline-none" />
        <span class="text-slate-500">–</span>
        <input type="date" bind:value={dateMax}
          min={dateRange.min} max={dateRange.max}
          class="bg-slate-800 border border-slate-600 rounded-md px-1.5 py-1 text-xs text-slate-200 w-[130px] focus:border-cyan-400 focus:outline-none" />
      </div>
    {/if}
    {#if lassoIds}
      <button
        class="px-3 py-1 rounded-md border border-cyan-400 bg-cyan-400/10 text-cyan-400 hover:bg-cyan-400/20 transition-colors"
        onclick={addLassoToCart}
      >
        Add {lassoIds.length} to Cart
      </button>
    {/if}
  </div>

  <!-- Main layout: sidebar (taxonomy drill-down in taxonomy mode, metadata
       filter otherwise) + scatter + detail panel. -->
  <div class="flex gap-6">
    {#if colorMode === 'taxonomy'}
      <TaxonomyDrillNav heightClass="h-[600px]" />
    {:else}
      <!-- Metadata filter sidebar: cycle header picks a categorical
           column, list below shows unique values sorted by descending
           frequency. Click a value to filter the scatter to matching
           samples; click the active value again to clear. -->
      <div class="w-56 shrink-0 h-[600px] flex flex-col bg-slate-800 rounded-lg border border-slate-700 p-3 text-xs">
        {#if categoricalMetaCols.length === 0}
          <p class="text-slate-500 italic text-center mt-4">
            Upload a metadata TSV with categorical columns to filter here.
          </p>
        {:else}
          <div class="flex items-center justify-between gap-2 mb-2">
            <button
              class="flex-1 px-2 py-1 rounded border border-slate-600 text-slate-200 hover:border-slate-500 transition-colors text-left font-semibold truncate"
              onclick={cycleFilterCol}
              title={`Cycle: ${categoricalMetaCols.join(' → ')}`}
            >{(filterCol || '').replace(/_/g, ' ') || '—'} &#x25BE;</button>
            {#if filterVal !== null}
              <button
                class="px-2 py-1 rounded border border-slate-600 text-slate-400 hover:text-rose-400 hover:border-rose-400/40 transition-colors"
                onclick={() => (filterVal = null)}
                title="Clear active filter"
              >✕</button>
            {/if}
          </div>
          {#if filterValueList.length === 0}
            <p class="text-slate-500 italic text-center mt-4">No values for this column.</p>
          {:else}
            <div class="text-[10px] uppercase tracking-wide text-slate-500 mb-1">
              {filterValueList.length} value{filterValueList.length === 1 ? '' : 's'}
            </div>
            <ul class="flex-1 overflow-y-auto divide-y divide-slate-700/50 -mx-1">
              {#each filterValueList as { val, count }}
                <li>
                  <button
                    class="w-full text-left px-2 py-1 rounded transition-colors flex items-baseline justify-between gap-2
                      {filterVal === val ? 'bg-cyan-400/15 text-cyan-300' : 'text-slate-300 hover:bg-slate-700/50'}"
                    onclick={() => toggleFilterValue(val)}
                  >
                    <span class="truncate" title={val}>{val}</span>
                    <span class="text-[10px] text-slate-500 font-mono shrink-0">{count}</span>
                  </button>
                </li>
              {/each}
            </ul>
          {/if}
        {/if}
      </div>
    {/if}
    <div class="flex-1 h-[600px] flex flex-col">
      {#if scatterData && scatterData.points.length > 0}
        <PlotlyScatter
          data={scatterData}
          colorBy={activeColorBy}
          sizeBy="raw"
          {sizeScale}
          mode="tsne"
          colorMap={colorMode === 'taxonomy' ? taxColorMap : {}}
          {coordExtents}
          {searchMatchIds}
          onselect={handleSelect}
          onclick={handleClick}
          exportName={`danaseq_sample_tsne_${$taxNav.level}${$taxNav.filter ? '_' + $taxNav.filter : ''}${colorMode !== 'taxonomy' ? '_' + activeColorBy : ''}`}
        />
      {:else if !$activeSubTaxa}
        <div class="flex items-center justify-center h-full text-slate-500">Loading taxonomy…</div>
      {:else}
        <div class="flex items-center justify-center h-full text-slate-500">No samples under the current drill-down filter.</div>
      {/if}
    </div>

    <!-- Detail panel -->
    <div class="w-72 shrink-0 bg-slate-800 rounded-lg border border-slate-700 p-4 space-y-3 h-fit max-h-[600px] overflow-y-auto">
      {#if selectedDetail}
        <div class="flex items-center justify-between gap-2">
          <h3 class="text-sm font-semibold text-cyan-400 font-mono truncate">{selectedDetail.id}</h3>
          <button
            class="text-xs px-2 py-1 rounded border transition-colors shrink-0
              {$cartItems.has(selectedDetail.id)
                ? 'bg-cyan-400/20 text-cyan-400 border-cyan-400/40'
                : 'text-slate-400 border-slate-600 hover:text-cyan-400 hover:border-cyan-400/40'}"
            onclick={() => toggleCart(selectedDetail.id)}
          >
            {$cartItems.has(selectedDetail.id) ? 'In Cart' : '+ Cart'}
          </button>
        </div>
        <div class="grid grid-cols-2 gap-2 text-xs">
          <div class="text-slate-400">Flowcell</div><div class="text-slate-200 font-mono">{selectedDetail.flowcell ?? '-'}</div>
          <div class="text-slate-400">Barcode</div><div class="text-slate-200 font-mono">{selectedDetail.barcode ?? '-'}</div>
          <div class="text-slate-400">Reads</div><div class="text-slate-200 font-mono">{selectedDetail.read_count?.toLocaleString() ?? '-'}</div>
          <div class="text-slate-400">Bases</div><div class="text-slate-200 font-mono">{selectedDetail.total_bases ? `${(selectedDetail.total_bases/1e6).toFixed(1)}M` : '-'}</div>
          <div class="text-slate-400">Avg Length</div><div class="text-slate-200 font-mono">{selectedDetail.avg_length ? Math.round(selectedDetail.avg_length).toLocaleString() : '-'}</div>
          <div class="text-slate-400">Diversity</div><div class="text-slate-200 font-mono">{selectedDetail.diversity?.toFixed(2) ?? '-'}</div>
          <div class="text-slate-400">Phylum</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.dominant_phylum ?? '-'}</div>
          <div class="text-slate-400">Class</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.dominant_class ?? '-'}</div>
          {#if selectedDetail.start_time}
            <div class="text-slate-400">Run Date</div><div class="text-slate-200 font-mono text-[11px]">{selectedDetail.start_time.slice(0, 10)}</div>
          {/if}
        </div>
        {#if selectedDetail.lat != null}
          <div class="border-t border-slate-700 pt-2">
            <div class="text-xs text-slate-400 mb-1">Location</div>
            <div class="grid grid-cols-2 gap-2 text-xs">
              <div class="text-slate-400">Lat</div><div class="text-slate-200 font-mono">{selectedDetail.lat}</div>
              <div class="text-slate-400">Lon</div><div class="text-slate-200 font-mono">{selectedDetail.lon}</div>
              {#if selectedDetail.station}<div class="text-slate-400">Station</div><div class="text-slate-200 font-mono">{selectedDetail.station}</div>{/if}
              {#if selectedDetail.depth_m != null}<div class="text-slate-400">Depth</div><div class="text-slate-200 font-mono">{selectedDetail.depth_m}m</div>{/if}
            </div>
          </div>
        {/if}
      {:else}
        <p class="text-slate-500 text-xs">Click a sample point or table row to view details.</p>
      {/if}
    </div>
  </div>

  <!-- Sample table -->
  {#if $samples}
    {@const filteredSamples = $samples.filter(passesFilters)}
    {@const tableRows = lassoIds
      ? filteredSamples.filter(s => lassoIds.includes(s.id))
      : $cartActive && $cartItems.size > 0
        ? filteredSamples.filter(s => $cartItems.has(s.id))
        : searchMatchIds
          ? filteredSamples.filter(s => searchMatchIds.has(s.id))
          : filteredSamples}
    <div>
      <h3 class="text-sm font-semibold text-slate-300 mb-2">
        {lassoIds ? `Selected (${lassoIds.length})`
          : $cartActive && $cartItems.size > 0 ? `Cart (${$cartItems.size})`
          : searchMatchIds ? `Matches (${tableRows.length}/${filteredSamples.length})`
          : `Barcodes (${filteredSamples.length}/${$samples.length})`}
      </h3>
      <DataTable
        columns={tableColumns}
        rows={tableRows}
        onRowClick={(row) => { selectedSample.set(row.id); }}
        selectedId={$selectedSample}
        idKey="id"
        maxHeight="300px"
        exportFilename="sample_stats"
        actionLabel={(row) => $cartItems.has(row.id) ? 'In Cart' : '+ Cart'}
        actionFn={(row) => toggleCart(row.id)}
        actionStyle={(row) => $cartItems.has(row.id)
          ? 'text-[10px] px-2 py-0.5 rounded border bg-cyan-400/20 text-cyan-400 border-cyan-400/40 transition-colors'
          : 'text-[10px] px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-cyan-400 hover:border-cyan-400/40 transition-colors'}
      />
    </div>
  {/if}
</div>

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
