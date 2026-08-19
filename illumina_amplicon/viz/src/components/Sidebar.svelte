<script>
  import { store, GROUP_HEX, buildTaxColorMap, getEffectiveColorLevel, findTaxonLevel, listAssays, assayHeading } from '../stores/data.svelte.js';
  import AutocompleteInput from './AutocompleteInput.svelte';

  let { activeTab = 'samples', filters = $bindable({}), open = false,
        onClose = null } = $props();

  function generateClusterColor(id, total) {
    const hue = ((id - 1) * 137.508) % 360;
    return `hsl(${hue}, 70%, 55%)`;
  }

  // When user picks from autocomplete, navigate to that taxon's level
  function navigateToTaxon(name) {
    if (!name || name === 'unclassified' || filters.colorMode !== 'taxonomy') return;
    const level = findTaxonLevel(name);
    if (level) {
      pushNav();
      filters.colorByLevel = level;
      filters.taxonFilter = name;
      filters.taxonFilterLevel = level;
    }
  }

  // A hand-typed filter is a substring of the whole lineage, so it has no rank.
  function clearFilterLevel() {
    filters.taxonFilterLevel = '';
  }

  // ── Point-scale sliders ───────────────────────────────────────────────────
  // The track runs 1x .. the scale at which the biggest marker hits the 100px
  // renderer ceiling, log-spaced so the small end (where the plot actually
  // changes) gets as much travel as the large end.

  /** Top of the track. Guarded above 1 so log(topScale) can't be 0. */
  function topScale(max) {
    return Math.max(1.0001, max || 1);
  }

  /** Scale actually in force — the stored value can outrun a newly-lowered top. */
  function effScale(value, max) {
    return Math.min(value ?? 1, topScale(max));
  }

  function scaleToSlider(value, max) {
    const top = topScale(max);
    return Math.round(Math.log(effScale(value, max)) / Math.log(top) * 100);
  }

  function sliderToScale(pos, max) {
    return Math.pow(topScale(max), +pos / 100);
  }

  function fmtScale(v) {
    return v >= 10 ? Math.round(v) : v.toFixed(1);
  }

  function pushNav() {
    filters.navStack = [
      ...(filters.navStack || []),
      { level: filters.colorByLevel, filter: filters.taxonFilter, filterLevel: filters.taxonFilterLevel },
    ];
  }

  let assays = $derived.by(() => listAssays());

  function toggleAssay(key) {
    const next = new Set(filters.assays || []);
    if (next.has(key)) next.delete(key); else next.add(key);
    filters.assays = next;
  }

  // Collapsible sections
  let sections = $state({
    assay: true,
    taxonomy: true,
    samples: true,
    network: true,
    phylogeny: true,
    heatmap: true,
  });

  function toggle(section) {
    sections[section] = !sections[section];
  }

  // Taxonomy levels from data
  let taxLevels = $derived(
    store.taxonomy[Object.keys(store.taxonomy)[0]]?.levels || []
  );

  // Autocomplete candidates
  let taxCandidates = $derived.by(() => {
    const db = Object.keys(store.taxonomy)[0];
    if (!db || !store.taxonomy[db]?.assignments) return [];
    const terms = new Set();
    for (const asvId in store.taxonomy[db].assignments) {
      const vals = store.taxonomy[db].assignments[asvId];
      for (const v of vals) {
        if (v && v !== 'unclassified') terms.add(v);
      }
    }
    return [...terms].sort();
  });

  let sampleCandidates = $derived(
    store.samples.map(s => s.id).filter(Boolean).sort()
  );
</script>

<!-- A bottom sheet below md and a left column at md and up.
     A sheet rather than a side drawer: filters are a thing you adjust while
     watching the plot change, and a side drawer covers the plot. Docked to the
     bottom it takes half the height, leaves the viz visible above it, and sits
     where a thumb already is. -->
<aside class="fixed inset-x-0 bottom-0 z-40 flex max-h-[55vh] transform flex-col
              overflow-y-auto overscroll-contain rounded-t-xl border-t border-slate-700
              bg-slate-900 shadow-2xl transition-transform duration-200
              md:static md:inset-auto md:z-auto md:max-h-none md:w-64 md:shrink-0
              md:translate-y-0 md:rounded-none md:border-r md:border-t-0
              md:border-slate-800 md:bg-slate-900/60 md:shadow-none
              {open ? 'translate-y-0' : 'translate-y-full'}">
  {#if onClose}
    <!-- Sticky so the handle stays reachable however far the filters scroll. -->
    <div class="sticky top-0 z-10 flex items-center justify-between border-b border-slate-800
                bg-slate-900 px-3 py-2 md:hidden">
      <span class="text-xs font-semibold uppercase tracking-wider text-slate-400">Filters</span>
      <button type="button" onclick={onClose} aria-label="Close filters"
              class="rounded px-3 py-1 text-lg leading-none text-slate-400 hover:text-white">&times;</button>
    </div>
  {/if}
  <div class="p-3 border-b border-slate-800">
    {#if store.runInfo?.title}
      <p class="truncate text-sm font-semibold text-slate-200" title={store.runInfo.title}>
        {store.runInfo.title}
      </p>
    {/if}
    {#if store.runInfo?.bioproject}
      <p class="text-[11px] text-slate-400">
        <a href="https://www.ncbi.nlm.nih.gov/bioproject/{store.runInfo.bioproject}"
           target="_blank" rel="noopener" class="hover:text-cyan-400">{store.runInfo.bioproject}</a>
        {#if store.runInfo.portal_url}
          <span class="text-slate-600">·</span>
          <a href={store.runInfo.portal_url} target="_blank" rel="noopener"
             class="hover:text-cyan-400">run {store.runInfo.slug}</a>
        {/if}
      </p>
    {/if}
    <p class="text-xs text-slate-500">{store.samples.length} samples | {store.asvs.length} ASVs</p>
  </div>

  <!-- ══ Taxonomy Filters (shared) ══ -->
  <div class="border-b border-slate-800">
    <button class="flex w-full items-center justify-between px-3 py-2 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200" onclick={() => toggle('taxonomy')}>
      Taxonomy
      <span class="text-[10px]">{sections.taxonomy ? '▾' : '▸'}</span>
    </button>
    {#if sections.taxonomy}
      <div class="space-y-3 px-3 pb-3">
        <AutocompleteInput
          bind:value={filters.taxonFilter}
          label="Filter"
          placeholder="e.g. Proteobacteria"
          candidates={taxCandidates}
          onPick={navigateToTaxon}
          onType={clearFilterLevel}
        />

        <label class="block">
          <span class="text-xs text-slate-400">Color by</span>
          <select bind:value={filters.colorMode} class="mt-1 w-full rounded border border-slate-700 bg-slate-800 px-2 py-1 text-sm text-slate-200"
            onchange={(e) => {
              if (e.target.value !== 'taxonomy') {
                // Save taxonomy nav state
                filters._savedTaxNav = {
                  colorByLevel: filters.colorByLevel,
                  navStack: filters.navStack,
                  taxonFilter: filters.taxonFilter,
                  taxonFilterLevel: filters.taxonFilterLevel,
                };
              } else {
                // Restore taxonomy nav state
                if (filters._savedTaxNav) {
                  filters.colorByLevel = filters._savedTaxNav.colorByLevel;
                  filters.navStack = filters._savedTaxNav.navStack;
                  filters.taxonFilter = filters._savedTaxNav.taxonFilter;
                  filters.taxonFilterLevel = filters._savedTaxNav.taxonFilterLevel ?? '';
                  filters._savedTaxNav = null;
                }
              }
            }}>
            <option value="taxonomy">Taxonomy</option>
            <option value="group">Group (broad)</option>
            {#if store.heatmap?.sampleClusters || store.heatmap?.asvClusters}
              <option value="cluster">Cluster</option>
            {/if}
          </select>
        </label>

        {#if filters.colorMode === 'cluster'}
          {#if store.heatmap?.sampleClusters && (activeTab === 'samples' || activeTab === 'heatmap')}
            <label class="block">
              <span class="text-xs text-slate-400">Sample clusters (k): {filters.sampleClusterK}</span>
              <input type="range" min="2" max={Math.max(2, ...Object.keys(store.heatmap.sampleClusters).map(Number))}
                step="1" bind:value={filters.sampleClusterK} class="mt-1 w-full accent-blue-500" />
            </label>
          {/if}
          {#if store.heatmap?.asvClusters && (activeTab === 'network' || activeTab === 'phylogeny' || activeTab === 'heatmap')}
            <label class="block">
              <span class="text-xs text-slate-400">ASV clusters (k): {filters.asvClusterK}</span>
              <input type="range" min="2" max={Math.max(2, ...Object.keys(store.heatmap.asvClusters).map(Number))}
                step="1" bind:value={filters.asvClusterK} class="mt-1 w-full accent-blue-500" />
            </label>
          {/if}
        {:else if filters.colorMode === 'group'}
          <fieldset class="space-y-1">
            <legend class="text-xs text-slate-400">Groups</legend>
            {#each Object.keys(filters.groupFlags || {}) as group}
              <label class="flex items-center gap-2 text-sm capitalize">
                <input type="checkbox" bind:checked={filters.groupFlags[group]} class="accent-blue-500" />
                <span class="inline-block h-2.5 w-2.5 rounded-full" style="background:{GROUP_HEX[group]}"></span>
                {group}
              </label>
            {/each}
          </fieldset>
        {:else}
          {@const effectiveLevel = getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter, filters.taxonFilterLevel)}
          {@const taxColors = buildTaxColorMap(effectiveLevel, filters.taxonFilter, filters.taxonFilterLevel)}
          <div class="space-y-0.5 max-h-48 overflow-y-auto">
            <!-- The number in each row means different things at different
                 levels — ASVs at a rank, reads for a single ASV — so the column
                 is labelled rather than left to be inferred. -->
            <div class="mb-1 flex items-baseline gap-2 text-[10px] text-slate-500">
              <span>{effectiveLevel === '_asv' ? 'ASV' : effectiveLevel} ({taxColors.ranked.length})</span>
              <span class="ml-auto">{effectiveLevel === '_asv' ? 'reads' : 'ASVs'}</span>
            </div>
            {#if filters.navStack?.length > 0}
              {@const prev = filters.navStack[filters.navStack.length - 1]}
              <button
                class="flex items-center gap-1.5 w-full text-left text-xs hover:bg-slate-800 rounded px-1 py-1 text-cyan-400 border-b border-slate-700 mb-1"
                onclick={() => {
                  const stack = [...filters.navStack];
                  const popped = stack.pop();
                  filters.navStack = stack;
                  filters.colorByLevel = popped.level;
                  filters.taxonFilter = popped.filter;
                  filters.taxonFilterLevel = popped.filterLevel ?? '';
                }}
              >
                &#x25B4; Up{prev.filter ? ` to ${prev.filter}` : ` to ${prev.level}`}
              </button>
            {:else if filters.taxonFilter}
              <button
                class="flex items-center gap-1.5 w-full text-left text-xs hover:bg-slate-800 rounded px-1 py-1 text-cyan-400 border-b border-slate-700 mb-1"
                onclick={() => {
                  filters.taxonFilter = '';
                  filters.taxonFilterLevel = '';
                }}
              >
                &#x25B4; Up to {filters.colorByLevel}
              </button>
            {/if}
            {#each taxColors.ranked as item}
              <button
                class="flex items-center gap-1.5 w-full text-left text-xs hover:bg-slate-800 rounded px-1 py-0.5"
                onclick={() => {
                  // Read once up front: {@const} is a reactive derived, so after
                  // the assignments below it re-evaluates to the level we just
                  // drilled INTO rather than the one this row belongs to.
                  const level = effectiveLevel;
                  if (level === '_asv') return;
                  if (item.name === 'unclassified') return;
                  pushNav();
                  filters.colorByLevel = level;
                  filters.taxonFilter = item.name;
                  filters.taxonFilterLevel = level;
                }}
              >
                <span class="inline-block h-2 w-2 rounded-full shrink-0" style="background:{item.color}"></span>
                <span class="truncate">{item.name}</span>
                <span class="ml-auto shrink-0 tabular-nums text-slate-500">{(item.count ?? 0).toLocaleString()}</span>
              </button>
            {/each}
          </div>
        {/if}
      </div>
    {/if}
  </div>

  <!-- ══ Assay / primer set (shared) ══
       One assay is worth stating but not worth a checkbox: selecting the only
       option cannot change anything. Filtering appears from two up. -->
  {#if assays.length === 1}
    <div class="border-b border-slate-800 px-3 py-2">
      <p class="text-[10px] uppercase tracking-wider text-slate-500">Assay</p>
      <p class="text-xs text-slate-300">{assayHeading(assays[0].gene, assays[0].region, assays[0].fwd, assays[0].rev, assays[0].place)}</p>
      {#if assays[0].span}
        <p class="text-[10px] text-slate-400">{assays[0].span}</p>
      {/if}
      {#if assays[0].fwd || assays[0].rev}
        <p class="break-all text-[10px] text-slate-500">{[assays[0].fwd, assays[0].rev].filter(Boolean).join(' / ')}</p>
      {:else if !assays[0].span}
        <p class="text-[10px] text-slate-500">no primers recorded</p>
      {/if}
      <p class="text-[10px] text-slate-500">
        all {assays[0].n} sample{assays[0].n === 1 ? '' : 's'}{assays[0].meanMatch != null ? ` · ${(assays[0].meanMatch * 100).toFixed(1)}% matched` : ''}
      </p>
    </div>
  {:else if assays.length > 1}
    <div class="border-b border-slate-800">
      <button class="flex w-full items-center justify-between px-3 py-2 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200" onclick={() => toggle('assay')}>
        Assay
        <span class="text-[10px]">{sections.assay ? '▾' : '▸'}</span>
      </button>
      {#if sections.assay}
        <div class="space-y-2 px-3 pb-3">
          <p class="text-[10px] text-slate-500">
            {assays.length} assay{assays.length === 1 ? '' : 's'}
            {#if filters.assays?.size}&mdash; {filters.assays.size} selected{:else}&mdash; all shown{/if}
          </p>
          {#each assays as a}
            <label class="flex items-start gap-2 text-xs cursor-pointer hover:bg-slate-800/50 rounded px-1 py-1">
              <input type="checkbox" class="mt-0.5 accent-blue-500"
                checked={filters.assays?.has(a.key) ?? false}
                onchange={() => toggleAssay(a.key)} />
              <span class="min-w-0 flex-1">
                <span class="block truncate text-slate-200">{assayHeading(a.gene, a.region, a.fwd, a.rev, a.place)}</span>
                {#if a.span}
                  <span class="block truncate text-[10px] text-slate-400">{a.span}</span>
                {/if}
                <!-- Where it sits on the gene and what was in the tube are
                     different facts, and a run can have one without the other:
                     a pre-trimmed end has a location and no primer, a
                     functional-gene assay has a primer and no location. -->
                {#if a.fwd || a.rev}
                  <span class="block truncate text-[10px] text-slate-500">
                    {[a.fwd, a.rev].filter(Boolean).join(' / ')}
                  </span>
                {:else if !a.span}
                  <span class="block text-[10px] text-slate-500">no primers recorded</span>
                {/if}
                <span class="block text-[10px] text-slate-500">
                  {a.n} sample{a.n === 1 ? '' : 's'}{a.meanMatch != null ? ` · ${(a.meanMatch * 100).toFixed(1)}% matched` : ''}
                </span>
              </span>
            </label>
          {/each}
          {#if filters.assays?.size}
            <button class="text-[11px] text-cyan-400 hover:text-cyan-300"
              onclick={() => { filters.assays = new Set(); }}>Clear assay filter</button>
          {/if}
        </div>
      {/if}
    </div>
  {/if}

  <!-- ══ Sample Controls ══ -->
  {#if activeTab === 'samples' || activeTab === 'network'}
    <div class="border-b border-slate-800">
      <button class="flex w-full items-center justify-between px-3 py-2 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200" onclick={() => toggle('samples')}>
        Sample Filters
        <span class="text-[10px]">{sections.samples ? '▾' : '▸'}</span>
      </button>
      {#if sections.samples}
        <div class="space-y-3 px-3 pb-3">
          <label class="block">
            <span class="text-xs text-slate-400">Min reads: {(filters.minReads || 0).toLocaleString()}</span>
            <input type="range" min="0" max="100" step="1"
              value={Math.round(Math.log2((filters.minReads || 0) + 1) / Math.log2(Math.max(...store.samples.map(s => s.total_reads || 0), 1) + 1) * 100)}
              oninput={(e) => { const mx = Math.max(...store.samples.map(s => s.total_reads || 0), 1); filters.minReads = Math.round(Math.pow(2, +e.target.value / 100 * Math.log2(mx + 1)) - 1); }}
              class="mt-1 w-full accent-blue-500" />
          </label>

          <AutocompleteInput
            bind:value={filters.sampleFilter}
            label="Sample filter"
            placeholder="e.g. Plate1"
            candidates={sampleCandidates}
          />

          {#if activeTab === 'samples'}
            <label class="block">
              <span class="text-xs text-slate-400">Point scale: {fmtScale(effScale(filters.pointScale, store.maxPointScale))}x</span>
              <input type="range" min="0" max="100" step="1"
                value={scaleToSlider(filters.pointScale, store.maxPointScale)}
                oninput={(e) => { filters.pointScale = sliderToScale(e.target.value, store.maxPointScale); }}
                class="mt-1 w-full accent-blue-500" />
            </label>
          {/if}
        </div>
      {/if}
    </div>
  {/if}

  <!-- ══ Network Controls ══ -->
  {#if activeTab === 'network' || activeTab === 'samples'}
    <div class="border-b border-slate-800">
      <button class="flex w-full items-center justify-between px-3 py-2 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200" onclick={() => toggle('network')}>
        Network Filters
        <span class="text-[10px]">{sections.network ? '▾' : '▸'}</span>
      </button>
      {#if sections.network}
        <div class="space-y-3 px-3 pb-3">
          <label class="block">
            <span class="text-xs text-slate-400">Min prevalence: {filters.minPrevalence || 0}</span>
            <input type="range" min="0" max="100" step="1"
              value={Math.round(Math.log2((filters.minPrevalence || 0) + 1) / Math.log2(Math.max(...store.asvs.map(a => a.n_samples || 0), 1) + 1) * 100)}
              oninput={(e) => { const mx = Math.max(...store.asvs.map(a => a.n_samples || 0), 1); filters.minPrevalence = Math.round(Math.pow(2, +e.target.value / 100 * Math.log2(mx + 1)) - 1); }}
              class="mt-1 w-full accent-blue-500" />
          </label>

          {#if store.network?.edges?.length > 0}
            <label class="block">
              <span class="text-xs text-slate-400">Correlation threshold: {(filters.corrThreshold || 0.3).toFixed(2)}</span>
              <input type="range" min="0" max="1" step="0.01" bind:value={filters.corrThreshold} class="mt-1 w-full accent-blue-500" />
            </label>

            <label class="flex items-center gap-2 text-sm">
              <input type="checkbox" bind:checked={filters.showEdges} class="accent-blue-500" />
              Show edges
            </label>
          {/if}

          {#if activeTab === 'network'}
            <label class="block">
              <span class="text-xs text-slate-400">ASV point scale: {fmtScale(effScale(filters.networkPointScale, store.maxNetworkPointScale))}x</span>
              <input type="range" min="0" max="100" step="1"
                value={scaleToSlider(filters.networkPointScale, store.maxNetworkPointScale)}
                oninput={(e) => { filters.networkPointScale = sliderToScale(e.target.value, store.maxNetworkPointScale); }}
                class="mt-1 w-full accent-blue-500" />
            </label>
          {/if}
        </div>
      {/if}
    </div>
  {/if}

  <!-- ══ Phylogeny Controls ══ -->
  {#if activeTab === 'phylogeny'}
    <div class="border-b border-slate-800">
      <button class="flex w-full items-center justify-between px-3 py-2 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200" onclick={() => toggle('phylogeny')}>
        Phylogeny
        <span class="text-[10px]">{sections.phylogeny ? '▾' : '▸'}</span>
      </button>
      {#if sections.phylogeny}
        <div class="space-y-3 px-3 pb-3">
          <label class="block">
            <span class="text-xs text-slate-400">Layout</span>
            <select bind:value={filters.treeLayout} class="mt-1 w-full rounded border border-slate-700 bg-slate-800 px-2 py-1 text-sm text-slate-200">
              <option value="rc">Rectangular</option>
              <option value="rd">Radial</option>
              <option value="cr">Circular</option>
              <option value="dg">Diagonal</option>
              <option value="hr">Hierarchical</option>
            </select>
          </label>

          <fieldset class="space-y-1">
            <legend class="text-xs text-slate-400">Tip labels</legend>
            <label class="flex items-center gap-2 text-xs">
              <input type="checkbox" checked={filters.treeLabelLevels?.includes('id')}
                onchange={(e) => {
                  const levels = [...(filters.treeLabelLevels || [])];
                  if (e.target.checked) levels.push('id');
                  else levels.splice(levels.indexOf('id'), 1);
                  filters.treeLabelLevels = levels;
                }} class="accent-blue-500" />
              ASV ID
            </label>
            {#each taxLevels as level}
              <label class="flex items-center gap-2 text-xs">
                <input type="checkbox" checked={filters.treeLabelLevels?.includes(level)}
                  onchange={(e) => {
                    const levels = [...(filters.treeLabelLevels || [])];
                    if (e.target.checked) levels.push(level);
                    else levels.splice(levels.indexOf(level), 1);
                    filters.treeLabelLevels = levels;
                  }} class="accent-blue-500" />
                {level}
              </label>
            {/each}
            <label class="flex items-center gap-2 text-xs">
              <input type="checkbox" checked={filters.treeLabelLevels?.includes('bootstrap')}
                onchange={(e) => {
                  const levels = [...(filters.treeLabelLevels || [])];
                  if (e.target.checked) levels.push('bootstrap');
                  else levels.splice(levels.indexOf('bootstrap'), 1);
                  filters.treeLabelLevels = levels;
                }} class="accent-blue-500" />
              Bootstrap
            </label>
          </fieldset>

          <label class="block">
            <span class="text-xs text-slate-400">Min bootstrap: {filters.treeMinBootstrap || 0}%</span>
            <input type="range" min="0" max="100" step="5" bind:value={filters.treeMinBootstrap} class="mt-1 w-full accent-blue-500" />
          </label>

          <label class="flex items-center gap-2 text-xs">
            <input type="checkbox" bind:checked={filters.treePrune} class="accent-blue-500" />
            Prune filtered tips
          </label>
        </div>
      {/if}
    </div>
  {/if}

  <!-- ══ Cluster Sliders (tables tab) ══ -->
  {#if activeTab === 'tables' && filters.colorMode === 'cluster'}
    <div class="border-b border-slate-800">
      <div class="space-y-3 px-3 py-3">
        {#if store.heatmap?.sampleClusters}
          <label class="block">
            <span class="text-xs text-slate-400">Sample clusters (k): {filters.sampleClusterK}</span>
            <input type="range" min="2" max={Math.max(2, ...Object.keys(store.heatmap.sampleClusters).map(Number))}
              step="1" bind:value={filters.sampleClusterK} class="mt-1 w-full accent-blue-500" />
          </label>
        {/if}
        {#if store.heatmap?.asvClusters}
          <label class="block">
            <span class="text-xs text-slate-400">ASV clusters (k): {filters.asvClusterK}</span>
            <input type="range" min="2" max={Math.max(2, ...Object.keys(store.heatmap.asvClusters).map(Number))}
              step="1" bind:value={filters.asvClusterK} class="mt-1 w-full accent-blue-500" />
          </label>
        {/if}
      </div>
    </div>
  {/if}

  <!-- ══ Heatmap Controls ══ -->
  {#if activeTab === 'heatmap'}
    <div class="border-b border-slate-800">
      <button class="flex w-full items-center justify-between px-3 py-2 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200" onclick={() => toggle('heatmap')}>
        Heatmap
        <span class="text-[10px]">{sections.heatmap ? '▾' : '▸'}</span>
      </button>
      {#if sections.heatmap}
        <div class="space-y-3 px-3 pb-3">
          <label class="block">
            <span class="text-xs text-slate-400">Min max(RA): {(filters.heatmapMinMaxRA ?? 1).toFixed(1)}%</span>
            <input type="range" min="0" max="100" step="1"
              value={Math.round(Math.log2((filters.heatmapMinMaxRA ?? 1) + 1) / Math.log2(101) * 100)}
              oninput={(e) => { filters.heatmapMinMaxRA = Math.round(10 * (Math.pow(2, +e.target.value / 100 * Math.log2(101)) - 1)) / 10; }}
              class="mt-1 w-full accent-blue-500" />
          </label>

          <label class="block">
            <span class="text-xs text-slate-400">Cell size: {filters.heatmapCellSize}px</span>
            <input type="range" min="1" max="20" step="1" bind:value={filters.heatmapCellSize} class="mt-1 w-full accent-blue-500" />
          </label>

          <label class="block">
            <span class="text-xs text-slate-400">ASV ordering</span>
            <select bind:value={filters.heatmapAsvTree} class="mt-1 w-full rounded border border-slate-700 bg-slate-800 px-2 py-1 text-sm text-slate-200">
              <option value="ward">Ward Clustering</option>
              {#if store.treeNewick}
                <option value="phylogeny">Phylogeny (NJ)</option>
              {/if}
            </select>
          </label>
        </div>
      {/if}
    </div>
  {/if}

  <!-- ══ Selection info ══ -->
  <div class="mt-auto border-t border-slate-800 p-3">
    {#if store.selectedSample != null}
      <p class="text-xs text-slate-400">Selected: {store.samples[store.selectedSample]?.id || ''}</p>
    {:else if store.selectedAsv != null}
      <p class="text-xs text-slate-400">Selected: {store.asvs[store.selectedAsv]?.id || ''}</p>
    {:else}
      <p class="text-xs text-slate-500">Click to select, Shift+drag to lasso</p>
      <p class="text-xs text-slate-500">Shift+double-click to deselect</p>
    {/if}
  </div>
</aside>
