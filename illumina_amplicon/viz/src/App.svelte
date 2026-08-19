<script>
  import { onMount } from 'svelte';
  import NavBar from './components/NavBar.svelte';
  import Sidebar from './components/Sidebar.svelte';
  import SampleView from './views/SampleView.svelte';
  import NetworkView from './views/NetworkView.svelte';
  import PhyloTreeView from './views/PhyloTreeView.svelte';
  import HeatmapView from './views/HeatmapView.svelte';
  import TablesView from './views/TablesView.svelte';
  import ProvenanceView from './views/ProvenanceView.svelte';
  import { store, loadData, asvIdsForAssays } from './stores/data.svelte.js';

  let activeTab = $state('samples');
  // Drawer state for narrow screens; ignored from md up, where the sidebar is a
  // column and always shown.
  let sidebarOpen = $state(false);

  // Shared filter state
  let filters = $state({
    // Taxonomy (shared across all views)
    taxonFilter: '',
    taxonFilterLevel: '',      // rank the filter names, when picked (not hand-typed)
    colorMode: 'taxonomy',     // 'taxonomy', 'group', 'cluster'
    sampleClusterK: 4,
    asvClusterK: 4,
    colorByLevel: 'Domain',    // current taxonomy nav level
    navStack: [],              // breadcrumb trail for drill-down navigation
    groupFlags: {
      prokaryote: true,
      eukaryote: true,
      chloroplast: true,
      mitochondria: true,
      unknown: true,
    },
    // Assay (primer set) — empty Set means "all assays", not "none"
    assays: new Set(),
    // Samples
    minReads: 0,
    sampleFilter: '',
    showOverlay: true,
    pointScale: 3,
    // Network
    minPrevalence: 0,
    corrThreshold: 0.3,
    showEdges: true,
    networkPointScale: 10,
    // Selections (shared across views)
    lassoSampleIds: new Set(),
    lassoAsvIds: new Set(),

    // Per-view zoom persistence
    sampleZoom: null,   // {xRange, yRange}
    networkZoom: null,

    // Phylogeny
    heatmapAsvTree: 'ward',    // 'ward' or 'phylogeny' — ASV ordering on heatmap
    heatmapCellSize: 3,
    heatmapMinMaxRA: 1.0,     // min(max(RA%)) to include an ASV in heatmap
    treeLayout: 'rc',
    treeMinBootstrap: 0,
    treePrune: false,
    treeLabelLevels: ['Genus', 'id'],
  });

  // Resolve the assay selection to an ASV set once, here, rather than in each
  // view — otherwise the sidebar counts and the plots can disagree.
  $effect(() => {
    store.assayAsvIds = asvIdsForAssays(filters.assays);
    store.assayCacheKey = [...(filters.assays || [])].sort().join(';');
  });

  function updateTab() {
    const hash = window.location.hash.replace('#', '') || 'samples';
    if (['samples', 'network', 'phylogeny', 'heatmap', 'tables', 'provenance'].includes(hash)) {
      activeTab = hash;
      // Picking a tab is a request to see it, not to keep looking at filters.
      sidebarOpen = false;
    }
  }

  onMount(() => {
    updateTab();
    window.addEventListener('hashchange', updateTab);
    loadData();
    return () => window.removeEventListener('hashchange', updateTab);
  });
</script>

<!-- 100vh is taller than the visible area on mobile browsers, where the address
     bar overlays the page; dvh tracks what is actually on screen. Declared twice
     rather than as utilities so the cascade does the fallback: a browser without
     dvh ignores the second and keeps 100vh. -->
<div class="flex flex-col" style="height:100vh;height:100dvh">
  <NavBar {activeTab} onToggleFilters={() => (sidebarOpen = !sidebarOpen)} />

  <div class="relative flex flex-1 overflow-hidden">
    {#if !store.loading && !store.error}
      <Sidebar {activeTab} bind:filters open={sidebarOpen}
               onClose={() => (sidebarOpen = false)} />
      <!-- Dim only the area above the sheet, and let a tap there dismiss it. No
           full-screen scrim: the point of a bottom sheet is that the plot stays
           visible and legible while the filters are open. -->
      {#if sidebarOpen}
        <button type="button" aria-label="Close filters" tabindex="-1"
                onclick={() => (sidebarOpen = false)}
                class="fixed inset-x-0 top-0 bottom-[55vh] z-30 bg-black/30 md:hidden"></button>
      {/if}
    {/if}

    <main class="min-w-0 flex-1 overflow-hidden">
      {#if store.loading}
        <div class="flex h-full items-center justify-center">
          <div class="text-center">
            <div class="mb-4 h-8 w-8 animate-spin rounded-full border-2 border-blue-500 border-t-transparent mx-auto"></div>
            <p class="text-sm text-muted">Loading data...</p>
          </div>
        </div>
      {:else if store.error}
        <div class="flex h-full items-center justify-center">
          <div class="rounded-lg border border-rose-300 bg-rose-50 p-6 text-center dark:border-red-800 dark:bg-red-950/50">
            <p class="text-rose-700 dark:text-red-400">{store.error}</p>
            <button
              class="mt-3 rounded bg-rose-600 px-4 py-1.5 text-sm text-white hover:bg-rose-700 dark:bg-red-800 dark:hover:bg-red-700"
              onclick={() => loadData()}
            >Retry</button>
          </div>
        </div>
      {:else}
        {#if activeTab === 'samples'}
          <SampleView {filters} />
        {:else if activeTab === 'network'}
          <NetworkView {filters} />
        {:else if activeTab === 'phylogeny'}
          <PhyloTreeView {filters} />
        {:else if activeTab === 'heatmap'}
          <HeatmapView {filters} />
        {:else if activeTab === 'tables'}
          <TablesView {filters} />
        {:else if activeTab === 'provenance'}
          <ProvenanceView />
        {/if}
      {/if}
    </main>
  </div>
</div>
