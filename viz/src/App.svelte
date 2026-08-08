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
  import { store, loadData } from './stores/data.svelte.js';

  let activeTab = $state('samples');

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
    // Samples
    minReads: 0,
    sampleFilter: '',
    showOverlay: true,
    pointScale: 20,
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

  function updateTab() {
    const hash = window.location.hash.replace('#', '') || 'samples';
    if (['samples', 'network', 'phylogeny', 'heatmap', 'tables', 'provenance'].includes(hash)) {
      activeTab = hash;
    }
  }

  onMount(() => {
    updateTab();
    window.addEventListener('hashchange', updateTab);
    loadData();
    return () => window.removeEventListener('hashchange', updateTab);
  });
</script>

<div class="flex h-screen flex-col">
  <NavBar {activeTab} />

  <div class="flex flex-1 overflow-hidden">
    {#if !store.loading && !store.error}
      <Sidebar {activeTab} bind:filters />
    {/if}

    <main class="flex-1 overflow-hidden">
      {#if store.loading}
        <div class="flex h-full items-center justify-center">
          <div class="text-center">
            <div class="mb-4 h-8 w-8 animate-spin rounded-full border-2 border-blue-500 border-t-transparent mx-auto"></div>
            <p class="text-sm text-slate-400">Loading data...</p>
          </div>
        </div>
      {:else if store.error}
        <div class="flex h-full items-center justify-center">
          <div class="rounded-lg border border-red-800 bg-red-950/50 p-6 text-center">
            <p class="text-red-400">{store.error}</p>
            <button
              class="mt-3 rounded bg-red-800 px-4 py-1.5 text-sm text-red-100 hover:bg-red-700"
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
