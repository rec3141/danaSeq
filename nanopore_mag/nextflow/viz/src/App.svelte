<script>
  import { onMount } from 'svelte';
  import NavBar from './components/layout/NavBar.svelte';
  import LoadingSpinner from './components/ui/LoadingSpinner.svelte';
  import { loading, error, loadAllData } from './stores/data.js';

  import OverviewView from './views/OverviewView.svelte';
  import QualityView from './views/QualityView.svelte';
  import TaxonomyView from './views/TaxonomyView.svelte';
  import KeggView from './views/KeggView.svelte';
  import MgeView from './views/MgeView.svelte';
  import CoverageView from './views/CoverageView.svelte';
  import EukaryoticView from './views/EukaryoticView.svelte';
  import ContigExplorerView from './views/ContigExplorerView.svelte';

  let activeTab = $state('overview');

  function updateHash() {
    const hash = window.location.hash.replace('#', '') || 'overview';
    activeTab = hash;
  }

  onMount(() => {
    updateHash();
    window.addEventListener('hashchange', updateHash);
    loadAllData();
    return () => window.removeEventListener('hashchange', updateHash);
  });
</script>

<div class="min-h-screen bg-slate-950">
  <NavBar {activeTab} />

  <main class="max-w-screen-2xl mx-auto px-4 py-6">
    {#if $loading}
      <LoadingSpinner />
    {:else if $error}
      <div class="bg-rose-900/30 border border-rose-700 rounded-lg p-4 text-rose-300">
        <h3 class="font-semibold mb-1">Data loading error</h3>
        <p class="text-sm">{$error}</p>
        <p class="text-xs mt-2 text-rose-400">Run <code class="bg-slate-800 px-1 rounded">npm run preprocess</code> to generate data files.</p>
      </div>
    {:else}
      {#if activeTab === 'overview'}
        <OverviewView />
      {:else if activeTab === 'quality'}
        <QualityView />
      {:else if activeTab === 'taxonomy'}
        <TaxonomyView />
      {:else if activeTab === 'kegg'}
        <KeggView />
      {:else if activeTab === 'mge'}
        <MgeView />
      {:else if activeTab === 'coverage'}
        <CoverageView />
      {:else if activeTab === 'eukaryotic'}
        <EukaryoticView />
      {:else if activeTab === 'contigs'}
        <ContigExplorerView />
      {/if}
    {/if}
  </main>
</div>
