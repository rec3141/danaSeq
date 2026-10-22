<script>
  import { onMount } from 'svelte';
  import NavBar from './components/layout/NavBar.svelte';
  import LoadingSpinner from './components/ui/LoadingSpinner.svelte';
  import { index, loading, error, loadIndex } from './stores/data.js';
  import { initMetadata } from './stores/metadata.js';

  import OverviewView from './views/OverviewView.svelte';
  import ExplorerView from './views/ExplorerView.svelte';
  import MagsView from './views/MagsView.svelte';
  import TaxonomyView from './views/TaxonomyView.svelte';
  import FunctionView from './views/FunctionView.svelte';
  import GenesView from './views/GenesView.svelte';
  import SamplesView from './views/SamplesView.svelte';

  let activeTab = $state('overview');

  function updateHash() {
    const hash = window.location.hash.replace('#', '') || 'overview';
    // Old tab ids from the first release of this dashboard
    activeTab = { tsne: 'explorer', functions: 'function' }[hash] ?? hash;
  }

  onMount(() => {
    updateHash();
    window.addEventListener('hashchange', updateHash);
    loadIndex().then(initMetadata);
    return () => window.removeEventListener('hashchange', updateHash);
  });

  $effect(() => { if ($index?.title) document.title = `${$index.title} · dānaSeq RNA`; });
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
        <p class="text-xs mt-2 text-rose-400">The viz reads <code class="bg-slate-800 px-1 rounded">index.json</code> and the per-reference tables from <code class="bg-slate-800 px-1 rounded">&lt;outdir&gt;/viz/</code>.</p>
      </div>
    {:else}
      {#if activeTab === 'overview'}
        <OverviewView />
      {:else if activeTab === 'explorer'}
        <ExplorerView />
      {:else if activeTab === 'mags'}
        <MagsView />
      {:else if activeTab === 'taxonomy'}
        <TaxonomyView />
      {:else if activeTab === 'function'}
        <FunctionView />
      {:else if activeTab === 'genes'}
        <GenesView />
      {:else if activeTab === 'samples'}
        <SamplesView />
      {/if}
    {/if}
  </main>
</div>
