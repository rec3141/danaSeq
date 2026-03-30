<script>
  import { onMount } from 'svelte';
  import NavBar from './components/layout/NavBar.svelte';
  import CartDrawer from './components/ui/CartDrawer.svelte';
  import LoadingSpinner from './components/ui/LoadingSpinner.svelte';
  import { loading, error, loadAllData } from './stores/data.js';

  import SampleExplorerView from './views/SampleExplorerView.svelte';
  import ReadExplorerView from './views/ReadExplorerView.svelte';
  import MapView from './views/MapView.svelte';
  import EnvironmentalView from './views/EnvironmentalView.svelte';
  import TaxonomyView from './views/TaxonomyView.svelte';
  import FunctionView from './views/FunctionView.svelte';

  let activeTab = $state('samples');
  let cartOpen = $state(false);

  function updateHash() {
    const hash = window.location.hash.replace('#', '') || 'samples';
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
  <NavBar {activeTab} onCartToggle={() => cartOpen = !cartOpen} />
  <CartDrawer open={cartOpen} />

  <main class="max-w-screen-2xl mx-auto px-4 py-6" class:mr-80={cartOpen}>
    {#if $loading}
      <LoadingSpinner />
    {:else if $error}
      <div class="bg-rose-900/30 border border-rose-700 rounded-lg p-4 text-rose-300">
        <h3 class="font-semibold mb-1">Data loading error</h3>
        <p class="text-sm">{$error}</p>
        <p class="text-xs mt-2 text-rose-400">Run <code class="bg-slate-800 px-1 rounded">npm run preprocess</code> to generate data files.</p>
      </div>
    {:else}
      {#if activeTab === 'samples'}
        <SampleExplorerView />
      {:else if activeTab === 'reads'}
        <ReadExplorerView />
      {:else if activeTab === 'map'}
        <MapView />
      {:else if activeTab === 'environmental'}
        <EnvironmentalView />
      {:else if activeTab === 'taxonomy'}
        <TaxonomyView />
      {:else if activeTab === 'function'}
        <FunctionView />
      {/if}
    {/if}
  </main>
</div>
