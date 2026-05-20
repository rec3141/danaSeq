<script>
  import { onMount } from 'svelte';
  import NavBar from './components/NavBar.svelte';
  import { loading, error, loadAllData } from './stores/data.js';

  import OverviewView   from './views/OverviewView.svelte';
  import SamplesView    from './views/SamplesView.svelte';
  import ExpressionView from './views/ExpressionView.svelte';
  import ReferencesView from './views/ReferencesView.svelte';

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
      <div class="text-slate-400 text-center py-12">Loading…</div>
    {:else if $error}
      <div class="bg-rose-900/30 border border-rose-700 rounded-lg p-4 text-rose-300">
        <h3 class="font-semibold mb-1">Data loading error</h3>
        <p class="text-sm">{$error}</p>
      </div>
    {:else}
      {#if activeTab === 'overview'}
        <OverviewView />
      {:else if activeTab === 'samples'}
        <SamplesView />
      {:else if activeTab === 'expression'}
        <ExpressionView />
      {:else if activeTab === 'references'}
        <ReferencesView />
      {/if}
    {/if}
  </main>
</div>
