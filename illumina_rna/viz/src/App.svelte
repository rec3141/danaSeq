<script>
  import { onMount } from 'svelte';
  import NavBar from './components/NavBar.svelte';
  import { index, loading, error, loadIndex } from './stores/data.js';

  import OverviewView  from './views/OverviewView.svelte';
  import MagsView      from './views/MagsView.svelte';
  import FunctionsView from './views/FunctionsView.svelte';
  import GenesView     from './views/GenesView.svelte';
  import SamplesView   from './views/SamplesView.svelte';

  let activeTab = $state('overview');

  function updateHash() {
    activeTab = window.location.hash.replace('#', '') || 'overview';
  }

  onMount(() => {
    updateHash();
    window.addEventListener('hashchange', updateHash);
    loadIndex();
    return () => window.removeEventListener('hashchange', updateHash);
  });

  $effect(() => { if ($index?.title) document.title = `${$index.title} · danaSeq RNA`; });
</script>

<div class="min-h-screen bg-slate-950">
  <NavBar {activeTab} />

  <main class="max-w-screen-2xl mx-auto px-4 py-6">
    {#if $error}
      <div class="bg-rose-900/30 border border-rose-700 rounded-lg p-4 text-rose-300">
        <h3 class="font-semibold mb-1">Data loading error</h3>
        <p class="text-sm">{$error}</p>
      </div>
    {:else if $loading}
      <div class="text-slate-400 text-center py-12">Loading…</div>
    {:else if activeTab === 'overview'}
      <OverviewView />
    {:else if activeTab === 'mags'}
      <MagsView />
    {:else if activeTab === 'functions'}
      <FunctionsView />
    {:else if activeTab === 'genes'}
      <GenesView />
    {:else if activeTab === 'samples'}
      <SamplesView />
    {/if}
  </main>
</div>
