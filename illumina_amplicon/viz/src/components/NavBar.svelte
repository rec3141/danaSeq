<script>
  const tabs = [
    { id: 'samples', label: 'Sample Explorer', hash: '#samples' },
    { id: 'network', label: 'ASV Network', hash: '#network' },
    { id: 'phylogeny', label: 'Phylogeny', hash: '#phylogeny' },
    { id: 'heatmap', label: 'Heatmap', hash: '#heatmap' },
    { id: 'tables', label: 'Data Tables', hash: '#tables' },
    { id: 'provenance', label: 'Provenance', hash: '#provenance' },
  ];

  let { activeTab = 'samples' } = $props();
</script>

<nav class="flex items-center gap-1 border-b border-slate-800 bg-slate-900/80 px-4 backdrop-blur">
  <a href="#" class="mr-4 flex items-center gap-2 py-3 text-sm font-bold tracking-wide text-slate-100">
    <span class="text-blue-400">&#x25C8;</span>
    Microscape
  </a>

  {#each tabs as tab}
    <a
      href={tab.hash}
      class="relative px-4 py-3 text-sm font-medium transition-colors
        {activeTab === tab.id
          ? 'text-blue-400'
          : 'text-slate-400 hover:text-slate-200'}"
    >
      {tab.label}
      {#if activeTab === tab.id}
        <span class="absolute inset-x-2 -bottom-px h-0.5 rounded-full bg-blue-500"></span>
      {/if}
    </a>
  {/each}

  <div class="ml-auto text-xs text-slate-500">Amplicon Pipeline Viz</div>
</nav>
