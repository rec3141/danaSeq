<script>
  // `short` is what a narrow screen gets. Wrapping to a second row beats
  // scrolling — a tab you cannot see is a tab you will not press — and six full
  // labels do not fit a phone on any number of rows worth having.
  const tabs = [
    { id: 'samples', label: 'Sample Explorer', short: 'Samples', hash: '#samples' },
    { id: 'network', label: 'ASV Network', short: 'Network', hash: '#network' },
    { id: 'phylogeny', label: 'Phylogeny', short: 'Tree', hash: '#phylogeny' },
    { id: 'heatmap', label: 'Heatmap', short: 'Heat', hash: '#heatmap' },
    { id: 'tables', label: 'Data Tables', short: 'Tables', hash: '#tables' },
    { id: 'provenance', label: 'Provenance', short: 'Prov', hash: '#provenance' },
  ];

  let { activeTab = 'samples', onToggleFilters = null } = $props();
</script>

<nav class="flex flex-wrap items-center gap-x-1 border-b border-slate-800 bg-slate-900/80 px-2 py-1 backdrop-blur sm:flex-nowrap sm:py-0 sm:px-4">
  {#if onToggleFilters}
    <!-- The sidebar is a drawer below md, so it needs a way back. -->
    <button type="button" onclick={onToggleFilters} aria-label="Filters"
            class="mr-1 shrink-0 rounded px-2 py-3 text-slate-300 hover:text-white md:hidden">
      <svg viewBox="0 0 20 20" fill="currentColor" class="h-5 w-5" aria-hidden="true">
        <path d="M3 5h14M3 10h14M3 15h14" stroke="currentColor" stroke-width="2" stroke-linecap="round"/>
      </svg>
    </button>
  {/if}

  <a href="#" class="mr-2 flex shrink-0 items-center gap-2 py-3 text-sm font-bold tracking-wide text-slate-100 sm:mr-4">
    <span class="text-blue-400">&#x25C8;</span>
    <!-- The full name is the first thing worth losing on a narrow screen. -->
    <span class="hidden sm:inline">dānaSeq Amplicon</span>
    <span class="sm:hidden">dānaSeq</span>
  </a>

  <!-- Wrap onto as many rows as it takes; every tab stays reachable without
       scrolling sideways. -->
  <div class="flex min-w-0 flex-1 flex-wrap items-center gap-x-1">

  {#each tabs as tab}
    <a
      href={tab.hash}
      class="relative shrink-0 px-2 py-2 text-sm font-medium transition-colors sm:px-4 sm:py-3
        {activeTab === tab.id
          ? 'text-blue-400'
          : 'text-slate-400 hover:text-slate-200'}"
    >
      <span class="hidden sm:inline">{tab.label}</span>
      <span class="sm:hidden">{tab.short}</span>
      {#if activeTab === tab.id}
        <span class="absolute inset-x-2 -bottom-px h-0.5 rounded-full bg-blue-500"></span>
      {/if}
    </a>
  {/each}
  </div>

  <div class="ml-auto hidden shrink-0 pl-2 text-xs text-slate-500 lg:block">Amplicon Pipeline Viz</div>
</nav>
