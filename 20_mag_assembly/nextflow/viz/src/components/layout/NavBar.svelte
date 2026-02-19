<script>
  const tabs = [
    { id: 'overview', label: 'Overview' },
    { id: 'contigs', label: 'Explorer' },
    { id: 'quality', label: 'Quality' },
    { id: 'taxonomy', label: 'Taxonomy' },
    { id: 'coverage', label: 'Abundance' },
    { id: 'kegg', label: 'Function' },
    { id: 'eukaryotic', label: 'Eukaryotic' },
    { id: 'mge', label: 'MGE' },
  ];

  let { activeTab = 'overview' } = $props();
</script>

<nav class="bg-slate-900 border-b border-slate-700 sticky top-0 z-50">
  <div class="max-w-screen-2xl mx-auto px-4">
    <div class="flex items-center h-14">
      <div class="flex items-center gap-2 mr-8 shrink-0">
        <span class="text-cyan-400 font-semibold text-lg tracking-tight">danaSeq</span>
        <span class="text-slate-400 text-sm">MAG Dashboard</span>
      </div>
      <div class="flex gap-1 overflow-x-auto">
        {#each tabs as tab}
          <a
            href="#{tab.id}"
            class="px-3 py-2 text-sm font-medium rounded-md transition-colors whitespace-nowrap
              {activeTab === tab.id
                ? 'bg-slate-800 text-cyan-400'
                : 'text-slate-400 hover:text-slate-200 hover:bg-slate-800/50'}"
          >
            {tab.label}
          </a>
        {/each}
      </div>
    </div>
  </div>
</nav>
