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
        <a href="https://github.com/rec3141/danaSeq" target="_blank" rel="noopener noreferrer"
           class="flex items-center gap-1.5 text-cyan-400 hover:text-cyan-300 transition-colors">
          <svg class="w-4 h-4" viewBox="0 0 16 16" fill="currentColor">
            <path d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38
              0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15
              -.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.21 1.87.87 2.33.66.07-.52.28-.87
              .51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12
              0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82
              2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65
              3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013
              0 0016 8c0-4.42-3.58-8-8-8z"/>
          </svg>
          <span class="font-semibold text-lg tracking-tight">d&#257;naSeq</span>
        </a>
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
