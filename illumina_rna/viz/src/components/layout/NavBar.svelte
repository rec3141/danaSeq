<script>
  import { index, refName, loadReference } from '../../stores/data.js';

  const tabs = [
    { id: 'overview', label: 'Overview' },
    { id: 'explorer', label: 'Explorer' },
    { id: 'mags', label: 'MAGs' },
    { id: 'function', label: 'Function' },
    { id: 'genes', label: 'Genes' },
    { id: 'samples', label: 'Samples' },
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
        <span class="text-slate-400 text-sm">RNA Dashboard</span>
        {#if $index?.title}
          <span class="text-slate-600 text-sm mx-1">/</span>
          <span class="text-slate-400 text-sm truncate max-w-[18rem]" title={$index.title}>{$index.title}</span>
        {/if}
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
      <div class="ml-auto flex items-center gap-3 pl-4 text-xs shrink-0">
        {#if ($index?.references?.length ?? 0) > 1}
          <select
            class="px-2 py-1 rounded-md border border-slate-600 bg-slate-800 text-slate-300 text-xs focus:border-cyan-400 focus:outline-none cursor-pointer"
            value={$refName}
            onchange={e => loadReference(e.currentTarget.value)}
          >
            {#each $index.references as r}<option value={r.name}>{r.name}</option>{/each}
          </select>
        {/if}
        {#if $index?.mag_url}
          <a href={$index.mag_url} target="_blank" rel="noopener noreferrer"
             class="px-3 py-1 rounded-md border border-slate-600 text-slate-400 hover:text-cyan-400 hover:border-cyan-400 transition-colors whitespace-nowrap">
            MAG Dashboard &#x2197;
          </a>
        {/if}
      </div>
    </div>
  </div>
</nav>
