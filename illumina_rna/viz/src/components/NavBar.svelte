<script>
  import { index, refName, loadReference } from '../stores/data.js';
  let { activeTab = 'overview' } = $props();
  const tabs = [
    { id: 'overview',  label: 'Overview' },
    { id: 'mags',      label: 'MAGs' },
    { id: 'functions', label: 'Functions' },
    { id: 'genes',     label: 'Genes' },
    { id: 'samples',   label: 'Samples' },
  ];
</script>

<nav class="bg-slate-900 border-b border-slate-800 sticky top-0 z-10">
  <div class="max-w-screen-2xl mx-auto px-4 flex items-center gap-6 flex-wrap">
    <div class="py-3">
      <div class="font-semibold text-cyan-400 leading-tight">danaSeq · RNA</div>
      {#if $index?.title}<div class="text-[11px] text-slate-400 leading-tight">{$index.title}</div>{/if}
    </div>
    <div class="flex">
      {#each tabs as tab}
        <a href="#{tab.id}"
           class="px-4 py-3 text-sm border-b-2 {activeTab === tab.id
             ? 'border-cyan-400 text-cyan-300'
             : 'border-transparent text-slate-400 hover:text-slate-200'}">
          {tab.label}
        </a>
      {/each}
    </div>
    <div class="ml-auto flex items-center gap-4 text-xs">
      {#if ($index?.references?.length ?? 0) > 1}
        <select class="bg-slate-950 border border-slate-700 rounded px-2 py-1 text-slate-200"
          value={$refName} onchange={e => loadReference(e.currentTarget.value)}>
          {#each $index.references as r}<option value={r.name}>{r.name}</option>{/each}
        </select>
      {/if}
      {#if $index?.mag_url}
        <a class="text-cyan-300 hover:underline" href={$index.mag_url} target="_blank" rel="noopener">MAG analysis ↗</a>
      {/if}
    </div>
  </div>
</nav>
