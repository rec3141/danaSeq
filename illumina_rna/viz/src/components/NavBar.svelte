<script>
  let { activeTab = 'overview' } = $props();
  const tabs = [
    { id: 'overview',   label: 'Overview' },
    { id: 'samples',    label: 'Samples' },
    { id: 'expression', label: 'Expression' },
    { id: 'references', label: 'References' },
  ];
</script>

<nav class="bg-slate-900 border-b border-slate-800 sticky top-0 z-10">
  <div class="max-w-screen-2xl mx-auto px-4 flex items-center">
    <div class="font-semibold text-cyan-400 py-3 mr-8">danaSeq · RNA</div>
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
  </div>
</nav>
