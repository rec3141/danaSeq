<script>
  import { index, refName, loadReference, samples } from '../../stores/data.js';
  import { sampleMeta, loadMetadataFile, clearUploadedMetadata, hasBundledMetadata,
           metadataTemplate, metadataTsv } from '../../stores/metadata.js';

  const tabs = [
    { id: 'overview', label: 'Overview' },
    { id: 'explorer', label: 'Explorer' },
    { id: 'mags', label: 'MAGs' },
    { id: 'taxonomy', label: 'Taxonomy' },
    { id: 'function', label: 'Function' },
    { id: 'genes', label: 'Genes' },
    { id: 'samples', label: 'Samples' },
  ];

  let { activeTab = 'overview' } = $props();

  // ── Sample metadata menu ──
  let metaOpen = $state(false);
  let metaInput;
  let metaStatus = $state('');
  let metaError = $state(false);
  let nSamples = $derived($samples?.samples?.length ?? 0);

  function save(name, text) {
    const url = URL.createObjectURL(new Blob([text], { type: 'text/tab-separated-values' }));
    const a = document.createElement('a');
    a.href = url; a.download = name; a.click();
    URL.revokeObjectURL(url);
  }
  async function upload(e) {
    const file = e.currentTarget.files?.[0];
    e.currentTarget.value = '';
    if (!file) return;
    const r = await loadMetadataFile(file);
    metaError = !!r.error;
    metaStatus = r.error ?? `${r.matched} of ${nSamples} samples matched, ${r.columns} column${r.columns === 1 ? '' : 's'}`;
    if (!r.error) setTimeout(() => { metaStatus = ''; }, 6000);
  }
  function onWindowClick(e) {
    if (metaOpen && !e.target.closest?.('[data-meta-menu]')) metaOpen = false;
  }
</script>

<svelte:window onclick={onWindowClick} />

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
        <div class="relative" data-meta-menu>
          <button
            class="px-3 py-1 rounded-md border transition-colors whitespace-nowrap
              {$sampleMeta ? 'border-cyan-400/60 text-cyan-400 hover:border-cyan-400' : 'border-slate-600 text-slate-400 hover:text-cyan-400 hover:border-cyan-400'}"
            onclick={() => metaOpen = !metaOpen}
            title="Sample metadata: colour, order and compare samples by it"
          >Metadata{$sampleMeta ? ` · ${$sampleMeta.columns.length}` : ''} &#x25BE;</button>
          {#if metaOpen}
            <div class="absolute right-0 mt-1 w-80 bg-slate-800 border border-slate-600 rounded-lg shadow-xl p-3 z-50 text-xs text-slate-300 space-y-2">
              {#if $sampleMeta}
                <div>
                  <div class="text-slate-200 font-medium truncate" title={$sampleMeta.source}>{$sampleMeta.source}</div>
                  <div class="text-slate-500">{$sampleMeta.origin === 'upload' ? 'Uploaded, kept in this browser only' : 'Shipped with this dashboard'}</div>
                  <div class="mt-1">{$sampleMeta.matched} of {nSamples} samples matched on <span class="font-mono">{$sampleMeta.idColumn}</span>{#if $sampleMeta.unmatchedRows.length}; <span class="text-amber-400">{$sampleMeta.unmatchedRows.length} row{$sampleMeta.unmatchedRows.length === 1 ? '' : 's'} unmatched</span>{/if}</div>
                  <div class="text-slate-400 mt-1">{$sampleMeta.columns.map(c => c.key).join(', ') || 'no columns with values'}</div>
                </div>
              {:else}
                <p class="text-slate-400">Add a table with one row per sample (first column, or a <span class="font-mono">sample</span> column, holds the sample name). Its columns become sample annotations: t-SNE colours, heatmap strips, sample order and expression-vs-metadata plots.</p>
              {/if}
              <div class="flex flex-wrap gap-2 pt-1">
                <button class="px-2 py-1 rounded border border-cyan-400 bg-cyan-400/10 text-cyan-400" onclick={() => metaInput.click()}>Upload TSV / CSV / XLSX</button>
                <button class="px-2 py-1 rounded border border-slate-600 text-slate-400 hover:text-slate-200" onclick={() => save('metadata_template.tsv', metadataTemplate())}>Template</button>
                {#if $sampleMeta}
                  <button class="px-2 py-1 rounded border border-slate-600 text-slate-400 hover:text-slate-200" onclick={() => save('metadata.tsv', metadataTsv())}>Download</button>
                {/if}
                {#if $sampleMeta?.origin === 'upload'}
                  <button class="px-2 py-1 rounded border border-slate-600 text-slate-400 hover:text-rose-300 hover:border-rose-400"
                    onclick={() => { clearUploadedMetadata(); metaStatus = ''; }}>{hasBundledMetadata() ? 'Use shipped table' : 'Remove'}</button>
                {/if}
              </div>
              {#if metaStatus}<p class={metaError ? 'text-rose-400' : 'text-emerald-400'}>{metaStatus}</p>{/if}
            </div>
          {/if}
          <input type="file" accept=".tsv,.csv,.txt,.xlsx,.xls" class="hidden" bind:this={metaInput} onchange={upload} />
        </div>
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
