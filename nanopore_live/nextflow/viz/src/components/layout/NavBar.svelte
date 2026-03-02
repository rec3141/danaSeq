<script>
  import { cartCount, cartActive } from '../../stores/cart.js';
  import { loadMetadataTsv, metadata, samples } from '../../stores/data.js';
  import { get } from 'svelte/store';

  const tabs = [
    { id: 'samples', label: 'Samples' },
    { id: 'reads', label: 'Reads' },
    { id: 'map', label: 'Map' },
    { id: 'environmental', label: 'Environmental' },
    { id: 'taxonomy', label: 'Taxonomy' },
    { id: 'function', label: 'Function' },
  ];

  let { activeTab = 'samples', onCartToggle = null } = $props();
  let metaInput;
  let metaStatus = $state('');

  let metaError = $state(false);

  function downloadMetadata() {
    const meta = get(metadata);
    const samps = get(samples);
    if (!meta || !Object.keys(meta).length) return;

    // Collect all column names across all entries
    const allCols = new Set();
    for (const v of Object.values(meta)) {
      for (const k of Object.keys(v)) allCols.add(k);
    }
    const cols = [...allCols].sort();
    const header = ['sample_id', ...cols];

    const rows = [header.join('\t')];
    // Use samples order if available, otherwise metadata keys
    const ids = samps ? samps.map(s => s.id) : Object.keys(meta);
    for (const id of ids) {
      const m = meta[id];
      if (!m) continue;
      const vals = [id, ...cols.map(c => m[c] != null ? String(m[c]) : '')];
      rows.push(vals.join('\t'));
    }

    const blob = new Blob([rows.join('\n')], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    const date = new Date().toISOString().slice(0, 10);
    a.href = url;
    a.download = `metadata_${date}.tsv`;
    a.click();
    URL.revokeObjectURL(url);
  }

  function handleMetaUpload(e) {
    const file = e.target.files?.[0];
    if (!file) return;
    const reader = new FileReader();
    reader.onload = () => {
      const result = loadMetadataTsv(reader.result);
      if (result.error) {
        metaStatus = result.error;
        metaError = true;
        // Errors stay visible until dismissed
      } else {
        metaStatus = `${result.matched}/${result.total} matched`;
        metaError = false;
        setTimeout(() => { metaStatus = ''; }, 5000);
      }
    };
    reader.readAsText(file);
    e.target.value = '';
  }
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
        <span class="text-slate-400 text-sm">Sample Dashboard</span>
      </div>
      <div class="flex gap-1 overflow-x-auto flex-1">
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
      <!-- Metadata upload -->
      <div class="ml-4 flex items-center gap-2 shrink-0">
        <input type="file" accept=".tsv,.txt,.csv" class="hidden" bind:this={metaInput} onchange={handleMetaUpload} />
        <button
          class="text-xs px-2 py-1.5 rounded border border-slate-600 text-slate-400 hover:text-slate-200 hover:border-slate-500 transition-colors"
          onclick={() => metaInput.click()}
          title="Upload metadata TSV (sample_id, barcode, lat, lon, ...)"
        >
          <svg class="w-3.5 h-3.5 inline -mt-0.5 mr-1" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4-8l-4-4m0 0L8 8m4-4v12"/>
          </svg>
          Metadata
        </button>
        {#if $metadata && Object.keys($metadata).length > 0}
          <button
            class="text-xs px-2 py-1.5 rounded border border-slate-600 text-slate-400 hover:text-slate-200 hover:border-slate-500 transition-colors"
            onclick={downloadMetadata}
            title="Download current metadata as TSV"
          >
            <svg class="w-3.5 h-3.5 inline -mt-0.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4 0l-4 4m0 0l-4-4m4 4V4"/>
            </svg>
          </button>
        {/if}
        {#if metaStatus}
          <span
            class="text-[10px] {metaError ? 'text-rose-400 cursor-pointer' : 'text-cyan-400'}"
            onclick={() => { metaStatus = ''; metaError = false; }}
          >{metaStatus}{#if metaError} ✕{/if}</span>
        {/if}
      </div>

      <!-- Cart controls -->
      <div class="ml-4 flex items-center gap-1">
        {#if $cartCount > 0}
          <button
            class="flex items-center gap-1 px-2 py-1.5 rounded-l-md text-xs font-medium transition-colors border
              {$cartActive
                ? 'bg-cyan-400/20 text-cyan-400 border-cyan-400/40'
                : 'text-slate-400 hover:text-slate-200 border-slate-600 hover:border-slate-500'}"
            onclick={() => cartActive.update(v => !v)}
            title={$cartActive ? 'Cart filter ON — click to disable' : 'Cart filter OFF — click to enable'}
          >
            <svg class="w-3.5 h-3.5" fill="none" stroke="currentColor" viewBox="0 0 24 24">
              <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2"
                d="M3 4a1 1 0 011-1h16a1 1 0 011 1v2.586a1 1 0 01-.293.707l-6.414 6.414a1 1 0 00-.293.707V17l-4 4v-6.586a1 1 0 00-.293-.707L3.293 7.293A1 1 0 013 6.586V4z"/>
            </svg>
            {$cartActive ? 'ON' : 'OFF'}
          </button>
        {/if}
        <button
          class="flex items-center gap-1.5 px-3 py-1.5 text-sm font-medium transition-colors border
            {$cartCount > 0
              ? ($cartActive ? 'rounded-r-md bg-cyan-400/10 text-cyan-400 border-cyan-400/40' : 'rounded-r-md text-slate-400 hover:text-slate-200 border-slate-600 hover:border-slate-500')
              : 'rounded-md text-slate-400 hover:text-slate-200 hover:bg-slate-800/50 border-transparent'}"
          onclick={() => { if (onCartToggle) onCartToggle(); }}
        >
          <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
            <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2"
              d="M3 3h2l.4 2M7 13h10l4-8H5.4M7 13L5.4 5M7 13l-2.293 2.293c-.63.63-.184 1.707.707 1.707H17m0 0a2 2 0 100 4 2 2 0 000-4zm-8 2a2 2 0 100 4 2 2 0 000-4z"/>
          </svg>
          Cart
          {#if $cartCount > 0}
            <span class="bg-cyan-400 text-slate-900 text-xs font-bold px-1.5 py-0.5 rounded-full min-w-[20px] text-center">
              {$cartCount}
            </span>
          {/if}
        </button>
      </div>
    </div>
  </div>
</nav>
