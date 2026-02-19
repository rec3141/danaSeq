<script>
  import D3Heatmap from '../components/charts/D3Heatmap.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { onMount } from 'svelte';
  import { keggHeatmap, checkm2All, loadCheckm2All } from '../stores/data.js';
  import { selectedMag } from '../stores/selection.js';

  let keggData = $derived($keggHeatmap);
  let allBins = $derived($checkm2All);
  let selected = $derived($selectedMag);

  onMount(() => { loadCheckm2All(); });

  // Taxonomy source for detail panel
  let taxSource = $state('kaiju');
  const taxSources = [
    { key: 'kaiju', label: 'Kaiju' },
    { key: 'kraken2', label: 'Kraken2' },
    { key: 'rrna', label: 'rRNA' },
  ];

  // Binner toggle buttons — filter heatmap rows by bin name prefix
  let activeBinners = $state(new Set(['dastool']));

  const binnerDefs = [
    { key: 'dastool',  label: 'DAS Tool', prefix: 'dastool-' },
    { key: 'semibin',  label: 'SemiBin2', prefix: 'semibin_' },
    { key: 'metabat',  label: 'MetaBAT2', prefix: 'metabat_' },
    { key: 'maxbin',   label: 'MaxBin2',  prefix: 'maxbin_' },
    { key: 'lorbin',   label: 'LorBin',   prefix: 'lorbin_' },
    { key: 'comebin',  label: 'COMEBin',  prefix: 'comebin_' },
  ];

  function toggleBinner(key) {
    const next = new Set(activeBinners);
    if (next.has(key)) {
      if (next.size > 1) next.delete(key);
    } else {
      next.add(key);
    }
    activeBinners = next;
  }

  // Detect which binner a MAG belongs to by name prefix
  function magBinner(id) {
    for (const { key, prefix } of binnerDefs) {
      if (key !== 'all' && id.startsWith(prefix)) return key;
    }
    return 'unknown';
  }

  // Count bins per binner from kegg heatmap MAG IDs
  let binnerCounts = $derived.by(() => {
    if (!keggData) return {};
    const counts = {};
    for (const id of keggData.mag_ids) {
      const b = magBinner(id);
      counts[b] = (counts[b] || 0) + 1;
    }
    return counts;
  });

  // Module category filters — derived from data, additive toggles
  let categoryDefs = $derived.by(() => {
    if (!keggData?.module_groups) return [];
    const catModules = {};
    keggData.module_ids.forEach((mid, i) => {
      const group = keggData.module_groups[i] || 'Other';
      if (!catModules[group]) catModules[group] = [];
      catModules[group].push(mid);
    });
    // Fixed display order
    const order = ['Carbon', 'Nitrogen', 'Sulfur', 'Energy', 'Biosynthesis', 'Other'];
    const result = [];
    for (const key of order) {
      if (catModules[key]) result.push({ key, modules: catModules[key] });
    }
    // Append any groups not in the fixed order
    for (const [key, modules] of Object.entries(catModules)) {
      if (!order.includes(key)) result.push({ key, modules });
    }
    return result;
  });

  let activeCategories = $state(null);  // null = not yet initialized

  // Initialize activeCategories when categoryDefs first loads
  $effect(() => {
    if (activeCategories === null && categoryDefs.length > 0) {
      activeCategories = new Set(categoryDefs.map(c => c.key));
    }
  });

  function toggleCategory(key) {
    const next = new Set(activeCategories);
    if (next.has(key)) {
      next.delete(key);
    } else {
      next.add(key);
    }
    activeCategories = next;
  }

  // Union of all selected category module IDs (all selected = no filtering)
  let allowedModules = $derived.by(() => {
    if (!activeCategories || activeCategories.size === 0 || activeCategories.size === categoryDefs.length) return null;
    const allowed = new Set();
    for (const { key, modules } of categoryDefs) {
      if (activeCategories.has(key)) {
        for (const m of modules) allowed.add(m);
      }
    }
    return allowed;
  });

  // Filter by both binner and module category
  let filteredData = $derived.by(() => {
    if (!keggData) return null;

    // Row filter: binner name prefix
    const rowIndices = keggData.mag_ids
      .map((id, i) => activeBinners.has(magBinner(id)) ? i : -1)
      .filter(i => i >= 0);

    // Column filter: module category (null = show all)
    let colIndices;
    if (!allowedModules) {
      colIndices = keggData.module_ids.map((_, i) => i);
    } else {
      colIndices = keggData.module_ids
        .map((id, i) => allowedModules.has(id) ? i : -1)
        .filter(i => i >= 0);
    }

    // Recompute row_order
    const rowSet = new Set(rowIndices);
    const filteredRowOrder = keggData.row_order.filter(i => rowSet.has(i));
    const oldToNewRow = {};
    rowIndices.forEach((oldIdx, newIdx) => { oldToNewRow[oldIdx] = newIdx; });

    // Recompute col_order
    const colSet = new Set(colIndices);
    const filteredColOrder = keggData.col_order.filter(i => colSet.has(i));
    const oldToNewCol = {};
    colIndices.forEach((oldIdx, newIdx) => { oldToNewCol[oldIdx] = newIdx; });

    return {
      mag_ids: rowIndices.map(i => keggData.mag_ids[i]),
      module_ids: colIndices.map(i => keggData.module_ids[i]),
      module_names: colIndices.map(i => keggData.module_names[i]),
      module_categories: colIndices.map(i => keggData.module_categories?.[i] || ''),
      matrix: rowIndices.map(ri => colIndices.map(ci => keggData.matrix[ri][ci])),
      row_order: filteredRowOrder.map(i => oldToNewRow[i]),
      col_order: filteredColOrder.map(i => oldToNewCol[i]),
    };
  });

  function handleRowClick(magId) {
    selectedMag.set(magId);
  }

  // Table: per-module completeness for selected MAG
  let selectedMagModules = $derived.by(() => {
    if (!selected || !keggData) return [];
    const magIdx = keggData.mag_ids.indexOf(selected);
    if (magIdx < 0) return [];
    return keggData.module_ids.map((id, j) => ({
      module: id,
      name: keggData.module_names[j],
      completeness: keggData.matrix[magIdx][j],
      group: keggData.module_groups?.[j] || '',
      category: keggData.module_categories?.[j] || '',
    })).filter(m => m.completeness > 0 && (!allowedModules || allowedModules.has(m.module)))
      .sort((a, b) => b.completeness - a.completeness);
  });

  function exportModuleTsv() {
    if (!selectedMagModules.length) return;
    const header = 'Module\tName\tCompleteness\tGroup\tCategory';
    const body = selectedMagModules.map(m =>
      `${m.module}\t${m.name}\t${(m.completeness * 100).toFixed(1)}%\t${m.group}\t${m.category}`
    ).join('\n');
    const blob = new Blob([header + '\n' + body], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `${selected}_modules.tsv`;
    a.click();
    URL.revokeObjectURL(url);
  }

  const moduleTableColumns = [
    { key: 'module', label: 'Module' },
    { key: 'completeness', label: '%', render: (v) => {
      const pct = (v * 100).toFixed(0);
      const cls = v >= 0.75 ? 'text-emerald-400' : v >= 0.5 ? 'text-amber-400' : 'text-slate-400';
      return `<span class="${cls}">${pct}%</span>`;
    }},
    { key: 'name', label: 'Name' },
    { key: 'group', label: 'Group' },
    { key: 'category', label: 'Category' },
  ];

  // Snapshot: serialize SVG to downloadable file
  let heatmapEl;
  function downloadSnapshot() {
    if (!heatmapEl) return;
    const svg = heatmapEl.querySelector('svg');
    if (!svg) return;
    const clone = svg.cloneNode(true);
    // Add white background for publication
    const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    bg.setAttribute('width', '100%');
    bg.setAttribute('height', '100%');
    bg.setAttribute('fill', '#0f172a');
    clone.insertBefore(bg, clone.firstChild);
    const blob = new Blob([new XMLSerializer().serializeToString(clone)], { type: 'image/svg+xml' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    const binners = [...activeBinners].join('+');
    a.download = `kegg_modules_${binners}.svg`;
    a.click();
    URL.revokeObjectURL(url);
  }

  // Summary stats
  let summaryStats = $derived.by(() => {
    if (!filteredData) return null;
    const nMags = filteredData.mag_ids.length;
    const nMods = filteredData.module_ids.length;
    let completeModules = 0;
    for (let j = 0; j < nMods; j++) {
      for (let i = 0; i < nMags; i++) {
        if (filteredData.matrix[i][j] >= 0.75) { completeModules++; break; }
      }
    }
    let sum = 0, count = 0;
    for (let i = 0; i < nMags; i++) {
      for (let j = 0; j < nMods; j++) {
        if (filteredData.matrix[i][j] > 0) { sum += filteredData.matrix[i][j]; count++; }
      }
    }
    return { nMags, nMods, completeModules, avgCompleteness: count > 0 ? sum / count : 0 };
  });

  // Unified detail lookup — checkm2All has taxonomy + MGE for all bins
  let selectedBinData = $derived.by(() => {
    if (!selected || !allBins) return null;
    return allBins.find(b => b.dastool_name === selected)
        || allBins.find(b => b.name === selected)
        || null;
  });
</script>

{#if keggData}
  <div class="mb-4 flex items-center gap-3 flex-wrap text-xs">
    <span class="text-slate-400">Binners:</span>
    {#each binnerDefs as { key, label }}
      {#if binnerCounts[key]}
        <button
          class="px-3 py-1 rounded-md border transition-colors text-center
            {activeBinners.has(key) ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
          onclick={() => toggleBinner(key)}
        >
          {label} <span class="text-slate-500">({binnerCounts[key]})</span>
        </button>
      {/if}
    {/each}
    <span class="text-slate-600">|</span>
    <span class="text-slate-400">Category:</span>
    {#each categoryDefs as { key, modules }}
      <button
        class="px-3 py-1 rounded-md border transition-colors
          {activeCategories?.has(key) ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        onclick={() => toggleCategory(key)}
      >
        {key} <span class="text-slate-500">({modules.length})</span>
      </button>
    {/each}
    {#if summaryStats}
      <span class="text-slate-500 ml-auto">
        {summaryStats.nMods} modules, {summaryStats.completeModules} detected (≥75%), avg {(summaryStats.avgCompleteness * 100).toFixed(0)}%
      </span>
    {/if}
  </div>

  <div class="grid grid-cols-1 lg:grid-cols-3 gap-6">
    <div bind:this={heatmapEl} class="lg:col-span-2 bg-slate-800 rounded-lg p-4 border border-slate-700 overflow-auto" style="max-height: 700px;">
      <div class="flex items-center justify-between mb-2">
        <h3 class="text-sm font-medium text-slate-400">
          KEGG Module Completeness ({filteredData?.mag_ids?.length || 0} bins × {filteredData?.module_ids?.length || 0} modules)
        </h3>
        <button
          class="px-2 py-1 text-xs rounded border border-slate-600 text-slate-400 hover:border-slate-500 hover:text-slate-300 transition-colors"
          onclick={downloadSnapshot}
        >SVG</button>
      </div>
      <D3Heatmap data={filteredData} onRowClick={handleRowClick} selectedRow={selected}
        tooltipFormat={(val, magId, modId, modName) => {
          const idx = filteredData?.module_ids?.indexOf(modId) ?? -1;
          const cat = idx >= 0 ? filteredData.module_categories[idx] : '';
          return `<strong>${modId}</strong><br>${modName}<br>${cat ? `<span style="color:#94a3b8">${cat}</span><br>` : ''}MAG: ${magId}<br>Completeness: ${(val * 100).toFixed(0)}%`;
        }}
      />
    </div>

    <div class="flex flex-col gap-4 overflow-auto" style="max-height: 700px;">
      {#if selectedBinData}
        <div class="bg-slate-800 rounded-lg p-4 border border-slate-700 shrink-0">
          <h3 class="text-cyan-400 font-semibold mb-3">{selectedBinData.name}</h3>
          <p class="text-xs text-slate-500 -mt-2 mb-3">{selectedBinData.binner}{selectedBinData.is_dastool ? ' (DAS Tool consensus)' : ''}</p>
          <dl class="space-y-2 text-sm">
            <div class="flex justify-between">
              <dt class="text-slate-400">Completeness</dt>
              <dd class="font-mono text-emerald-400">{selectedBinData.completeness}%</dd>
            </div>
            <div class="flex justify-between">
              <dt class="text-slate-400">Contamination</dt>
              <dd class="font-mono text-amber-400">{selectedBinData.contamination}%</dd>
            </div>
            <div class="flex justify-between">
              <dt class="text-slate-400">GC Content</dt>
              <dd class="font-mono">{(selectedBinData.gc * 100).toFixed(1)}%</dd>
            </div>
            <div class="flex justify-between">
              <dt class="text-slate-400">Genome Size</dt>
              <dd class="font-mono">{(selectedBinData.genome_size / 1e6).toFixed(2)} Mbp</dd>
            </div>
            <div class="flex justify-between">
              <dt class="text-slate-400">N50</dt>
              <dd class="font-mono">{(selectedBinData.n50 / 1e3).toFixed(1)} Kb</dd>
            </div>
            {#if selectedBinData.taxonomy && Object.keys(selectedBinData.taxonomy).length > 0}
              <div class="border-t border-slate-700 pt-2 mt-2">
                <dt class="text-slate-400 mb-1 flex items-center gap-2">
                  Taxonomy
                  {#each taxSources as { key, label }}
                    {#if selectedBinData.taxonomy[key]}
                      <button
                        class="text-[10px] px-1.5 py-0.5 rounded border transition-colors
                          {taxSource === key ? 'border-cyan-400 text-cyan-400' : 'border-slate-600 text-slate-500 hover:text-slate-300'}"
                        onclick={() => taxSource = key}
                      >{label}</button>
                    {/if}
                  {/each}
                </dt>
                {#if selectedBinData.taxonomy[taxSource]}
                  <dd class="font-mono text-xs">
                    {#each Object.entries(selectedBinData.taxonomy[taxSource]).filter(([,v]) => v) as [rank, val]}
                      <div><span class="text-slate-500">{rank}:</span> {val}</div>
                    {/each}
                  </dd>
                {:else}
                  <dd class="text-xs text-slate-500 italic">No {taxSource} data for this bin</dd>
                {/if}
              </div>
            {/if}
            <div class="border-t border-slate-700 pt-2 mt-2">
              <dt class="text-slate-400 mb-1">MGE Summary</dt>
              <dd class="font-mono text-xs grid grid-cols-2 gap-1">
                <span>Virus: {selectedBinData.n_virus}</span>
                <span>Plasmid: {selectedBinData.n_plasmid}</span>
                <span>Defense: {selectedBinData.n_defense}</span>
                <span>Integron: {selectedBinData.n_integron}</span>
              </dd>
            </div>
          </dl>
        </div>
      {/if}

      <div class="bg-slate-800 rounded-lg p-4 border border-slate-700 flex flex-col min-h-0 flex-1">
        {#if selected && selectedMagModules.length > 0}
          <div class="flex items-center justify-between mb-2">
            <h3 class="text-sm font-medium text-slate-400">
              Modules in <span class="text-cyan-400">{selected}</span> ({selectedMagModules.length} detected)
            </h3>
            <button
              class="px-2 py-1 text-xs rounded border border-slate-600 text-slate-400 hover:border-slate-500 hover:text-slate-300 transition-colors"
              onclick={exportModuleTsv}
            >TSV</button>
          </div>
          <div class="flex-1 min-h-0">
            <DataTable
              columns={moduleTableColumns}
              rows={selectedMagModules}
              maxHeight="100%"
              hideExport={true}
            />
          </div>
        {:else}
          <p class="text-slate-500 text-sm py-8 text-center">Click a bin row to view its modules</p>
        {/if}
      </div>
    </div>
  </div>

  <div class="mt-4 text-xs text-slate-500">
    Ward hierarchical clustering on both axes. Color: dark = 0%, bright = 100% completeness.
  </div>
{/if}
