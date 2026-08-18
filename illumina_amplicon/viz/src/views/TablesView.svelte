<script>
  import { store, GROUP_HEX, buildTaxColorMap, getAsvColor, getEffectiveColorLevel, getClusterColor, makeTaxonMatcher } from '../stores/data.svelte.js';

  let { filters = {} } = $props();

  let activeTable = $state('asvs');
  let search = $state('');
  let sortCol = $state(null);
  let sortAsc = $state(true);

  // ── Taxonomy helpers ──
  let primaryDb = $derived(Object.keys(store.taxonomy)[0] || null);
  let taxLevels = $derived(primaryDb ? (store.taxonomy[primaryDb]?.levels || []) : []);
  let taxAssignments = $derived(primaryDb ? (store.taxonomy[primaryDb]?.assignments || {}) : {});
  let taxBootstraps = $derived(primaryDb ? (store.taxonomy[primaryDb]?.bootstraps || {}) : {});

  let effectiveColorLevel = $derived(
    filters.colorMode === 'taxonomy'
      ? getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter, filters.taxonFilterLevel)
      : null
  );
  let taxCmap = $derived(
    effectiveColorLevel && effectiveColorLevel !== 'group'
      ? buildTaxColorMap(effectiveColorLevel, filters.taxonFilter, filters.taxonFilterLevel).colorMap
      : null
  );

  // ── Dynamic columns ──
  let asvCols = $derived.by(() => {
    const cols = [
      { key: '_color', label: '', width: '30px' },
      { key: 'id', label: 'ASV' },
      { key: 'total_reads', label: 'Reads', numeric: true },
      { key: 'n_samples', label: 'Prevalence', numeric: true },
    ];

    if (store.heatmap?.asvClusters) {
      const k = filters.asvClusterK || 4;
      cols.push({ key: '_asvCluster', label: `Cluster (k=${k})`, numeric: true });
    }

    cols.push({ key: 'group', label: 'Group' });

    for (const level of taxLevels) {
      cols.push({ key: `_tax_${level}`, label: level });
    }

    return cols;
  });

  let sampleCols = $derived.by(() => {
    const cols = [
      { key: '_color', label: '', width: '30px' },
      { key: 'id', label: 'Sample ID' },
      { key: 'total_reads', label: 'Reads', numeric: true },
      { key: 'n_asvs', label: 'Richness', numeric: true },
    ];

    if (store.heatmap?.sampleClusters) {
      const k = filters.sampleClusterK || 4;
      cols.push({ key: '_sampleCluster', label: `Sample Cluster (k=${k})`, numeric: true });
    }

    return cols;
  });

  // ASV x sample counts. One row per ASV, one column per sample, with the
  // taxonomy carried alongside so the exported file stands on its own — a matrix
  // of bare counts is not something you can read six months later.
  const MATRIX_ROW_CAP = 250;

  let sampleIds = $derived(store.counts?.samples || []);

  let matrixCols = $derived.by(() => {
    const cols = [
      { key: '_color', label: '', width: '30px' },
      { key: 'id', label: 'ASV' },
    ];
    for (const level of taxLevels) cols.push({ key: `_tax_${level}`, label: level });
    cols.push({ key: 'total_reads', label: 'Total', numeric: true });
    for (const sid of sampleIds) {
      cols.push({ key: `_s_${sid}`, label: sid, numeric: true, sample: true });
    }
    return cols;
  });

  // sequence id -> {sample id -> reads}, built once from the packed counts.
  let countsByAsv = $derived.by(() => {
    const c = store.counts;
    if (!c?.data) return {};
    const out = {};
    for (const [si, ai, reads] of c.data) {
      const aid = c.asvs[ai];
      (out[aid] ||= {})[c.samples[si]] = reads;
    }
    return out;
  });

  let cols = $derived(
    activeTable === 'asvs' ? asvCols
    : activeTable === 'samples' ? sampleCols
    : matrixCols
  );

  // ── Build enriched rows ──
  let rawRows = $derived.by(() => {
    if (activeTable === 'matrix') {
      const per = countsByAsv;
      return store.asvs.map(a => {
        const row = { id: a.id, total_reads: a.total_reads, group: a.group };
        const tax = taxAssignments[a.id];
        if (tax) {
          for (let i = 0; i < taxLevels.length; i++) row[`_tax_${taxLevels[i]}`] = tax[i] || '';
        }
        const counts = per[a.id] || {};
        for (const sid of sampleIds) row[`_s_${sid}`] = counts[sid] ?? 0;
        row._color = GROUP_HEX[a.group] || GROUP_HEX.unknown;
        return row;
      });
    }
    if (activeTable === 'asvs') {
      const asvClusters = store.heatmap?.asvClusters?.[String(filters.asvClusterK || 4)] || {};
      return store.asvs.map(a => {
        const row = { ...a };
        // Add taxonomy levels
        const tax = taxAssignments[a.id];
        if (tax) {
          for (let i = 0; i < taxLevels.length; i++) {
            row[`_tax_${taxLevels[i]}`] = tax[i] || '';
          }
        }
        // Add cluster
        row._asvCluster = asvClusters[a.id] || '';
        // Add color
        if (filters.colorMode === 'cluster') {
          row._color = getClusterColor(a.id, 'asvCluster', filters.asvClusterK);
        } else if (filters.colorMode === 'group') {
          row._color = GROUP_HEX[a.group] || GROUP_HEX.unknown;
        } else if (taxCmap) {
          row._color = getAsvColor(a.id, effectiveColorLevel, taxCmap);
        } else {
          row._color = GROUP_HEX[a.group] || GROUP_HEX.unknown;
        }
        return row;
      });
    } else {
      const sampleClusters = store.heatmap?.sampleClusters?.[String(filters.sampleClusterK || 4)] || {};
      return store.samples.map(s => {
        const row = { ...s };
        row._sampleCluster = sampleClusters[s.id] || '';
        if (filters.colorMode === 'cluster') {
          row._color = getClusterColor(s.id, 'sampleCluster', filters.sampleClusterK);
        } else {
          row._color = '#475569';
        }
        return row;
      });
    }
  });

  // ── Filter by search and taxonomy filter ──
  let filteredRows = $derived.by(() => {
    let rows = rawRows;

    // Taxonomy filter
    if (filters.taxonFilter && activeTable !== 'samples') {
      const matches = makeTaxonMatcher(filters.taxonFilter, filters.taxonFilterLevel);
      rows = rows.filter(r => matches(r.id ?? '', r.taxonomy ?? ''));
    }

    // Search filter
    if (search) {
      const q = search.toLowerCase();
      rows = rows.filter(r =>
        Object.entries(r).some(([k, v]) =>
          k !== '_color' && v != null && String(v).toLowerCase().includes(q)
        )
      );
    }

    // Sort
    if (sortCol && sortCol !== '_color') {
      const col = cols.find(c => c.key === sortCol);
      const numeric = col?.numeric ?? false;
      rows = [...rows].sort((a, b) => {
        let va = a[sortCol] ?? '';
        let vb = b[sortCol] ?? '';
        if (numeric) {
          va = Number(va) || 0;
          vb = Number(vb) || 0;
          return sortAsc ? va - vb : vb - va;
        }
        va = String(va).toLowerCase();
        vb = String(vb).toLowerCase();
        return sortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
      });
    }

    return rows;
  });

  // A 3,000 x 66 grid is 200,000 cells and the browser will not thank you for
  // it. Draw the head of the table and export the whole thing.
  let pageRows = $derived(
    activeTable === 'matrix' ? filteredRows.slice(0, MATRIX_ROW_CAP) : filteredRows
  );
  let matrixTruncated = $derived(
    activeTable === 'matrix' && filteredRows.length > MATRIX_ROW_CAP
  );

  function toggleSort(key) {
    if (key === '_color') return;
    if (sortCol === key) { sortAsc = !sortAsc; }
    else { sortCol = key; sortAsc = true; }
  }

  function exportCsv() {
    // filteredRows, not pageRows: the drawn table is capped for the browser's
    // sake, and an export that silently stopped at 250 rows would be a trap.
    const exportCols = cols.filter(c => c.key !== '_color');
    const header = exportCols.map(c => c.label).join(',');
    const body = filteredRows.map(r =>
      exportCols.map(c => {
        const v = r[c.key] ?? '';
        return String(v).includes(',') ? `"${v}"` : v;
      }).join(',')
    ).join('\n');
    const blob = new Blob([header + '\n' + body], { type: 'text/csv' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    const name = activeTable === 'matrix' ? 'asv_by_sample' : activeTable;
    a.href = url; a.download = `${name}.csv`; a.click();
    URL.revokeObjectURL(url);
  }
</script>

<div class="flex h-full flex-col overflow-hidden p-4">
  <div class="mb-4 flex flex-wrap items-center gap-3">
    <div class="flex rounded-lg border border-slate-700 bg-slate-800">
      <button
        class="px-4 py-1.5 text-sm font-medium transition-colors rounded-l-lg
          {activeTable === 'asvs' ? 'bg-blue-600 text-white' : 'text-slate-400 hover:text-slate-200'}"
        onclick={() => { activeTable = 'asvs'; sortCol = null; }}
      >ASVs</button>
      <button
        class="px-4 py-1.5 text-sm font-medium transition-colors rounded-r-lg
          {activeTable === 'samples' ? 'bg-blue-600 text-white' : 'text-slate-400 hover:text-slate-200'}"
        onclick={() => { activeTable = 'samples'; sortCol = null; }}
      >Samples</button>
      <button
        class="rounded px-3 py-1.5 text-sm font-medium transition-colors
          {activeTable === 'matrix' ? 'bg-blue-600 text-white' : 'text-slate-400 hover:text-slate-200'}"
        onclick={() => { activeTable = 'matrix'; sortCol = null; }}
      >ASV &times; Sample</button>
    </div>

    <input
      type="text" bind:value={search} placeholder="Search..."
      class="rounded border border-slate-700 bg-slate-800 px-3 py-1.5 text-sm text-slate-200 placeholder-slate-500 focus:border-blue-500 focus:outline-none w-64"
    />

    <span class="text-xs text-slate-500">{filteredRows.length} rows</span>

    <button
      class="ml-auto rounded border border-slate-700 bg-slate-800 px-3 py-1.5 text-xs font-medium text-slate-300 hover:bg-slate-700 hover:text-slate-100 transition-colors"
      onclick={exportCsv}
    >Export CSV</button>
  </div>

  {#if matrixTruncated}
    <p class="mb-2 text-xs text-slate-500">
      Showing the first {MATRIX_ROW_CAP} of {filteredRows.length.toLocaleString()} ASVs —
      export includes all of them, with taxonomy.
    </p>
  {/if}

  <div class="flex-1 overflow-auto rounded-lg border border-slate-800">
    <table class="w-full text-sm">
      <thead class="sticky top-0 z-10 bg-slate-900 text-left">
        <tr>
          {#each cols as col}
            <th
              class="cursor-pointer select-none border-b border-slate-800 px-3 py-2 text-xs font-semibold uppercase tracking-wider text-slate-400 hover:text-slate-200 transition-colors"
              style={col.width ? `width:${col.width}` : ''}
              onclick={() => toggleSort(col.key)}
            >
              {#if col.key !== '_color'}
                <span class="inline-flex items-center gap-1">
                  {col.label}
                  {#if sortCol === col.key}
                    <span class="text-blue-400">{sortAsc ? '▲' : '▼'}</span>
                  {/if}
                </span>
              {/if}
            </th>
          {/each}
        </tr>
      </thead>
      <tbody>
        {#each pageRows as row}
          <tr class="border-t border-slate-800/50 hover:bg-slate-800/30 transition-colors">
            {#each cols as col}
              <td class="px-3 py-1.5 text-slate-300 {col.numeric ? 'text-right font-mono' : ''}">
                {#if col.key === '_color'}
                  <span class="inline-block h-3 w-3 rounded-full" style="background:{row._color}"></span>
                {:else if col.key === 'group' && row[col.key]}
                  <span class="inline-flex items-center gap-1.5 capitalize">
                    <span class="inline-block h-2 w-2 rounded-full" style="background:{GROUP_HEX[row[col.key]] ?? '#888'}"></span>
                    {row[col.key]}
                  </span>
                {:else if col.numeric}
                  {(row[col.key] ?? 0).toLocaleString()}
                {:else}
                  <span class="max-w-xs truncate block text-xs">{row[col.key] ?? ''}</span>
                {/if}
              </td>
            {/each}
          </tr>
        {/each}

        {#if pageRows.length === 0}
          <tr>
            <td colspan={cols.length} class="px-4 py-8 text-center text-slate-500">
              No data available.
            </td>
          </tr>
        {/if}
      </tbody>
    </table>
  </div>

</div>
