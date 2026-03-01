<script>
  let { columns = [], rows = [], onRowClick = null, selectedId = null, idKey = 'id',
        maxHeight = '400px', hideExport = false, actionLabel = null, actionFn = null, actionStyle = null } = $props();

  let sortCol = $state(null);
  let sortAsc = $state(true);

  function handleSort(col) {
    if (sortCol === col) { sortAsc = !sortAsc; }
    else { sortCol = col; sortAsc = true; }
  }

  let sortedRows = $derived.by(() => {
    if (!sortCol) return rows;
    const key = sortCol;
    return [...rows].sort((a, b) => {
      const va = a[key], vb = b[key];
      if (va == null) return 1;
      if (vb == null) return -1;
      if (typeof va === 'number') return sortAsc ? va - vb : vb - va;
      return sortAsc ? String(va).localeCompare(String(vb)) : String(vb).localeCompare(String(va));
    });
  });

  function exportTsv() {
    const header = columns.map(c => c.label).join('\t');
    const body = sortedRows.map(row =>
      columns.map(c => row[c.key] != null ? row[c.key] : '').join('\t')
    ).join('\n');
    const blob = new Blob([header + '\n' + body], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = 'table.tsv'; a.click();
    URL.revokeObjectURL(url);
  }
</script>

<div class="overflow-auto rounded-lg border border-slate-700" style="max-height: {maxHeight}">
  <table class="w-full text-sm">
    <thead class="sticky top-0 bg-slate-800 z-10">
      <tr>
        {#each columns as col}
          <th
            class="px-3 py-2 text-left text-xs font-medium text-slate-400 uppercase tracking-wider cursor-pointer hover:text-cyan-400 whitespace-nowrap"
            onclick={() => handleSort(col.key)}
          >
            {col.label}
            {#if sortCol === col.key}
              <span class="ml-1">{sortAsc ? '\u2191' : '\u2193'}</span>
            {/if}
          </th>
        {/each}
        {#if actionLabel}
          <th class="px-2 py-2 text-xs font-medium text-slate-400 uppercase tracking-wider"></th>
        {/if}
        {#if !hideExport}
          <th class="px-2 py-2 text-right">
            <button
              class="text-[10px] px-1.5 py-0.5 rounded border border-slate-600 text-slate-500 hover:text-slate-300 hover:border-slate-500 transition-colors"
              onclick={exportTsv} title="Export as TSV"
            >TSV</button>
          </th>
        {/if}
      </tr>
    </thead>
    <tbody class="divide-y divide-slate-700/50">
      {#each sortedRows as row}
        <tr
          class="hover:bg-slate-800/80 transition-colors {row[idKey] === selectedId ? 'bg-cyan-400/10' : ''}"
          class:cursor-pointer={!!onRowClick}
          onclick={() => onRowClick && onRowClick(row)}
        >
          {#each columns as col}
            <td class="px-3 py-2 whitespace-nowrap text-slate-300 font-mono text-xs">
              {#if col.render}
                {@html col.render(row[col.key], row)}
              {:else}
                {row[col.key] != null ? row[col.key] : '-'}
              {/if}
            </td>
          {/each}
          {#if actionLabel && actionFn}
            <td class="px-2 py-1 whitespace-nowrap">
              <button
                class={actionStyle ? actionStyle(row) : 'text-[10px] px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-cyan-400 hover:border-cyan-400/40 transition-colors'}
                onclick={(e) => { e.stopPropagation(); actionFn(row); }}
              >{typeof actionLabel === 'function' ? actionLabel(row) : actionLabel}</button>
            </td>
          {/if}
        </tr>
      {/each}
    </tbody>
  </table>
</div>
