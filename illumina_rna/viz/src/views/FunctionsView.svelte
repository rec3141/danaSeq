<script>
  import { functions_ } from '../stores/data.js';

  let refKeys = $derived(Object.keys($functions_ || {}).sort());
  let selectedRef = $state(null);
  let query = $state('');
  let sortKey = $state('total');
  let sortDir = $state('desc');

  $effect(() => {
    if (!selectedRef && refKeys.length) selectedRef = refKeys[0];
  });

  let rows = $derived.by(() => {
    const d = $functions_?.[selectedRef];
    if (!d) return { samples: [], rows: [] };
    const out = [];
    for (let i = 0; i < d.genes.length; i++) {
      const counts = d.matrix[i];
      const total = counts.reduce((s, x) => s + x, 0);
      const a = d.annotations[i] || {};
      out.push({
        gene_id: d.genes[i],
        product: a.product || '—',
        gene:    a.gene    || '',
        seqid:   a.seqid   || '',
        total,
        counts,
      });
    }
    return { samples: d.samples, rows: out };
  });

  let filtered = $derived.by(() => {
    let r = rows.rows;
    if (query.trim()) {
      const q = query.toLowerCase();
      r = r.filter(x =>
        x.gene_id.toLowerCase().includes(q) ||
        (x.product || '').toLowerCase().includes(q) ||
        (x.gene    || '').toLowerCase().includes(q)
      );
    }
    r = [...r].sort((a, b) => {
      const av = a[sortKey], bv = b[sortKey];
      if (typeof av === 'number') return sortDir === 'desc' ? bv - av : av - bv;
      return sortDir === 'desc' ? String(bv).localeCompare(String(av)) : String(av).localeCompare(String(bv));
    });
    return r.slice(0, 500);
  });

  function setSort(k) {
    if (sortKey === k) sortDir = sortDir === 'desc' ? 'asc' : 'desc';
    else { sortKey = k; sortDir = 'desc'; }
  }
  function chev(k) { return sortKey === k ? (sortDir === 'desc' ? ' ▼' : ' ▲') : ''; }
  function fmt(n) { return n == null ? '—' : n.toLocaleString(); }

  // Color cell by count: log-scale, slate (0) → cyan (high)
  function cellBg(count, max) {
    if (!count || max <= 0) return 'transparent';
    const t = Math.log2(count + 1) / Math.log2(max + 1);
    const alpha = (0.05 + 0.6 * t).toFixed(2);
    return `rgba(34, 211, 238, ${alpha})`;
  }
</script>

<div class="space-y-3">
  <div class="bg-slate-900 border border-slate-800 rounded-lg p-3 flex flex-wrap items-center gap-4">
    <label class="text-sm text-slate-400">
      MAG:
      <select bind:value={selectedRef} class="ml-2 bg-slate-800 border border-slate-700 rounded px-2 py-1 text-slate-200 text-sm">
        {#each refKeys as k}<option value={k}>{k}</option>{/each}
      </select>
    </label>
    <label class="text-sm text-slate-400 flex-1 min-w-[200px]">
      Search:
      <input type="text" bind:value={query} placeholder="locus / gene / product"
             class="ml-2 bg-slate-800 border border-slate-700 rounded px-2 py-1 text-slate-200 text-sm w-64"/>
    </label>
    <div class="text-xs text-slate-500">{filtered.length} of {rows.rows.length} genes (capped at 500)</div>
  </div>

  {#if refKeys.length === 0}
    <div class="text-slate-500 text-sm py-8 text-center bg-slate-900 border border-slate-800 rounded-lg">
      No function data. Run with annotated references (Bakta GFFs).
    </div>
  {:else}
    {@const maxc = Math.max(...filtered.flatMap(r => r.counts))}
    <div class="bg-slate-900 border border-slate-800 rounded-lg overflow-x-auto max-h-[75vh] overflow-y-auto">
      <table class="text-xs">
        <thead class="text-slate-500 bg-slate-900 sticky top-0 z-10">
          <tr>
            <th class="text-left px-3 py-2 cursor-pointer hover:text-slate-300" onclick={() => setSort('gene_id')}>Locus{chev('gene_id')}</th>
            <th class="text-left px-3 cursor-pointer hover:text-slate-300" onclick={() => setSort('gene')}>Gene{chev('gene')}</th>
            <th class="text-left px-3 cursor-pointer hover:text-slate-300" onclick={() => setSort('product')}>Product{chev('product')}</th>
            <th class="text-right px-3 cursor-pointer hover:text-slate-300" onclick={() => setSort('total')}>Σ counts{chev('total')}</th>
            {#each rows.samples as s}
              <th class="text-right px-2 font-normal text-slate-600" title={s}>{s.replace(/^V?\d{4,}_/, '').slice(0, 7)}</th>
            {/each}
          </tr>
        </thead>
        <tbody>
          {#each filtered as r}
            <tr class="border-t border-slate-800/60 hover:bg-slate-800/40">
              <td class="px-3 py-1 text-slate-300 font-mono">{r.gene_id}</td>
              <td class="px-3 text-cyan-400">{r.gene}</td>
              <td class="px-3 text-slate-300 max-w-[280px] truncate" title={r.product}>{r.product}</td>
              <td class="px-3 text-right text-cyan-300 font-medium">{fmt(r.total)}</td>
              {#each r.counts as c}
                <td class="text-right px-2 text-slate-200" style:background-color={cellBg(c, maxc)}>
                  {c || ''}
                </td>
              {/each}
            </tr>
          {/each}
        </tbody>
      </table>
    </div>
  {/if}
</div>
