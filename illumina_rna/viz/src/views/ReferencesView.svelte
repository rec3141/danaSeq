<script>
  import { references, taxonomy } from '../stores/data.js';

  let selectedName = $state(null);
  let sortKey = $state('total_length');
  let sortDir = $state('desc');

  let rows = $derived(($references || []).map(r => ({
    name: r.name,
    n_contigs: r.n_contigs,
    total_length: r.total_length,
    completeness: r.quality?.completeness ?? null,
    contamination: r.quality?.contamination ?? null,
    gc: r.quality?.gc != null ? r.quality.gc * 100 : null,
    n50: r.quality?.contig_n50 ?? null,
    n_cds: r.quality?.n_cds ?? null,
    total_reads: (r.contigs || []).reduce((s, c) => s + (c.reads || 0), 0),
    contigs: r.contigs || [],
    taxonomy: $taxonomy?.[r.name]?.lineage ?? null,
    ani: $taxonomy?.[r.name]?.ani ?? null,
  })));

  let sortedRows = $derived([...rows].sort((a, b) => {
    const av = a[sortKey], bv = b[sortKey];
    if (av == null && bv == null) return 0;
    if (av == null) return 1;
    if (bv == null) return -1;
    return sortDir === 'desc' ? bv - av : av - bv;
  }));

  let selected = $derived(rows.find(r => r.name === selectedName) || null);

  function fmt(n) { return n == null ? '—' : typeof n === 'number' ? n.toLocaleString(undefined, { maximumFractionDigits: 2 }) : n; }
  function pct(n) { return n == null ? '—' : n.toFixed(1) + '%'; }
  function setSort(k) {
    if (sortKey === k) sortDir = sortDir === 'desc' ? 'asc' : 'desc';
    else { sortKey = k; sortDir = 'desc'; }
  }
  function chev(k) { return sortKey === k ? (sortDir === 'desc' ? ' ▼' : ' ▲') : ''; }

  function qualBadge(comp, cont) {
    if (comp == null) return null;
    // MIMAG quality tiers
    if (comp >= 90 && cont < 5)  return { label: 'high',   cls: 'bg-emerald-900/40 text-emerald-300' };
    if (comp >= 50 && cont < 10) return { label: 'medium', cls: 'bg-amber-900/40 text-amber-300' };
    return { label: 'low', cls: 'bg-slate-800 text-slate-400' };
  }
</script>

<div class="space-y-4">
  <div class="bg-slate-900 border border-slate-800 rounded-lg overflow-x-auto">
    <table class="text-xs w-full">
      <thead class="text-slate-500 bg-slate-900/80 sticky top-0">
        <tr>
          <th class="text-left px-3 py-2">MAG</th>
          <th class="text-left px-3">Taxonomy (sendsketch)</th>
          <th class="text-right px-3 cursor-pointer hover:text-slate-300" onclick={() => setSort('ani')}>ANI{chev('ani')}</th>
          <th class="text-right px-3 cursor-pointer hover:text-slate-300" onclick={() => setSort('total_length')}>Genome size{chev('total_length')}</th>
          <th class="text-right px-3 cursor-pointer hover:text-slate-300" onclick={() => setSort('completeness')}>Completeness{chev('completeness')}</th>
          <th class="text-right px-3 cursor-pointer hover:text-slate-300" onclick={() => setSort('contamination')}>Contam.{chev('contamination')}</th>
          <th class="text-right px-3 cursor-pointer hover:text-slate-300" onclick={() => setSort('gc')}>GC{chev('gc')}</th>
          <th class="text-right px-3 cursor-pointer hover:text-slate-300" onclick={() => setSort('n50')}>N50{chev('n50')}</th>
          <th class="text-right px-3 cursor-pointer hover:text-slate-300" onclick={() => setSort('n_cds')}>CDS{chev('n_cds')}</th>
          <th class="text-right px-3 cursor-pointer hover:text-slate-300" onclick={() => setSort('total_reads')}>RNA reads{chev('total_reads')}</th>
          <th class="px-3"></th>
        </tr>
      </thead>
      <tbody>
        {#each sortedRows as r}
          {@const q = qualBadge(r.completeness, r.contamination)}
          <tr class="border-t border-slate-800 hover:bg-slate-800/40 {selectedName === r.name ? 'bg-cyan-900/20' : ''}">
            <td class="px-3 py-1 text-slate-200 font-mono">
              {r.name}
              {#if q}<span class="ml-2 px-1.5 py-0.5 rounded text-[10px] uppercase {q.cls}">{q.label}</span>{/if}
            </td>
            <td class="px-3 text-slate-300 italic max-w-[260px] truncate" title={r.taxonomy}>{r.taxonomy || '—'}</td>
            <td class="px-3 text-right text-slate-300">{r.ani != null ? r.ani.toFixed(1) : '—'}</td>
            <td class="px-3 text-right text-slate-300">{fmt(r.total_length)}</td>
            <td class="px-3 text-right text-slate-300">{pct(r.completeness)}</td>
            <td class="px-3 text-right text-slate-300">{pct(r.contamination)}</td>
            <td class="px-3 text-right text-slate-300">{pct(r.gc)}</td>
            <td class="px-3 text-right text-slate-300">{fmt(r.n50)}</td>
            <td class="px-3 text-right text-slate-300">{fmt(r.n_cds)}</td>
            <td class="px-3 text-right text-cyan-300">{fmt(r.total_reads)}</td>
            <td class="px-3 text-right">
              <button class="text-cyan-400 hover:text-cyan-300 text-xs"
                onclick={() => selectedName = (selectedName === r.name ? null : r.name)}>
                {selectedName === r.name ? 'close' : 'contigs'}
              </button>
            </td>
          </tr>
        {/each}
      </tbody>
    </table>
  </div>

  {#if selected}
    <div class="bg-slate-900 border border-slate-800 rounded-lg p-4 max-h-[60vh] overflow-y-auto">
      <div class="text-sm font-semibold text-slate-200 mb-2">{selected.name} · per-contig RNA reads ({selected.contigs.length} contigs)</div>
      <table class="text-xs w-full">
        <thead class="text-slate-500 sticky top-0 bg-slate-900">
          <tr>
            <th class="text-left py-1">Contig</th>
            <th class="text-right">Length</th>
            <th class="text-right">Mapped reads</th>
          </tr>
        </thead>
        <tbody>
          {#each selected.contigs.toSorted((a, b) => (b.reads || 0) - (a.reads || 0)).slice(0, 200) as c}
            <tr class="border-t border-slate-800">
              <td class="py-1 text-slate-300 font-mono">{c.name}</td>
              <td class="py-1 text-right text-slate-200">{fmt(c.length)}</td>
              <td class="py-1 text-right text-cyan-300">{fmt(c.reads)}</td>
            </tr>
          {/each}
        </tbody>
      </table>
    </div>
  {/if}
</div>
