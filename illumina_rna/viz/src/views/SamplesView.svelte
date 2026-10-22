<script>
  import { samples, fmtInt, fmtPct } from '../stores/data.js';

  const cols = [
    ['reads_raw', 'Raw reads', fmtInt], ['reads_filtered', 'After QC', fmtInt],
    ['pct_rrna', 'rRNA', fmtPct], ['reads_nohuman', 'Non-rRNA, non-human', fmtInt],
    ['reads_mapped', 'Mapped', fmtInt], ['pct_mapped', 'Mapped %', fmtPct],
    ['pairs_assigned', 'Pairs on genes', fmtInt], ['pct_assigned', 'Assigned %', fmtPct],
    ['genes_detected', 'Genes detected', fmtInt],
  ];
  let sortKey = $state('sample');
  let dir = $state(1);
  let rows = $derived(($samples?.rows ?? []).map((r, i) => ({ ...r, genes_detected: $samples.genes_detected?.[i] })));
  let sorted = $derived([...rows].sort((a, b) => {
    const x = a[sortKey], y = b[sortKey];
    return (typeof x === 'number' ? x - y : String(x).localeCompare(String(y))) * dir;
  }));
  const by = k => { if (sortKey === k) dir = -dir; else { sortKey = k; dir = k === 'sample' ? 1 : -1; } };
</script>

<div class="bg-slate-900 border border-slate-800 rounded-lg p-4 overflow-auto">
  <table class="w-full text-xs tabular-nums">
    <thead class="text-slate-400 sticky top-0 bg-slate-900">
      <tr>
        <th class="text-left p-1 cursor-pointer" onclick={() => by('sample')}>Sample</th>
        {#each cols as [k, label]}<th class="text-right p-1 cursor-pointer" onclick={() => by(k)}>{label}</th>{/each}
      </tr>
    </thead>
    <tbody class="text-slate-200">
      {#each sorted as r}
        <tr class="border-t border-slate-800">
          <td class="p-1 font-mono">{r.sample}</td>
          {#each cols as [k, , f]}<td class="text-right p-1">{f(r[k])}</td>{/each}
        </tr>
      {/each}
    </tbody>
  </table>
</div>
