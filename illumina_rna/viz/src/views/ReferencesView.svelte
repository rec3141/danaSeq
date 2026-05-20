<script>
  import { references } from '../stores/data.js';

  let selectedName = $state(null);
  let selected = $derived($references.find(r => r.name === selectedName) || null);

  function fmt(n) { return n == null ? '—' : n.toLocaleString(); }
</script>

<div class="grid grid-cols-3 gap-4">
  <div class="col-span-1 bg-slate-900 border border-slate-800 rounded-lg p-3 max-h-[80vh] overflow-y-auto">
    <div class="text-xs uppercase text-slate-500 mb-2 px-1">References ({$references.length})</div>
    {#each $references as r}
      <button
        class="block w-full text-left px-3 py-2 rounded text-sm
               {selectedName === r.name ? 'bg-cyan-900/40 text-cyan-200' : 'hover:bg-slate-800 text-slate-300'}"
        onclick={() => selectedName = r.name}>
        <div class="font-mono">{r.name}</div>
        <div class="text-xs text-slate-500">{r.n_contigs} contigs · {fmt(r.total_length)} bp</div>
      </button>
    {/each}
  </div>

  <div class="col-span-2 space-y-4">
    {#if !selected}
      <div class="bg-slate-900 border border-slate-800 rounded-lg p-6 text-slate-500 text-sm">
        Select a reference on the left.
      </div>
    {:else}
      <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
        <div class="text-sm font-semibold text-slate-200 mb-2">{selected.name}</div>
        <div class="grid grid-cols-2 gap-4 text-xs">
          <div>
            <div class="text-slate-500">Contigs</div>
            <div class="text-slate-200 text-lg font-mono">{selected.n_contigs}</div>
          </div>
          <div>
            <div class="text-slate-500">Total length</div>
            <div class="text-slate-200 text-lg font-mono">{fmt(selected.total_length)} bp</div>
          </div>
        </div>
      </div>

      <div class="bg-slate-900 border border-slate-800 rounded-lg p-4 max-h-[60vh] overflow-y-auto">
        <div class="text-sm font-semibold text-slate-200 mb-2">Per-contig read totals</div>
        <table class="text-xs w-full">
          <thead class="text-slate-500 sticky top-0 bg-slate-900">
            <tr>
              <th class="text-left py-1">Contig</th>
              <th class="text-right">Length</th>
              <th class="text-right">Mapped reads</th>
            </tr>
          </thead>
          <tbody>
            {#each selected.contigs as c}
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
</div>
