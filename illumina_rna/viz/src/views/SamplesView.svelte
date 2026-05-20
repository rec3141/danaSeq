<script>
  import { samples, overview } from '../stores/data.js';

  let selectedId = $state(null);
  let selected = $derived($samples.find(s => s.id === selectedId) || null);

  function pct(x) { return x == null ? '—' : (x * 100).toFixed(1) + '%'; }
  function fmt(n) { return n == null ? '—' : n.toLocaleString(); }
</script>

<div class="grid grid-cols-3 gap-4">
  <div class="col-span-1 bg-slate-900 border border-slate-800 rounded-lg p-3 max-h-[80vh] overflow-y-auto">
    <div class="text-xs uppercase text-slate-500 mb-2 px-1">Samples ({$samples.length})</div>
    {#each $samples as s}
      <button
        class="block w-full text-left px-3 py-2 rounded text-sm
               {selectedId === s.id ? 'bg-cyan-900/40 text-cyan-200' : 'hover:bg-slate-800 text-slate-300'}"
        onclick={() => selectedId = s.id}>
        {s.id}
      </button>
    {/each}
  </div>

  <div class="col-span-2 space-y-4">
    {#if !selected}
      <div class="bg-slate-900 border border-slate-800 rounded-lg p-6 text-slate-500 text-sm">
        Select a sample on the left.
      </div>
    {:else}
      <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
        <div class="text-sm font-semibold text-slate-200 mb-2">{selected.id} · read flow</div>
        <table class="text-xs w-full">
          <thead class="text-slate-500">
            <tr><th class="text-left py-1">Stage</th><th class="text-right">Reads</th></tr>
          </thead>
          <tbody>
            {#each Object.entries(selected.read_flow || {}) as [stage, n]}
              <tr class="border-t border-slate-800">
                <td class="py-1 text-slate-300">{stage}</td>
                <td class="py-1 text-right font-mono text-slate-200">{fmt(n)}</td>
              </tr>
            {/each}
          </tbody>
        </table>
      </div>

      <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
        <div class="text-sm font-semibold text-slate-200 mb-2">{selected.id} · per-reference mapping</div>
        <table class="text-xs w-full">
          <thead class="text-slate-500">
            <tr>
              <th class="text-left py-1">Reference</th>
              <th class="text-right">Total reads</th>
              <th class="text-right">Mapped</th>
              <th class="text-right">Rate</th>
            </tr>
          </thead>
          <tbody>
            {#each Object.entries(selected.per_reference || {}) as [ref, st]}
              <tr class="border-t border-slate-800">
                <td class="py-1 text-slate-300 font-mono">{ref}</td>
                <td class="py-1 text-right text-slate-200">{fmt(st.total_reads)}</td>
                <td class="py-1 text-right text-slate-200">{fmt(st.mapped_reads)}</td>
                <td class="py-1 text-right text-cyan-300">{pct(st.mapping_rate)}</td>
              </tr>
            {/each}
          </tbody>
        </table>
      </div>
    {/if}
  </div>
</div>
