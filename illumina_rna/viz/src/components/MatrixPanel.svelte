<script>
  /**
   * Entities × samples expression: heatmap of the top N by mean TPM
   * (log10(TPM+1), sequential blue) above a sortable, filterable table.
   * Clicking a row (or a heatmap cell) shows that entity across samples.
   *
   * rows:    [{ id, tpm: number[], ...extra }]
   * columns: [{ key, label, fmt?, align?, link? }] extra table columns
   */
  import PlotlyChart from './PlotlyChart.svelte';
  import { SEQ_BLUE, mean, fmtTpm } from '../stores/data.js';

  let { rows = [], samples = [], columns = [], title = '', subtitle = '', topN = 40,
        idLabel = 'ID', labelOf = r => r.id } = $props();

  let filter = $state('');
  let sortKey = $state('mean');
  let sortDir = $state(-1);
  let selected = $state(null);

  let enriched = $derived(rows.map(r => ({ ...r, mean: mean(r.tpm), max: Math.max(...r.tpm) })));
  let top = $derived([...enriched].sort((a, b) => b.mean - a.mean).slice(0, topN));

  let heat = $derived([{
    type: 'heatmap',
    z: top.map(r => r.tpm.map(v => Math.log10(v + 1))),
    customdata: top.map(r => r.tpm),
    x: samples,
    y: top.map(r => labelOf(r)),
    colorscale: SEQ_BLUE,
    zmin: 0,
    colorbar: { title: { text: 'log₁₀(TPM+1)', side: 'right' }, thickness: 10, outlinewidth: 0 },
    xgap: 1, ygap: 1,
    hovertemplate: '%{y}<br>%{x}<br>TPM %{customdata:.1f}<extra></extra>',
  }]);
  let heatLayout = $derived({
    xaxis: { tickangle: -60, tickfont: { size: 8 }, automargin: true, showgrid: false },
    yaxis: { autorange: 'reversed', tickfont: { size: 9 }, automargin: true, showgrid: false },
    margin: { l: 10, r: 10, t: 10, b: 100 },
    dragmode: false,
  });

  let shown = $derived.by(() => {
    const q = filter.trim().toLowerCase();
    const f = q ? enriched.filter(r => [r.id, ...columns.map(c => r[c.key])].some(v => v != null && String(v).toLowerCase().includes(q))) : enriched;
    const k = sortKey;
    return [...f].sort((a, b) => {
      const x = a[k], y = b[k];
      if (x == null) return 1;
      if (y == null) return -1;
      return (typeof x === 'number' ? x - y : String(x).localeCompare(String(y))) * sortDir;
    }).slice(0, 500);
  });

  function sortBy(k) {
    if (sortKey === k) sortDir = -sortDir; else { sortKey = k; sortDir = -1; }
  }

  let profile = $derived(selected ? [{
    type: 'bar', x: samples, y: selected.tpm,
    marker: { color: '#3987e5' },
    hovertemplate: '%{x}<br>TPM %{y:.1f}<extra></extra>',
  }] : []);
  let profileLayout = {
    xaxis: { tickangle: -60, tickfont: { size: 8 }, automargin: true },
    yaxis: { title: 'TPM', gridcolor: '#1e293b' },
    margin: { l: 60, r: 10, t: 10, b: 100 }, bargap: 0.2,
  };

  function onHeatClick(ev) {
    const y = ev?.points?.[0]?.pointIndex?.[0];
    if (y != null) selected = top[y];
  }
</script>

<div class="space-y-4">
  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
    <div class="text-sm font-semibold text-slate-200">{title}</div>
    {#if subtitle}<div class="text-xs text-slate-400 mb-2">{subtitle}</div>{/if}
    {#if top.length}
      <PlotlyChart data={heat} layout={heatLayout} config={{ displayModeBar: false }} height="{Math.max(360, 18 * top.length + 140)}px" onclick={onHeatClick} />
    {:else}
      <div class="text-slate-500 text-sm py-8 text-center">No data for this reference.</div>
    {/if}
  </div>

  {#if selected}
    <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
      <div class="flex justify-between items-start">
        <div>
          <div class="text-sm font-semibold text-slate-200 font-mono">{labelOf(selected)}</div>
          <div class="text-xs text-slate-400">mean TPM {fmtTpm(selected.mean)} · max {fmtTpm(selected.max)}</div>
        </div>
        <button class="text-xs text-slate-400 hover:text-slate-200" onclick={() => (selected = null)}>close</button>
      </div>
      <PlotlyChart data={profile} layout={profileLayout} height="300px" />
    </div>
  {/if}

  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4">
    <div class="flex items-center justify-between mb-2 gap-4">
      <div class="text-xs text-slate-400">{enriched.length.toLocaleString()} rows{shown.length < enriched.length ? `, showing ${shown.length}` : ''}</div>
      <input class="bg-slate-950 border border-slate-700 rounded px-2 py-1 text-xs text-slate-200 w-64" placeholder="Filter…" bind:value={filter} />
    </div>
    <div class="overflow-auto max-h-[600px]">
      <table class="w-full text-xs tabular-nums">
        <thead class="text-slate-400 sticky top-0 bg-slate-900">
          <tr>
            <th class="text-left p-1 cursor-pointer" onclick={() => sortBy('id')}>{idLabel}</th>
            {#each columns as c}
              <th class="p-1 cursor-pointer {c.align === 'right' ? 'text-right' : 'text-left'}" onclick={() => sortBy(c.key)}>{c.label}</th>
            {/each}
            <th class="text-right p-1 cursor-pointer" onclick={() => sortBy('mean')}>Mean TPM</th>
            <th class="text-right p-1 cursor-pointer" onclick={() => sortBy('max')}>Max TPM</th>
          </tr>
        </thead>
        <tbody class="text-slate-200">
          {#each shown as r (r.id)}
            <tr class="border-t border-slate-800 hover:bg-slate-800/50 cursor-pointer {selected?.id === r.id ? 'bg-slate-800' : ''}" onclick={() => (selected = r)}>
              <td class="p-1 font-mono whitespace-nowrap">{r.id}</td>
              {#each columns as c}
                <td class="p-1 {c.align === 'right' ? 'text-right' : ''} {c.wrap ? '' : 'whitespace-nowrap'}">
                  {#if c.link && c.link(r)}
                    <a class="text-cyan-300 hover:underline" href={c.link(r)} target="_blank" rel="noopener" onclick={e => e.stopPropagation()}>{c.fmt ? c.fmt(r[c.key], r) : r[c.key] ?? '—'}</a>
                  {:else}
                    {c.fmt ? c.fmt(r[c.key], r) : r[c.key] ?? '—'}
                  {/if}
                </td>
              {/each}
              <td class="text-right p-1">{fmtTpm(r.mean)}</td>
              <td class="text-right p-1">{fmtTpm(r.max)}</td>
            </tr>
          {/each}
        </tbody>
      </table>
    </div>
  </div>
</div>
