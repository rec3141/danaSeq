<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import { store } from '../stores/data.svelte.js';

  let plotDiv = $state(null);
  let heatDiv = $state(null);
  let selected = $state('__all__');

  const prov = $derived(store.provenance);
  const stages = $derived(prov?.stages ?? []);
  const sampleIds = $derived(Object.keys(prov?.samples ?? {}).sort());

  // Counts for whatever is selected: the combined total, or one sample.
  const counts = $derived.by(() => {
    if (!prov) return null;
    return selected === '__all__' ? prov.total : prov.samples?.[selected];
  });

  // Retention relative to the raw FASTQ — the number that matters, since a
  // sample can look "shallow" when in fact primer removal discarded it.
  const rows = $derived.by(() => {
    if (!counts || !stages.length) return [];
    const raw = counts[stages[0].id] || 0;
    return stages.map(s => {
      const n = counts[s.id] ?? 0;
      return { id: s.id, label: s.label, n, pct: raw ? (100 * n) / raw : 0 };
    });
  });

  // Samples whose reads mostly died at primer removal — worth surfacing, it
  // usually means that run used a different primer pair than the one applied.
  const lostSamples = $derived.by(() => {
    if (!prov?.samples || stages.length < 2) return [];
    const [rawKey, primerKey] = [stages[0].id, stages[1].id];
    return Object.entries(prov.samples)
      .map(([id, c]) => ({ id, raw: c[rawKey] ?? 0, kept: c[primerKey] ?? 0 }))
      .filter(s => s.raw >= 1000 && s.kept / s.raw < 0.5)
      .map(s => ({ ...s, pct: (100 * s.kept) / s.raw }))
      .sort((a, b) => a.pct - b.pct);
  });

  function draw() {
    if (!plotDiv || !rows.length) return;
    Plotly.react(
      plotDiv,
      [{
        type: 'bar',
        x: rows.map(r => r.label),
        y: rows.map(r => r.n),
        marker: { color: rows.map((_, i) => ['#38bdf8', '#818cf8', '#34d399'][i % 3]) },
        text: rows.map(r => `${r.n.toLocaleString()}<br>${r.pct.toFixed(1)}%`),
        textposition: 'auto',
        hovertemplate: '%{x}<br>%{y:,} reads<extra></extra>',
      }],
      {
        margin: { l: 70, r: 20, t: 30, b: 60 },
        plot_bgcolor: 'rgba(2, 6, 15, 1)',
        paper_bgcolor: 'rgba(2, 6, 15, 1)',
        font: { color: '#94a3b8', size: 12 },
        xaxis: { gridcolor: '#1e293b' },
        yaxis: { title: 'reads', gridcolor: '#1e293b', rangemode: 'tozero' },
        showlegend: false,
      },
      { responsive: true, displaylogo: false },
    );
  }

  // All-samples heatmap: rows = samples, columns = stages, colour = log10(reads).
  // Gives a one-glance view of the whole dataset's funnel — a sample that
  // collapses at a stage (e.g. an outlier stripped by the prevalence filter)
  // shows up as a dark cell in that column. Rows are ordered by overall
  // retention so the worst-retained samples cluster together.
  function drawHeat() {
    if (!heatDiv || selected !== '__all__' || !stages.length || !sampleIds.length) return;
    const rawId = stages[0].id, finalId = stages[stages.length - 1].id;
    const retention = (id) => {
      const r = prov.samples[id][rawId] || 0;
      return r ? (prov.samples[id][finalId] || 0) / r : 0;
    };
    const order = [...sampleIds].sort((a, b) => retention(a) - retention(b));
    const x = stages.map(s => s.label);
    const z = order.map(id => stages.map(s => {
      const n = prov.samples[id][s.id] ?? 0;
      return n > 0 ? Math.log10(n) : null;   // null → blank cell for zero reads
    }));
    const raw = order.map(id => stages.map(s => prov.samples[id][s.id] ?? 0));
    Plotly.react(
      heatDiv,
      [{
        type: 'heatmap',
        x, y: order, z,
        customdata: raw,
        colorscale: 'Viridis',
        hoverongaps: false,
        hovertemplate: '%{y}<br>%{x}<br>%{customdata:,} reads<extra></extra>',
        colorbar: { title: { text: 'log₁₀(reads)', side: 'right' }, thickness: 12,
                    tickfont: { color: '#94a3b8', size: 10 } },
      }],
      {
        margin: { l: 110, r: 10, t: 10, b: 70 },
        plot_bgcolor: 'rgba(2, 6, 15, 1)',
        paper_bgcolor: 'rgba(2, 6, 15, 1)',
        font: { color: '#94a3b8', size: 11 },
        xaxis: { tickangle: -30, automargin: true },
        yaxis: { automargin: true, autorange: 'reversed' },
      },
      { responsive: true, displaylogo: false },
    );
  }

  onMount(() => { draw(); drawHeat(); });
  $effect(() => { rows; draw(); });
  $effect(() => { selected; stages; sampleIds; drawHeat(); });
</script>

<div class="flex h-full flex-col gap-3 overflow-y-auto p-3">
  {#if !prov}
    <div class="flex flex-1 items-center justify-center text-sm text-slate-500">
      No provenance data — <code class="mx-1">data/provenance.json</code> was not produced for this run.
    </div>
  {:else}
    <div class="flex items-center gap-3">
      <label for="prov-sample" class="text-xs text-slate-400">Sample</label>
      <select
        id="prov-sample"
        bind:value={selected}
        class="rounded border border-slate-700 bg-slate-900 px-2 py-1 text-xs text-slate-200"
      >
        <option value="__all__">All samples combined ({sampleIds.length})</option>
        {#each sampleIds as id}
          <option value={id}>{id}</option>
        {/each}
      </select>
      {#if rows.length}
        <span class="text-xs text-slate-500">
          {rows[0].n.toLocaleString()} raw &rarr; {rows[rows.length - 1].n.toLocaleString()} retained
          ({rows[rows.length - 1].pct.toFixed(1)}%)
        </span>
      {/if}
    </div>

    <div bind:this={plotDiv} class="h-72 shrink-0"></div>

    {#if selected === '__all__' && sampleIds.length}
      <div class="shrink-0">
        <p class="mb-1 text-xs text-slate-400">
          All samples &middot; log<sub>10</sub>(reads) at each stage
          <span class="text-slate-500">— rows sorted by retention (lowest first)</span>
        </p>
        <div
          bind:this={heatDiv}
          style="height: {Math.min(700, 70 + sampleIds.length * 16)}px"
        ></div>
      </div>
    {/if}

    <div class="max-h-40 overflow-y-auto">
      <table class="w-full text-xs">
        <thead class="text-slate-400">
          <tr><th class="p-1 text-left">Step</th><th class="p-1 text-right">Reads</th><th class="p-1 text-right">% of raw</th></tr>
        </thead>
        <tbody class="text-slate-300">
          {#each rows as r}
            <tr class="border-t border-slate-800">
              <td class="p-1">{r.label}</td>
              <td class="p-1 text-right">{r.n.toLocaleString()}</td>
              <td class="p-1 text-right" class:text-rose-400={r.pct < 50}>{r.pct.toFixed(1)}%</td>
            </tr>
          {/each}
        </tbody>
      </table>

      {#if selected === '__all__' && lostSamples.length}
        <p class="mt-2 text-xs text-amber-400">
          {lostSamples.length} sample{lostSamples.length === 1 ? '' : 's'} lost over half their reads at
          primer removal — usually a different primer pair than the one applied:
        </p>
        <table class="w-full text-xs">
          <tbody class="text-slate-300">
            {#each lostSamples.slice(0, 20) as s}
              <tr class="border-t border-slate-800">
                <td class="p-1">{s.id}</td>
                <td class="p-1 text-right">{s.raw.toLocaleString()} &rarr; {s.kept.toLocaleString()}</td>
                <td class="p-1 text-right text-rose-400">{s.pct.toFixed(1)}%</td>
              </tr>
            {/each}
          </tbody>
        </table>
      {/if}
    </div>
  {/if}
</div>
