<script>
  /**
   * Side panel for one MAG / function / gene: key facts as a definition list
   * (same layout as the MAG dashboard's bin panel) and its expression across
   * samples.
   *
   * facts: [{ label, value, cls? }]   cls colours the value (e.g. 'text-emerald-400')
   *
   * The profile is coloured by the shared sample annotation and can plot the
   * expression against it (scatter for numbers and dates, boxes per category).
   */
  import PlotlyChart from './charts/PlotlyChart.svelte';
  import { fmtTpm, mean } from '../stores/data.js';
  import { sampleAnnot } from '../stores/metadata.js';

  let { title = '', subtitle = '', facts = [], tpm = [], samples = [], onclose = null, links = [], horizontal = false, unit = 'TPM' } = $props();

  let detected = $derived(tpm.filter(v => v > 0).length);
  let mode = $state('samples');   // 'samples' | 'annot'
  let an = $derived($sampleAnnot);
  let canPlotVs = $derived(!!an);
  let vs = $derived(mode === 'annot' && canPlotVs);

  let profile = $derived.by(() => {
    if (!vs) return [{
      type: 'bar', x: samples, y: tpm,
      marker: { color: an ? samples.map(an.colorOf) : '#22d3ee', opacity: 0.9 },
      hovertext: samples.map(n => (an ? `${an.label}: ${an.display(n)}` : '')),
      hovertemplate: `%{x}<br>${unit} %{y:.2f}<br>%{hovertext}<extra></extra>`,
    }];
    const idx = samples.map((_, j) => j).filter(j => an.value(samples[j]) != null);
    if (an.type === 'category') {
      return an.legend.levels.map(l => {
        const js = idx.filter(j => String(an.value(samples[j])) === l);
        return {
          type: 'box', name: l, y: js.map(j => tpm[j]), text: js.map(j => samples[j]),
          boxpoints: 'all', jitter: 0.4, pointpos: 0, marker: { color: an.legend.colorFor(l), size: 5 },
          line: { color: an.legend.colorFor(l), width: 1.5 }, fillcolor: 'rgba(0,0,0,0)',
          hovertemplate: `%{text}<br>${unit} %{y:.2f}<extra>${l}</extra>`,
        };
      });
    }
    return [{
      type: 'scatter', mode: 'markers',
      x: idx.map(j => an.value(samples[j])), y: idx.map(j => tpm[j]), text: idx.map(j => samples[j]),
      marker: { color: '#22d3ee', size: 7, line: { color: '#0f172a', width: 1 } },
      hovertemplate: `%{text}<br>${an.label}: %{x}<br>${unit} %{y:.2f}<extra></extra>`,
    }];
  });
  // Spearman correlation for a numeric/date metadata column
  let rho = $derived.by(() => {
    if (!vs || an.type === 'category') return null;
    const toNum = v => (an.type === 'date' ? Date.parse(v) : v);
    const pts = samples.map((n, j) => [toNum(an.value(n)), tpm[j]]).filter(([x]) => Number.isFinite(x));
    if (pts.length < 4) return null;
    const rank = a => { const o = a.map((v, i) => [v, i]).sort((p, q) => p[0] - q[0]); const r = new Array(a.length);
      for (let i = 0; i < o.length;) { let k = i; while (k + 1 < o.length && o[k + 1][0] === o[i][0]) k++;
        for (let m = i; m <= k; m++) r[o[m][1]] = (i + k) / 2; i = k + 1; } return r; };
    const rx = rank(pts.map(p => p[0])), ry = rank(pts.map(p => p[1]));
    const mx = mean(rx), my = mean(ry);
    let sxy = 0, sxx = 0, syy = 0;
    rx.forEach((x, i) => { sxy += (x - mx) * (ry[i] - my); sxx += (x - mx) ** 2; syy += (ry[i] - my) ** 2; });
    return sxx && syy ? { r: sxy / Math.sqrt(sxx * syy), n: pts.length } : null;
  });
</script>

<div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
  <div class="flex items-start justify-between mb-1 gap-2">
    <h3 class="text-cyan-400 font-semibold font-mono break-all">{title}</h3>
    {#if onclose}
      <button class="text-slate-500 hover:text-slate-300 text-lg leading-none" onclick={onclose} title="Close">&times;</button>
    {/if}
  </div>
  {#if subtitle}<p class="text-xs text-slate-500 mb-3">{subtitle}</p>{/if}
  <div class={horizontal ? 'grid grid-cols-1 lg:grid-cols-3 gap-6' : ''}>
  <div>
  <dl class="space-y-2 text-sm">
    {#each facts.filter(f => f.value != null && f.value !== '') as f}
      <div class="flex justify-between gap-4">
        <dt class="text-slate-400 shrink-0">{f.label}</dt>
        <dd class="font-mono text-right {f.cls ?? ''}">{f.value}</dd>
      </div>
    {/each}
    <div class="flex justify-between">
      <dt class="text-slate-400">Mean {unit}</dt>
      <dd class="font-mono text-cyan-400">{fmtTpm(mean(tpm))}</dd>
    </div>
    <div class="flex justify-between">
      <dt class="text-slate-400">Detected in</dt>
      <dd class="font-mono">{detected} / {samples.length} samples</dd>
    </div>
  </dl>
  {#if links.length}
    <div class="flex gap-3 mt-3 text-xs">
      {#each links as l}
        <a href={l.href} target="_blank" rel="noopener noreferrer" class="text-slate-500 hover:text-cyan-400 transition-colors">{l.label} &#x2197;</a>
      {/each}
    </div>
  {/if}
  </div>
  <div class={horizontal ? 'lg:col-span-2' : 'border-t border-slate-700 pt-2 mt-3'}>
    {#if canPlotVs}
      <div class="flex justify-end gap-1 text-[11px] mb-1">
        {#each [['samples', 'Per sample'], ['annot', `vs ${an.label}`]] as [k, l]}
          <button class="px-2 py-0.5 rounded border transition-colors
            {mode === k ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
            onclick={() => mode = k}>{l}</button>
        {/each}
      </div>
    {/if}
    <PlotlyChart minHeight="260px" data={profile} layout={{
      title: { text: vs ? `${unit} vs ${an.label}${rho ? ` (Spearman ρ = ${rho.r.toFixed(2)}, n = ${rho.n})` : ''}` : `${unit} per sample${an ? ` · coloured by ${an.label}` : ''}`, font: { size: 12 } },
      height: 260,
      margin: { t: 30, r: 10, b: 60, l: 50 },
      xaxis: vs
        ? { title: an.type === 'category' ? '' : an.label, type: an.type === 'date' ? 'date' : undefined, tickfont: { size: 9 } }
        : { tickangle: -60, tickfont: { size: 7 }, showticklabels: samples.length <= 60 },
      yaxis: { title: unit },
      showlegend: false,
      bargap: 0.15,
    }} />
  </div>
  </div>
</div>
