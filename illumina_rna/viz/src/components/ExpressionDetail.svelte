<script>
  /**
   * Side panel for one MAG / function / gene: key facts as a definition list
   * (same layout as the MAG dashboard's bin panel) and its expression across
   * samples.
   *
   * facts: [{ label, value, cls? }]   cls colours the value (e.g. 'text-emerald-400')
   */
  import PlotlyChart from './charts/PlotlyChart.svelte';
  import { fmtTpm, mean } from '../stores/data.js';

  let { title = '', subtitle = '', facts = [], tpm = [], samples = [], onclose = null, links = [], horizontal = false, unit = 'TPM' } = $props();

  let detected = $derived(tpm.filter(v => v > 0).length);
  let profile = $derived([{
    type: 'bar', x: samples, y: tpm,
    marker: { color: '#22d3ee', opacity: 0.8 },
    hovertemplate: `%{x}<br>${unit} %{y:.2f}<extra></extra>`,
  }]);
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
    <PlotlyChart minHeight="260px" data={profile} layout={{
      title: { text: `${unit} per sample`, font: { size: 12 } },
      height: 260,
      margin: { t: 30, r: 10, b: 60, l: 50 },
      xaxis: { tickangle: -60, tickfont: { size: 7 }, showticklabels: samples.length <= 60 },
      yaxis: { title: unit },
      bargap: 0.15,
    }} />
  </div>
  </div>
</div>
