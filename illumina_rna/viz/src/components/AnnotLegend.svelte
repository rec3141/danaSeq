<script>
  /** HTML legend for the shared sample annotation: swatches for categories,
   *  a gradient for numbers and dates. */
  import { sampleAnnot, MISSING_COLOR } from '../stores/metadata.js';

  let { max = 14 } = $props();
  let a = $derived($sampleAnnot);
  let grad = $derived(a?.legend?.colorAt
    ? `linear-gradient(to right, ${[0, 0.25, 0.5, 0.75, 1].map(t => a.legend.colorAt(t)).join(', ')})` : '');
</script>

{#if a?.legend}
  <div class="flex items-center gap-x-3 gap-y-1 flex-wrap text-[11px] text-slate-400">
    <span class="text-slate-500">{a.label}:</span>
    {#if a.legend.kind === 'category'}
      {#each a.legend.items.slice(0, max) as it}
        <span class="flex items-center gap-1"><span class="inline-block w-2.5 h-2.5 rounded-sm" style="background:{it.color}"></span>{it.label}</span>
      {/each}
      {#if a.legend.items.length > max}<span class="text-slate-500">+{a.legend.items.length - max} more</span>{/if}
    {:else if grad}
      <span class="font-mono">{a.legend.minLabel}</span>
      <span class="inline-block w-24 h-2.5 rounded-sm" style="background:{grad}"></span>
      <span class="font-mono">{a.legend.maxLabel}</span>
    {/if}
    {#if a.missing}
      <span class="flex items-center gap-1"><span class="inline-block w-2.5 h-2.5 rounded-sm" style="background:{MISSING_COLOR}"></span>no value ({a.missing})</span>
    {/if}
  </div>
{/if}
