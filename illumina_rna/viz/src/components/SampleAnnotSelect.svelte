<script>
  /** Picks the sample annotation shared by every view (t-SNE colour, heatmap
   *  strip, composition bars, expression panels). */
  import { sampleAnnotations, sampleAnnotKey } from '../stores/metadata.js';

  let { label = 'Samples by:' } = $props();
  let groups = $derived([...new Set($sampleAnnotations.map(a => a.group))]);
</script>

<label class="flex items-center gap-2 text-xs text-slate-400">
  {label}
  <select class="bg-slate-800 border border-slate-600 text-slate-200 text-xs rounded px-2 py-1 max-w-[12rem]" bind:value={$sampleAnnotKey}>
    {#each groups as g}
      <optgroup label={g}>
        {#each $sampleAnnotations.filter(a => a.group === g) as a}<option value={a.key}>{a.label}</option>{/each}
      </optgroup>
    {/each}
  </select>
</label>
