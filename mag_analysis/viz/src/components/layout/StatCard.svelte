<script>
  let { label = '', value = '', sub = '', color = 'cyan' } = $props();

  const colorMap = {
    cyan: 'text-cyan-400',
    emerald: 'text-emerald-400',
    amber: 'text-amber-400',
    rose: 'text-rose-400',
    purple: 'text-purple-400',
    slate: 'text-slate-300',
  };
</script>

<div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
  <div class="text-slate-400 text-xs font-medium uppercase tracking-wider mb-1">{label}</div>
  <div class="{colorMap[color] || colorMap.cyan} text-2xl font-semibold font-mono">{value}</div>
  {#if sub}
    <div class="text-slate-500 text-xs mt-1">{sub}</div>
  {/if}
</div>
