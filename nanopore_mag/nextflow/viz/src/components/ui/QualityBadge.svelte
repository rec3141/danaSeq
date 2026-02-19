<script>
  let { quality = 'LQ' } = $props();

  const styles = {
    HQ: 'bg-emerald-400/20 text-emerald-400 border-emerald-400/30',
    MQ: 'bg-amber-400/20 text-amber-400 border-amber-400/30',
    LQ: 'bg-slate-400/20 text-slate-400 border-slate-400/30',
  };
</script>

<span class="inline-block px-2 py-0.5 text-xs font-semibold rounded border {styles[quality] || styles.LQ}">
  {quality}
</span>
