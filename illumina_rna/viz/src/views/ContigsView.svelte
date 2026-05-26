<script>
  import PlotlyChart from '../components/PlotlyChart.svelte';
  import { contigs, taxonomy } from '../stores/data.js';

  let colorBy = $state('mag');      // 'mag' | 'reads' | 'gc' | 'genus'
  let sizeBy  = $state('reads');    // 'reads' | 'length' | 'none'
  let minLen  = $state(1000);
  let magFilter = $state('all');

  let magList = $derived.by(() => {
    const set = new Set(($contigs?.contigs || []).map(c => c.mag));
    return Array.from(set).sort();
  });

  const COLORS = [
    '#22d3ee', '#34d399', '#a78bfa', '#f472b6', '#fbbf24',
    '#fb7185', '#60a5fa', '#facc15', '#84cc16', '#06b6d4',
    '#c084fc', '#f87171', '#2dd4bf', '#a3e635', '#fb923c',
    '#818cf8', '#ec4899',
  ];
  function palette(keys) {
    const m = {};
    keys.forEach((k, i) => { m[k] = COLORS[i % COLORS.length]; });
    return m;
  }
  function genusOf(mag) {
    const lin = $taxonomy?.[mag]?.lineage;
    return lin ? lin.split(/\s/)[0] : null;
  }

  let filtered = $derived(($contigs?.contigs || []).filter(c =>
    c.length >= minLen && (magFilter === 'all' || c.mag === magFilter)
  ));

  function hover(c) {
    const tax = $taxonomy?.[c.mag];
    const taxLine = tax?.lineage ? `<br>↳ ${tax.lineage} (ANI ${tax.ani?.toFixed(1)})` : '';
    return `<b>${c.id}</b><br>${c.mag}${taxLine}` +
           `<br>${c.length.toLocaleString()} bp · GC ${(c.gc * 100).toFixed(1)}%` +
           `<br>Σ reads: ${c.total_reads.toLocaleString()}`;
  }

  function sizeArr(arr) {
    if (sizeBy === 'none') return 6;
    const key = sizeBy === 'reads' ? 'total_reads' : 'length';
    const vals = arr.map(c => Math.log10((c[key] || 0) + 1));
    const max = Math.max(1e-3, ...vals);
    return vals.map(v => 4 + 16 * (v / max));
  }

  let traces = $derived.by(() => {
    const out = [];
    if (colorBy === 'mag') {
      const grouped = {};
      filtered.forEach(c => { (grouped[c.mag] ||= []).push(c); });
      const pal = palette(magList);
      for (const mag of Object.keys(grouped).sort()) {
        const arr = grouped[mag];
        const tax = $taxonomy?.[mag];
        out.push({
          x: arr.map(c => c.x), y: arr.map(c => c.y),
          mode: 'markers', type: 'scattergl',
          name: tax?.lineage ? `${mag} · ${tax.lineage.split(' ').slice(0, 2).join(' ')}` : mag,
          marker: { color: pal[mag], opacity: 0.75, line: { width: 0 }, size: sizeArr(arr) },
          text: arr.map(hover), hovertemplate: '%{text}<extra></extra>',
        });
      }
    } else if (colorBy === 'genus') {
      const byGenus = {};
      filtered.forEach(c => {
        const g = genusOf(c.mag) || '(unknown)';
        (byGenus[g] ||= []).push(c);
      });
      const genera = Object.keys(byGenus).sort();
      const pal = palette(genera);
      for (const g of genera) {
        const arr = byGenus[g];
        out.push({
          x: arr.map(c => c.x), y: arr.map(c => c.y),
          mode: 'markers', type: 'scattergl', name: g,
          marker: { color: pal[g], opacity: 0.75, line: { width: 0 }, size: sizeArr(arr) },
          text: arr.map(hover), hovertemplate: '%{text}<extra></extra>',
        });
      }
    } else {
      const colorVals = filtered.map(c =>
        colorBy === 'reads' ? Math.log10((c.total_reads || 0) + 1) : c.gc
      );
      out.push({
        x: filtered.map(c => c.x), y: filtered.map(c => c.y),
        mode: 'markers', type: 'scattergl', name: 'contigs',
        marker: {
          color: colorVals,
          colorscale: colorBy === 'reads' ? 'Viridis' : 'RdBu',
          showscale: true,
          colorbar: { title: colorBy === 'reads' ? 'log₁₀(Σ+1)' : 'GC' },
          opacity: 0.8, line: { width: 0 }, size: sizeArr(filtered),
        },
        text: filtered.map(hover), hovertemplate: '%{text}<extra></extra>',
      });
    }
    return out;
  });

  let layout = $derived({
    xaxis: { title: 't-SNE 1 (TNF)' },
    yaxis: { title: 't-SNE 2 (TNF)' },
    showlegend: (colorBy === 'mag' || colorBy === 'genus'),
    legend: { font: { size: 10 } },
  });
</script>

<div class="space-y-4">
  <div class="bg-slate-900 border border-slate-800 rounded-lg p-4 flex flex-wrap items-center gap-6">
    <label class="text-sm text-slate-400">
      MAG:
      <select bind:value={magFilter} class="ml-2 bg-slate-800 border border-slate-700 rounded px-2 py-1 text-slate-200 text-sm">
        <option value="all">all MAGs</option>
        {#each magList as m}<option value={m}>{m}</option>{/each}
      </select>
    </label>
    <label class="text-sm text-slate-400">
      Color by:
      <select bind:value={colorBy} class="ml-2 bg-slate-800 border border-slate-700 rounded px-2 py-1 text-slate-200 text-sm">
        <option value="mag">MAG</option>
        <option value="genus">Genus (sendsketch)</option>
        <option value="reads">RNA expression</option>
        <option value="gc">GC content</option>
      </select>
    </label>
    <label class="text-sm text-slate-400">
      Size by:
      <select bind:value={sizeBy} class="ml-2 bg-slate-800 border border-slate-700 rounded px-2 py-1 text-slate-200 text-sm">
        <option value="reads">RNA expression</option>
        <option value="length">Contig length</option>
        <option value="none">Constant</option>
      </select>
    </label>
    <label class="text-sm text-slate-400">
      Min length:
      <input type="number" min="0" step="500" bind:value={minLen}
             class="ml-2 bg-slate-800 border border-slate-700 rounded px-2 py-1 text-slate-200 text-sm w-24"/>
    </label>
    <div class="text-xs text-slate-500">
      {filtered.length} of {($contigs?.contigs || []).length} contigs · scroll to zoom, drag to pan
    </div>
  </div>

  <div class="bg-slate-900 border border-slate-800 rounded-lg p-2">
    {#if !$contigs}
      <div class="text-slate-500 text-sm py-12 text-center">No contig TNF data — run the augmentation step.</div>
    {:else}
      <PlotlyChart data={traces} {layout} height="720px" />
    {/if}
  </div>
</div>
