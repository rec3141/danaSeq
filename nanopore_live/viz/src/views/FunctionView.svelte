<script>
  import { onMount } from 'svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, sampleFunction, readExplorer, loadReadExplorer } from '../stores/data.js';
  import { cartItems, cartActive } from '../stores/cart.js';

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c'];

  // EC top-class names — first digit of an EC number maps to a coarse class.
  const EC_CLASS_NAMES = {
    '1': 'Oxidoreductases', '2': 'Transferases', '3': 'Hydrolases',
    '4': 'Lyases', '5': 'Isomerases', '6': 'Ligases', '7': 'Translocases',
  };

  // Lazy-load read_explorer so the "Dominant product" column can populate.
  // The store is shared so /reads or /map visits already cached it.
  onMount(() => { loadReadExplorer(); });

  let activeSamples = $derived.by(() => {
    if (!$samples) return [];
    if ($cartActive && $cartItems.size > 0) return $samples.filter(s => $cartItems.has(s.id));
    return $samples;
  });

  // --- Gene type bar (kept) -------------------------------------------------
  let geneBarTraces = $derived.by(() => {
    if (!$sampleFunction || !activeSamples.length) return [];
    const ids = activeSamples.map(s => s.id);
    const types = ['cds_count', 'rrna_count', 'trna_count'];
    const labels = ['CDS', 'rRNA', 'tRNA'];
    return types.map((t, i) => ({
      type: 'bar', name: labels[i],
      x: ids,
      y: ids.map(id => $sampleFunction[id]?.[t] ?? 0),
      marker: { color: PALETTE[i] },
    }));
  });
  let geneBarLayout = $derived({
    barmode: 'stack', height: 400,
    xaxis: { title: { text: 'Sample', font: { color: '#94a3b8' } }, tickangle: -45, tickfont: { size: 9 } },
    yaxis: { title: { text: 'Feature Count', font: { color: '#94a3b8' } } },
    legend: { font: { size: 10 } },
    margin: { b: 80 },
  });

  // --- Hypothetical % bar (kept) -------------------------------------------
  let hypotheticalTraces = $derived.by(() => {
    if (!$sampleFunction || !activeSamples.length) return [];
    const ids = activeSamples.map(s => s.id);
    return [{
      type: 'bar',
      x: ids,
      y: ids.map(id => {
        const f = $sampleFunction[id];
        return f?.hypothetical_pct != null ? f.hypothetical_pct : 0;
      }),
      marker: { color: '#fb923c' },
    }];
  });
  let hypotheticalLayout = $derived({
    height: 300,
    xaxis: { title: { text: 'Sample', font: { color: '#94a3b8' } }, tickangle: -45, tickfont: { size: 9 } },
    yaxis: { title: { text: '% Hypothetical', font: { color: '#94a3b8' } }, range: [0, 100] },
    margin: { b: 80 },
    showlegend: false,
  });

  // --- EC abundance heatmap (NEW) ------------------------------------------
  // Collapse depth: 1 = top class (e.g. "2"), 2 = subclass ("2.5"),
  // 3 = sub-subclass ("2.5.1"), 4 = full EC ("2.5.1.61").
  let ecDepth = $state(3);
  let ecTopN = $state(30);
  let ecValueMode = $state('pct'); // 'pct' = %-of-sample-annotations, 'log' = log10 raw count

  /** Collapse "1.2.3.4" → first `depth` segments joined by '.'. */
  function collapseEc(ec, depth) {
    if (!ec) return null;
    const parts = String(ec).split('.');
    if (parts.length === 0) return null;
    return parts.slice(0, depth).join('.');
  }

  /** Pretty label: append top-class name when depth=1; else just the code. */
  function ecLabel(ec, depth) {
    if (depth === 1 && EC_CLASS_NAMES[ec]) return `${ec}. ${EC_CLASS_NAMES[ec]}`;
    return ec;
  }

  // Per-sample collapsed EC counts: { sid: { collapsedEc: count } }.
  let ecBySample = $derived.by(() => {
    if (!$sampleFunction) return {};
    const out = {};
    for (const s of activeSamples) {
      const raw = $sampleFunction[s.id]?.ec_counts;
      if (!raw) continue;
      const collapsed = {};
      for (const [ec, count] of Object.entries(raw)) {
        const c = collapseEc(ec, ecDepth);
        if (!c) continue;
        collapsed[c] = (collapsed[c] || 0) + count;
      }
      out[s.id] = collapsed;
    }
    return out;
  });

  // Top-N collapsed EC ids by total count across active samples.
  let topEcIds = $derived.by(() => {
    const totals = {};
    for (const sid in ecBySample) {
      for (const [ec, n] of Object.entries(ecBySample[sid])) {
        totals[ec] = (totals[ec] || 0) + n;
      }
    }
    return Object.entries(totals)
      .sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0]))
      .slice(0, ecTopN)
      .map(([ec]) => ec);
  });

  let ecHeatmap = $derived.by(() => {
    if (!topEcIds.length || !activeSamples.length) return null;
    const x = activeSamples.map(s => s.id);
    const y = topEcIds;
    const z = y.map(ec => x.map(sid => {
      const cnt = ecBySample[sid]?.[ec] ?? 0;
      if (ecValueMode === 'log') return cnt > 0 ? Math.log10(cnt + 1) : null;
      // pct mode: count / n_genes * 100
      const ng = $sampleFunction[sid]?.n_genes ?? 0;
      if (!cnt) return null;
      return ng > 0 ? (cnt / ng) * 100 : null;
    }));
    return {
      traces: [{
        type: 'heatmap', x, y, z,
        colorscale: 'Viridis',
        colorbar: {
          title: { text: ecValueMode === 'log' ? 'log₁₀(count+1)' : '% of annotations', font: { size: 10 } },
          tickfont: { size: 10 },
        },
        hovertemplate: '%{y}<br>%{x}<br>' + (ecValueMode === 'log' ? 'log₁₀: %{z:.2f}' : '%{z:.2f}%') + '<extra></extra>',
      }],
      layout: {
        height: Math.max(280, 18 * y.length + 120),
        xaxis: { tickangle: -45, tickfont: { size: 9 }, automargin: true },
        yaxis: {
          tickfont: { size: 10 },
          // When depth=1, show class names alongside the digit.
          ticktext: y.map(ec => ecLabel(ec, ecDepth)),
          tickvals: y,
          automargin: true,
        },
        margin: { b: 80, l: 80, t: 8 },
      },
    };
  });

  // --- Per-sample dominant EC + dominant product ---------------------------
  /** For one sample: the EC subclass (depth=3 by default) with the highest count. */
  function dominantEc(sid) {
    const ec = $sampleFunction?.[sid]?.ec_counts;
    if (!ec) return null;
    const collapsed = {};
    for (const [k, c] of Object.entries(ec)) {
      const ck = collapseEc(k, 3);
      if (ck) collapsed[ck] = (collapsed[ck] || 0) + c;
    }
    let best = null, bestN = 0;
    for (const [k, n] of Object.entries(collapsed)) {
      if (n > bestN) { best = k; bestN = n; }
    }
    return best;
  }

  /** EC top class digit → name for one sample. */
  function dominantEcClass(sid) {
    const ec = $sampleFunction?.[sid]?.ec_counts;
    if (!ec) return null;
    const totals = {};
    for (const [k, c] of Object.entries(ec)) {
      const cls = String(k).split('.')[0];
      if (cls) totals[cls] = (totals[cls] || 0) + c;
    }
    let best = null, bestN = 0;
    for (const [k, n] of Object.entries(totals)) {
      if (n > bestN) { best = k; bestN = n; }
    }
    return best ? (EC_CLASS_NAMES[best] ?? best) : null;
  }

  /** Build per-sample top-product index from $readExplorer (lazy). */
  let dominantProductBySample = $derived.by(() => {
    if (!$readExplorer?.reads) return null;
    const counts = {};
    for (const r of $readExplorer.reads) {
      const sid = r.sample;
      const products = r.products;
      if (!sid || !products) continue;
      // products is "Foo; Bar; Baz" — count each
      const seen = (counts[sid] ||= {});
      for (const p of String(products).split('; ')) {
        const t = p.trim();
        if (t) seen[t] = (seen[t] || 0) + 1;
      }
    }
    const out = {};
    for (const [sid, byProd] of Object.entries(counts)) {
      let best = null, bestN = 0;
      for (const [p, n] of Object.entries(byProd)) {
        if (n > bestN) { best = p; bestN = n; }
      }
      if (best) out[sid] = { product: best, count: bestN };
    }
    return out;
  });

  // --- Functional summary table (augmented) --------------------------------
  let funcTableRows = $derived.by(() => {
    if (!$sampleFunction) return [];
    return activeSamples.map(s => {
      const f = $sampleFunction[s.id] || {};
      const dp = dominantProductBySample?.[s.id];
      return {
        id: s.id,
        n_genes: f.n_genes ?? 0,
        cds_count: f.cds_count ?? 0,
        rrna_count: f.rrna_count ?? 0,
        trna_count: f.trna_count ?? 0,
        hypothetical_pct: f.hypothetical_pct ?? 0,
        n_ec: f.ec_counts ? Object.keys(f.ec_counts).length : 0,
        dominant_ec_class: dominantEcClass(s.id) ?? '—',
        dominant_ec: dominantEc(s.id) ?? '—',
        dominant_product: dp ? dp.product : ($readExplorer ? '—' : '…'),
      };
    });
  });

  const funcColumns = [
    { key: 'id', label: 'Sample' },
    { key: 'n_genes', label: 'Total Genes', render: v => v.toLocaleString() },
    { key: 'cds_count', label: 'CDS', render: v => v.toLocaleString() },
    { key: 'rrna_count', label: 'rRNA' },
    { key: 'trna_count', label: 'tRNA' },
    { key: 'hypothetical_pct', label: '% Hyp', render: v => v.toFixed(1) },
    { key: 'n_ec', label: 'EC Numbers' },
    { key: 'dominant_ec_class', label: 'Dominant EC class' },
    { key: 'dominant_ec', label: 'Dominant EC' },
    { key: 'dominant_product', label: 'Dominant product' },
  ];

  // --- Top EC numbers table (kept) -----------------------------------------
  let topEc = $derived.by(() => {
    if (!$sampleFunction) return [];
    const counts = {};
    for (const s of activeSamples) {
      const ec = $sampleFunction[s.id]?.ec_counts;
      if (!ec) continue;
      for (const [ec_num, count] of Object.entries(ec)) {
        counts[ec_num] = (counts[ec_num] || 0) + count;
      }
    }
    return Object.entries(counts)
      .map(([ec, count]) => ({ ec, count }))
      .sort((a, b) => b.count - a.count)
      .slice(0, 30);
  });

  const ecColumns = [
    { key: 'ec', label: 'EC Number' },
    { key: 'count', label: 'Occurrences', render: v => v.toLocaleString() },
  ];
</script>

<div class="space-y-6">
  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6">
    <!-- Gene type breakdown -->
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">Gene Features per Sample</h3>
      {#if geneBarTraces.length}
        <PlotlyChart traces={geneBarTraces} layout={geneBarLayout} exportName="danaseq_gene_features_per_sample" />
      {:else}
        <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">No functional data</div>
      {/if}
    </div>

    <!-- Hypothetical protein ratio -->
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">Hypothetical Protein %</h3>
      {#if hypotheticalTraces.length}
        <PlotlyChart traces={hypotheticalTraces} layout={hypotheticalLayout} exportName="danaseq_hypothetical_protein_pct" />
      {:else}
        <div class="h-[300px] flex items-center justify-center text-slate-500 text-sm">No data</div>
      {/if}
    </div>
  </div>

  <!-- EC abundance heatmap (NEW) -->
  <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
    <div class="flex flex-wrap items-center gap-3 mb-3">
      <h3 class="text-sm font-semibold text-slate-300 mr-2">EC abundance heatmap</h3>
      <label class="flex items-center gap-2 text-slate-400 text-xs">
        <span>Depth</span>
        <select bind:value={ecDepth}
          class="px-2 py-1 rounded-md border border-cyan-400 bg-slate-900 text-cyan-400 hover:bg-cyan-400/10 transition-colors focus:outline-none focus:ring-1 focus:ring-cyan-400 text-xs">
          <option value={1}>Class (1)</option>
          <option value={2}>Subclass (1.2)</option>
          <option value={3}>Sub-subclass (1.2.3)</option>
          <option value={4}>Full EC (1.2.3.4)</option>
        </select>
      </label>
      <label class="flex items-center gap-2 text-slate-400 text-xs">
        <span>Top N</span>
        <input type="range" min="5" max="100" step="5" bind:value={ecTopN} class="accent-cyan-400" />
        <span class="font-mono text-slate-500 w-8">{ecTopN}</span>
      </label>
      <label class="flex items-center gap-2 text-slate-400 text-xs">
        <span>Value</span>
        <select bind:value={ecValueMode}
          class="px-2 py-1 rounded-md border border-cyan-400 bg-slate-900 text-cyan-400 hover:bg-cyan-400/10 transition-colors focus:outline-none focus:ring-1 focus:ring-cyan-400 text-xs">
          <option value="pct">% of annotations</option>
          <option value="log">log₁₀(count+1)</option>
        </select>
      </label>
      <span class="ml-auto text-[11px] text-slate-500">
        {topEcIds.length}/{ecTopN} EC × {activeSamples.length} samples
      </span>
    </div>
    {#if ecHeatmap}
      <PlotlyChart traces={ecHeatmap.traces} layout={ecHeatmap.layout} exportName="danaseq_ec_heatmap" />
    {:else}
      <div class="h-[300px] flex items-center justify-center text-slate-500 text-sm">No EC data</div>
    {/if}
  </div>

  <!-- Functional summary table -->
  <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
    <h3 class="text-sm font-semibold text-slate-300 mb-2">Functional Summary per sample</h3>
    {#if !$readExplorer}
      <p class="text-[11px] text-slate-500 italic mb-2">
        Loading per-read product annotations… (the Dominant-product column will fill in when ready)
      </p>
    {/if}
    <DataTable columns={funcColumns} rows={funcTableRows} maxHeight="360px" exportFilename="functional_summary" />
  </div>

  <!-- Top EC numbers -->
  {#if topEc.length}
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">Top EC Numbers (uncollapsed)</h3>
      <DataTable columns={ecColumns} rows={topEc} maxHeight="300px" idKey="ec" exportFilename="top_ec_numbers" />
    </div>
  {/if}
</div>
