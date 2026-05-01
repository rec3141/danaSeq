<script>
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import { samples, sampleTaxonomy, metadata } from '../stores/data.js';
  import { cartItems, cartActive } from '../stores/cart.js';

  const RANK_LEVELS = [
    { code: 'P', label: 'Phylum' },
    { code: 'C', label: 'Class' },
    { code: 'O', label: 'Order' },
    { code: 'F', label: 'Family' },
    { code: 'G', label: 'Genus' },
    { code: 'S', label: 'Species' },
  ];

  // Start at Genus (index 4) — most interpretable granularity for site comparison.
  let rankIdx = $state(4);
  let topN = $state(30);
  let valueMode = $state('pct');  // 'pct' | 'log'

  let currentRank = $derived(RANK_LEVELS[rankIdx]);

  function cycleRank() { rankIdx = (rankIdx + 1) % RANK_LEVELS.length; }
  function cycleValue() { valueMode = valueMode === 'pct' ? 'log' : 'pct'; }

  let activeSamples = $derived.by(() => {
    if (!$samples) return [];
    if ($cartActive && $cartItems.size > 0) return $samples.filter(s => $cartItems.has(s.id));
    return $samples;
  });

  // Children map (inverted parents) for fast descendant walk at any rank.
  let childrenMap = $derived.by(() => {
    const parents = $sampleTaxonomy?.parents;
    if (!parents) return {};
    const out = {};
    for (const [child, parent] of Object.entries(parents)) {
      (out[parent] ||= []).push(child);
    }
    return out;
  });

  // Taxa at the current rank keyed by name.
  let taxaAtRank = $derived.by(() => {
    const ranks = $sampleTaxonomy?.ranks;
    if (!ranks) return [];
    const out = [];
    for (const [name, r] of Object.entries(ranks)) if (r === currentRank.code) out.push(name);
    return out;
  });

  // Descendant set per taxon at the current rank, used to roll up per-sample
  // direct counts into this rank's buckets.
  let descendantSets = $derived.by(() => {
    const kids = childrenMap;
    const out = {};
    function walk(taxon, acc) {
      acc.push(taxon);
      const cs = kids[taxon];
      if (cs) for (const c of cs) walk(c, acc);
    }
    for (const t of taxaAtRank) {
      const acc = [];
      walk(t, acc);
      out[t] = acc;
    }
    return out;
  });

  // Per-sample rollup: { sampleId: { taxon_at_rank: aggregated_count } }.
  let perSampleRankCounts = $derived.by(() => {
    const out = {};
    const samplesData = $sampleTaxonomy?.samples || {};
    for (const sid in samplesData) {
      const direct = samplesData[sid]?.direct || {};
      if (currentRank.code === 'P') { out[sid] = { ...(samplesData[sid].phylum || {}) }; continue; }
      if (currentRank.code === 'C') { out[sid] = { ...(samplesData[sid].class || {}) }; continue; }
      const row = {};
      for (const [taxon, descendants] of Object.entries(descendantSets)) {
        let n = 0;
        for (const d of descendants) if (direct[d]) n += direct[d];
        if (n > 0) row[taxon] = n;
      }
      out[sid] = row;
    }
    return out;
  });

  // Top-N taxa across active samples, ranked by summed counts.
  let topTaxa = $derived.by(() => {
    const totals = {};
    for (const s of activeSamples) {
      const row = perSampleRankCounts[s.id] || {};
      for (const [taxon, count] of Object.entries(row)) totals[taxon] = (totals[taxon] || 0) + count;
    }
    return Object.entries(totals)
      .sort((a, b) => b[1] - a[1])
      .slice(0, topN)
      .map(([name, count]) => ({ name, count }));
  });

  // Per-sample denominator: sum of classified reads at current rank (so the
  // heatmap cells are "% of reads classified into this rank" — a within-sample
  // fraction that's directly comparable across sites without depth bias).
  function sampleDenom(sid) {
    const row = perSampleRankCounts[sid] || {};
    let s = 0;
    for (const v of Object.values(row)) s += v;
    return s || 1;
  }

  // ---- Column ordering ---------------------------------------------------
  //
  // The order button cycles through: default (sample id) → each metadata
  // column → Ward clustering. Ward runs on squared Euclidean distances over
  // 4th-root-transformed relative abundances — a standard downweighting of
  // dominant taxa used in community ecology (Anderson & Willis 2003; Chao
  // et al., equivalent to the γ=¼ Tukey ladder). Leaves come out in the
  // dendrogram's traversal order so columns that cluster together sit
  // side-by-side in the heatmap.

  // Detect metadata columns present across active samples.
  let metaColumns = $derived.by(() => {
    if (!$metadata) return [];
    const keys = new Set();
    for (const s of activeSamples) {
      const m = $metadata[s.id];
      if (!m) continue;
      for (const k of Object.keys(m)) {
        if (k === 'lat' || k === 'lon') continue;
        keys.add(k);
      }
    }
    return [...keys].sort();
  });

  // Order modes: 'default' (sample id) | 'meta:<col>' | 'cluster'
  let orderModes = $derived([
    { id: 'default', label: 'Default' },
    ...metaColumns.map(c => ({ id: `meta:${c}`, label: c })),
    { id: 'cluster', label: 'Ward' },
  ]);
  let orderIdx = $state(0);
  let currentOrder = $derived(orderModes[orderIdx % Math.max(1, orderModes.length)]);

  function cycleOrder() { orderIdx = (orderIdx + 1) % orderModes.length; }

  // Label modes are independent from Order so sorting by a metadata column
  // doesn't force that column onto the x-axis tick (and vice versa). Default
  // label is FLOWCELL:barcodeNN.
  let labelModes = $derived([
    { id: 'default', label: 'FC:bc' },
    ...metaColumns.map(c => ({ id: `meta:${c}`, label: c })),
  ]);
  let labelIdx = $state(0);
  let currentLabel = $derived(labelModes[labelIdx % Math.max(1, labelModes.length)]);

  function cycleLabel() { labelIdx = (labelIdx + 1) % labelModes.length; }

  // Ward linkage on 4th-root-transformed relative abundances. The 4th-root
  // transform (x^¼) compresses the contribution of dominant taxa so rare
  // taxa inform the geometry more strongly — standard practice for
  // clustering community composition. Ward's Lance-Williams update
  // operates on squared Euclidean distances between cluster centroids
  // (in sample space); we carry those through the merges directly.
  function wardOrder(sampleIds) {
    const n = sampleIds.length;
    if (n <= 1) return sampleIds.slice();

    // Build per-sample vectors in a shared taxon order. Count-of-zero taxa
    // still contribute (x^¼ = 0), just without cost — no need to sparsify.
    const taxaSet = new Set();
    for (const sid of sampleIds) {
      for (const t of Object.keys(perSampleRankCounts[sid] || {})) taxaSet.add(t);
    }
    const taxa = [...taxaSet];
    const vecs = sampleIds.map(sid => {
      const counts = perSampleRankCounts[sid] || {};
      const denom = sampleDenom(sid);
      return taxa.map(t => Math.pow((counts[t] || 0) / denom, 0.25));
    });

    const sqEuc = (a, b) => {
      let s = 0;
      for (let k = 0; k < a.length; k++) { const d = a[k] - b[k]; s += d * d; }
      return s;
    };

    let clusters = sampleIds.map(sid => ({ size: 1, leaves: [sid] }));
    const D = new Map();
    const key = (i, j) => i < j ? `${i},${j}` : `${j},${i}`;
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) D.set(key(i, j), sqEuc(vecs[i], vecs[j]));
    }
    const live = new Set([...Array(n).keys()]);
    while (live.size > 1) {
      let minD = Infinity, pi = -1, pj = -1;
      const arr = [...live];
      for (let a = 0; a < arr.length; a++) {
        for (let b = a + 1; b < arr.length; b++) {
          const d = D.get(key(arr[a], arr[b])) ?? Infinity;
          if (d < minD) { minD = d; pi = arr[a]; pj = arr[b]; }
        }
      }
      const ni = clusters[pi].size, nj = clusters[pj].size;
      const dij = D.get(key(pi, pj)) ?? 0;
      // Lance-Williams coefficients for Ward on squared Euclidean:
      //   d(uv, k) = ((n_u + n_k) d(u,k) + (n_v + n_k) d(v,k) - n_k d(u,v))
      //             / (n_u + n_v + n_k)
      for (const k of live) {
        if (k === pi || k === pj) continue;
        const nk = clusters[k].size;
        const dik = D.get(key(pi, k)) ?? 0;
        const djk = D.get(key(pj, k)) ?? 0;
        const dnew = ((ni + nk) * dik + (nj + nk) * djk - nk * dij) / (ni + nj + nk);
        D.set(key(pi, k), dnew);
      }
      clusters[pi] = { size: ni + nj, leaves: [...clusters[pi].leaves, ...clusters[pj].leaves] };
      live.delete(pj);
    }
    return clusters[[...live][0]].leaves;
  }

  let orderedSamples = $derived.by(() => {
    if (!activeSamples.length) return [];
    const meta = $metadata || {};
    const mode = currentOrder.id;
    if (mode === 'default') {
      return [...activeSamples].sort((a, b) => (a.id || '').localeCompare(b.id || ''));
    }
    if (mode === 'cluster') {
      const ids = activeSamples.map(s => s.id);
      const order = wardOrder(ids);
      const byId = new Map(activeSamples.map(s => [s.id, s]));
      return order.map(id => byId.get(id)).filter(Boolean);
    }
    // meta:<col> — alphanumeric sort with stable (id) tiebreak; numeric
    // values compared as numbers so temperature/depth sort naturally.
    const col = mode.slice(5);
    const getVal = (s) => {
      const v = meta[s.id]?.[col];
      if (v == null) return '';
      return typeof v === 'number' ? v : String(v).replace(/^"|"$/g, '');
    };
    return [...activeSamples].sort((a, b) => {
      const va = getVal(a), vb = getVal(b);
      if (typeof va === 'number' && typeof vb === 'number') {
        if (va !== vb) return va - vb;
      } else {
        const sa = String(va), sb = String(vb);
        const c = sa.localeCompare(sb);
        if (c !== 0) return c;
      }
      return (a.id || '').localeCompare(b.id || '');
    });
  });

  // ---- Heatmap traces + layout -------------------------------------------

  // Single-line x-tick label controlled by `currentLabel` — never stacks two
  // strings so Plotly doesn't render overlapping primary/subline pairs. Falls
  // back to FLOWCELL:barcodeNN when a chosen metadata column is missing for
  // a given sample so no tick is ever blank.
  function xDisplayLabel(s) {
    const fc = s.flowcell && s.barcode ? `${s.flowcell}:${s.barcode}` : s.id;
    const mode = currentLabel.id;
    if (mode.startsWith('meta:')) {
      const col = mode.slice(5);
      const v = $metadata?.[s.id]?.[col];
      if (v != null && v !== '') return String(v).replace(/^"|"$/g, '');
    }
    return fc;
  }

  let heatmapTrace = $derived.by(() => {
    if (!topTaxa.length || !orderedSamples.length) return null;
    // Numeric indices for x — avoids Plotly's category-axis quirk where
    // both the auto category label AND our tickmode-array ticktext render
    // as overlapping (white + grey) labels. Ticks get their display name
    // via tickvals/ticktext only.
    const xVals = orderedSamples.map((_, i) => i);
    const xDisplay = orderedSamples.map(xDisplayLabel);
    const xHoverText = orderedSamples.map(s => xDisplayLabel(s).replace(/<[^>]*>/g, ''));
    const yLabels = topTaxa.map(t => t.name);
    const z = topTaxa.map(taxon =>
      orderedSamples.map(s => {
        const count = (perSampleRankCounts[s.id] || {})[taxon.name] ?? 0;
        if (valueMode === 'log') return count > 0 ? Math.log10(count) : null;
        return (count / sampleDenom(s.id)) * 100;
      })
    );
    // customdata: one value per column, broadcast per row by Plotly.
    const customdata = topTaxa.map(() => xHoverText);
    return {
      trace: {
        type: 'heatmap',
        x: xVals,
        y: yLabels,
        z,
        xgap: 0,
        ygap: 0,
        colorscale: 'Viridis',
        hoverongaps: false,
        colorbar: {
          title: { text: valueMode === 'log' ? 'log10(reads)' : '% of classified', font: { color: '#94a3b8', size: 11 } },
          tickfont: { color: '#94a3b8', size: 10 },
          thickness: 12,
        },
        customdata,
        hovertemplate: (valueMode === 'log'
          ? '<b>%{y}</b><br>%{customdata}<br>log10(reads): %{z:.2f}<extra></extra>'
          : '<b>%{y}</b><br>%{customdata}<br>%{z:.2f}%%<extra></extra>'),
      },
      xVals, xDisplay, yLabels,
    };
  });

  let heatmapLayout = $derived.by(() => {
    if (!heatmapTrace) return {};
    const { xVals, xDisplay, yLabels } = heatmapTrace;
    return {
      height: Math.max(400, 40 + yLabels.length * 18),
      // Bigger top margin now that sample labels live there, smaller bottom
      // since only the colorbar tail sits below.
      margin: { t: 110, r: 20, b: 30, l: 220 },
      xaxis: {
        // Linear numeric axis so Plotly doesn't also emit category labels
        // alongside our custom ticktext. Ticks are placed explicitly per
        // sample column via tickvals.
        side: 'top',
        tickmode: 'array',
        tickvals: xVals,
        ticktext: xDisplay,
        tickangle: -30,
        tickfont: { size: 10, color: '#cbd5e1' },
        automargin: true,
        showgrid: false,
        zeroline: false,
        range: [-0.5, xVals.length - 0.5],
      },
      yaxis: {
        type: 'category',
        tickmode: 'array',
        tickvals: yLabels,
        ticktext: yLabels,
        tickfont: { size: 10, color: '#cbd5e1' },
        automargin: true,
        autorange: 'reversed',   // most abundant at top
      },
    };
  });
</script>

<div class="space-y-4">
  <!-- Controls -->
  <div class="flex items-center gap-3 flex-wrap text-xs">
    <button
      class="px-3 py-1 rounded-md border border-cyan-400 bg-cyan-400/10 text-cyan-400 hover:bg-cyan-400/20 transition-colors"
      style="min-width: 5.5rem"
      onclick={cycleRank}
      title={`Click to cycle rank: ${RANK_LEVELS.map(r => r.label).join(' → ')}`}
    >
      {currentRank.label} &#x25BE;
    </button>

    <button
      class="px-3 py-1 rounded-md border border-cyan-400 bg-cyan-400/10 text-cyan-400 hover:bg-cyan-400/20 transition-colors"
      style="min-width: 6rem"
      onclick={cycleValue}
      title="Click to cycle: % of classified ↔ log10(reads)"
    >
      {valueMode === 'pct' ? '% classified' : 'log10(reads)'} &#x25BE;
    </button>

    <label class="flex items-center gap-2 text-slate-400">
      <span>Order:</span>
      <select
        bind:value={orderIdx}
        class="px-2 py-1 rounded-md border border-cyan-400 bg-slate-900 text-cyan-400 hover:bg-cyan-400/10 transition-colors focus:outline-none focus:ring-1 focus:ring-cyan-400"
      >
        {#each orderModes as m, i}
          <option value={i}>{m.label}</option>
        {/each}
      </select>
    </label>

    <label class="flex items-center gap-2 text-slate-400">
      <span>Label:</span>
      <select
        bind:value={labelIdx}
        class="px-2 py-1 rounded-md border border-cyan-400 bg-slate-900 text-cyan-400 hover:bg-cyan-400/10 transition-colors focus:outline-none focus:ring-1 focus:ring-cyan-400"
      >
        {#each labelModes as m, i}
          <option value={i}>{m.label}</option>
        {/each}
      </select>
    </label>

    <div class="flex items-center gap-2 text-slate-400">
      <span>Top N</span>
      <input type="range" min="5" max="100" step="1" bind:value={topN}
        class="w-28 accent-cyan-400" />
      <span class="text-slate-500 w-10 font-mono tabular-nums">{topN}</span>
    </div>

    <span class="text-slate-500">
      {activeSamples.length} samples
      {#if $cartActive && $cartItems.size > 0}<span class="text-cyan-400">(cart filtered)</span>{/if}
    </span>
  </div>

  <!-- Heatmap -->
  <div class="rounded-lg border border-slate-700 bg-slate-900/40 p-2">
    {#if heatmapTrace}
      <PlotlyChart
        traces={[heatmapTrace.trace]}
        layout={heatmapLayout}
        exportName={`danaseq_heatmap_${currentRank.label.toLowerCase()}_${currentOrder.label.toLowerCase().replace(/[^a-z0-9]+/g, '-')}_by-${currentLabel.label.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`}
      />
    {:else}
      <div class="p-8 text-center text-slate-500 text-sm">
        {#if !$sampleTaxonomy}
          Loading taxonomy data…
        {:else if !activeSamples.length}
          No samples in the current view.
        {:else}
          No classifications at {currentRank.label} rank for the active samples.
        {/if}
      </div>
    {/if}
  </div>

  <p class="text-[11px] text-slate-500 italic">
    Cells show each sample's {currentRank.label}-level composition
    {#if valueMode === 'log'}(log10 raw read counts){:else}(% of classified reads){/if}.
    Columns ordered by <code class="bg-slate-800 px-1 rounded">{currentOrder.label}</code>; taxa ranked by total count across the current sample set.
  </p>
</div>
