<script>
  import D3Sunburst from '../components/charts/D3Sunburst.svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { samples, sampleTaxonomy, taxonomySunburst } from '../stores/data.js';
  import { cartItems, cartActive } from '../stores/cart.js';

  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8',
                   '#94a3b8','#d4d4d8','#78716c'];

  const RANK_LEVELS = [
    { code: 'P', label: 'Phylum' },
    { code: 'C', label: 'Class' },
    { code: 'O', label: 'Order' },
    { code: 'F', label: 'Family' },
    { code: 'G', label: 'Genus' },
    { code: 'S', label: 'Species' },
  ];

  let rankIdx = $state(0);
  let showPct = $state(false);

  let currentRank = $derived(RANK_LEVELS[rankIdx]);

  function cycleRank() {
    rankIdx = (rankIdx + 1) % RANK_LEVELS.length;
  }

  // Filtered sample list
  let activeSamples = $derived.by(() => {
    if (!$samples) return [];
    if ($cartActive && $cartItems.size > 0) return $samples.filter(s => $cartItems.has(s.id));
    return $samples;
  });

  let cartIsFiltering = $derived($cartActive && $cartItems.size > 0);

  // Build children map (inverted parent map) for clade count computation
  let childrenMap = $derived.by(() => {
    const parents = $sampleTaxonomy?.parents;
    if (!parents) return {};
    const children = {};
    for (const [child, parent] of Object.entries(parents)) {
      if (!children[parent]) children[parent] = [];
      children[parent].push(child);
    }
    return children;
  });

  // Set of taxa at each rank
  let taxaByRank = $derived.by(() => {
    const ranks = $sampleTaxonomy?.ranks;
    if (!ranks) return {};
    const byRank = {};
    for (const [name, rank] of Object.entries(ranks)) {
      if (!byRank[rank]) byRank[rank] = new Set();
      byRank[rank].add(name);
    }
    return byRank;
  });

  // Precompute descendant sets for current rank (tree walk happens once per rank change)
  let descendantSets = $derived.by(() => {
    const taxaAtRank = taxaByRank[currentRank.code];
    if (!taxaAtRank) return {};
    // For each taxon at the target rank, collect all descendant names (including self)
    const sets = {};
    function collectDescendants(taxon, acc) {
      acc.push(taxon);
      const kids = childrenMap[taxon];
      if (kids) for (const kid of kids) collectDescendants(kid, acc);
    }
    for (const taxon of taxaAtRank) {
      const acc = [];
      collectDescendants(taxon, acc);
      sets[taxon] = acc;
    }
    return sets;
  });

  // Get taxonomy at current rank for a sample
  function getRankTax(sampleId) {
    const tax = $sampleTaxonomy?.samples?.[sampleId];
    if (!tax) return {};
    // Use precomputed clade counts for phylum/class (faster, more accurate)
    if (currentRank.code === 'P') return tax.phylum || {};
    if (currentRank.code === 'C') return tax.class || {};
    // Sum descendant direct counts using precomputed sets
    const direct = tax.direct || {};
    const result = {};
    for (const [taxon, descendants] of Object.entries(descendantSets)) {
      let count = 0;
      for (const d of descendants) {
        if (direct[d]) count += direct[d];
      }
      if (count > 0) result[taxon] = count;
    }
    return result;
  }

  // Per-sample totals for percentage mode
  function sampleTotal(id) {
    const tax = getRankTax(id);
    return Object.values(tax).reduce((s, v) => s + v, 0) || 1;
  }

  // Build cart-filtered sunburst from per-sample direct counts + parent map
  let sunburstData = $derived.by(() => {
    if (!$sampleTaxonomy) return $taxonomySunburst;
    if (!cartIsFiltering) return $taxonomySunburst;

    const parents = $sampleTaxonomy.parents || {};
    const sampleData = $sampleTaxonomy.samples || {};

    // Aggregate direct counts across carted samples
    const merged = {};
    for (const s of activeSamples) {
      const direct = sampleData[s.id]?.direct;
      if (!direct) continue;
      for (const [taxon, count] of Object.entries(direct)) {
        merged[taxon] = (merged[taxon] || 0) + count;
      }
    }

    if (Object.keys(merged).length === 0) return $taxonomySunburst;

    // Build tree from parent map + merged counts
    const tree = {};

    function ensurePath(taxon) {
      const path = [taxon];
      let cur = taxon;
      while (parents[cur]) {
        path.unshift(parents[cur]);
        cur = parents[cur];
      }
      let node = tree;
      for (const name of path) {
        if (!node[name]) node[name] = { _count: 0, _children: {} };
        if (name === taxon) node[name]._count += (merged[taxon] || 0);
        node = node[name]._children;
      }
    }

    for (const taxon of Object.keys(merged)) {
      ensurePath(taxon);
    }

    function toSunburst(children) {
      const result = [];
      for (const [name, data] of Object.entries(children)) {
        const node = { name };
        const childList = toSunburst(data._children);
        if (childList.length) node.children = childList;
        if (data._count > 0) node.value = data._count;
        if ('children' in node || 'value' in node) result.push(node);
      }
      function subtreeTotal(n) {
        let v = n.value || 0;
        for (const c of n.children || []) v += subtreeTotal(c);
        return v;
      }
      result.sort((a, b) => subtreeTotal(b) - subtreeTotal(a));
      return result;
    }

    const children = toSunburst(tree);
    return { name: 'Community', children };
  });

  // Stacked bar chart: per-sample composition at current rank
  let stackedTraces = $derived.by(() => {
    if (!$sampleTaxonomy?.samples || !activeSamples.length) return [];
    const allTaxa = new Set();
    for (const s of activeSamples) {
      const tax = getRankTax(s.id);
      Object.keys(tax).forEach(k => allTaxa.add(k));
    }
    const taxaList = [...allTaxa].sort();
    const sampleIds = activeSamples.map(s => s.id);

    return taxaList.map((taxon, i) => ({
      type: 'bar',
      name: taxon.length > 20 ? taxon.slice(0, 18) + '..' : taxon,
      x: sampleIds,
      y: sampleIds.map(id => {
        const count = getRankTax(id)[taxon] ?? 0;
        return showPct ? (count / sampleTotal(id)) * 100 : count;
      }),
      marker: { color: PALETTE[i % PALETTE.length] },
    }));
  });

  let stackedLayout = $derived({
    barmode: 'stack',
    height: 400,
    xaxis: { title: { text: 'Sample', font: { color: '#94a3b8' } }, tickangle: -45, tickfont: { size: 9 } },
    yaxis: { title: { text: showPct ? '% of Classified' : 'Read Count', font: { color: '#94a3b8' } },
             range: showPct ? [0, 100] : undefined },
    legend: { font: { size: 9 }, bgcolor: 'rgba(15,23,42,0.7)' },
    margin: { b: 80 },
  });

  // Top taxa table (aggregated across active samples)
  let topTaxa = $derived.by(() => {
    if (!$sampleTaxonomy?.samples) return [];
    const counts = {};
    for (const s of activeSamples) {
      const tax = getRankTax(s.id);
      for (const [taxon, count] of Object.entries(tax)) {
        counts[taxon] = (counts[taxon] || 0) + count;
      }
    }
    const total = Object.values(counts).reduce((s, v) => s + v, 0);
    return Object.entries(counts)
      .map(([taxon, count]) => ({
        taxon,
        count,
        pct: total > 0 ? Math.round((count / total) * 1000) / 10 : 0,
      }))
      .sort((a, b) => b.count - a.count);
  });

  let taxColumns = $derived([
    { key: 'taxon', label: currentRank.label },
    { key: 'count', label: 'Reads', render: v => v?.toLocaleString() ?? '-' },
    { key: 'pct', label: '%', render: v => `${v.toFixed(1)}%` },
  ]);
</script>

<div class="space-y-6">
  <!-- Rank toggle (shared across table + bar) -->
  <div class="flex items-center gap-3">
    <span class="text-xs text-slate-500">Taxonomic Level:</span>
    <button
      class="text-xs px-2.5 py-1 rounded border font-medium transition-colors
        bg-cyan-400/20 text-cyan-400 border-cyan-400/40 hover:bg-cyan-400/30"
      onclick={cycleRank}
    >
      {currentRank.label}
    </button>
    <span class="text-[10px] text-slate-600">
      {RANK_LEVELS.map(r => r.label === currentRank.label ? `[${r.label}]` : r.label).join(' > ')}
    </span>
  </div>

  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6">
    <!-- Sunburst -->
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">
        Taxonomy Overview
        {#if cartIsFiltering}
          <span class="text-xs text-cyan-400 ml-2">(cart filtered)</span>
        {/if}
        <span class="text-xs text-slate-500 font-normal ml-2">click to zoom, click center to zoom out</span>
      </h3>
      {#if sunburstData}
        <D3Sunburst data={sunburstData} colorDepth={3} />
      {:else}
        <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">No taxonomy data</div>
      {/if}
    </div>

    <!-- Top taxa table -->
    <div class="bg-slate-800 rounded-lg border border-slate-700 p-4">
      <h3 class="text-sm font-semibold text-slate-300 mb-2">
        Top {currentRank.label === 'Family' ? 'Families' : currentRank.label === 'Class' ? 'Classes' : currentRank.label === 'Genus' ? 'Genera' : currentRank.label + 's'}
        {#if cartIsFiltering}
          <span class="text-xs text-cyan-400 ml-2">(cart filtered)</span>
        {/if}
      </h3>
      <DataTable columns={taxColumns} rows={topTaxa} maxHeight="450px" idKey="taxon" exportFilename={`top_${currentRank.label.toLowerCase()}`} />
    </div>
  </div>

  <!-- Stacked bar chart -->
  <div class="bg-slate-800 rounded-lg border border-slate-700 p-4 overflow-hidden">
    <div class="flex items-center justify-between mb-2">
      <h3 class="text-sm font-semibold text-slate-300">
        Per-Sample {currentRank.label} Composition
        {#if cartIsFiltering}
          <span class="text-xs text-cyan-400 ml-2">(cart filtered)</span>
        {/if}
      </h3>
      <button
        class="text-xs px-2 py-1 rounded border transition-colors
          {showPct
            ? 'bg-cyan-400/20 text-cyan-400 border-cyan-400/40'
            : 'text-slate-400 border-slate-600 hover:text-slate-200 hover:border-slate-500'}"
        onclick={() => showPct = !showPct}
      >
        {showPct ? '%' : '#'}
      </button>
    </div>
    {#if stackedTraces.length}
      <PlotlyChart traces={stackedTraces} layout={stackedLayout} />
    {:else}
      <div class="h-[400px] flex items-center justify-center text-slate-500 text-sm">No per-sample taxonomy data</div>
    {/if}
  </div>
</div>
