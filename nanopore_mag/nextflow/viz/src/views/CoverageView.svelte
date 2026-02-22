<script>
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import D3Heatmap from '../components/charts/D3Heatmap.svelte';
  import { onMount } from 'svelte';
  import { coverage, checkm2All, loadCheckm2All } from '../stores/data.js';
  import { selectedMag } from '../stores/selection.js';

  let covData = $derived($coverage);
  let allBins = $derived($checkm2All);
  let selected = $derived($selectedMag);

  onMount(() => { loadCheckm2All(); });

  // Binner toggle state — DAS Tool on by default
  let activeBinners = $state(new Set(['dastool']));

  // Taxonomy source
  let taxSource = $state('kaiju');
  const taxSources = [
    { key: 'kaiju', label: 'Kaiju' },
    { key: 'kraken2', label: 'Kraken2' },
    { key: 'rrna', label: 'rRNA' },
  ];

  function cycleSource() {
    const keys = taxSources.map(s => s.key);
    taxSource = keys[(keys.indexOf(taxSource) + 1) % keys.length];
  }
  function taxSourceLabel() {
    return taxSources.find(s => s.key === taxSource)?.label || taxSource;
  }

  const binnerDefs = [
    { key: 'dastool', label: 'DAS Tool' },
    { key: 'semibin', label: 'SemiBin2' },
    { key: 'metabat', label: 'MetaBAT2' },
    { key: 'maxbin',  label: 'MaxBin2' },
    { key: 'lorbin',  label: 'LorBin' },
    { key: 'comebin', label: 'COMEBin' },
  ];

  function toggleBinner(key) {
    const next = new Set(activeBinners);
    if (next.has(key)) {
      if (next.size > 1) next.delete(key);
    } else {
      next.add(key);
    }
    activeBinners = next;
  }

  // Binner counts for button labels
  let binnerCounts = $derived.by(() => {
    if (!covData) return {};
    const counts = {};
    for (const b of covData.bins) {
      counts[b.binner] = (counts[b.binner] || 0) + 1;
    }
    return counts;
  });

  // Quality filter sliders
  let minCompleteness = $state(0);
  let maxContamination = $state(100);

  // Quality tier helper
  function qualityTier(b) {
    if (b.completeness >= 90 && b.contamination < 5) return 'HQ';
    if (b.completeness >= 50 && b.contamination < 10) return 'MQ';
    return 'LQ';
  }

  // Lookup map: bin name -> checkm2All record
  let binMetaMap = $derived.by(() => {
    if (!allBins) return {};
    const map = {};
    for (const b of allBins) {
      map[b.name] = b;
      if (b.dastool_name) map[b.dastool_name] = b;
    }
    return map;
  });

  // Sort-by cycling
  const sortDefs = [
    { key: 'clustered', label: 'Clustered' },
    { key: 'name', label: 'A-Z' },
    { key: 'completeness', label: 'Completeness' },
    { key: 'contamination', label: 'Contamination' },
    { key: 'quality', label: 'Quality' },
    { key: 'size', label: 'Size' },
    { key: 'taxonomy', label: 'Taxonomy' },
    { key: 'total_depth', label: 'Total Depth' },
  ];
  let sortBy = $state('clustered');

  function cycleSort() {
    const keys = sortDefs.map(d => d.key);
    sortBy = keys[(keys.indexOf(sortBy) + 1) % keys.length];
  }

  // Filter bins by active binners + quality thresholds
  let filteredIndices = $derived.by(() => {
    if (!covData) return [];
    const hasQualityFilter = minCompleteness > 0 || maxContamination < 100;
    return covData.bins
      .map((b, i) => {
        if (!activeBinners.has(b.binner)) return -1;
        const meta = binMetaMap[b.id];
        if (meta) {
          if (meta.completeness < minCompleteness) return -1;
          if (meta.contamination > maxContamination) return -1;
        } else if (hasQualityFilter) {
          return -1;
        }
        return i;
      })
      .filter(i => i >= 0);
  });

  // Sort filtered indices
  let sortedIndices = $derived.by(() => {
    if (!covData || !filteredIndices.length) return filteredIndices;
    if (sortBy === 'clustered') return filteredIndices;

    return [...filteredIndices].sort((a, b) => {
      const idA = covData.bins[a].id;
      const idB = covData.bins[b].id;
      const metaA = binMetaMap[idA];
      const metaB = binMetaMap[idB];
      switch (sortBy) {
        case 'name':
          return idA.localeCompare(idB);
        case 'completeness':
          return (metaB?.completeness ?? 0) - (metaA?.completeness ?? 0);
        case 'contamination':
          return (metaA?.contamination ?? 0) - (metaB?.contamination ?? 0);
        case 'quality': {
          const tierOrder = { HQ: 0, MQ: 1, LQ: 2 };
          const tA = metaA ? qualityTier(metaA) : 'LQ';
          const tB = metaB ? qualityTier(metaB) : 'LQ';
          const d = (tierOrder[tA] ?? 3) - (tierOrder[tB] ?? 3);
          return d !== 0 ? d : (metaB?.completeness ?? 0) - (metaA?.completeness ?? 0);
        }
        case 'size':
          return (metaB?.genome_size ?? 0) - (metaA?.genome_size ?? 0);
        case 'taxonomy': {
          const taxA = metaA?.taxonomy?.[taxSource];
          const taxB = metaB?.taxonomy?.[taxSource];
          const linA = taxA ? ['domain', 'phylum', 'class', 'order', 'family', 'genus'].map(r => taxA[r]).filter(Boolean).join(';') : '';
          const linB = taxB ? ['domain', 'phylum', 'class', 'order', 'family', 'genus'].map(r => taxB[r]).filter(Boolean).join(';') : '';
          return linA.localeCompare(linB);
        }
        case 'total_depth': {
          const totalA = covData.matrix[a].reduce((s, v) => s + v, 0);
          const totalB = covData.matrix[b].reduce((s, v) => s + v, 0);
          return totalB - totalA;
        }
        default: return 0;
      }
    });
  });

  // Find the max depth for color scaling (across filtered bins)
  let maxDepth = $derived.by(() => {
    if (!covData || !sortedIndices.length) return 1;
    let mx = 0;
    for (const i of sortedIndices) {
      for (const v of covData.matrix[i]) if (v > mx) mx = v;
    }
    return mx || 1;
  });

  // Log-scaled Viridis color function for coverage
  import * as d3 from 'd3';

  let covColorScale = $derived.by(() => {
    const logMax = Math.log10(maxDepth + 1);
    const viridis = d3.scaleSequential(d3.interpolateViridis).domain([0, logMax]);
    return (val) => {
      if (val <= 0) return '#0f172a';
      return viridis(Math.log10(val + 1));
    };
  });

  // Heatmap data: sorted+filtered bins with raw depths
  let heatmapData = $derived.by(() => {
    if (!covData || !sortedIndices.length) return null;
    return {
      mag_ids: sortedIndices.map(i => covData.bins[i].id),
      module_ids: covData.sample_names.map(s => s),
      module_names: covData.sample_names,
      matrix: sortedIndices.map(i =>
        covData.matrix[i].map(v => Math.round(v * 1000) / 1000)
      ),
      row_order: sortedIndices.map((_, i) => i),
      col_order: covData.sample_names.map((_, i) => i),
    };
  });

  // Row annotations for coverage heatmap
  let covRowAnnotations = $derived.by(() => {
    if (!heatmapData) return null;
    const ids = heatmapData.mag_ids;

    const qualityColor = (v) => v === 'HQ' ? '#34d399' : v === 'MQ' ? '#fbbf24' : '#64748b';
    const compColor = (v) => {
      const n = parseFloat(v);
      return n >= 90 ? '#34d399' : n >= 50 ? '#fbbf24' : '#64748b';
    };
    const contColor = (v) => {
      const n = parseFloat(v);
      return n < 5 ? '#34d399' : n < 10 ? '#fbbf24' : '#ef4444';
    };

    return [
      {
        label: 'Quality',
        values: ids.map(id => {
          const m = binMetaMap[id];
          return m ? qualityTier(m) : '-';
        }),
        colorFn: qualityColor,
      },
      {
        label: 'Comp%',
        values: ids.map(id => {
          const m = binMetaMap[id];
          return m ? m.completeness.toFixed(0) : '-';
        }),
        colorFn: compColor,
      },
      {
        label: 'Cont%',
        values: ids.map(id => {
          const m = binMetaMap[id];
          return m ? m.contamination.toFixed(1) : '-';
        }),
        colorFn: contColor,
      },
      {
        label: 'Size',
        values: ids.map(id => {
          const m = binMetaMap[id];
          return m ? (m.genome_size / 1e6).toFixed(1) + 'M' : '-';
        }),
      },
      {
        label: `Tax (${taxSourceLabel()})`,
        values: ids.map(id => {
          const m = binMetaMap[id];
          const tax = m?.taxonomy?.[taxSource];
          if (!tax) return '-';
          return ['domain', 'phylum', 'class', 'order', 'family', 'genus']
            .map(r => tax[r]).filter(Boolean).join('; ') || '-';
        }),
        align: 'start',
        width: 800,
      },
    ];
  });

  // Tooltip formatter for coverage
  function covTooltip(val, magId, modId, modName) {
    return `<strong>${modId}</strong><br>Bin: ${magId}<br>Depth: ${val.toFixed(2)}x`;
  }

  // Box plot: depth distribution per sample (filtered bins only)
  let boxData = $derived.by(() => {
    if (!covData || !sortedIndices.length) return [];
    const palette = ['#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa', '#fb923c', '#818cf8', '#f472b6'];
    return covData.sample_names.map((sname, j) => ({
      type: 'box',
      name: sname,
      y: sortedIndices.map(i => covData.matrix[i][j]).filter(v => v > 0),
      marker: { color: palette[j % palette.length] },
      boxpoints: 'outliers',
    }));
  });

  // Bar chart: total depth per bin (sorted, filtered)
  let barData = $derived.by(() => {
    if (!covData || !sortedIndices.length) return [];
    const totals = sortedIndices.map(i => ({
      name: covData.bins[i].id,
      total: covData.matrix[i].reduce((a, b) => a + b, 0),
    })).sort((a, b) => b.total - a.total);

    return [{
      type: 'bar',
      x: totals.map(t => t.name),
      y: totals.map(t => t.total),
      marker: { color: '#22d3ee', opacity: 0.8 },
      hovertemplate: '%{x}: %{y:.1f}x<extra></extra>',
    }];
  });

  let showDetail = $state(false);

  function handleRowClick(magId) {
    selectedMag.set(magId);
    showDetail = true;
  }

  // Unified detail lookup — checkm2All has taxonomy + MGE for all bins
  let selectedBinData = $derived.by(() => {
    if (!selected || !allBins) return null;
    return allBins.find(b => b.dastool_name === selected)
        || allBins.find(b => b.name === selected)
        || null;
  });
</script>

{#if covData}
  <div class="sticky top-14 z-10 bg-slate-950 pb-2 -mx-4 px-4 pt-1">
  <div class="flex items-center gap-3 flex-wrap text-xs">
    <span class="text-slate-400">Binners:</span>
    {#each binnerDefs as { key, label }}
      {#if binnerCounts[key]}
        <button
          class="px-3 py-1 rounded-md border transition-colors text-center
            {activeBinners.has(key) ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
          onclick={() => toggleBinner(key)}
        >
          {label} <span class="text-slate-500">({binnerCounts[key]})</span>
        </button>
      {/if}
    {/each}
    <span class="text-slate-500">{sortedIndices.length} bins</span>
    <span class="text-slate-600">|</span>
    <span class="text-slate-400">Classifier:</span>
    <button
      class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
      style="min-width: 5rem"
      onclick={cycleSource}
      title={`Click to cycle: ${taxSources.map(s => s.label).join(' \u2192 ')}`}
    >{taxSourceLabel()} &#x25BE;</button>
    <span class="text-slate-600">|</span>
    <label class="text-slate-400 flex items-center gap-1">
      Comp &ge;
      <input type="range" min="0" max="100" step="5" bind:value={minCompleteness}
        class="w-28 h-1 accent-cyan-400" />
      <span class="text-cyan-400 font-mono w-8">{minCompleteness}%</span>
    </label>
    <label class="text-slate-400 flex items-center gap-1">
      Cont &le;
      <input type="range" min="0" max="100" step="1" bind:value={maxContamination}
        class="w-28 h-1 accent-cyan-400" />
      <span class="text-cyan-400 font-mono w-8">{maxContamination}%</span>
    </label>
  </div>
  </div>

  <div class="grid grid-cols-1 {selectedBinData && showDetail ? 'lg:grid-cols-3' : ''} gap-6 mb-6">
    <div class="{selectedBinData && showDetail ? 'lg:col-span-2' : ''}">
      <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
        <div class="flex items-center justify-between mb-2">
          <h3 class="text-sm font-medium text-slate-400">
            Coverage Heatmap ({sortedIndices.length} bins x {covData.sample_names.length} samples)
          </h3>
          <button
            class="text-xs px-2 py-0.5 rounded border border-cyan-400 bg-cyan-400/10 text-cyan-400"
            onclick={cycleSort}
            title={`Sort rows: ${sortDefs.map(d => d.label).join(' \u2192 ')}`}
          >Sort: {sortDefs.find(d => d.key === sortBy)?.label} &#x25BE;</button>
        </div>
        {#if heatmapData}
          <div class="overflow-auto" style="max-height: 700px;">
          <D3Heatmap
            data={heatmapData}
            onRowClick={handleRowClick}
            tooltipFormat={covTooltip}
            colorScale={covColorScale}
            legendLabel={['0x', `${maxDepth.toFixed(0)}x`]}
            selectedRow={selected}
            rowAnnotations={covRowAnnotations}
          />
          </div>
        {/if}
      </div>
    </div>

    {#if selectedBinData && showDetail}
      <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
        <div class="flex items-center justify-between mb-3">
          <h3 class="text-cyan-400 font-semibold">{selectedBinData.name}</h3>
          <button
            class="text-slate-500 hover:text-slate-300 text-lg leading-none"
            onclick={() => showDetail = false}
            title="Close"
          >&times;</button>
        </div>
        <p class="text-xs text-slate-500 mb-3">{selectedBinData.binner}{selectedBinData.is_dastool ? ' (DAS Tool consensus)' : ''}</p>
        <dl class="space-y-2 text-sm">
          <div class="flex justify-between">
            <dt class="text-slate-400">Completeness</dt>
            <dd class="font-mono text-emerald-400">{selectedBinData.completeness}%</dd>
          </div>
          <div class="flex justify-between">
            <dt class="text-slate-400">Contamination</dt>
            <dd class="font-mono text-amber-400">{selectedBinData.contamination}%</dd>
          </div>
          <div class="flex justify-between">
            <dt class="text-slate-400">GC Content</dt>
            <dd class="font-mono">{(selectedBinData.gc * 100).toFixed(1)}%</dd>
          </div>
          <div class="flex justify-between">
            <dt class="text-slate-400">Genome Size</dt>
            <dd class="font-mono">{(selectedBinData.genome_size / 1e6).toFixed(2)} Mbp</dd>
          </div>
          <div class="flex justify-between">
            <dt class="text-slate-400">N50</dt>
            <dd class="font-mono">{(selectedBinData.n50 / 1e3).toFixed(1)} Kb</dd>
          </div>
          {#if selectedBinData.taxonomy && Object.keys(selectedBinData.taxonomy).length > 0}
            <div class="border-t border-slate-700 pt-2 mt-2">
              <dt class="text-slate-400 mb-1 flex items-center gap-2">
                Taxonomy
                {#each taxSources as { key, label }}
                  {#if selectedBinData.taxonomy[key]}
                    <button
                      class="text-[10px] px-1.5 py-0.5 rounded border transition-colors
                        {taxSource === key ? 'border-cyan-400 text-cyan-400' : 'border-slate-600 text-slate-500 hover:text-slate-300'}"
                      onclick={() => taxSource = key}
                    >{label}</button>
                  {/if}
                {/each}
              </dt>
              {#if selectedBinData.taxonomy[taxSource]}
                <dd class="font-mono text-xs">
                  {#each Object.entries(selectedBinData.taxonomy[taxSource]).filter(([,v]) => v) as [rank, val]}
                    <div><span class="text-slate-500">{rank}:</span> {val}</div>
                  {/each}
                </dd>
              {:else}
                <dd class="text-xs text-slate-500 italic">No {taxSource} data for this bin</dd>
              {/if}
            </div>
          {/if}
          <div class="border-t border-slate-700 pt-2 mt-2">
            <dt class="text-slate-400 mb-1">MGE Summary</dt>
            <dd class="font-mono text-xs grid grid-cols-2 gap-1">
              <span>Virus: {selectedBinData.n_virus}</span>
              <span>Plasmid: {selectedBinData.n_plasmid}</span>
              <span>Defense: {selectedBinData.n_defense}</span>
              <span>Integron: {selectedBinData.n_integron}</span>
            </dd>
          </div>
        </dl>
      </div>
    {/if}
  </div>

  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6">
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart
        data={boxData}
        layout={{
          title: { text: 'Depth Distribution per Sample', font: { size: 14 } },
          yaxis: { title: 'Depth (x)', type: 'log' },
          showlegend: false,
        }}
      />
    </div>
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart
        data={barData}
        layout={{
          title: { text: 'Total Depth per Bin', font: { size: 14 } },
          xaxis: { tickangle: -45, tickfont: { size: 8 } },
          yaxis: { title: 'Total depth (x)' },
        }}
      />
    </div>
  </div>
{/if}
