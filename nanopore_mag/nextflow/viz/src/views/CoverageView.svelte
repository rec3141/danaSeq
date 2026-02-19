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

  // Taxonomy source for detail panel
  let taxSource = $state('kaiju');
  const taxSources = [
    { key: 'kaiju', label: 'Kaiju' },
    { key: 'kraken2', label: 'Kraken2' },
    { key: 'rrna', label: 'rRNA' },
  ];

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
      if (next.size > 1) next.delete(key);  // keep at least one active
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

  // Filter bins by active binners, preserving pre-computed Bray-Curtis order
  let filteredIndices = $derived.by(() => {
    if (!covData) return [];
    return covData.bins
      .map((b, i) => ({ b, i }))
      .filter(({ b }) => activeBinners.has(b.binner))
      .map(({ i }) => i);
  });

  // Find the max depth for color scaling (across filtered bins)
  let maxDepth = $derived.by(() => {
    if (!covData || !filteredIndices.length) return 1;
    let mx = 0;
    for (const i of filteredIndices) {
      for (const v of covData.matrix[i]) if (v > mx) mx = v;
    }
    return mx || 1;
  });

  // Heatmap data: filtered bins with raw depths
  let heatmapData = $derived.by(() => {
    if (!covData || !filteredIndices.length) return null;
    return {
      mag_ids: filteredIndices.map(i => covData.bins[i].id),
      module_ids: covData.sample_names.map(s => s),
      module_names: covData.sample_names,
      matrix: filteredIndices.map(i =>
        covData.matrix[i].map(v => Math.round(v * 1000) / 1000)
      ),
    };
  });

  // Tooltip formatter for coverage
  function covTooltip(val, magId, modId, modName) {
    return `<strong>${modId}</strong><br>Bin: ${magId}<br>Depth: ${val.toFixed(2)}x`;
  }

  // Box plot: depth distribution per sample (filtered bins only)
  let boxData = $derived.by(() => {
    if (!covData || !filteredIndices.length) return [];
    const palette = ['#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa', '#fb923c', '#818cf8', '#f472b6'];
    return covData.sample_names.map((sname, j) => ({
      type: 'box',
      name: sname,
      y: filteredIndices.map(i => covData.matrix[i][j]).filter(v => v > 0),
      marker: { color: palette[j % palette.length] },
      boxpoints: 'outliers',
    }));
  });

  // Bar chart: total depth per bin (sorted, filtered)
  let barData = $derived.by(() => {
    if (!covData || !filteredIndices.length) return [];
    const totals = filteredIndices.map(i => ({
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

  function handleRowClick(magId) {
    selectedMag.set(magId);
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
  <div class="mb-4 flex items-center gap-3 flex-wrap text-xs">
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
    <span class="text-slate-500">{filteredIndices.length} bins shown</span>
  </div>

  <div class="grid grid-cols-1 {selectedBinData ? 'lg:grid-cols-3' : ''} gap-6 mb-6">
    <div class="{selectedBinData ? 'lg:col-span-2' : ''}">
      <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
        <h3 class="text-sm font-medium text-slate-400 mb-2">
          Coverage Heatmap ({filteredIndices.length} bins x {covData.sample_names.length} samples, Bray-Curtis clustered)
        </h3>
        {#if heatmapData}
          <div class="overflow-y-auto" style="max-height: 700px;">
          <D3Heatmap
            data={heatmapData}
            onRowClick={handleRowClick}
            tooltipFormat={covTooltip}
            colorDomain={[0, maxDepth]}
            logScale={true}
            legendLabel={['0x', `${maxDepth.toFixed(0)}x`]}
          />
          </div>
        {/if}
      </div>
    </div>

    {#if selectedBinData}
      <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
        <h3 class="text-cyan-400 font-semibold mb-3">{selectedBinData.name}</h3>
        <p class="text-xs text-slate-500 -mt-2 mb-3">{selectedBinData.binner}{selectedBinData.is_dastool ? ' (DAS Tool consensus)' : ''}</p>
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
