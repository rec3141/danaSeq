<script>
  import { onMount } from 'svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import QualityBadge from '../components/ui/QualityBadge.svelte';
  import { checkm2All, loadCheckm2All, contigExplorer, loadContigExplorer } from '../stores/data.js';
  import { selectedMag } from '../stores/selection.js';

  let allBins = $derived($checkm2All);
  let contigs = $derived($contigExplorer);
  let selected = $derived($selectedMag);

  onMount(() => { loadCheckm2All(); loadContigExplorer(); });

  // Binner toggles — DAS Tool on by default
  let activeBinners = $state(new Set(['dastool']));

  const binnerDefs = [
    { key: 'dastool', label: 'DAS Tool' },
    { key: 'semibin', label: 'SemiBin2' },
    { key: 'metabat', label: 'MetaBAT2' },
    { key: 'maxbin',  label: 'MaxBin2' },
    { key: 'lorbin',  label: 'LorBin' },
    { key: 'comebin', label: 'COMEBin' },
  ];

  const binnerColors = {
    dastool: '#22d3ee', semibin: '#22d3ee', metabat: '#34d399', maxbin: '#fbbf24',
    lorbin: '#a78bfa', comebin: '#fb923c',
  };
  const binnerSymbols = {
    semibin: 'circle', metabat: 'square', maxbin: 'diamond',
    lorbin: 'triangle-up', comebin: 'pentagon',
  };

  function toggleBinner(key) {
    const next = new Set(activeBinners);
    if (next.has(key)) {
      if (next.size > 1) next.delete(key);
    } else {
      next.add(key);
    }
    activeBinners = next;
  }

  let binnerCounts = $derived.by(() => {
    if (!allBins) return {};
    const counts = {};
    for (const b of allBins) {
      const bk = b.is_dastool ? 'dastool' : b.binner;
      counts[bk] = (counts[bk] || 0) + 1;
    }
    return counts;
  });

  // Filtered bins
  let filteredBins = $derived.by(() => {
    if (!allBins) return [];
    return allBins.filter(b => {
      const bk = b.is_dastool ? 'dastool' : b.binner;
      return activeBinners.has(bk);
    });
  });

  // Quality tier helper
  function qualityTier(b) {
    if (b.completeness >= 90 && b.contamination < 5) return 'HQ';
    if (b.completeness >= 50 && b.contamination < 10) return 'MQ';
    return 'LQ';
  }

  let scatterData = $derived.by(() => {
    if (!filteredBins.length) return [];

    const tiers = [
      { quality: 'HQ', color: '#34d399', label: 'HQ' },
      { quality: 'MQ', color: '#fbbf24', label: 'MQ' },
      { quality: 'LQ', color: '#64748b', label: 'LQ' },
    ];

    return tiers.map(tier => {
      const bins = filteredBins.filter(b => qualityTier(b) === tier.quality);
      return {
        type: 'scatter',
        mode: 'markers',
        name: tier.label,
        x: bins.map(b => b.contamination),
        y: bins.map(b => b.completeness),
        text: bins.map(b => {
          const name = b.is_dastool ? (b.dastool_name || b.name) : b.name;
          return `${name}<br>Comp: ${b.completeness}%<br>Cont: ${b.contamination}%<br>Binner: ${b.binner}<br>Size: ${(b.genome_size/1e6).toFixed(1)} Mbp`;
        }),
        customdata: bins.map(b => b.is_dastool ? (b.dastool_name || b.name) : b.name),
        marker: {
          color: tier.color,
          size: bins.map(b => Math.max(9, Math.sqrt(b.genome_size / 1e5) * 1.5)),
          symbol: bins.map(b => binnerSymbols[b.binner] || 'circle'),
          line: { color: '#0f172a', width: 1 },
        },
        hoverinfo: 'text',
      };
    }).filter(t => t.x.length > 0);
  });

  let scatterLayout = $derived({
    title: { text: 'CheckM2 Quality Assessment', font: { size: 14 } },
    xaxis: { title: 'Contamination (%)', range: [0, 15], constrain: 'domain', constraintoward: 'left' },
    yaxis: { title: 'Completeness (%)', range: [0, 105], fixedrange: true },
    shapes: [
      // MIMAG thresholds
      { type: 'line', x0: 5, x1: 5, y0: 0, y1: 105, line: { color: '#34d399', width: 1, dash: 'dash' } },
      { type: 'line', x0: 10, x1: 10, y0: 0, y1: 105, line: { color: '#fbbf24', width: 1, dash: 'dash' } },
      { type: 'line', x0: -1, x1: 50, y0: 90, y1: 90, line: { color: '#34d399', width: 1, dash: 'dash' } },
      { type: 'line', x0: -1, x1: 50, y0: 50, y1: 50, line: { color: '#fbbf24', width: 1, dash: 'dash' } },
    ],
    annotations: [
      { x: 2.5, y: 103, text: 'HQ', showarrow: false, font: { color: '#34d399', size: 11 } },
      { x: 7.5, y: 103, text: 'MQ', showarrow: false, font: { color: '#fbbf24', size: 11 } },
    ],
  });

  function handleClick(event) {
    const cd = event?.points?.[0]?.customdata;
    if (cd) selectedMag.set(cd);
  }

  // Taxonomy source for detail panel
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

  let tableColumns = $derived([
    { key: 'name', label: 'Bin' },
    { key: 'quality', label: 'Quality', render: (v) => {
      const cls = v === 'HQ' ? 'text-emerald-400' : v === 'MQ' ? 'text-amber-400' : 'text-slate-400';
      return `<span class="${cls} font-semibold">${v}</span>`;
    }},
    { key: 'completeness', label: 'Comp %' },
    { key: 'contamination', label: 'Cont %' },
    { key: 'binner', label: 'Binner' },
    { key: 'size', label: 'Size', render: (v) => (v/1e6).toFixed(2) + ' Mb' },
    { key: 'n50', label: 'N50', render: (v) => (v/1e3).toFixed(1) + ' Kb' },
    { key: 'contigs', label: 'Contigs' },
    { key: 'taxonomy', label: `Taxonomy (${taxSourceLabel()})` },
  ]);

  let tableRows = $derived.by(() => {
    return filteredBins.map(b => {
      const tax = b.taxonomy?.[taxSource] || {};
      const lineage = ['domain', 'phylum', 'class', 'order', 'family', 'genus']
        .map(r => tax[r]).filter(Boolean).join('; ');
      return {
        name: b.is_dastool ? (b.dastool_name || b.name) : b.name,
        quality: qualityTier(b),
        completeness: b.completeness,
        contamination: b.contamination,
        binner: b.binner,
        size: b.genome_size,
        n50: b.n50,
        contigs: b.n_contigs ?? 0,
        taxonomy: lineage || '-',
      };
    });
  });

  let selectedBinData = $derived.by(() => {
    if (!selected || !allBins) return null;
    return allBins.find(b => b.dastool_name === selected)
        || allBins.find(b => b.name === selected)
        || null;
  });

  // Diagnostic histogram cycling
  let diagMetric = $state('gc');
  const diagDefs = [
    { key: 'gc', label: 'GC%', field: 'gc', unit: '%' },
    { key: 'depth', label: 'Coverage', field: 'depth', unit: 'x', log: true },
    { key: 'length', label: 'Length', field: 'length', unit: 'bp', log: true },
  ];
  function cycleDiag() {
    const keys = diagDefs.map(d => d.key);
    diagMetric = keys[(keys.indexOf(diagMetric) + 1) % keys.length];
  }
  function diagLabel() {
    return diagDefs.find(d => d.key === diagMetric)?.label || diagMetric;
  }

  // Find contigs belonging to selected bin
  let selectedBinContigs = $derived.by(() => {
    if (!selected || !contigs?.contigs) return [];
    const allCtgs = contigs.contigs;
    // Try DAS Tool bin field first
    let binCtgs = allCtgs.filter(c => c.bin === selected);
    if (binCtgs.length > 0) return binCtgs;
    // Try per-binner fields
    const binnerFields = { semibin: 'semibin_bin', metabat: 'metabat_bin', maxbin: 'maxbin_bin',
                           lorbin: 'lorbin_bin', comebin: 'comebin_bin' };
    for (const [, field] of Object.entries(binnerFields)) {
      binCtgs = allCtgs.filter(c => c[field] === selected);
      if (binCtgs.length > 0) return binCtgs;
    }
    return [];
  });

  // Histogram data: background (all contigs) + foreground (selected bin)
  let histogramData = $derived.by(() => {
    if (!contigs?.contigs) return { traces: [], xaxisConfig: {} };
    const def = diagDefs.find(d => d.key === diagMetric);
    if (!def) return { traces: [], xaxisConfig: {} };

    const allVals = contigs.contigs
      .map(c => c[def.field])
      .filter(v => v != null && v > 0);
    const binVals = selectedBinContigs
      .map(c => c[def.field])
      .filter(v => v != null && v > 0);

    const xvals = def.log ? allVals.map(v => Math.log10(v)) : allVals;
    const xbin = def.log ? binVals.map(v => Math.log10(v)) : binVals;

    // Compute shared bin edges
    const min = Math.min(...xvals);
    const max = Math.max(...xvals);
    const nBins = 50;
    const step = (max - min) / nBins;
    const edges = Array.from({ length: nBins + 1 }, (_, i) => min + i * step);

    function histogram(vals) {
      const counts = new Array(nBins).fill(0);
      for (const v of vals) {
        let idx = Math.floor((v - min) / step);
        if (idx >= nBins) idx = nBins - 1;
        if (idx < 0) idx = 0;
        counts[idx]++;
      }
      const total = vals.length || 1;
      return counts.map(c => c / total * 100);
    }

    const bgPct = histogram(xvals);
    const fgPct = xbin.length ? histogram(xbin) : [];

    // Bin centers for x-axis (keep in log10 space for log metrics)
    const centers = edges.slice(0, -1).map(e => e + step / 2);

    // For log metrics, build tick labels showing real values
    let xaxisConfig;
    if (def.log) {
      // Nice tick values in real space
      const tickVals = [];
      const tickText = [];
      const logMin = Math.floor(min);
      const logMax = Math.ceil(max);
      for (let e = logMin; e <= logMax; e++) {
        const realVal = Math.pow(10, e);
        tickVals.push(e);
        if (realVal >= 1e6) tickText.push((realVal / 1e6).toFixed(0) + 'M');
        else if (realVal >= 1e3) tickText.push((realVal / 1e3).toFixed(0) + 'K');
        else tickText.push(realVal.toFixed(realVal < 1 ? 2 : 0));
      }
      xaxisConfig = { title: diagLabel(), tickvals: tickVals, ticktext: tickText };
    } else {
      xaxisConfig = { title: diagLabel() };
    }

    // Hover: show real values for log metrics
    const hoverX = def.log
      ? centers.map(c => { const v = Math.pow(10, c); return v >= 1e3 ? (v/1e3).toFixed(1)+'K' : v.toFixed(1); })
      : centers.map(c => c.toFixed(1));

    const traces = [{
      type: 'bar',
      name: 'All contigs',
      x: centers,
      y: bgPct,
      text: hoverX,
      textposition: 'none',
      marker: { color: 'rgba(100,116,139,0.35)' },
      hovertemplate: '%{text} ' + (def.unit || '') + ': %{y:.1f}%<extra>All</extra>',
    }];

    if (fgPct.length) {
      traces.push({
        type: 'bar',
        name: selectedBinData?.name || selected,
        x: centers,
        y: fgPct,
        text: hoverX,
        textposition: 'none',
        marker: { color: '#22d3ee', opacity: 0.6 },
        hovertemplate: '%{text} ' + (def.unit || '') + ': %{y:.1f}%<extra>' + (selectedBinData?.name || selected) + '</extra>',
      });
    }

    return { traces, xaxisConfig };
  });

  let histogramLayout = $derived.by(() => {
    return {
      barmode: 'overlay',
      xaxis: histogramData.xaxisConfig || { title: diagLabel() },
      yaxis: { title: '% of contigs', fixedrange: true },
      showlegend: false,
      margin: { t: 10, r: 10, b: 50, l: 50 },
    };
  });

</script>

{#if allBins}
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
    <span class="text-slate-500">{filteredBins.length} bins</span>
    <span class="text-slate-600">|</span>
    <span class="text-slate-400">Classifier:</span>
    <button
      class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
      style="min-width: 5rem"
      onclick={cycleSource}
      title={`Click to cycle: ${taxSources.map(s => s.label).join(' \u2192 ')}`}
    >{taxSourceLabel()} &#x25BE;</button>
  </div>

  <div class="grid grid-cols-1 lg:grid-cols-3 gap-6 mb-6">
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700 modebar-left">
      <PlotlyChart
        data={scatterData}
        layout={scatterLayout}
        onclick={handleClick}
        config={{
          modeBarButtonsToRemove: [
            'select2d', 'lasso2d', 'autoScale2d',
            'toggleSpikelines', 'hoverCompareCartesian', 'hoverClosestCartesian',
            'zoomIn2d', 'zoomOut2d', 'zoom2d',
          ],
          modeBarOrientation: 'v',
        }}
      />
      <div class="flex gap-4 mt-2 text-xs text-slate-400 justify-center">
        {#each Object.entries(binnerSymbols) as [name, sym]}
          <span class="flex items-center gap-1">
            <span class="text-slate-400">{sym === 'circle' ? '\u25CF' : sym === 'square' ? '\u25A0' : sym === 'diamond' ? '\u25C6' : sym === 'triangle-up' ? '\u25B2' : '\u2B1F'}</span>
            {name}
          </span>
        {/each}
      </div>
    </div>

    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      {#if selectedBinData}
        <h3 class="text-cyan-400 font-semibold mb-1">{selectedBinData.name}</h3>
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
          <div class="flex justify-between">
            <dt class="text-slate-400">Contigs</dt>
            <dd class="font-mono">{selectedBinData.n_contigs ?? '-'}</dd>
          </div>
          {#if selectedBinData.taxonomy && Object.keys(selectedBinData.taxonomy).length > 0}
            <div class="border-t border-slate-700 pt-2 mt-2">
              <dt class="text-slate-400 mb-1">Taxonomy ({taxSourceLabel()})</dt>
              {#if selectedBinData.taxonomy[taxSource]}
                <dd class="font-mono text-xs">
                  {#each Object.entries(selectedBinData.taxonomy[taxSource]).filter(([,v]) => v) as [rank, val]}
                    <div><span class="text-slate-500">{rank}:</span> {val}</div>
                  {/each}
                </dd>
              {:else}
                <dd class="text-xs text-slate-500 italic">No {taxSourceLabel()} data for this bin</dd>
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
      {:else}
        <p class="text-slate-500 text-sm py-8 text-center">Click a point or table row to view bin details</p>
      {/if}
    </div>

    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700 flex flex-col">
      <div class="flex items-center justify-between mb-2">
        <h3 class="text-sm font-medium text-slate-400">
          {diagLabel()} Distribution
        </h3>
        <button
          class="text-xs px-2 py-0.5 rounded border border-cyan-400 bg-cyan-400/10 text-cyan-400"
          onclick={cycleDiag}
          title={`Click to cycle: ${diagDefs.map(d => d.label).join(' \u2192 ')}`}
        >{diagLabel()} &#x25BE;</button>
      </div>
      {#if selectedBinData}
        <p class="text-xs text-slate-500 mb-1">
          <span class="text-cyan-400">{selectedBinData.name}</span>
          ({selectedBinContigs.length} contigs)
        </p>
      {/if}
      <div class="flex-1 min-h-0">
        {#if histogramData.traces?.length}
          <PlotlyChart
            data={histogramData.traces}
            layout={histogramLayout}
            config={{ displayModeBar: false }}
          />
        {:else}
          <p class="text-slate-500 text-sm py-8 text-center">Loading contig data...</p>
        {/if}
      </div>
    </div>
  </div>

  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
    <h3 class="text-sm font-medium text-slate-400 mb-2">Bin Quality ({filteredBins.length} bins)</h3>
    <DataTable
      columns={tableColumns}
      rows={tableRows}
      onRowClick={(row) => selectedMag.set(row.name)}
      selectedId={selected}
      maxHeight="500px"
    />
  </div>
{/if}

<style>
  .modebar-left :global(.modebar) {
    left: 0 !important;
    right: auto !important;
  }
</style>
