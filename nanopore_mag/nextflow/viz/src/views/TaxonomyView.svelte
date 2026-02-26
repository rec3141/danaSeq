<script>
  import D3Sunburst from '../components/charts/D3Sunburst.svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { onMount } from 'svelte';
  import { taxonomySunburst, checkm2All, loadCheckm2All } from '../stores/data.js';
  import { selectedMag } from '../stores/selection.js';

  let taxData = $derived($taxonomySunburst);
  let allBins = $derived($checkm2All);

  onMount(() => { loadCheckm2All(); });

  // Classifier source (cycles)
  let source = $state('kaiju');
  const sourceDefs = [
    { key: 'kaiju', label: 'Kaiju' },
    { key: 'kraken2', label: 'Kraken2' },
    { key: 'sendsketch', label: 'BBSketch' },
    { key: 'rrna', label: 'rRNA' },
  ];

  function cycleSource() {
    const keys = sourceDefs.map(s => s.key);
    source = keys[(keys.indexOf(source) + 1) % keys.length];
  }
  function sourceLabel() {
    return sourceDefs.find(s => s.key === source)?.label || source;
  }

  // Rank for composition bar chart (cycles)
  let rank = $state('phylum');
  const rankDefs = [
    { key: 'domain', label: 'Domain' },
    { key: 'phylum', label: 'Phylum' },
    { key: 'class', label: 'Class' },
    { key: 'order', label: 'Order' },
    { key: 'family', label: 'Family' },
    { key: 'genus', label: 'Genus' },
  ];

  function cycleRank() {
    const keys = rankDefs.map(r => r.key);
    rank = keys[(keys.indexOf(rank) + 1) % keys.length];
  }
  function rankLabel() {
    return rankDefs.find(r => r.key === rank)?.label || rank;
  }

  // Per-bin bar mode: raw Mb vs percentage
  let pctMode = $state(false);

  // Binner toggles
  let activeBinners = $state(new Set(['dastool']));

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

  let binnerCounts = $derived.by(() => {
    if (!allBins) return {};
    const counts = {};
    for (const b of allBins) {
      const bk = b.is_dastool ? 'dastool' : b.binner;
      counts[bk] = (counts[bk] || 0) + 1;
    }
    return counts;
  });

  // Filtered bins for table
  let filteredBins = $derived.by(() => {
    if (!allBins) return [];
    return allBins.filter(b => {
      const bk = b.is_dastool ? 'dastool' : b.binner;
      return activeBinners.has(bk);
    });
  });

  const palette = ['#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa',
                   '#fb923c', '#2dd4bf', '#818cf8', '#f472b6', '#4ade80',
                   '#e879f9', '#38bdf8', '#94a3b8', '#d4d4d8', '#78716c',
                   '#64748b', '#475569'];

  // Sunburst for current source
  let sunburstData = $derived.by(() => {
    if (!taxData?.sunbursts) return null;
    return taxData.sunbursts[source] || null;
  });

  // Globally stable color assignments: merge ALL sources × ALL binners × each rank
  // Colors never change regardless of which source/binner/rank is selected
  let allRankColors = $derived.by(() => {
    if (!taxData?.composition) return {};
    const colors = {};
    for (const rk of rankDefs.map(r => r.key)) {
      const merged = {};
      for (const binner of Object.keys(taxData.composition)) {
        for (const src of sourceDefs.map(s => s.key)) {
          const rankData = taxData.composition[binner]?.[src]?.[rk];
          if (!rankData) continue;
          for (const [taxon, bp] of Object.entries(rankData)) {
            merged[taxon] = (merged[taxon] || 0) + bp;
          }
        }
      }
      const sorted = Object.entries(merged)
        .filter(([t]) => t !== 'Unclassified')
        .sort((a, b) => b[1] - a[1]);
      const map = {};
      sorted.forEach(([t], i) => { map[t] = palette[i % palette.length]; });
      map['Unclassified'] = '#475569';
      map['Other'] = '#64748b';
      colors[rk] = map;
    }
    return colors;
  });

  // Current rank's color map (for bar charts)
  let taxonColorMap = $derived(allRankColors[rank] || {});

  // Sunburst gets a per-depth color map: {1: {Bacteria: '#...', ...}, 2: {Actinomycetota: '#...', ...}, ...}
  const sunburstRanks = ['domain', 'phylum', 'class', 'order'];
  let sunburstColorMaps = $derived.by(() => {
    const maps = {};
    sunburstRanks.forEach((rk, i) => {
      maps[i + 1] = allRankColors[rk] || {};
    });
    return maps;
  });

  // Composition bar chart — rebuild top list from composition at selected rank
  let barData = $derived.by(() => {
    if (!taxData?.composition) return [];
    const merged = {};
    for (const binner of activeBinners) {
      const rankData = taxData.composition[binner]?.[source]?.[rank];
      if (!rankData) continue;
      for (const [taxon, bp] of Object.entries(rankData)) {
        merged[taxon] = (merged[taxon] || 0) + bp;
      }
    }
    const entries = Object.entries(merged)
      .filter(([name]) => name !== 'Unclassified')
      .sort((a, b) => b[1] - a[1]);
    const unclassified = merged['Unclassified'] || 0;
    const top = entries.slice(0, 15);
    const otherBp = entries.slice(15).reduce((s, [, v]) => s + v, 0);
    if (otherBp > 0) top.push(['Other', otherBp]);
    if (unclassified > 0) top.push(['Unclassified', unclassified]);

    if (!top.length) return [];
    return [{
      type: 'bar',
      x: top.map(([name]) => name),
      y: top.map(([, bp]) => bp / 1e6),
      marker: {
        color: top.map(([name]) => taxonColorMap[name] || '#64748b'),
      },
      hovertemplate: '%{x}: %{y:.2f} Mb<extra></extra>',
    }];
  });

  // Per-bin stacked bar: bins on x-axis, length-weighted composition on y, stacked by taxon
  let stackedBarData = $derived.by(() => {
    if (!filteredBins.length) return [];

    // Extract per-bin composition at the selected source+rank
    const bins = filteredBins.map(b => {
      const name = b.is_dastool ? (b.dastool_name || b.name) : b.name;
      const comp = b.composition?.[source]?.[rank] || {};
      const totalBp = Object.values(comp).reduce((s, v) => s + v, 0);
      // Dominant taxon = taxon with most bp in this bin
      let dominant = 'Unclassified';
      let maxBp = 0;
      for (const [t, bp] of Object.entries(comp)) {
        if (bp > maxBp) { maxBp = bp; dominant = t; }
      }
      return { name, comp, totalBp, dominant };
    });

    // Sort: group by dominant taxon (ordered by global total), then by size within group
    // First pass: rank taxa by global total across all bins
    const globalTotals = {};
    for (const b of bins) {
      globalTotals[b.dominant] = (globalTotals[b.dominant] || 0) + b.totalBp;
    }
    const taxonRank = {};
    Object.entries(globalTotals)
      .filter(([t]) => t !== 'Unclassified')
      .sort((a, b) => b[1] - a[1])
      .forEach(([t], i) => { taxonRank[t] = i; });
    taxonRank['Unclassified'] = 9999;

    bins.sort((a, b) =>
      (taxonRank[a.dominant] ?? 999) - (taxonRank[b.dominant] ?? 999)
      || b.totalBp - a.totalBp
    );

    // Accumulate global taxon totals to determine top taxa
    const taxonTotals = {};
    for (const b of bins) {
      for (const [taxon, bp] of Object.entries(b.comp)) {
        taxonTotals[taxon] = (taxonTotals[taxon] || 0) + bp;
      }
    }
    const taxa = Object.entries(taxonTotals)
      .filter(([t]) => t !== 'Unclassified')
      .sort((a, b) => b[1] - a[1])
      .map(([t]) => t);

    const topTaxa = taxa.slice(0, 15);
    const otherTaxa = new Set(taxa.slice(15));
    const allTaxa = [...topTaxa];
    if (otherTaxa.size > 0) allTaxa.push('Other');
    if (taxonTotals['Unclassified']) allTaxa.push('Unclassified');

    const binNames = bins.map(b => b.name);

    // Build one trace per taxon, each bin contributes its bp for that taxon
    return allTaxa.map((taxon, i) => {
      const color = taxonColorMap[taxon] || palette[i % palette.length];
      return {
        type: 'bar',
        name: taxon,
        x: binNames,
        y: bins.map(b => {
          let raw;
          if (taxon === 'Other') {
            raw = 0;
            for (const [t, bp] of Object.entries(b.comp)) {
              if (otherTaxa.has(t)) raw += bp;
            }
          } else {
            raw = b.comp[taxon] || 0;
          }
          return pctMode ? (b.totalBp > 0 ? raw / b.totalBp * 100 : 0) : raw / 1e6;
        }),
        marker: { color },
        hovertemplate: pctMode
          ? `${taxon}<br>%{x}: %{y:.1f}%<extra></extra>`
          : `${taxon}<br>%{x}: %{y:.2f} Mb<extra></extra>`,
      };
    });
  });

  // Table
  let tableColumns = $derived([
    { key: 'name', label: 'Bin' },
    { key: 'domain', label: 'Domain' },
    { key: 'phylum', label: 'Phylum' },
    { key: 'class', label: 'Class' },
    { key: 'order', label: 'Order' },
    { key: 'family', label: 'Family' },
    { key: 'genus', label: 'Genus' },
  ]);

  let tableRows = $derived.by(() => {
    return filteredBins.map(b => {
      const tax = b.taxonomy?.[source] || {};
      return {
        name: b.is_dastool ? (b.dastool_name || b.name) : b.name,
        domain: tax.domain || '-',
        phylum: tax.phylum || '-',
        class: tax.class || '-',
        order: tax.order || '-',
        family: tax.family || '-',
        genus: tax.genus || '-',
      };
    });
  });
</script>

<div class="sticky top-14 z-10 bg-slate-950 pb-2 -mx-4 px-4 pt-1">
<div class="flex items-center gap-3 flex-wrap text-xs">
  <span class="text-slate-400">Classifier:</span>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
    style="min-width: 5rem"
    onclick={cycleSource}
    title={`Click to cycle: ${sourceDefs.map(s => s.label).join(' → ')}`}
  >
    {sourceLabel()} &#x25BE;
  </button>
  <span class="text-slate-400">Rank:</span>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
    style="min-width: 4.5rem"
    onclick={cycleRank}
    title={`Click to cycle: ${rankDefs.map(r => r.label).join(' → ')}`}
  >
    {rankLabel()} &#x25BE;
  </button>
  <span class="text-slate-600">|</span>
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
</div>
</div>

{#if taxData || filteredBins.length}
  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-6">
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <h3 class="text-sm font-medium text-slate-400 mb-2">Community Composition ({sourceLabel()}, all contigs)</h3>
      {#if sunburstData}
        <D3Sunburst data={sunburstData} colorMaps={sunburstColorMaps} />
      {:else}
        <p class="text-slate-500 text-sm py-20 text-center">No sunburst data for {sourceLabel()}</p>
      {/if}
    </div>
    <div class="flex flex-col gap-6">
      <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
        <h3 class="text-sm font-medium text-slate-400 mb-2">{rankLabel()} Composition ({sourceLabel()}, length-weighted)</h3>
        {#if barData.length}
          <PlotlyChart
            data={barData}
            layout={{
              xaxis: { tickangle: -45, tickfont: { size: 9 } },
              yaxis: { title: 'Total length (Mb)', fixedrange: true },
              showlegend: false,
              margin: { b: 120 },
              height: 300,
            }}
          />
        {/if}
      </div>
      <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
        <div class="flex items-center justify-between mb-2">
          <h3 class="text-sm font-medium text-slate-400">Per-Bin {rankLabel()} ({sourceLabel()}, {filteredBins.length} bins)</h3>
          <button
            class="text-xs px-2 py-0.5 rounded border transition-colors
              {pctMode ? 'border-cyan-400 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
            onclick={() => pctMode = !pctMode}
          >{pctMode ? '%' : 'Mb'}</button>
        </div>
        {#if stackedBarData.length}
          <PlotlyChart
            data={stackedBarData}
            layout={{
              barmode: 'stack',
              xaxis: { tickangle: -45, tickfont: { size: 8 } },
              yaxis: { title: pctMode ? 'Composition (%)' : 'Genome size (Mb)', fixedrange: true },
              showlegend: true,
              legend: { font: { size: 9 }, orientation: 'h', y: -0.55 },
              margin: { b: 160 },
              height: 420,
            }}
          />
        {/if}
      </div>
    </div>
  </div>

  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
    <h3 class="text-sm font-medium text-slate-400 mb-2">Per-Bin Taxonomy ({sourceLabel()}, {filteredBins.length} bins)</h3>
    <DataTable
      columns={tableColumns}
      rows={tableRows}
      onRowClick={(row) => selectedMag.set(row.name)}
      maxHeight="400px"
    />
  </div>
{/if}
