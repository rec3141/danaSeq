<script>
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import D3Heatmap from '../components/charts/D3Heatmap.svelte';
  import { onMount } from 'svelte';
  import { coverage, binQuality, loadBinQuality, samples, sampleMetadata, ensureGtdbTaxonomy } from '../stores/data.js';
  import { selectedMag } from '../stores/selection.js';
  import { tsne, brayCurtis } from '../lib/tsne.js';

  let covData = $derived($coverage);
  let allBins = $derived($binQuality);
  let selected = $derived($selectedMag);

  onMount(() => { loadBinQuality(); });

  // Binner toggle state — DAS Tool on by default
  let activeBinners = $state(new Set(['dastool']));

  // Taxonomy source
  let taxSource = $state('kaiju');
  const taxSources = [
    { key: 'kaiju', label: 'Kaiju' },
    { key: 'kraken2', label: 'Kraken2' },
    { key: 'sendsketch', label: 'BBSketch' },
    { key: 'rrna', label: 'rRNA' },
    { key: 'gtdbtk', label: 'GTDB-Tk' },
  ];
  // GTDB-Tk lineages may need filling in from phylotree.json (ensureGtdbTaxonomy)
  $effect(() => { if (taxSource === 'gtdbtk') ensureGtdbTaxonomy(); });

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

  // Lookup map: bin name -> binQuality record
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

  // Sample t-SNE: Bray-Curtis on fourth-root relative abundance of the filtered
  // bins (a bin's relabund = summed share of the sample's mapped contig bases).
  // Samples with nothing in the filtered bins are left out. Colouring is kept
  // separate so changing it doesn't rerun the embedding.
  let tsneEmbedding = $derived.by(() => {
    if (!covData?.relabund || !sortedIndices.length) return null;
    const rows = [];
    covData.sample_names.forEach((name, j) => {
      const ra = sortedIndices.map(i => covData.relabund[i][j]);
      const binned = ra.reduce((a, b) => a + b, 0);
      if (binned > 0) rows.push({ name, ra, binned });
    });
    if (rows.length < 4) return { tooFew: true, n: rows.length };
    const coords = tsne(brayCurtis(rows.map(s => s.ra.map(v => Math.pow(v, 0.25)))));
    const ids = sortedIndices.map(i => covData.bins[i].id);
    for (const s of rows) {
      s.top = s.ra.map((v, k) => [v, ids[k]]).sort((a, b) => b[0] - a[0]).slice(0, 3)
        .map(([v, id]) => `${id}: ${(v * 100).toFixed(1)}%`).join('<br>');
    }
    return { n: rows.length, samples: rows, coords };
  });

  // Sample metadata by coverage sample name (samples.json maps the name to the
  // FLOWCELL:barcodeNN id metadata.json is keyed on).
  let metaBySample = $derived.by(() => {
    if (!$samples || !$sampleMetadata) return {};
    const out = {};
    for (const s of $samples) {
      const m = $sampleMetadata[s.id];
      if (m) out[s.name] = m;
    }
    return out;
  });

  // Colour-by choices: the binned fraction, then every metadata column that
  // has a value for at least one sample in the embedding.
  let tsneColorOptions = $derived.by(() => {
    const cols = new Set();
    for (const s of tsneEmbedding?.samples ?? []) {
      for (const k of Object.keys(metaBySample[s.name] ?? {})) cols.add(k);
    }
    return [...cols].sort();
  });
  let tsneColorBy = $state('');
  let tsneColorCol = $derived(tsneColorOptions.includes(tsneColorBy) ? tsneColorBy : '');

  const categoryPalette = ['#22d3ee', '#f472b6', '#fbbf24', '#34d399', '#a78bfa', '#fb923c',
    '#60a5fa', '#f87171', '#a3e635', '#e879f9', '#2dd4bf', '#facc15', '#c084fc', '#4ade80'];

  let sampleTsne = $derived.by(() => {
    const emb = tsneEmbedding;
    if (!emb || emb.tooFew) return emb;
    const col = tsneColorCol;
    const marker = { size: 11, line: { color: '#0f172a', width: 1 } };
    const point = (idx) => ({
      type: 'scatter',
      mode: 'markers',
      x: idx.map(i => emb.coords[i][0]),
      y: idx.map(i => emb.coords[i][1]),
      hovertext: idx.map(i => {
        const s = emb.samples[i];
        const v = col ? metaBySample[s.name]?.[col] : null;
        return `<b>${s.name}</b>${col ? `<br>${col}: ${v ?? '\u2013'}` : ''}` +
          `<br>In shown bins: ${(s.binned * 100).toFixed(1)}%<br>${s.top}`;
      }),
      hovertemplate: '%{hovertext}<extra></extra>',
    });
    const all = emb.samples.map((_, i) => i);
    if (!col) {
      return { n: emb.n, showlegend: false, traces: [{ ...point(all), marker: {
        ...marker,
        color: emb.samples.map(s => s.binned * 100),
        colorscale: 'Viridis',
        colorbar: { title: { text: '% in bins', side: 'right' }, thickness: 10 },
      } }] };
    }
    const value = (i) => metaBySample[emb.samples[i].name]?.[col];
    const has = all.filter(i => value(i) != null && value(i) !== '');
    const missing = all.filter(i => !has.includes(i));
    const missingTrace = missing.length
      ? [{ ...point(missing), name: 'No value', marker: { ...marker, color: '#475569' } }]
      : [];
    if (has.every(i => typeof value(i) === 'number')) {
      return { n: emb.n, showlegend: false, traces: [...missingTrace, { ...point(has), marker: {
        ...marker,
        color: has.map(value),
        colorscale: 'Viridis',
        colorbar: { title: { text: col, side: 'right' }, thickness: 10 },
      } }] };
    }
    const groups = new Map();
    for (const i of has) {
      const k = String(value(i));
      if (!groups.has(k)) groups.set(k, []);
      groups.get(k).push(i);
    }
    const traces = [...groups.keys()].sort((a, b) => a.localeCompare(b, undefined, { numeric: true }))
      .map((k, g) => ({ ...point(groups.get(k)), name: k,
        marker: { ...marker, color: categoryPalette[g % categoryPalette.length] } }));
    return { n: emb.n, showlegend: true, traces: [...traces, ...missingTrace] };
  });

  // Bin t-SNE (the transpose): one embedding of every binner's bins over samples,
  // computed in preprocess (Bray-Curtis on fourth-root relabund); the binner and
  // quality filters only choose which bins are drawn.
  let tsneMode = $state('samples'); // 'samples' | 'bins'
  let hasBinTsne = $derived(!!covData?.bin_tsne?.some(Boolean));
  let showBins = $derived(tsneMode === 'bins' && hasBinTsne);

  const binColorDefs = [
    { key: 'phylum', label: 'Phylum' },
    { key: 'quality', label: 'Quality' },
    { key: 'binner', label: 'Binner' },
    { key: 'relabund', label: 'Total relabund' },
  ];
  let binColorBy = $state('phylum');
  const binnerLabel = (k) => binnerDefs.find(d => d.key === k)?.label ?? k;

  let binTsne = $derived.by(() => {
    if (!showBins) return null;
    const idx = sortedIndices.filter(i => covData.bin_tsne[i]);
    if (!idx.length) return { n: 0, traces: [] };
    const lineage = (tax) => tax
      ? ['domain', 'phylum', 'class', 'order', 'family', 'genus'].map(r => tax[r]).filter(Boolean).join('; ')
      : '';
    const total = (i) => covData.relabund[i].reduce((a, b) => a + b, 0);
    const point = (sub) => ({
      type: 'scattergl',
      mode: 'markers',
      x: sub.map(i => covData.bin_tsne[i][0]),
      y: sub.map(i => covData.bin_tsne[i][1]),
      customdata: sub.map(i => covData.bins[i].id),
      hovertext: sub.map(i => {
        const b = covData.bins[i], m = binMetaMap[b.id];
        const q = m ? `<br>${qualityTier(m)} · ${m.completeness.toFixed(0)}% comp · ${m.contamination.toFixed(1)}% cont` : '';
        const tax = lineage(m?.taxonomy?.[taxSource]);
        return `<b>${b.id}</b> (${binnerLabel(b.binner)})${q}` +
          `${tax ? `<br>${tax}` : ''}<br>Total relabund: ${(total(i) * 100).toFixed(2)}%`;
      }),
      hovertemplate: '%{hovertext}<extra></extra>',
    });
    const marker = { size: idx.length > 5000 ? 4 : 6, opacity: 0.85 };

    if (binColorBy === 'relabund') {
      return { n: idx.length, showlegend: false, traces: [{ ...point(idx), marker: {
        ...marker,
        color: idx.map(i => Math.log10(Math.max(total(i), 1e-6))),
        colorscale: 'Viridis',
        colorbar: { title: { text: 'log10 relabund', side: 'right' }, thickness: 10 },
      } }] };
    }
    const groupOf = {
      phylum: (i) => binMetaMap[covData.bins[i].id]?.taxonomy?.[taxSource]?.phylum || 'Unclassified',
      quality: (i) => { const m = binMetaMap[covData.bins[i].id]; return m ? qualityTier(m) : 'No CheckM2'; },
      binner: (i) => binnerLabel(covData.bins[i].binner),
    }[binColorBy];
    const groups = new Map();
    for (const i of idx) {
      const k = groupOf(i);
      if (!groups.has(k)) groups.set(k, []);
      groups.get(k).push(i);
    }
    // Largest groups get the palette; the tail and the unassigned are grey
    const tierOrder = ['HQ', 'MQ', 'LQ', 'No CheckM2'];
    let keys = [...groups.keys()].sort((a, b) => binColorBy === 'quality'
      ? tierOrder.indexOf(a) - tierOrder.indexOf(b)
      : groups.get(b).length - groups.get(a).length);
    const grey = new Set(['Unclassified', 'No CheckM2', 'Other']);
    if (binColorBy === 'phylum' && keys.length > categoryPalette.length) {
      const top = keys.filter(k => !grey.has(k)).slice(0, categoryPalette.length);
      const other = keys.filter(k => !top.includes(k) && k !== 'Unclassified').flatMap(k => groups.get(k));
      for (const k of keys) if (!top.includes(k) && k !== 'Unclassified') groups.delete(k);
      if (other.length) groups.set('Other', other);
      keys = [...top, ...['Other', 'Unclassified'].filter(k => groups.has(k))];
    }
    const fixed = { HQ: '#34d399', MQ: '#fbbf24', LQ: '#64748b' };
    let g = 0;
    const traces = keys.map((k, rank) => ({ ...point(groups.get(k)), name: `${k} (${groups.get(k).length})`,
      legendrank: rank,
      marker: { ...marker, color: fixed[k] && binColorBy === 'quality' ? fixed[k]
        : grey.has(k) ? '#475569' : categoryPalette[g++ % categoryPalette.length] } }));
    // Draw grey groups first so coloured points sit on top (legendrank keeps the legend order)
    traces.sort((a, b) => (a.marker.color === '#475569' ? 0 : 1) - (b.marker.color === '#475569' ? 0 : 1));
    return { n: idx.length, showlegend: true, traces };
  });

  function handleBinTsneClick(e) {
    const id = e?.points?.[0]?.customdata;
    if (id) handleRowClick(id);
  }

  let showDetail = $state(false);

  function handleRowClick(magId) {
    selectedMag.set(magId);
    showDetail = true;
  }

  // Unified detail lookup — binQuality has taxonomy + MGE for all bins
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

  {#if sortedIndices.length}
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
    {#if sampleTsne || hasBinTsne}
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <div class="flex items-center gap-3 flex-wrap text-xs text-slate-400 mb-1">
        {#if hasBinTsne}
          <div class="flex rounded-md border border-slate-600 overflow-hidden">
            {#each [['samples', 'Samples'], ['bins', 'Bins']] as [key, label]}
              <button
                class="px-3 py-1 transition-colors {tsneMode === key ? 'bg-cyan-600 text-white' : 'text-slate-400 hover:text-slate-200'}"
                onclick={() => tsneMode = key}
              >{label}</button>
            {/each}
          </div>
        {/if}
        {#if showBins}
          <label class="flex items-center gap-2">
            Color by:
            <select
              class="bg-slate-800 border border-slate-600 text-slate-200 text-xs rounded px-2 py-1"
              bind:value={binColorBy}
            >
              {#each binColorDefs as { key, label }}
                <option value={key}>{label}{key === 'phylum' ? ` (${taxSourceLabel()})` : ''}</option>
              {/each}
            </select>
          </label>
        {:else if sampleTsne && !sampleTsne.tooFew && tsneColorOptions.length}
          <label class="flex items-center gap-2">
            Color by:
            <select
              class="bg-slate-800 border border-slate-600 text-slate-200 text-xs rounded px-2 py-1"
              bind:value={tsneColorBy}
            >
              <option value="">% in shown bins</option>
              {#each tsneColorOptions as c}
                <option value={c}>{c}</option>
              {/each}
            </select>
          </label>
        {/if}
      </div>
      {#if showBins}
        {#if binTsne?.n}
          <PlotlyChart
            data={binTsne.traces}
            onclick={handleBinTsneClick}
            layout={{
              title: { text: `Bin t-SNE (${binTsne.n} bins, Bray-Curtis on \u2074\u221A relabund)`, font: { size: 14 } },
              xaxis: { showticklabels: false, zeroline: false, title: 't-SNE 1' },
              yaxis: { showticklabels: false, zeroline: false, title: 't-SNE 2' },
              showlegend: binTsne.showlegend,
              legend: { font: { size: 10 }, itemsizing: 'constant' },
              hovermode: 'closest',
            }}
          />
        {:else}
          <p class="text-slate-500 text-sm py-12 text-center">No shown bins have reads in any sample</p>
        {/if}
      {:else if !sampleTsne || sampleTsne.tooFew}
        <h3 class="text-sm font-medium text-slate-400 mb-2">Sample t-SNE</h3>
        <p class="text-slate-500 text-sm py-12 text-center">
          Needs at least 4 samples with reads in the shown bins ({sampleTsne?.n ?? 0} now)
        </p>
      {:else}
      <PlotlyChart
        data={sampleTsne.traces}
        layout={{
          title: { text: `Sample t-SNE (${sampleTsne.n} samples, Bray-Curtis on \u2074\u221A relabund)`, font: { size: 14 } },
          xaxis: { showticklabels: false, zeroline: false, title: 't-SNE 1' },
          yaxis: { showticklabels: false, zeroline: false, title: 't-SNE 2' },
          showlegend: sampleTsne.showlegend,
          legend: { font: { size: 10 } },
        }}
      />
      {/if}
    </div>
    {/if}
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
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700 {sampleTsne || hasBinTsne ? 'lg:col-span-2' : ''}">
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
  {:else}
  <div class="flex flex-col items-center justify-center py-24 text-center">
    <p class="text-slate-400 text-lg mb-2">No bins match the current filters</p>
    <p class="text-slate-500 text-sm">Try adjusting the binner selection or quality thresholds above, or run binning to generate bins.</p>
  </div>
  {/if}
{/if}
