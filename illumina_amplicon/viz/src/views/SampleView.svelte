<script>
  import { onMount, untrack } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import {
    store, countsBySample,
    GROUP_HEX,
    buildTaxColorMap, getAsvColor, getEffectiveColorLevel, getClusterColor,
    taxaMatchingFilter, makeTaxonMatcher, scaleMarkerSizes, maxUsefulScale,
    sampleInAssays, assayKey, listAssays, assayHeading, assaySpan, assayPlace,
    lineageOf,
    markerRadius,
  } from '../stores/data.svelte.js';
  import { chrome } from '../stores/theme.svelte.js';

  let { filters = {} } = $props();

  let plotDiv;
  let hasPlot = false;
  let lastPoints = null;   // point sizes for the scale-only restyle path

  // The table recomputed these for EVERY row. They depend only on the filters, so
  // derive them once — buildTaxColorMap walks every ASV assignment per call.
  let tableColorLevel = $derived(
    filters.colorMode === 'group'
      ? 'group'
      : getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter, filters.taxonFilterLevel)
  );
  let tableColorMap = $derived.by(() =>
    tableColorLevel !== 'group'
      ? buildTaxColorMap(tableColorLevel, filters.taxonFilter, filters.taxonFilterLevel).colorMap
      : null
  );

  // ── Derived data ──────────────────────────────────────────────────────────

  let filteredSamples = $derived.by(() => {
    let s = store.samples.filter(s =>
      (s.total_reads ?? 0) >= (filters.minReads || 0) && sampleInAssays(s, filters.assays)
    );
    if (filters.sampleFilter) {
      // Substring, matching the taxonomy box. A sample id is not a pattern
      // language, and an invalid one used to silently match everything.
      const q = filters.sampleFilter.toLowerCase();
      const matched = s.filter(sample => (sample.id ?? '').toLowerCase().includes(q));
      if (matched.length > 0) s = matched;
    }
    return s;
  });

  let cMap = $derived.by(() => countsBySample());

  let taxonMatches = $derived.by(() =>
    makeTaxonMatcher(filters.taxonFilter, filters.taxonFilterLevel)
  );

  let topTaxa = $derived.by(() => {
    const matches = taxonMatches;
    const gf = filters.groupFlags || {};

    // Merge counts: single click > lasso > all filtered
    const merged = new Map();
    let sampleIds;
    if (store.selectedSample != null) {
      sampleIds = [store.samples[store.selectedSample]?.id];
    } else if (filters.lassoSampleIds.size > 0) {
      sampleIds = [...filters.lassoSampleIds];
    } else {
      sampleIds = filteredSamples.map(s => s.id);
    }

    for (const sampleId of sampleIds) {
      for (const e of (cMap.get(sampleId) ?? [])) {
        const prev = merged.get(e.asv_idx) || 0;
        merged.set(e.asv_idx, prev + e.count);
      }
    }

    const total = [...merged.values()].reduce((s, c) => s + c, 0) || 1;

    return [...merged.entries()]
      .map(([asv_idx, count]) => ({
        asv: store.asvs[asv_idx],
        count,
        pct: ((count / total) * 100).toFixed(1),
      }))
      .filter(e => {
        if (!e.asv) return false;
        const group = e.asv.group ?? 'prokaryote';
        if (gf[group] === false) return false;
        if (matches && !matches(e.asv.id ?? '', lineageOf(e.asv.id, e.asv.taxonomy))) return false;
        return true;
      })
      .sort((a, b) => b.count - a.count)
      .slice(0, 50);
  });

  let selectedSampleObj = $derived(
    store.selectedSample != null ? store.samples[store.selectedSample] : null
  );

  let totalAssayCount = $derived.by(() => listAssays().length);

  // ── Build plotly traces ───────────────────────────────────────────────────

  // Point scale alone does not change the data, only marker size — restyle instead of
  // rebuilding every trace.
  $effect(() => {
    const scale = filters.pointScale ?? 1;
    if (!plotDiv || !hasPlot || !lastPoints) return;
    Plotly.restyle(plotDiv, { 'marker.size': [scaleMarkerSizes(lastPoints, scale)] }, [0]);
  });

  $effect(() => {
    if (!plotDiv) return;
    // Reading the palette here is what redraws the plot when the theme flips.
    const c = chrome();
    if (filteredSamples.length === 0) {
      if (hasPlot) Plotly.react(plotDiv, [{ x: [], y: [], type: 'scattergl' }], {
        plot_bgcolor: c.plotBg, paper_bgcolor: c.plotBg,
        xaxis: { visible: false }, yaxis: { visible: false },
        annotations: [{ text: 'No samples match filters', showarrow: false,
          font: { color: c.faint, size: 14 }, xref: 'paper', yref: 'paper', x: 0.5, y: 0.5 }],
      });
      lastPoints = null;
      return;
    }

    const colorLevel = filters.colorMode === 'group' ? 'group' : getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter, filters.taxonFilterLevel);
    const cmap = colorLevel !== 'group' ? buildTaxColorMap(colorLevel, filters.taxonFilter, filters.taxonFilterLevel).colorMap : null;
    const matches = taxonMatches;
    const gf = filters.groupFlags || {};
    // untracked: a scale change is handled by the restyle effect below, which does not
    // rebuild the plot. Reading it reactively here would defeat that.
    const scale = untrack(() => filters.pointScale ?? 1);

    // Use pre-aggregated counts if available, otherwise fall back to on-the-fly
    const isAsvLevel = colorLevel === '_asv';
    const aggLevel = isAsvLevel ? null : colorLevel;
    const aggData = aggLevel && store.aggCounts?.[aggLevel];
    // Once drilled down, filters.taxonFilter is an ancestor of the aggregated
    // taxa, so it must be resolved to the matching taxa at this level rather
    // than tested against each child name (which matched nothing).
    const allowedTaxa = taxaMatchingFilter(aggLevel, filters.taxonFilter, filters.taxonFilterLevel);

    const aggNameQuery = filters.taxonFilter && !allowedTaxa
      ? filters.taxonFilter.toLowerCase()
      : null;

    const allPoints = [];

    if (aggData && !isAsvLevel && filters.colorMode !== 'cluster') {
      // Fast path: use pre-aggregated data
      const { data: aggEntries, samples: aggSamples, taxa: aggTaxa } = aggData;
      const sampleSet = new Set(filteredSamples.map(s => s.id));
      const sampleLookup = {};
      for (const s of filteredSamples) sampleLookup[s.id] = s;

      // Build per-sample totals from agg data
      const sampleTotals = {};
      for (const [si, ti, count, prop] of aggEntries) {
        const sid = aggSamples[si];
        if (!sampleSet.has(sid)) continue;
        sampleTotals[sid] = (sampleTotals[sid] || 0) + count;
      }

      // Build color map for taxa
      const taxonColors = {};
      if (cmap) {
        for (const taxon of aggTaxa) {
          taxonColors[taxon] = cmap[taxon] || '#475569';
        }
      }

      for (const [si, ti, count, prop] of aggEntries) {
        const sid = aggSamples[si];
        if (!sampleSet.has(sid)) continue;
        const sample = sampleLookup[sid];
        const taxon = aggTaxa[ti];

        // Apply taxonomy filter (lineage-aware, see allowedTaxa above). The
        // fallback only fires when the level could not be resolved, so all it
        // can do is test the aggregated taxon name itself.
        if (allowedTaxa) {
          if (!allowedTaxa.has(taxon)) continue;
        } else if (aggNameQuery && !String(taxon).toLowerCase().includes(aggNameQuery)) continue;

        const proportion = prop || (count / (sampleTotals[sid] || 1));
        const color = cmap ? (taxonColors[taxon] || '#475569') : (GROUP_HEX[taxon] || '#475569');

        allPoints.push({
          x: sample.x,
          y: sample.y,
          size: markerRadius(proportion, 17),
          color,
          proportion,
          text: `${sid}<br>${taxon}: ${(proportion * 100).toFixed(1)}%<br>${(sample.total_reads ?? 0).toLocaleString()} reads`,
        });
      }
    } else {
      // Slow path: aggregate on the fly (for ASV level, cluster mode, or missing agg data)
      for (const sample of filteredSamples) {
        const entries = cMap.get(sample.id) ?? [];
        const totalCount = entries.reduce((s, e) => s + e.count, 0) || 1;

        for (const { asv_idx, count } of entries) {
          const asv = store.asvs[asv_idx];
          if (!asv) continue;

          const group = asv.group ?? 'prokaryote';
          if (gf[group] === false) continue;
          if ((asv.n_samples ?? 0) < (filters.minPrevalence || 0)) continue;
          if (matches && !matches(asv.id ?? '', lineageOf(asv.id, asv.taxonomy))) continue;

          const proportion = count / totalCount;
          let color;
          if (filters.colorMode === 'cluster') {
            color = getClusterColor(sample.id, 'sampleCluster', filters.sampleClusterK);
          } else if (cmap) {
            color = getAsvColor(asv.id, colorLevel, cmap);
          } else {
            color = GROUP_HEX[group] ?? GROUP_HEX.unknown;
          }

          allPoints.push({
            x: sample.x,
            y: sample.y,
            size: markerRadius(proportion, 17),
            color,
            proportion,
            text: `${sample.id}<br>${asv.id}: ${(proportion * 100).toFixed(1)}%<br>${(sample.total_reads ?? 0).toLocaleString()} reads`,
          });
        }
      }
    }

    // Sort largest first
    allPoints.sort((a, b) => b.proportion - a.proportion);
    const totalPoints = allPoints.length;

    // Base sizes, unscaled. Both this path and the restyle effect run them through
    // scaleMarkerSizes so the two cannot drift.
    lastPoints = allPoints.map(p => p.size * 0.3);
    store.maxPointScale = maxUsefulScale(lastPoints);

    const overlayTraces = [{
      x: allPoints.map(p => p.x),
      y: allPoints.map(p => p.y),
      mode: 'markers',
      type: 'scattergl',
      marker: {
        size: scaleMarkerSizes(lastPoints, scale),
        color: allPoints.map(p => p.color),
        opacity: 0.7,
        line: { width: 0 },
      },
      text: allPoints.map(p => p.text),
      hoverinfo: 'text',
      showlegend: false,
    }];

    const layout = {
      dragmode: 'pan',
      xaxis: { title: '', zeroline: false, showgrid: false, showticklabels: false },
      yaxis: { title: '', zeroline: false, showgrid: false, showticklabels: false, scaleanchor: 'x' },
      plot_bgcolor: c.plotBg,
      paper_bgcolor: c.plotBg,
      font: { color: c.axis },
      legend: {
        bgcolor: 'rgba(15, 23, 42, 0.8)',
        font: { size: 10 },
      },
      margin: { l: 20, r: 20, t: 10, b: 20 },
      title: { text: `${filteredSamples.length} samples, ${totalPoints.toLocaleString()} points (${isAsvLevel ? 'per ASV' : 'per ' + (colorLevel === 'group' ? 'group' : colorLevel)})`,
               font: { size: 12, color: c.faint }, x: 0.01, y: 0.99 },
    };

    const config = { scrollZoom: true, displayModeBar: false, doubleClick: 'reset+autosize' };

    if (!hasPlot) {
      Plotly.newPlot(plotDiv, overlayTraces, layout, config);
      hasPlot = true;

      plotDiv.on('plotly_click', (data) => {
        if (data.points?.[0]) {
          const pt = data.points[0];
          let bestIdx = -1, bestDist = Infinity;
          filteredSamples.forEach((s, i) => {
            const d = (s.x - pt.x) ** 2 + (s.y - pt.y) ** 2;
            if (d < bestDist) { bestDist = d; bestIdx = i; }
          });
          const sId = bestIdx >= 0 ? filteredSamples[bestIdx]?.id : null;
          const sIdx = sId ? store.samples.findIndex(s => s.id === sId) : -1;
          store.selectedSample = sIdx >= 0 ? sIdx : null;
          filters.lassoSampleIds = new Set();
        }
      });

      plotDiv.on('plotly_selected', (data) => {
        if (data?.points?.length > 0) {
          const ids = new Set();
          for (const pt of data.points) {
            let bestIdx = -1, bestDist = Infinity;
            filteredSamples.forEach((s, i) => {
              const d = (s.x - pt.x) ** 2 + (s.y - pt.y) ** 2;
              if (d < bestDist) { bestDist = d; bestIdx = i; }
            });
            if (bestIdx >= 0) ids.add(filteredSamples[bestIdx].id);
          }
          filters.lassoSampleIds = ids;
          store.selectedSample = null;
        }
      });

      plotDiv.on('plotly_deselect', () => {
        filters.lassoSampleIds = new Set();
      });

    } else {
      // Preserve current zoom on data updates
      const curLayout = plotDiv.layout;
      if (curLayout?.xaxis?.range) layout.xaxis.range = curLayout.xaxis.range;
      if (curLayout?.yaxis?.range) layout.yaxis.range = curLayout.yaxis.range;
      Plotly.react(plotDiv, overlayTraces, layout, config);
    }
  });

  function handleKey(e) {
    if (!plotDiv || !hasPlot) return;
    if (e.key === 'Shift') {
      Plotly.relayout(plotDiv, { dragmode: e.type === 'keydown' ? 'lasso' : 'pan' });
    } else if (e.key === 'Escape' && e.type === 'keydown') {
      store.selectedSample = null;
      filters.lassoSampleIds = new Set();
      Plotly.restyle(plotDiv, { selectedpoints: [null] });
    }
  }

  onMount(() => {
    document.addEventListener('keydown', handleKey);
    document.addEventListener('keyup', handleKey);
    // The table below appears and changes height as the selection changes, which
    // shrinks this frame under a canvas that was sized for the old one.
    const ro = new ResizeObserver(() => {
      if (plotDiv && hasPlot) Plotly.Plots.resize(plotDiv);
    });
    if (plotDiv?.parentElement) ro.observe(plotDiv.parentElement);
    return () => {
      ro.disconnect();
      document.removeEventListener('keydown', handleKey);
      document.removeEventListener('keyup', handleKey);
      if (plotDiv && hasPlot) Plotly.purge(plotDiv);
    };
  });
</script>

<div class="flex h-full flex-col">
  <!-- The plot is absolutely positioned, so it paints above the static pane
       below it: any moment its canvas is taller than this frame, it covers the
       table. Plotly sizes the canvas when it draws and does not re-size it when
       the frame does, so the frame clips as well as the observer resizing. -->
  <div class="relative min-h-0 flex-1 overflow-hidden">
    <div bind:this={plotDiv} class="absolute inset-0"></div>
  </div>

  {#if topTaxa.length > 0}
    <!-- The plot above keeps flex-1; this pane takes what it needs up to its own
         cap, and min-h-0 lets it actually shrink when the window is short. -->
    <div class="flex min-h-0 shrink-0 flex-col border-t border-line bg-surface/80 p-4">
      <div class="mb-2 flex items-center justify-between">
        <h3 class="text-sm font-semibold text-fg">
          {#if selectedSampleObj}
            {selectedSampleObj.id} &mdash; {(selectedSampleObj.total_reads ?? 0).toLocaleString()} reads
          {:else if filters.lassoSampleIds.size > 0}
            {filters.lassoSampleIds.size} selected samples &mdash; top {topTaxa.length} ASVs
          {:else}
            All samples ({filteredSamples.length}) &mdash; top {topTaxa.length} ASVs
          {/if}
        </h3>
        {#if selectedSampleObj}
          <button
            class="text-xs text-faint hover:text-fg2"
            onclick={() => store.selectedSample = null}
          >Clear selection</button>
        {/if}
      </div>

      <div class="mb-2">
        {#if selectedSampleObj && assayKey(selectedSampleObj)}
          {@const s = selectedSampleObj}
          <p class="text-xs text-muted">
            <span class="text-faint">Assay:</span>
            {assayHeading(s.assay_gene, s.assay_region, s.assay_primer_fwd, s.assay_primer_rev, assayPlace(s), s.assay_conflict)}
            {#if assayPlace(s)}
              <span class="text-faint">·</span> {assaySpan(assayPlace(s))}
            {:else if s.assay_primer_fwd || s.assay_primer_rev}
              <span class="text-faint">·</span> {[s.assay_primer_fwd, s.assay_primer_rev].filter(Boolean).join(' / ')}
            {/if}
            {#if Number.isFinite(Number(s.assay_match_fraction))}
              <span class="text-faint">·</span> {(Number(s.assay_match_fraction) * 100).toFixed(1)}% of reads matched
            {/if}
            {#if s.assay_source}<span class="text-faint"> (via {s.assay_source})</span>{/if}
          </p>
        {:else if filters.assays?.size}
          <p class="text-xs text-faint">
            Assay filter active: {filters.assays.size} of {totalAssayCount} assay{totalAssayCount === 1 ? '' : 's'}
          </p>
        {/if}
      </div>

        <!-- Sized against the window rather than pinned at 12rem, which showed
             six rows on any screen and cut the seventh through the middle. Kept
             to a fifth of it: this table is the plot's companion, and the plot
             is what the tab is for. The header needs its own stacking context or
             a row scrolling beneath it draws over it instead of under. -->
        <div class="max-h-[20vh] min-h-0 overflow-y-auto">
          <table class="w-full text-xs">
            <!-- Sticky on the cells rather than on <thead>: Safari does not
                 make a table section a sticky container, so a header set there
                 scrolls away and the rows arrive under nothing. -->
            <thead class="text-left text-muted">
              <tr>
                <th class="sticky top-0 z-10 bg-surface py-1 pr-4 shadow-[0_1px_0_0] shadow-line">ASV</th>
                <th class="sticky top-0 z-10 bg-surface py-1 pr-4 shadow-[0_1px_0_0] shadow-line">Taxonomy</th>
                <th class="sticky top-0 z-10 bg-surface py-1 pr-4 text-right shadow-[0_1px_0_0] shadow-line">Reads</th>
                <th class="sticky top-0 z-10 bg-surface py-1 text-right shadow-[0_1px_0_0] shadow-line">%</th>
              </tr>
            </thead>
            <tbody class="text-fg2">
              {#each topTaxa as row}
                {@const colorLevel = tableColorLevel}
                {@const cmap = tableColorMap}
                {@const rowColor = cmap ? getAsvColor(row.asv.id, colorLevel, cmap) : (GROUP_HEX[row.asv.group] ?? GROUP_HEX.prokaryote)}
                <tr class="border-t border-line/50 hover:bg-raised/30">
                  <td class="py-1 pr-4 font-mono">
                    <span class="inline-block h-2.5 w-2.5 rounded-full mr-1.5" style="background:{rowColor}"></span>
                    {row.asv.id ?? ''}
                  </td>
                  <td class="py-1 pr-4">{lineageOf(row.asv.id, row.asv.taxonomy)}</td>
                  <td class="py-1 pr-4 text-right font-mono">{row.count.toLocaleString()}</td>
                  <td class="py-1 text-right font-mono">{row.pct}</td>
                </tr>
              {/each}
            </tbody>
          </table>
        </div>
    </div>
  {/if}
</div>
