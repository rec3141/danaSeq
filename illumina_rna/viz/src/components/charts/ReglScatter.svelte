<script>
  import { onMount } from 'svelte';
  import createScatterplot from 'regl-scatterplot';
  import * as d3 from 'd3';

  let { data = null, colorBy = 'bin', sizeBy = 'length', sizeScale = 1.0, mode = 'pca', colorMap = {}, sizeRange = null, coordExtents = null, onselect = null, onclick = null, searchMatchIds = null, sampleDepthData = null, selectedSample = '', tipFormat = null, sizeExponent = 1,
        sizeOf = null, sizeMax = null } = $props();

  let canvasEl;
  let scatterplot;
  let indexMap = [];   // regl point index → contig object

  // Tooltip state
  let tipText = $state('');
  let tipX = $state(0);
  let tipY = $state(0);
  let tipShow = $state(false);
  let mouseX = 0, mouseY = 0;
  let hoveredIdx = -1;  // track hovered point for click
  let prevMode = null;  // track embedding mode to detect coordinate changes

  const BIN_LABELS = { bin: 'DAS Tool', semibin_bin: 'SemiBin2', metabat_bin: 'MetaBAT2', maxbin_bin: 'MaxBin2', lorbin_bin: 'LorBin', comebin_bin: 'COMEBin', vamb_bin: 'VAMB', vamb_tax_bin: 'VAMB-tax', binette_bin: 'Binette', magscot_bin: 'MAGScoT' };
  const BG = [15/255, 23/255, 42/255, 1]; // #0f172a
  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8',
                   '#94a3b8','#d4d4d8','#78716c'];
  const VIRIDIS = Array.from({length: 64}, (_, i) =>
    d3.color(d3.interpolateViridis(i / 63)).formatHex()
  );

  const CONTINUOUS = {
    depth:  c => Math.log10(c.depth + 0.01),
    length: c => Math.log10(c.length + 1),
    gc:     c => c.gc ?? 50,
  };

  function init() {
    if (!canvasEl || scatterplot) return;

    scatterplot = createScatterplot({
      canvas: canvasEl,
      width: 'auto',
      height: 'auto',
      backgroundColor: BG,
      pointSize: 4,
      deselectOnDblClick: true,
      deselectOnEscape: true,
      // Above a million points the library builds its spatial index in a worker
      // made from Function.toString(), which the minified bundle breaks (the
      // worker dies with a SyntaxError); build it on the main thread instead.
      spatialIndexUseWorker: false,
      opacityInactiveMax: 0.3,
      opacityInactiveScale: 0.5,
    });

    // Track mouse for tooltip positioning
    canvasEl.addEventListener('mousemove', (e) => {
      const rect = canvasEl.getBoundingClientRect();
      mouseX = e.clientX - rect.left;
      mouseY = e.clientY - rect.top;
    });

    scatterplot.subscribe('select', ({ points: indices }) => {
      if (onselect && indices.length) {
        onselect(indices.map(i => indexMap[i]?.id).filter(Boolean));
      }
    });

    scatterplot.subscribe('deselect', () => {
      if (onselect) onselect(null);
    });

    scatterplot.subscribe('pointOver', (idx) => {
      hoveredIdx = idx;
      const c = indexMap[idx];
      if (!c) return;
      // Pages whose points are not DNA contigs (e.g. the RNA dashboard) pass
      // their own tooltip text.
      if (tipFormat) {
        tipText = tipFormat(c);
        tipX = mouseX; tipY = mouseY - 16; tipShow = true;
        return;
      }
      const tax = c.kaiju_phylum || c.kraken2_phylum || c.sendsketch_phylum || c.rrna_phylum || '?';
      const isBin = colorBy === 'bin' || colorBy.endsWith('_bin');
      const binPart = isBin ? ` | ${BIN_LABELS[colorBy] || colorBy}: ${c[colorBy] || 'none'}` : '';
      let depthPart = `depth: ${c.depth}x`;
      if (colorBy === 'sample_depth' && sampleDepthData?.depthOf && selectedSample) {
        const sIdx = sampleDepthData.samples.indexOf(selectedSample);
        const sDepth = sIdx >= 0 ? sampleDepthData.depthOf(c, sIdx) : 0;
        depthPart = `${selectedSample}: ${sDepth.toFixed(2)}x | total: ${c.depth}x`;
      }
      tipText = `${c.id} | ${c.length.toLocaleString()} bp | ${depthPart} | GC: ${c.gc ?? '?'}%${binPart} | ${tax}`;
      tipX = mouseX;
      tipY = mouseY - 16;
      tipShow = true;
    });

    scatterplot.subscribe('pointOut', () => {
      hoveredIdx = -1;
      tipShow = false;
    });

    // Click on a hovered point → open detail panel
    canvasEl.addEventListener('click', () => {
      if (onclick && hoveredIdx >= 0) {
        const c = indexMap[hoveredIdx];
        if (c) onclick(c.id);
      }
    });
  }

  function applySearchFilter() {
    if (!scatterplot || !indexMap.length) return;
    if (searchMatchIds) {
      const indices = [];
      for (let i = 0; i < indexMap.length; i++) {
        if (searchMatchIds.has(indexMap[i].id)) indices.push(i);
      }
      scatterplot.filter(indices, { preventEvent: true });
    } else {
      scatterplot.unfilter({ preventEvent: true });
    }
  }

  // Drawing is split by what a control changes, so dragging a slider does
  // only its own work: colour/data/mode → sort + colours + upload; size field
  // or exponent → sizes + upload; size scale/max → point-size lookup only.
  let prep = null;      // { n, x, y, valueA, colorCfg, isCont }
  let valueB = null;

  function prepare() {
    const contigs = data.contigs;
    let xKey = 'pca_x', yKey = 'pca_y';
    if (mode === 'tsne' && data.has_tsne) { xKey = 'tsne_x'; yKey = 'tsne_y'; }
    if (mode === 'umap' && data.has_umap) { xKey = 'umap_x'; yKey = 'umap_y'; }

    // Build effective continuous map (add sample_depth dynamically)
    let effectiveContinuous = CONTINUOUS;
    if (colorBy === 'sample_depth' && sampleDepthData?.depthOf && selectedSample) {
      const sIdx = sampleDepthData.samples.indexOf(selectedSample);
      if (sIdx >= 0) {
        effectiveContinuous = {
          ...CONTINUOUS,
          sample_depth: c => Math.log10(sampleDepthData.depthOf(c, sIdx) + 0.01),
        };
      }
    }
    const isCont = !!effectiveContinuous[colorBy];
    const isBin = colorBy === 'bin' || colorBy.endsWith('_bin');
    const bgLabel = isBin ? 'unbinned' : 'Unknown';
    const n0 = contigs.length;

    // Draw order: low values / background first, high values / foreground on top.
    // Keys are computed once into a typed array; the sort compares numbers only.
    let sorted, raw = null;
    if (isCont) {
      const fn = effectiveContinuous[colorBy];
      const keys = new Float64Array(n0);
      for (let i = 0; i < n0; i++) keys[i] = fn(contigs[i]);
      const order = new Uint32Array(n0);
      for (let i = 0; i < n0; i++) order[i] = i;
      order.sort((a, b) => keys[a] - keys[b]);
      sorted = new Array(n0);
      raw = new Float64Array(n0);
      for (let i = 0; i < n0; i++) { sorted[i] = contigs[order[i]]; raw[i] = keys[order[i]]; }
    } else {
      // stable partition: background first, then the rest in original order
      sorted = [];
      const fg = [];
      for (const c of contigs) ((!c[colorBy] || c[colorBy] === bgLabel) ? sorted : fg).push(c);
      for (const c of fg) sorted.push(c);
    }
    indexMap = sorted;

    const n = sorted.length;
    const x = new Float32Array(n);
    const y = new Float32Array(n);
    const valueA = new Float32Array(n);

    // Extract and normalize coordinates to [-1, 1] with uniform scale (1:1 aspect)
    // Use full-dataset extents (from coordExtents) so filtering doesn't re-zoom
    for (let i = 0; i < n; i++) {
      x[i] = sorted[i][xKey] || 0;
      y[i] = sorted[i][yKey] || 0;
    }
    const modeKey = mode === 'tsne' && data.has_tsne ? 'tsne' : mode === 'umap' && data.has_umap ? 'umap' : 'pca';
    const ext = coordExtents?.[modeKey];
    const xMin = ext ? ext.xMin : d3.min(x);
    const xMax = ext ? ext.xMax : d3.max(x);
    const yMin = ext ? ext.yMin : d3.min(y);
    const yMax = ext ? ext.yMax : d3.max(y);
    const xSpan = xMax - xMin || 1;
    const ySpan = yMax - yMin || 1;
    const maxSpan = Math.max(xSpan, ySpan);
    const xMid = (xMin + xMax) / 2;
    const yMid = (yMin + yMax) / 2;
    const half = maxSpan * 0.55;
    for (let i = 0; i < n; i++) {
      x[i] = (x[i] - xMid) / half;
      y[i] = (y[i] - yMid) / half;
    }

    // ---- Color encoding ----
    let colorCfg;
    if (isCont) {
      // Use global fixed range for sample_depth (consistent across samples),
      // rank-based normalization for other continuous modes (spreads colors
      // evenly by percentile so narrow distributions show full color range).
      // Points are already sorted by value, so a point's rank is its position.
      if (colorBy === 'sample_depth' && sampleDepthData?.maxDepth != null) {
        const vMin = Math.log10(0.01);
        const vMax = Math.log10(sampleDepthData.maxDepth + 0.01);
        const vRange = vMax - vMin || 1;
        for (let i = 0; i < n; i++) valueA[i] = Math.max(0, Math.min(1, (raw[i] - vMin) / vRange));
      } else {
        for (let i = 0; i < n; i++) valueA[i] = i / (n - 1 || 1);
      }
      colorCfg = { colorBy: 'valueA', pointColor: VIRIDIS, opacityBy: null, opacity: 0.75 };
    } else {
      // Use stable colorMap from parent (derived from full unfiltered dataset)
      // so colors don't shift when filters change
      const useMap = colorMap && Object.keys(colorMap).length > 0;
      const catMap = new Map();
      catMap.set(bgLabel, 0);
      if (useMap) {
        // Build index from the stable colorMap keys (preserves full-dataset ordering)
        const sortedNames = Object.keys(colorMap).filter(k => k !== bgLabel).sort();
        for (const nm of sortedNames) catMap.set(nm, catMap.size);
      } else {
        const names = new Set();
        for (const c of sorted) names.add(c[colorBy] || bgLabel);
        const sortedNames = [...names].filter(k => k !== bgLabel).sort();
        for (const nm of sortedNames) catMap.set(nm, catMap.size);
      }

      for (let i = 0; i < n; i++) {
        const key = sorted[i][colorBy] || bgLabel;
        valueA[i] = catMap.has(key) ? catMap.get(key) : 0;
      }

      const colors = [colorMap[bgLabel] || '#475569'];
      const opacities = [0.25];
      const allNames = [...catMap.keys()].filter(k => k !== bgLabel);
      for (const nm of allNames) {
        colors.push(useMap ? colorMap[nm] : PALETTE[(catMap.get(nm) - 1) % PALETTE.length]);
        opacities.push(0.75);
      }
      colorCfg = { colorBy: 'valueA', pointColor: colors, opacityBy: 'valueA', opacity: opacities };
    }
    prep = { n, x, y, valueA, colorCfg, isCont };
  }

  // ---- Size encoding: px = sizeScale + ((1 + value)^sizeExponent − 1) ----
  // sizeScale is the size of a zero-value point (the minimum). 1 + value ≥ 1,
  // so points only grow as the value or the exponent rises. Raw value, no
  // normalization: filtering never resizes a point. Units: length in kb,
  // depth/expression as given (RPKM); sizeOf overrides the value.
  const MAX_PX = 80;
  function computeSizes() {
    const { n } = prep;
    valueB = new Float32Array(n);
    if (sizeBy === 'fixed') return;
    const sorted = indexMap;
    const get = sizeOf ?? (sizeBy === 'depth' ? (c => c.depth) : (c => c.length / 1e3));
    for (let i = 0; i < n; i++) {
      const px = sizeScale + Math.pow(1 + Math.max(0, get(sorted[i])), sizeExponent) - 1;
      valueB[i] = Math.min(1, px / MAX_PX);
    }
  }

  function sizeConfig() {
    if (sizeBy === 'fixed') return { sizeBy: null, pointSize: sizeScale };
    // valueB = px / MAX_PX; a fine lookup so sizes are near-continuous
    return { sizeBy: 'valueB', pointSize: Array.from({length: 256}, (_, i) => Math.max(1, (i / 255) * MAX_PX)) };
  }

  function upload() {
    // Save camera view before redraw (preserve zoom/pan across settings changes)
    const modeChanged = mode !== prevMode;
    const savedView = (!modeChanged && prevMode !== null) ? scatterplot.get('cameraView') : null;
    prevMode = mode;

    scatterplot.set({ ...prep.colorCfg, ...sizeConfig() });

    scatterplot.draw({ x: prep.x, y: prep.y, valueA: prep.valueA, valueB }, {
      zDataType: prep.isCont ? 'continuous' : 'categorical',
      wDataType: 'continuous',
    }).then(() => {
      // Restore camera after redraw (only if coordinates didn't change)
      if (savedView) {
        scatterplot.view(savedView, { preventEvent: true });
      }
      applySearchFilter();
    });
  }

  let lastColorKey = [], lastSizeKey = [];
  let sizeFrame = 0;
  const same = (a, b) => a.length === b.length && a.every((v, i) => v === b[i]);

  onMount(() => {
    return () => {
      if (scatterplot) { scatterplot.destroy(); scatterplot = null; }
    };
  });

  $effect(() => {
    const colorKey = [data, colorBy, mode, colorMap, coordExtents, sampleDepthData, selectedSample];
    const sizeKey = [sizeBy, sizeExponent, sizeScale, sizeOf];
    if (!canvasEl || !data?.contigs?.length) return;
    if (!scatterplot) init();
    if (!prep || !same(colorKey, lastColorKey)) {
      prepare(); computeSizes(); upload();
    } else if (!same(sizeKey, lastSizeKey)) {
      // sliders: at most one resize + upload per frame while dragging
      if (!sizeFrame) sizeFrame = requestAnimationFrame(() => { sizeFrame = 0; if (scatterplot) { computeSizes(); upload(); } });
    }
    lastColorKey = colorKey; lastSizeKey = sizeKey;
  });

  // Re-apply search filter when only searchMatchIds changes (no full redraw needed)
  $effect(() => {
    const _deps = [searchMatchIds];
    applySearchFilter();
  });
</script>

<div class="w-full relative flex-1 min-h-0">
  <canvas bind:this={canvasEl} class="absolute inset-0 w-full h-full rounded" style="background: #0f172a;"></canvas>
  {#if tipShow}
    <div
      class="absolute pointer-events-none bg-slate-900/95 text-slate-200 text-xs px-2 py-1 rounded shadow-lg border border-slate-600 whitespace-nowrap z-10"
      style="left: {tipX}px; top: {tipY}px; transform: translate(-50%, -100%)"
    >
      {tipText}
    </div>
  {/if}
</div>
