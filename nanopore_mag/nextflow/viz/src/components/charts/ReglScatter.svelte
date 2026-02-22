<script>
  import { onMount } from 'svelte';
  import createScatterplot from 'regl-scatterplot';
  import * as d3 from 'd3';

  let { data = null, colorBy = 'bin', sizeBy = 'length', sizeScale = 1.0, mode = 'pca', colorMap = {}, sizeRange = null, onselect = null, onclick = null } = $props();

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

  const BIN_LABELS = { bin: 'DAS Tool', semibin_bin: 'SemiBin2', metabat_bin: 'MetaBAT2', maxbin_bin: 'MaxBin2', lorbin_bin: 'LorBin', comebin_bin: 'COMEBin' };
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
      const tax = c.kaiju_phylum || c.kraken2_phylum || c.rrna_phylum || '?';
      const isBin = colorBy === 'bin' || colorBy.endsWith('_bin');
      const binPart = isBin ? ` | ${BIN_LABELS[colorBy] || colorBy}: ${c[colorBy] || 'none'}` : '';
      tipText = `${c.id} | ${c.length.toLocaleString()} bp | depth: ${c.depth} | GC: ${c.gc ?? '?'}%${binPart} | ${tax}`;
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

  function draw() {
    if (!scatterplot || !data?.contigs?.length) return;

    const contigs = data.contigs;
    let xKey = 'pca_x', yKey = 'pca_y';
    if (mode === 'tsne' && data.has_tsne) { xKey = 'tsne_x'; yKey = 'tsne_y'; }
    if (mode === 'umap' && data.has_umap) { xKey = 'umap_x'; yKey = 'umap_y'; }

    const isCont = !!CONTINUOUS[colorBy];
    const isBin = colorBy === 'bin' || colorBy.endsWith('_bin');
    const bgLabel = isBin ? 'unbinned' : 'Unknown';

    // Sort: background points first (drawn underneath)
    let sorted;
    if (!isCont) {
      sorted = [...contigs];
      sorted.sort((a, b) => {
        const aIsBg = !(a[colorBy]) || a[colorBy] === bgLabel;
        const bIsBg = !(b[colorBy]) || b[colorBy] === bgLabel;
        if (aIsBg && !bIsBg) return -1;
        if (!aIsBg && bIsBg) return 1;
        return 0;
      });
    } else {
      sorted = contigs;
    }
    indexMap = sorted;

    const n = sorted.length;
    const x = new Float32Array(n);
    const y = new Float32Array(n);
    const valueA = new Float32Array(n);
    const valueB = new Float32Array(n);

    // Extract and normalize coordinates to [-1, 1] with uniform scale (1:1 aspect)
    for (let i = 0; i < n; i++) {
      x[i] = sorted[i][xKey] || 0;
      y[i] = sorted[i][yKey] || 0;
    }
    const xExt = d3.extent(x);
    const yExt = d3.extent(y);
    const xSpan = xExt[1] - xExt[0] || 1;
    const ySpan = yExt[1] - yExt[0] || 1;
    const maxSpan = Math.max(xSpan, ySpan);
    const xMid = (xExt[0] + xExt[1]) / 2;
    const yMid = (yExt[0] + yExt[1]) / 2;
    const half = maxSpan * 0.55;
    for (let i = 0; i < n; i++) {
      x[i] = (x[i] - xMid) / half;
      y[i] = (y[i] - yMid) / half;
    }

    // ---- Color encoding ----
    let colorCfg;
    if (isCont) {
      const fn = CONTINUOUS[colorBy];
      const raw = new Float64Array(n);
      for (let i = 0; i < n; i++) raw[i] = fn(sorted[i]);
      const [vMin, vMax] = d3.extent(raw);
      const vRange = vMax - vMin || 1;
      for (let i = 0; i < n; i++) valueA[i] = (raw[i] - vMin) / vRange;
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

    // ---- Size encoding (use full-dataset extents from sizeRange to stay stable across filters) ----
    let sizeCfg;
    if (sizeBy === 'fixed') {
      sizeCfg = { sizeBy: null, pointSize: 3 * sizeScale };
    } else {
      // Match Plotly's absolute sizing: sqrt(depth)*0.8 or length^0.3*0.25, clamped at 1.5 min
      const sfn = sizeBy === 'depth'
        ? c => Math.max(1.5, Math.sqrt(c.depth) * 0.8)
        : c => Math.max(1.5, Math.pow(c.length, 0.3) * 0.25);
      // Normalize to [0,1] for regl but using absolute sizes as the basis
      let sMin = Infinity, sMax = -Infinity;
      for (let i = 0; i < n; i++) {
        const s = sfn(sorted[i]);
        if (s < sMin) sMin = s;
        if (s > sMax) sMax = s;
      }
      const sRange = sMax - sMin || 1;
      for (let i = 0; i < n; i++) {
        valueB[i] = Math.max(0, Math.min(1, (sfn(sorted[i]) - sMin) / sRange));
      }
      // Scale pointSize range to match the absolute pixel sizes from Plotly
      sizeCfg = { sizeBy: 'valueB', pointSize: [sMin * sizeScale, sMax * sizeScale] };
    }

    scatterplot.set({ ...colorCfg, ...sizeCfg });

    scatterplot.draw({ x, y, valueA, valueB }, {
      zDataType: isCont ? 'continuous' : 'categorical',
      wDataType: 'continuous',
    });
  }

  onMount(() => {
    return () => {
      if (scatterplot) { scatterplot.destroy(); scatterplot = null; }
    };
  });

  $effect(() => {
    const _deps = [data, colorBy, sizeBy, sizeScale, mode, colorMap];
    if (!canvasEl || !data) return;
    if (!scatterplot) init();
    draw();
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
