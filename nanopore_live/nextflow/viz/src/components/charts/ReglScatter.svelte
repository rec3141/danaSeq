<script>
  import { onMount } from 'svelte';
  import createScatterplot from 'regl-scatterplot';
  import * as d3 from 'd3';

  let { data = null, colorBy = 'sample', sizeBy = 'fixed', sizeScale = 1.0, mode = 'tsne', colorMap = {}, coordExtents = null, onselect = null, onclick = null, searchMatchIds = null, idField = 'id', exportName = 'danaseq_scatter' } = $props();

  let canvasEl;
  let scatterplot;
  let indexMap = [];

  let tipText = $state('');
  let tipX = $state(0);
  let tipY = $state(0);
  let tipShow = $state(false);
  let mouseX = 0, mouseY = 0;
  let hoveredIdx = -1;
  let prevMode = null;

  const BG = [15/255, 23/255, 42/255, 1];
  const PALETTE = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8',
                   '#94a3b8','#d4d4d8','#78716c'];
  const VIRIDIS = Array.from({length: 64}, (_, i) =>
    d3.color(d3.interpolateViridis(i / 63)).formatHex()
  );

  const CONTINUOUS = {
    gc:     c => c.gc ?? 50,
    length: c => Math.log10((c.length || 1) + 1),
    read_count: c => Math.log10((c.read_count || 1) + 1),
    total_bases: c => Math.log10((c.total_bases || 1) + 1),
    diversity: c => c.diversity ?? 0,
  };

  function exportPng() {
    if (!canvasEl) return;
    // regl-scatterplot preserves drawing buffer, capture directly
    const url = canvasEl.toDataURL('image/png');
    const a = document.createElement('a');
    a.href = url; a.download = `${exportName}.png`; a.click();
  }

  function init() {
    if (!canvasEl || scatterplot) return;
    scatterplot = createScatterplot({
      canvas: canvasEl, width: 'auto', height: 'auto',
      backgroundColor: BG, pointSize: 4,
      deselectOnDblClick: true, deselectOnEscape: true,
      opacityInactiveMax: 0.3, opacityInactiveScale: 0.5,
    });

    canvasEl.addEventListener('mousemove', (e) => {
      const rect = canvasEl.getBoundingClientRect();
      mouseX = e.clientX - rect.left;
      mouseY = e.clientY - rect.top;
    });

    scatterplot.subscribe('select', ({ points: indices }) => {
      if (onselect && indices.length) onselect(indices.map(i => indexMap[i]?.[idField]).filter(Boolean));
    });
    scatterplot.subscribe('deselect', () => { if (onselect) onselect(null); });

    scatterplot.subscribe('pointOver', (idx) => {
      hoveredIdx = idx;
      const c = indexMap[idx];
      if (!c) return;
      const parts = [];
      if (c.sample) parts.push(c.sample);
      if (c.length != null) parts.push(`${c.length.toLocaleString()} bp`);
      if (c.read_count != null) parts.push(`${c.read_count.toLocaleString()} reads`);
      if (c.gc != null) parts.push(`GC: ${typeof c.gc === 'number' ? c.gc.toFixed(1) : c.gc}%`);
      if (c[colorBy] != null && !CONTINUOUS[colorBy]) parts.push(`${colorBy}: ${c[colorBy]}`);
      if (!parts.length) parts.push(c[idField] || c.id);
      tipText = parts.join(' | ');
      tipX = mouseX; tipY = mouseY - 16; tipShow = true;
    });
    scatterplot.subscribe('pointOut', () => { hoveredIdx = -1; tipShow = false; });

    canvasEl.addEventListener('click', () => {
      if (onclick && hoveredIdx >= 0) {
        const c = indexMap[hoveredIdx];
        if (c) onclick(c[idField] || c.id);
      }
    });
  }

  function applySearchFilter() {
    if (!scatterplot || !indexMap.length) return;
    if (searchMatchIds) {
      const indices = [];
      for (let i = 0; i < indexMap.length; i++) {
        if (searchMatchIds.has(indexMap[i][idField] || indexMap[i].id)) indices.push(i);
      }
      scatterplot.filter(indices, { preventEvent: true });
    } else {
      scatterplot.unfilter({ preventEvent: true });
    }
  }

  function draw() {
    if (!scatterplot || !data?.points?.length) return;

    const points = data.points;
    const xKey = `${mode}_x`, yKey = `${mode}_y`;

    const isCont = !!CONTINUOUS[colorBy];

    let sorted;
    if (isCont) {
      const fn = CONTINUOUS[colorBy];
      sorted = [...points].sort((a, b) => fn(a) - fn(b));
    } else {
      sorted = [...points].sort((a, b) => {
        const aIsBg = !(a[colorBy]) || a[colorBy] === 'Unknown';
        const bIsBg = !(b[colorBy]) || b[colorBy] === 'Unknown';
        if (aIsBg && !bIsBg) return -1;
        if (!aIsBg && bIsBg) return 1;
        return 0;
      });
    }
    indexMap = sorted;

    const n = sorted.length;
    const x = new Float32Array(n);
    const y = new Float32Array(n);
    const valueA = new Float32Array(n);
    const valueB = new Float32Array(n);

    for (let i = 0; i < n; i++) {
      x[i] = sorted[i][xKey] || 0;
      y[i] = sorted[i][yKey] || 0;
    }

    const ext = coordExtents;
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

    let colorCfg;
    if (isCont) {
      const fn = CONTINUOUS[colorBy];
      const raw = new Float64Array(n);
      for (let i = 0; i < n; i++) raw[i] = fn(sorted[i]);
      // Rank-based normalization: spread colors evenly by percentile
      // so narrow distributions (e.g. GC 40-55%) still show full color range
      const indices = Array.from({length: n}, (_, i) => i);
      indices.sort((a, b) => raw[a] - raw[b]);
      for (let rank = 0; rank < n; rank++) valueA[indices[rank]] = rank / (n - 1 || 1);
      colorCfg = { colorBy: 'valueA', pointColor: VIRIDIS, opacityBy: null, opacity: 0.75 };
    } else {
      const useMap = colorMap && Object.keys(colorMap).length > 0;
      const catMap = new Map();
      catMap.set('Unknown', 0);
      if (useMap) {
        const sortedNames = Object.keys(colorMap).filter(k => k !== 'Unknown').sort();
        for (const nm of sortedNames) catMap.set(nm, catMap.size);
      } else {
        const names = new Set();
        for (const c of sorted) names.add(c[colorBy] || 'Unknown');
        const sortedNames = [...names].filter(k => k !== 'Unknown').sort();
        for (const nm of sortedNames) catMap.set(nm, catMap.size);
      }

      for (let i = 0; i < n; i++) {
        const key = sorted[i][colorBy] || 'Unknown';
        valueA[i] = catMap.has(key) ? catMap.get(key) : 0;
      }

      const colors = [colorMap['Unknown'] || '#475569'];
      const opacities = [0.25];
      const allNames = [...catMap.keys()].filter(k => k !== 'Unknown');
      for (const nm of allNames) {
        colors.push(useMap ? colorMap[nm] : PALETTE[(catMap.get(nm) - 1) % PALETTE.length]);
        opacities.push(0.75);
      }
      colorCfg = { colorBy: 'valueA', pointColor: colors, opacityBy: 'valueA', opacity: opacities };
    }

    let sizeCfg;
    if (sizeBy === 'fixed') {
      sizeCfg = { sizeBy: null, pointSize: 5 * sizeScale };
    } else {
      const sfn = CONTINUOUS[sizeBy] || (c => 1);
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
      const pMin = 2 * sizeScale, pMax = 12 * sizeScale;
      const sizeMap = Array.from({length: 64}, (_, i) => pMin + (i / 63) * (pMax - pMin));
      sizeCfg = { sizeBy: 'valueB', pointSize: sizeMap };
    }

    const modeChanged = mode !== prevMode;
    const savedView = (!modeChanged && prevMode !== null) ? scatterplot.get('cameraView') : null;
    prevMode = mode;

    scatterplot.set({ ...colorCfg, ...sizeCfg });

    scatterplot.draw({ x, y, valueA, valueB }, {
      zDataType: isCont ? 'continuous' : 'categorical',
      wDataType: 'continuous',
    }).then(() => {
      if (savedView) scatterplot.lookAt(savedView, { preventEvent: true });
      applySearchFilter();
    });
  }

  onMount(() => {
    return () => { if (scatterplot) { scatterplot.destroy(); scatterplot = null; } };
  });

  $effect(() => {
    const _deps = [data, colorBy, sizeBy, sizeScale, mode, colorMap, coordExtents];
    if (!canvasEl || !data) return;
    if (!scatterplot) init();
    draw();
  });

  $effect(() => {
    const _deps = [searchMatchIds];
    applySearchFilter();
  });
</script>

<div class="w-full relative flex-1 min-h-0">
  <canvas bind:this={canvasEl} class="absolute inset-0 w-full h-full rounded" style="background: #0f172a;"></canvas>
  <button
    class="absolute top-2 left-2 z-10 p-1 rounded text-slate-500 hover:text-cyan-400 hover:bg-slate-800/80 transition-colors"
    onclick={exportPng}
    title="Download PNG"
  >
    <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4 0l-4 4m0 0l-4-4m4 4V4"/>
    </svg>
  </button>
  {#if tipShow}
    <div
      class="absolute pointer-events-none bg-slate-900/95 text-slate-200 text-xs px-2 py-1 rounded shadow-lg border border-slate-600 whitespace-nowrap z-10"
      style="left: {tipX}px; top: {tipY}px; transform: translate(-50%, -100%)"
    >
      {tipText}
    </div>
  {/if}
</div>
