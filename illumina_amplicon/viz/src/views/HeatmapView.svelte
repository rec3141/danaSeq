<script>
  import { onMount } from 'svelte';
  import { store, getClusterColor, GROUP_HEX, buildTaxColorMap, getAsvColor, getEffectiveColorLevel, makeTaxonMatcher, sampleInAssays, activeDb} from '../stores/data.svelte.js';
  import { chrome } from '../stores/theme.svelte.js';

  let { filters = {} } = $props();

  let rawHeatmapData = null;  // full data, never reactive
  let heatmapData = $state(null);
  let visibleAsvCount = $state(0);
  let container;
  let gridEl;
  let canvas;
  let rowDendroSvg;
  let colDendroSvg;
  let colColorBar;
  let rowColorBar;
  let tooltip = $state({ show: false, x: 0, y: 0, text: '' });

  // Layout constants
  const ROW_DENDRO_W = 80;
  const COL_DENDRO_H = 120;
  const COLOR_BAR_W = 8;
  const COLOR_BAR_H = 8;
  let cellSize = $derived(filters.heatmapCellSize ?? 3);

  // Pre-compute max z per column once (expensive, do only at load)
  let colMaxZ = null;

  function computeColMax(data) {
    if (!data?.z?.length) return [];
    const nCols = data.z[0]?.length || 0;
    const maxes = new Float64Array(nCols);
    for (let i = 0; i < data.z.length; i++) {
      for (let j = 0; j < nCols; j++) {
        if (data.z[i][j] > maxes[j]) maxes[j] = data.z[i][j];
      }
    }
    return maxes;
  }

  // An assay filter cuts both axes: its samples are the rows, and the ASVs it
  // observed are the columns. Either dendrogram is dropped once its axis is
  // subset, since the stored linkage no longer describes what is drawn.
  function filterHeatmap(data, minMaxRA, assayAsvs, assaySamples) {
    if (!data?.z?.length) return data;
    const nCols = data.z[0]?.length || 0;
    const nRows = data.sampleIds?.length ?? data.z.length;

    if (!colMaxZ || colMaxZ.length !== nCols) {
      colMaxZ = computeColMax(data);
    }

    const zThreshold = Math.pow(minMaxRA / 100, 0.25);
    const keepCols = [];
    for (let j = 0; j < nCols; j++) {
      if (colMaxZ[j] < zThreshold) continue;
      if (assayAsvs && !assayAsvs.has(data.asvIds[j])) continue;
      keepCols.push(j);
    }

    const keepRows = [];
    for (let i = 0; i < nRows; i++) {
      if (assaySamples && !assaySamples.has(data.sampleIds[i])) continue;
      keepRows.push(i);
    }

    const colsCut = keepCols.length !== nCols;
    const rowsCut = keepRows.length !== nRows;
    if (!colsCut && !rowsCut) return data;

    return {
      ...data,
      z: keepRows.map(i => keepCols.map(j => data.z[i][j])),
      sampleIds: keepRows.map(i => data.sampleIds[i]),
      nSamples: keepRows.length,
      asvIds: keepCols.map(j => data.asvIds[j]),
      asvLabels: keepCols.map(j => data.asvLabels[j]),
      colDendro: colsCut ? { icoord: [], dcoord: [] } : data.colDendro,
      rowDendro: rowsCut ? { icoord: [], dcoord: [] } : data.rowDendro,
    };
  }

  let assayAsvs = $derived(store.assayAsvIds);
  let assaySamples = $derived.by(() => {
    if (!filters.assays?.size) return null;
    const ids = new Set();
    for (const s of store.samples) if (sampleInAssays(s, filters.assays)) ids.add(s.id);
    return ids;
  });

  onMount(async () => {
    try {
      let res = await fetch('./data/heatmap.json.gz');
      if (!res.ok) res = await fetch('./data/heatmap.json');
      if (res.ok) {
        rawHeatmapData = await res.json();
        const filtered = filterHeatmap(rawHeatmapData, filters.heatmapMinMaxRA ?? 1.0, assayAsvs, assaySamples);
        heatmapData = filtered;
        visibleAsvCount = filtered?.asvIds?.length ?? 0;
      }
    } catch (e) {
      console.error('Failed to load heatmap data:', e);
    }
  });

  // Viridis-ish colormap (256 entries)
  function viridis(t) {
    // Simplified viridis: dark purple -> teal -> yellow
    const r = Math.round(Math.max(0, Math.min(255, -510 * t * t + 713 * t + 68)));
    const g = Math.round(Math.max(0, Math.min(255, -270 * t * t + 442 * t + 24)));
    const b = Math.round(Math.max(0, Math.min(255, 320 * t * t - 570 * t + 280)));
    return `rgb(${r},${g},${b})`;
  }

  // Viridis sampled at 9 stops, linearly interpolated to 256
  const VIRIDIS_STOPS = [
    [68, 1, 84], [72, 35, 116], [64, 67, 135], [52, 94, 141],
    [33, 145, 140], [53, 183, 121], [109, 205, 89], [180, 222, 44], [253, 231, 37],
  ];

  function buildViridisLUT() {
    const lut = [];
    const nStops = VIRIDIS_STOPS.length;
    for (let i = 0; i < 256; i++) {
      const t = i / 255 * (nStops - 1);
      const lo = Math.floor(t);
      const hi = Math.min(lo + 1, nStops - 1);
      const f = t - lo;
      lut.push([
        Math.round(VIRIDIS_STOPS[lo][0] * (1 - f) + VIRIDIS_STOPS[hi][0] * f),
        Math.round(VIRIDIS_STOPS[lo][1] * (1 - f) + VIRIDIS_STOPS[hi][1] * f),
        Math.round(VIRIDIS_STOPS[lo][2] * (1 - f) + VIRIDIS_STOPS[hi][2] * f),
      ]);
    }
    return lut;
  }

  const viridisLUT = buildViridisLUT();

  // ── Newick parser → leaf order + dendrogram coords ──
  function parseNewick(s) {
    let pos = 0;
    function readNode() {
      const node = { children: [], label: '', branchLength: 0 };
      if (s[pos] === '(') {
        pos++;
        node.children.push(readNode());
        while (s[pos] === ',') { pos++; node.children.push(readNode()); }
        pos++; // skip )
      }
      let label = '';
      while (pos < s.length && s[pos] !== ':' && s[pos] !== ',' && s[pos] !== ')' && s[pos] !== ';') {
        label += s[pos++];
      }
      node.label = label;
      if (pos < s.length && s[pos] === ':') {
        pos++;
        let bl = '';
        while (pos < s.length && s[pos] !== ',' && s[pos] !== ')' && s[pos] !== ';') bl += s[pos++];
        node.branchLength = parseFloat(bl) || 0;
      }
      return node;
    }
    return readNode();
  }

  function getLeaves(node) {
    if (node.children.length === 0) return [node.label];
    return node.children.flatMap(getLeaves);
  }

  /**
   * Drop every leaf outside `keep` and collapse the internal nodes left with a
   * single child, folding their branch lengths into it.
   *
   * The heatmap columns are already filtered down by min max(RA), so drawing the
   * dendrogram from the unpruned tree gave one leaf per ASV in the whole dataset
   * above a handful of visible columns. Returns null if nothing is kept.
   */
  function pruneTree(node, keep) {
    if (node.children.length === 0) {
      return keep.has(node.label) ? { ...node } : null;
    }
    const kids = node.children.map(c => pruneTree(c, keep)).filter(Boolean);
    if (kids.length === 0) return null;
    if (kids.length === 1) {
      const only = kids[0];
      return { ...only, branchLength: only.branchLength + node.branchLength };
    }
    return { ...node, children: kids };
  }

  // Convert parsed tree to line segments for SVG rendering
  // Returns { lines: [{x1,y1,x2,y2}, ...], nLeaves }
  function treeToDendroLines(tree) {
    const lines = [];
    let leafIdx = 0;

    // First pass: compute max depth from root to leaves to set total tree height
    function maxDepth(node) {
      if (node.children.length === 0) return 0;
      return Math.max(...node.children.map(c => c.branchLength + maxDepth(c)));
    }

    const totalH = maxDepth(tree);

    // Layout: each node gets a position (leaf index) and absolute height from root
    function layout(node, depthFromRoot) {
      if (node.children.length === 0) {
        const pos = leafIdx;
        leafIdx++;
        return { pos, depth: depthFromRoot };
      }

      const childResults = node.children.map(c =>
        layout(c, depthFromRoot + c.branchLength)
      );

      const leftPos = childResults[0].pos;
      const rightPos = childResults[childResults.length - 1].pos;

      // Horizontal bar at this node's depth
      lines.push({ x1: leftPos, y1: depthFromRoot, x2: rightPos, y2: depthFromRoot });

      // Vertical drop from this node down to each child
      for (const r of childResults) {
        lines.push({ x1: r.pos, y1: depthFromRoot, x2: r.pos, y2: r.depth });
      }

      const centerPos = (leftPos + rightPos) / 2;
      return { pos: centerPos, depth: depthFromRoot };
    }

    layout(tree, 0);
    return { lines, nLeaves: leafIdx, maxDepth: totalH };
  }

  // Re-filter when slider changes (debounced — 65M cell scan is expensive)
  let filterTimer = null;
  $effect(() => {
    const threshold = filters.heatmapMinMaxRA;
    const aAsvs = assayAsvs, aSamples = assaySamples;
    if (rawHeatmapData) {
      if (filterTimer) clearTimeout(filterTimer);
      filterTimer = setTimeout(() => {
        const filtered = filterHeatmap(rawHeatmapData, threshold ?? 1.0, aAsvs, aSamples);
        heatmapData = filtered;
        visibleAsvCount = filtered?.asvIds?.length ?? 0;
      }, 300);
    }
  });

  let usePhyloOrder = $derived(filters.heatmapAsvTree === 'phylogeny' && !!store.treeNewick);

  $effect(() => {
    const _phyloOrder = usePhyloOrder;
    if (!heatmapData || !canvas || !rowDendroSvg || !colDendroSvg) return;

    let { z, sampleIds, asvIds, asvLabels, rowDendro, colDendro } = heatmapData;

    // If phylogeny ordering selected, reorder columns by NJ tree leaf order
    if (usePhyloOrder) {
      try {
        const tree = parseNewick(store.treeNewick);
        const phyloLeaves = getLeaves(tree);

        // Build mapping: asvId → column index in heatmap data
        const heatmapColIdx = {};
        asvIds.forEach((id, i) => { heatmapColIdx[id] = i; });

        // Filter to only ASVs present in the heatmap
        const newOrder = phyloLeaves.filter(id => id in heatmapColIdx);

        if (newOrder.length > 0) {
          const colMap = newOrder.map(id => heatmapColIdx[id]);

          // Reorder columns
          z = z.map(row => colMap.map(ci => row[ci]));
          asvIds = newOrder;
          asvLabels = colMap.map(ci => heatmapData.asvLabels[ci]);

          // Generate dendrogram lines from the phylogeny, pruned to the columns
          // actually drawn so its leaves line up with the heatmap.
          const pruned = pruneTree(tree, new Set(newOrder));
          if (pruned) colDendro = { _phyloLines: treeToDendroLines(pruned) };
        }
      } catch (e) {
        console.warn('Phylogeny ordering failed:', e);
      }
    }

    const nRows = z.length;
    const nCols = z[0]?.length || 0;
    if (nRows === 0 || nCols === 0) return;

    // Find max value for color scaling
    let zMax = 0;
    for (const row of z) for (const v of row) if (v > zMax) zMax = v;
    if (zMax === 0) zMax = 1;

    // Get container size
    const rect = container.getBoundingClientRect();
    const viewW = rect.width - ROW_DENDRO_W - COLOR_BAR_W;
    const viewH = rect.height - COL_DENDRO_H - COLOR_BAR_H;

    // Ensure minimum cell size — expand beyond viewport if needed (scrollable)
    const heatW = Math.max(viewW, nCols * cellSize);
    const heatH = Math.max(viewH, nRows * cellSize);

    // Size the grid to fit content (may exceed viewport → scrollable)
    if (gridEl) {
      gridEl.style.width = (ROW_DENDRO_W + COLOR_BAR_W + heatW) + 'px';
      gridEl.style.height = (COL_DENDRO_H + COLOR_BAR_H + heatH) + 'px';
    }

    const cellW = heatW / nCols;
    const cellH = heatH / nRows;

    // ── Draw heatmap on canvas ──
    canvas.width = heatW * window.devicePixelRatio;
    canvas.height = heatH * window.devicePixelRatio;
    canvas.style.width = heatW + 'px';
    canvas.style.height = heatH + 'px';

    const ctx = canvas.getContext('2d');
    ctx.scale(window.devicePixelRatio, window.devicePixelRatio);

    const imgData = ctx.createImageData(Math.ceil(heatW * window.devicePixelRatio), Math.ceil(heatH * window.devicePixelRatio));
    const pixels = imgData.data;
    const dpr = window.devicePixelRatio;

    // Direct pixel painting for speed
    for (let row = 0; row < nRows; row++) {
      const y0 = Math.floor(row * cellH * dpr);
      const y1 = Math.floor((row + 1) * cellH * dpr);
      for (let col = 0; col < nCols; col++) {
        const t = Math.min(z[row][col] / zMax, 1);
        const ci = Math.round(t * 255);
        const [r, g, b] = viridisLUT[ci];
        const x0 = Math.floor(col * cellW * dpr);
        const x1 = Math.floor((col + 1) * cellW * dpr);
        for (let py = y0; py < y1; py++) {
          for (let px = x0; px < x1; px++) {
            const idx = (py * imgData.width + px) * 4;
            pixels[idx] = r;
            pixels[idx + 1] = g;
            pixels[idx + 2] = b;
            pixels[idx + 3] = 255;
          }
        }
      }
    }
    ctx.putImageData(imgData, 0, 0);

    // ── Color bar helpers ──
    const colorLevel = (filters.colorMode === 'taxonomy' || filters.colorMode === undefined)
      ? getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter, filters.taxonFilterLevel)
      : null;
    const taxCmap = colorLevel && colorLevel !== 'group'
      ? buildTaxColorMap(colorLevel, filters.taxonFilter, filters.taxonFilterLevel).colorMap
      : null;

    // Build filtered set for dimming
    let filteredAsvSet = null;
    const taxonMatches = makeTaxonMatcher(filters.taxonFilter, filters.taxonFilterLevel);
    if (taxonMatches) {
      filteredAsvSet = new Set();
      const db = activeDb();
      const assigns = db ? store.taxonomy[db]?.assignments : {};
      for (const asvId in assigns) {
        if (taxonMatches(asvId)) filteredAsvSet.add(asvId);
      }
    }

    function getColor(id, type) {
      if (filters.colorMode === 'cluster') {
        const mode = type === 'sample' ? 'sampleCluster' : 'asvCluster';
        const k = type === 'sample' ? filters.sampleClusterK : filters.asvClusterK;
        return getClusterColor(id, mode, k);
      }
      if (type === 'asv') {
        if (filteredAsvSet && !filteredAsvSet.has(id)) return chrome().masked;
        if (filters.colorMode === 'group') {
          const asv = store.asvs.find(a => a.id === id);
          return GROUP_HEX[asv?.group ?? 'prokaryote'] ?? GROUP_HEX.unknown;
        }
        if (taxCmap) return getAsvColor(id, colorLevel, taxCmap);
        if (filteredAsvSet) return '#22d3ee';
      }
      if (type === 'sample' && filters.colorMode === 'cluster') {
        return getClusterColor(id, 'sampleCluster', filters.sampleClusterK);
      }
      return chrome().cell;
    }

    // ── Draw column dendrogram (SVG, above heatmap) ──
    colDendroSvg.innerHTML = '';
    colDendroSvg.setAttribute('width', heatW);
    colDendroSvg.setAttribute('height', COL_DENDRO_H);

    if (colDendro._phyloLines) {
      // NJ tree: line segments format (depth from root: 0=root, maxDepth=leaves)
      const { lines: dLines, nLeaves, maxDepth: treeMaxD } = colDendro._phyloLines;
      const nL = Math.max(nLeaves - 1, 1);
      const halfCell = cellW / 2;
      const pad = 2;
      for (const seg of dLines) {
        if (seg.x1 === seg.x2 && seg.y1 === seg.y2) continue;
        const x1 = (seg.x1 / nL) * (heatW - cellW) + halfCell;
        const x2 = (seg.x2 / nL) * (heatW - cellW) + halfCell;
        // Flip: root at top (y=0), leaves at bottom (y=COL_DENDRO_H)
        const y1 = pad + (seg.y1 / treeMaxD) * (COL_DENDRO_H - 2 * pad);
        const y2 = pad + (seg.y2 / treeMaxD) * (COL_DENDRO_H - 2 * pad);
        const line = document.createElementNS('http://www.w3.org/2000/svg', 'line');
        line.setAttribute('x1', x1); line.setAttribute('y1', y1);
        line.setAttribute('x2', x2); line.setAttribute('y2', y2);
        line.setAttribute('stroke', chrome().faint);
        line.setAttribute('stroke-width', '0.5');
        colDendroSvg.appendChild(line);
      }
    } else if (colDendro.icoord?.length > 0) {
      // Scipy format: U-shaped polylines
      const maxDist = Math.max(...colDendro.dcoord.flat());
      for (let i = 0; i < colDendro.icoord.length; i++) {
        const ix = colDendro.icoord[i];
        const iy = colDendro.dcoord[i];
        const points = [];
        for (let j = 0; j < 4; j++) {
          const px = ((ix[j] - 5) / (nCols * 10 - 10)) * (heatW - cellW) + cellW / 2;
          const py = COL_DENDRO_H - (iy[j] / maxDist) * (COL_DENDRO_H - 4);
          points.push(`${px},${py}`);
        }
        const line = document.createElementNS('http://www.w3.org/2000/svg', 'polyline');
        line.setAttribute('points', points.join(' '));
        line.setAttribute('fill', 'none');
        line.setAttribute('stroke', chrome().rule);
        line.setAttribute('stroke-width', '0.8');
        colDendroSvg.appendChild(line);
      }
    }

    // ── Draw row dendrogram (SVG, left of heatmap) ──
    rowDendroSvg.innerHTML = '';
    rowDendroSvg.setAttribute('width', ROW_DENDRO_W);
    rowDendroSvg.setAttribute('height', heatH);

    if (rowDendro.icoord.length > 0) {
      const maxDist = Math.max(...rowDendro.dcoord.flat());
      for (let i = 0; i < rowDendro.icoord.length; i++) {
        const ix = rowDendro.icoord[i];
        const iy = rowDendro.dcoord[i];
        const points = [];
        for (let j = 0; j < 4; j++) {
          const py = ((ix[j] - 5) / (nRows * 10 - 10)) * (heatH - cellH) + cellH / 2;
          const px = ROW_DENDRO_W - (iy[j] / maxDist) * (ROW_DENDRO_W - 4);
          points.push(`${px},${py}`);
        }
        const line = document.createElementNS('http://www.w3.org/2000/svg', 'polyline');
        line.setAttribute('points', points.join(' '));
        line.setAttribute('fill', 'none');
        line.setAttribute('stroke', chrome().rule);
        line.setAttribute('stroke-width', '0.8');
        rowDendroSvg.appendChild(line);
      }
    }

    // ── Column color bar (below col dendrogram, above heatmap) ──
    if (colColorBar) {
      colColorBar.width = heatW * dpr;
      colColorBar.height = COLOR_BAR_H * dpr;
      colColorBar.style.width = heatW + 'px';
      colColorBar.style.height = COLOR_BAR_H + 'px';
      const cbCtx = colColorBar.getContext('2d');
      cbCtx.scale(dpr, dpr);
      for (let col = 0; col < nCols; col++) {
        cbCtx.fillStyle = getColor(asvIds[col], 'asv');
        cbCtx.fillRect(col * cellW, 0, Math.ceil(cellW), COLOR_BAR_H);
      }
    }

    // ── Row color bar (right of row dendrogram, left of heatmap) ──
    if (rowColorBar) {
      rowColorBar.width = COLOR_BAR_W * dpr;
      rowColorBar.height = heatH * dpr;
      rowColorBar.style.width = COLOR_BAR_W + 'px';
      rowColorBar.style.height = heatH + 'px';
      const rbCtx = rowColorBar.getContext('2d');
      rbCtx.scale(dpr, dpr);
      for (let row = 0; row < nRows; row++) {
        rbCtx.fillStyle = getColor(sampleIds[row], 'sample');
        rbCtx.fillRect(0, row * cellH, COLOR_BAR_W, Math.ceil(cellH));
      }
    }

    // ── Hover handler ──
    canvas.onmousemove = (e) => {
      const cr = canvas.getBoundingClientRect();
      const mx = e.clientX - cr.left;
      const my = e.clientY - cr.top;
      const col = Math.floor(mx / cellW);
      const row = Math.floor(my / cellH);
      if (row >= 0 && row < nRows && col >= 0 && col < nCols) {
        const val = z[row][col];
        tooltip = {
          show: true,
          x: e.clientX,
          y: e.clientY,
          text: `${sampleIds[row]} | ${asvIds[col]} ${asvLabels[col]} | ${(val ** 4 * 100).toFixed(2)}%`,
        };
      } else {
        tooltip = { show: false, x: 0, y: 0, text: '' };
      }
    };
    canvas.onmouseleave = () => {
      tooltip = { show: false, x: 0, y: 0, text: '' };
    };
  });
</script>

<div class="flex h-full flex-col">
  {#if !heatmapData}
    <div class="flex-1 flex items-center justify-center">
      <div class="text-center">
        <div class="mb-4 h-8 w-8 animate-spin rounded-full border-2 border-blue-500 border-t-transparent mx-auto"></div>
        <p class="text-sm text-muted">Loading heatmap data...</p>
      </div>
    </div>
  {:else}
    <div class="flex-1 relative" bind:this={container}>
      <!-- Grid: dendro | colorbar | heatmap, with col dendro and colorbar on top -->
      <div class="absolute inset-0 overflow-auto" style="scrollbar-color: {chrome().cell} {chrome().masked};">
      <div bind:this={gridEl} class="grid" style="grid-template-columns: {ROW_DENDRO_W}px {COLOR_BAR_W}px auto; grid-template-rows: {COL_DENDRO_H}px {COLOR_BAR_H}px auto;">
        <!-- Row 1: spacer | spacer | col dendrogram -->
        <div></div>
        <div></div>
        <div class="overflow-hidden">
          <svg bind:this={colDendroSvg}></svg>
        </div>
        <!-- Row 2: spacer | spacer | col color bar -->
        <div></div>
        <div></div>
        <div class="overflow-hidden">
          <canvas bind:this={colColorBar} class="block"></canvas>
        </div>
        <!-- Row 3: row dendrogram | row color bar | heatmap -->
        <div class="overflow-hidden">
          <svg bind:this={rowDendroSvg}></svg>
        </div>
        <div class="overflow-hidden">
          <canvas bind:this={rowColorBar} class="block"></canvas>
        </div>
        <div class="overflow-hidden">
          <canvas bind:this={canvas} class="block"></canvas>
        </div>
      </div>
      </div>

      <!-- Title -->
      <div class="absolute top-1 left-1/2 -translate-x-1/2 text-xs text-faint pointer-events-none">
        {heatmapData.nSamples} samples × {visibleAsvCount}/{heatmapData.nAsvs} ASVs (max RA ≥{(filters.heatmapMinMaxRA ?? 1).toFixed(1)}%) — {usePhyloOrder ? 'Phylogeny (NJ)' : 'Ward'}
      </div>

      <!-- Tooltip -->
      {#if tooltip.show}
        <div class="fixed z-50 rounded bg-raised/95 px-3 py-1.5 text-xs text-fg shadow-lg pointer-events-none"
          style="left: {tooltip.x + 12}px; top: {tooltip.y - 20}px;">
          {tooltip.text}
        </div>
      {/if}
    </div>
  {/if}
</div>
