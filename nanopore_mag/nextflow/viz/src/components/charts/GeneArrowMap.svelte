<script>
  import { onMount } from 'svelte';
  import * as d3 from 'd3';

  let { features = [], contigLength = 0, height = 160, highlightId = null } = $props();
  let svgEl;
  let containerEl;

  const COLORS = {
    cds:    '#22d3ee',   // cyan
    trna:   '#34d399',   // emerald
    rrna:   '#fbbf24',   // amber
    tmrna:  '#a78bfa',   // violet
    ncrna:  '#f472b6',   // pink
    crispr: '#f87171',   // red
  };

  // Tooltip state
  let tipText = $state('');
  let tipX = $state(0);
  let tipY = $state(0);
  let tipShow = $state(false);

  // Zoom state
  let currentTransform = d3.zoomIdentity;
  let zoomBehavior = null;

  function arrowPath(x1, x2, y, h, strand) {
    const w = Math.abs(x2 - x1);
    const tipLen = Math.min(w * 0.3, h * 0.7);
    if (w < 1) return ''; // too small to draw

    if (strand === 1) {
      const tx = x2 - tipLen;
      return `M${x1},${y-h/2} L${tx},${y-h/2} L${x2},${y} L${tx},${y+h/2} L${x1},${y+h/2} Z`;
    } else {
      const tx = x1 + tipLen;
      return `M${x2},${y-h/2} L${tx},${y-h/2} L${x1},${y} L${tx},${y+h/2} L${x2},${y+h/2} Z`;
    }
  }

  function render(transform) {
    if (!svgEl) return;
    const svg = d3.select(svgEl);
    svg.selectAll('*').remove();

    if (!features.length || !contigLength) return;

    const rect = svgEl.getBoundingClientRect();
    const width = rect.width || 600;
    const margin = { left: 10, right: 10, top: 30, bottom: 28 };
    const plotW = width - margin.left - margin.right;
    const plotH = height - margin.top - margin.bottom;

    svg.attr('width', width).attr('height', height);

    // Clip path for gene area
    const clipId = 'gene-clip-' + Math.random().toString(36).slice(2, 8);
    svg.append('defs').append('clipPath').attr('id', clipId)
      .append('rect')
      .attr('x', 0).attr('y', -margin.top)
      .attr('width', plotW).attr('height', plotH + margin.top);

    const baseXScale = d3.scaleLinear()
      .domain([0, contigLength])
      .range([0, plotW]);

    // Apply zoom transform to x scale
    const xScale = transform ? transform.rescaleX(baseXScale) : baseXScale;

    const g = svg.append('g')
      .attr('transform', `translate(${margin.left},${margin.top})`);

    const geneG = g.append('g').attr('clip-path', `url(#${clipId})`);

    // Contig backbone
    geneG.append('line')
      .attr('x1', xScale(0)).attr('y1', plotH / 2)
      .attr('x2', xScale(contigLength)).attr('y2', plotH / 2)
      .attr('stroke', '#334155').attr('stroke-width', 1.5);

    // Track layout: forward strand above center, reverse below
    const trackH = 18;
    const gap = 4;
    const fwdY = plotH / 2 - gap - trackH / 2;
    const revY = plotH / 2 + gap + trackH / 2;

    // Only draw features visible in the current viewport
    const [domainMin, domainMax] = xScale.domain();

    for (const f of features) {
      // Skip features entirely outside viewport
      if (f.e < domainMin || f.s > domainMax) continue;

      const x1 = xScale(f.s);
      const x2 = xScale(f.e);
      const y = f.d === 1 ? fwdY : revY;
      const color = COLORS[f.t] || COLORS.cds;
      const path = arrowPath(x1, x2, y, trackH, f.d);
      if (!path) continue;

      geneG.append('path')
        .attr('d', path)
        .attr('fill', color)
        .attr('fill-opacity', 0.8)
        .attr('stroke', color)
        .attr('stroke-opacity', 1)
        .attr('stroke-width', 0.5)
        .attr('class', 'gene-arrow')
        .on('mouseenter', (event) => {
          const gene = f.g || '';
          const product = f.p || '';
          const label = gene || product || 'hypothetical protein';
          const coords = `${f.s.toLocaleString()}..${f.e.toLocaleString()} (${f.d === 1 ? '+' : '-'})`;
          const lines = [label];
          if (product && product !== label) lines.push(product);
          lines.push(`${f.t.toUpperCase()} | ${coords} | ${(f.e - f.s + 1).toLocaleString()} bp`);
          tipText = lines.join('\n');
          const r = svgEl.getBoundingClientRect();
          let tx = event.clientX - r.left;
          // Clamp tooltip to stay within container bounds
          const containerW = containerEl?.offsetWidth || width;
          tx = Math.max(100, Math.min(tx, containerW - 100));
          tipX = tx;
          tipY = y + margin.top - trackH / 2 - 4;
          tipShow = true;
          d3.select(event.target).attr('fill-opacity', 1).attr('stroke', '#fff').attr('stroke-width', 1.5);
        })
        .on('mouseleave', (event) => {
          tipShow = false;
          d3.select(event.target).attr('fill-opacity', 0.8).attr('stroke', color).attr('stroke-width', 0.5);
        });
    }

    // Axis
    const axisG = g.append('g')
      .attr('transform', `translate(0,${plotH})`);

    const axis = d3.axisBottom(xScale)
      .ticks(6)
      .tickSize(4)
      .tickFormat(d => {
        if (d >= 1e6) return (d / 1e6).toFixed(1) + ' Mb';
        if (d >= 1e3) return (d / 1e3).toFixed(0) + ' Kb';
        return d + '';
      });

    axisG.call(axis);
    axisG.selectAll('line, path').attr('stroke', '#475569');
    axisG.selectAll('text').attr('fill', '#64748b').attr('font-size', '10px');

    // Legend
    const legendTypes = [...new Set(features.map(f => f.t))].sort();
    const legendG = g.append('g')
      .attr('transform', `translate(${plotW - legendTypes.length * 55}, -2)`);
    legendTypes.forEach((t, i) => {
      const lg = legendG.append('g').attr('transform', `translate(${i * 55}, 0)`);
      lg.append('rect')
        .attr('width', 8).attr('height', 8).attr('rx', 1)
        .attr('fill', COLORS[t] || COLORS.cds).attr('opacity', 0.85);
      lg.append('text')
        .attr('x', 11).attr('y', 7)
        .attr('fill', '#94a3b8').attr('font-size', '9px')
        .text(t.toUpperCase());
    });

    // Attach zoom behavior to SVG
    // Allow panning beyond the contig extent so edge features can be centered
    const padding = plotW * 0.3;
    if (!zoomBehavior) {
      zoomBehavior = d3.zoom()
        .scaleExtent([1, 50])
        .translateExtent([[-padding, 0], [width + padding, height]])
        .extent([[0, 0], [width, height]])
        .filter(event => {
          if (event.type === 'dblclick') return false;
          return true;
        })
        .on('zoom', (event) => {
          currentTransform = event.transform;
          tipShow = false;
          render(currentTransform);
        });

      svg.call(zoomBehavior);
      svg.on('dblclick.zoom', null);
    } else {
      // Update translate extent for current width
      zoomBehavior.translateExtent([[-padding, 0], [width + padding, height]]);
      svg.call(zoomBehavior);
      svg.on('dblclick.zoom', null);
    }
  }

  // Reset zoom when contig changes
  $effect(() => {
    const _deps = [features, contigLength];
    currentTransform = d3.zoomIdentity;
    zoomBehavior = null;
  });

  onMount(() => { render(currentTransform); });

  $effect(() => {
    const _deps = [features, contigLength, height, highlightId];
    render(currentTransform);
  });
</script>

<div class="relative w-full overflow-x-clip overflow-y-visible" bind:this={containerEl}>
  <svg bind:this={svgEl} class="w-full cursor-grab active:cursor-grabbing" style="height: {height}px;"></svg>
  {#if tipShow}
    <div
      class="absolute pointer-events-none bg-slate-900/95 text-slate-200 text-xs px-2 py-1 rounded shadow-lg border border-slate-600 z-10 whitespace-pre-line max-w-sm"
      style="left: {tipX}px; top: {tipY}px; transform: translate(-50%, -100%)"
    >
      {tipText}
    </div>
  {/if}
  <div class="text-[9px] text-slate-600 text-right pr-1 -mt-1">scroll to zoom, drag to pan</div>
</div>
