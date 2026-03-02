<script>
  import { onMount } from 'svelte';
  import * as d3 from 'd3';

  let { data = null, colorDepth = 2, exportName = 'danaseq_sunburst' } = $props();
  let container;

  function exportSvg() {
    const svg = container?.querySelector('svg');
    if (!svg) return;
    const clone = svg.cloneNode(true);
    clone.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
    // Add dark background rect so export looks correct
    const bg = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
    const vb = clone.getAttribute('viewBox')?.split(' ').map(Number) || [0, 0, 500, 500];
    bg.setAttribute('x', vb[0]); bg.setAttribute('y', vb[1]);
    bg.setAttribute('width', vb[2]); bg.setAttribute('height', vb[3]);
    bg.setAttribute('fill', '#1e293b');
    clone.insertBefore(bg, clone.firstChild);
    const blob = new Blob([new XMLSerializer().serializeToString(clone)], { type: 'image/svg+xml' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url; a.download = `${exportName}.svg`; a.click();
    URL.revokeObjectURL(url);
  }

  onMount(() => { if (data) render(data); });

  $effect(() => {
    if (container && data) { container.innerHTML = ''; render(data); }
  });

  const palette = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
            '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8',
            '#94a3b8','#d4d4d8','#78716c','#64748b','#475569'];

  function render(data) {
    const width = container.clientWidth || 500;
    const height = width;
    const radius = width / 2;

    const colorScale = d3.scaleOrdinal().range(palette);
    const hierarchy = d3.hierarchy(data).sum(d => d.value || 0).sort((a, b) => b.value - a.value);
    const root = d3.partition().size([2 * Math.PI, radius])(hierarchy);

    // Store initial angles for zoom transitions
    root.each(d => { d._x0 = d.x0; d._x1 = d.x1; d._y0 = d.y0; d._y1 = d.y1; });

    let currentRoot = root;

    const arc = d3.arc()
      .startAngle(d => d._cx0)
      .endAngle(d => d._cx1)
      .padAngle(d => Math.min((d._cx1 - d._cx0) / 2, 0.003))
      .padRadius(radius / 3)
      .innerRadius(d => d._cy0)
      .outerRadius(d => Math.max(d._cy0, d._cy1 - 1));

    // Current view coordinates
    root.each(d => { d._cx0 = d.x0; d._cx1 = d.x1; d._cy0 = d.y0; d._cy1 = d.y1; });

    const svg = d3.select(container).append('svg')
      .attr('viewBox', `${-width/2} ${-height/2} ${width} ${height}`)
      .attr('width', '100%').attr('height', height);

    const tooltip = d3.select(container).append('div')
      .attr('class', 'absolute bg-slate-800 text-slate-200 text-xs px-2 py-1 rounded border border-slate-600 pointer-events-none opacity-0 z-10')
      .style('transition', 'opacity 0.15s');

    function getColor(d) {
      let node = d;
      while (node.depth > colorDepth && node.parent) node = node.parent;
      return colorScale(node.data.name);
    }

    function isVisible(d) {
      return d._cx1 > d._cx0 + 0.001 && d._cy1 > d._cy0 + 1;
    }

    const paths = svg.selectAll('path')
      .data(root.descendants().filter(d => d.depth > 0))
      .join('path')
      .attr('fill', getColor)
      .attr('fill-opacity', d => Math.max(0.3, 1 - (d.depth - 1) * 0.08))
      .attr('d', arc)
      .attr('display', d => isVisible(d) ? null : 'none')
      .style('cursor', 'pointer')
      .on('mouseover', function (event, d) {
        d3.select(this).attr('fill-opacity', 1).attr('stroke', '#e2e8f0').attr('stroke-width', 1);
        const pct = ((d.value / root.value) * 100).toFixed(1);
        const pathStr = d.ancestors().map(a => a.data.name).reverse().slice(1).join(' > ');
        tooltip.style('opacity', 1)
          .html(`<strong>${d.data.name}</strong><br><span style="color:#64748b">${pathStr}</span><br>${pct}%`);
      })
      .on('mousemove', function (event) {
        const rect = container.getBoundingClientRect();
        tooltip.style('left', (event.clientX - rect.left + 12) + 'px')
          .style('top', (event.clientY - rect.top - 10) + 'px');
      })
      .on('mouseout', function (event, d) {
        d3.select(this).attr('fill-opacity', Math.max(0.3, 1 - (d.depth - 1) * 0.08))
          .attr('stroke', null);
        tooltip.style('opacity', 0);
      })
      .on('click', function (event, d) {
        event.stopPropagation();
        zoomTo(d);
      });

    // Center label
    const centerText = svg.append('text')
      .attr('text-anchor', 'middle').attr('dy', '0.35em')
      .attr('fill', '#94a3b8').attr('font-size', '13px')
      .attr('font-family', 'Inter, sans-serif')
      .style('cursor', 'pointer')
      .text('Community')
      .on('click', () => zoomTo(root));

    function isDescendant(node, ancestor) {
      // In partition layout, descendants have angular range within ancestor's range
      return node.x0 >= ancestor.x0 && node.x1 <= ancestor.x1 && node.depth > ancestor.depth;
    }

    function zoomTo(target) {
      currentRoot = target;

      // Compute new angular and radial extents
      const xScale = d3.scaleLinear()
        .domain([target.x0, target.x1])
        .range([0, 2 * Math.PI]);
      const yScale = d3.scaleLinear()
        .domain([target.y0, radius])
        .range([0, radius]);

      const t = svg.transition().duration(600);

      paths.transition(t)
        .tween('data', d => {
          const show = target === root ? d.depth > 0 : isDescendant(d, target);
          const tx0 = show ? xScale(d.x0) : 0;
          const tx1 = show ? xScale(d.x1) : 0;
          const ty0 = show ? Math.max(0, yScale(d.y0)) : 0;
          const ty1 = show ? Math.max(0, yScale(d.y1)) : 0;
          const ix0 = d3.interpolate(d._cx0, tx0);
          const ix1 = d3.interpolate(d._cx1, tx1);
          const iy0 = d3.interpolate(d._cy0, ty0);
          const iy1 = d3.interpolate(d._cy1, ty1);
          return t => {
            d._cx0 = ix0(t); d._cx1 = ix1(t);
            d._cy0 = iy0(t); d._cy1 = iy1(t);
          };
        })
        .attrTween('d', d => () => arc(d))
        .attr('fill-opacity', d => {
          const relDepth = d.depth - target.depth;
          return relDepth <= 0 ? 0 : Math.max(0.3, 1 - (relDepth - 1) * 0.08);
        })
        .on('end', function(d) {
          const show = target === root ? d.depth > 0 : isDescendant(d, target);
          d3.select(this).attr('display', show ? null : 'none');
        });

      // Update center label
      const path = target.ancestors().map(a => a.data.name).reverse().slice(1).join(' > ');
      centerText.text(target === root ? 'Community' : path || target.data.name);
    }
  }
</script>

<div bind:this={container} class="relative w-full">
  <button
    class="absolute top-1 right-1 z-10 p-1 rounded text-slate-500 hover:text-cyan-400 hover:bg-slate-800/80 transition-colors"
    onclick={exportSvg}
    title="Download SVG"
  >
    <svg class="w-4 h-4" fill="none" stroke="currentColor" viewBox="0 0 24 24">
      <path stroke-linecap="round" stroke-linejoin="round" stroke-width="2" d="M4 16v1a3 3 0 003 3h10a3 3 0 003-3v-1m-4 0l-4 4m0 0l-4-4m4 4V4"/>
    </svg>
  </button>
</div>
