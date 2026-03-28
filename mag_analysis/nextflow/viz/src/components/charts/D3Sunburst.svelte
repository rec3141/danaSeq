<script>
  import { onMount } from 'svelte';
  import * as d3 from 'd3';

  let { data = null, colorMap = null, colorMaps = null, colorDepth = 1 } = $props();
  let container;

  onMount(() => {
    if (!data) return;
    render(data);
  });

  $effect(() => {
    if (container && data) {
      container.innerHTML = '';
      render(data);
    }
  });

  const defaultPalette = ['#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa', '#fb923c',
              '#2dd4bf', '#818cf8', '#f472b6', '#4ade80', '#e879f9', '#38bdf8',
              '#94a3b8', '#d4d4d8', '#78716c', '#64748b', '#475569'];

  function render(data) {
    const width = container.clientWidth || 500;
    const height = width;
    const radius = width / 2;

    const fallbackColor = d3.scaleOrdinal().range(defaultPalette);

    const hierarchy = d3.hierarchy(data)
      .sum(d => d.value || 0)
      .sort((a, b) => b.value - a.value);

    const partition = d3.partition()
      .size([2 * Math.PI, radius]);

    const root = partition(hierarchy);

    const arc = d3.arc()
      .startAngle(d => d.x0)
      .endAngle(d => d.x1)
      .padAngle(d => Math.min((d.x1 - d.x0) / 2, 0.005))
      .padRadius(radius / 2)
      .innerRadius(d => d.y0)
      .outerRadius(d => d.y1 - 1);

    const svg = d3.select(container)
      .append('svg')
      .attr('viewBox', `${-width / 2} ${-height / 2} ${width} ${height}`)
      .attr('width', '100%')
      .attr('height', height);

    // Tooltip
    const tooltip = d3.select(container)
      .append('div')
      .attr('class', 'absolute bg-slate-800 text-slate-200 text-xs px-2 py-1 rounded border border-slate-600 pointer-events-none opacity-0 z-10')
      .style('transition', 'opacity 0.15s');

    const paths = svg.selectAll('path')
      .data(root.descendants().filter(d => d.depth > 0))
      .join('path')
      .attr('fill', d => {
        // Per-depth color maps: each ring gets its own stable colors
        if (colorMaps && colorMaps[d.depth]) {
          const c = colorMaps[d.depth][d.data.name];
          if (c) return c;
        }
        // Fallback: walk up to colorDepth ancestor
        let node = d;
        while (node.depth > colorDepth && node.parent) node = node.parent;
        const name = node.data.name;
        if (colorMap && colorMap[name]) return colorMap[name];
        return fallbackColor(name);
      })
      .attr('fill-opacity', d => 1 - d.depth * 0.15)
      .attr('d', arc)
      .style('cursor', 'pointer')
      .on('mouseover', function (event, d) {
        d3.select(this).attr('fill-opacity', 1);
        const pct = ((d.value / root.value) * 100).toFixed(1);
        const path = d.ancestors().map(a => a.data.name).reverse().slice(1).join(' > ');
        tooltip
          .style('opacity', 1)
          .html(`<strong>${d.data.name}</strong><br>${path}<br>${pct}% (${(d.value / 1e6).toFixed(1)} Mbp)`);
      })
      .on('mousemove', function (event) {
        const rect = container.getBoundingClientRect();
        tooltip
          .style('left', (event.clientX - rect.left + 12) + 'px')
          .style('top', (event.clientY - rect.top - 10) + 'px');
      })
      .on('mouseout', function (event, d) {
        d3.select(this).attr('fill-opacity', 1 - d.depth * 0.15);
        tooltip.style('opacity', 0);
      });

    // Center label
    svg.append('text')
      .attr('text-anchor', 'middle')
      .attr('dy', '0.35em')
      .attr('fill', '#94a3b8')
      .attr('font-size', '14px')
      .attr('font-family', 'Inter, sans-serif')
      .text('Community');
  }
</script>

<div bind:this={container} class="relative w-full"></div>
