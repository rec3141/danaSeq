<script>
  import { onMount } from 'svelte';
  import * as d3 from 'd3';

  let { data = null, layout = 'radial', colorBy = 'binner', bins = [] } = $props();
  let container;

  const palette = ['#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa', '#fb923c',
                   '#2dd4bf', '#818cf8', '#f472b6', '#4ade80', '#e879f9', '#38bdf8',
                   '#94a3b8', '#d4d4d8', '#78716c', '#64748b', '#475569'];

  onMount(() => {
    if (data) render();
  });

  $effect(() => {
    if (container && data) {
      container.innerHTML = '';
      render();
    }
  });

  function getColorScale() {
    const binMap = {};
    for (const b of bins) {
      binMap[b.name] = b;
    }

    if (colorBy === 'binner') {
      const binners = [...new Set(bins.map(b => b.binner))];
      const scale = d3.scaleOrdinal().domain(binners).range(palette);
      return (d) => {
        const b = binMap[d.data.name];
        return b ? scale(b.binner) : '#475569';
      };
    }
    if (colorBy === 'domain') {
      const domains = [...new Set(bins.map(b => b.lineage?.domain).filter(Boolean))];
      const scale = d3.scaleOrdinal().domain(domains).range(palette);
      return (d) => {
        const b = binMap[d.data.name];
        return b?.lineage?.domain ? scale(b.lineage.domain) : '#475569';
      };
    }
    if (colorBy === 'phylum') {
      const phyla = [...new Set(bins.map(b => b.lineage?.phylum).filter(Boolean))];
      const scale = d3.scaleOrdinal().domain(phyla).range(palette);
      return (d) => {
        const b = binMap[d.data.name];
        return b?.lineage?.phylum ? scale(b.lineage.phylum) : '#475569';
      };
    }
    if (colorBy === 'quality') {
      return (d) => {
        const b = binMap[d.data.name];
        if (!b || b.completeness == null) return '#475569';
        const comp = b.completeness;
        const cont = b.contamination || 0;
        if (comp >= 90 && cont < 5) return '#34d399';  // high quality
        if (comp >= 50 && cont < 10) return '#fbbf24'; // medium quality
        return '#f87171'; // low quality
      };
    }
    return () => '#22d3ee';
  }

  function render() {
    if (!data || !container) return;

    const width = container.clientWidth || 600;
    const height = layout === 'radial' ? width : Math.max(400, bins.length * 16 + 80);
    const margin = { top: 20, right: 120, bottom: 20, left: 40 };

    const root = d3.hierarchy(data)
      .sort((a, b) => d3.ascending(a.data.name, b.data.name));

    const colorFn = getColorScale();

    const svg = d3.select(container)
      .append('svg')
      .attr('width', '100%')
      .attr('height', height);

    // Tooltip
    const tooltip = d3.select(container)
      .append('div')
      .attr('class', 'absolute bg-slate-800 text-slate-200 text-xs px-2 py-1 rounded border border-slate-600 pointer-events-none opacity-0 z-10')
      .style('transition', 'opacity 0.15s')
      .style('max-width', '300px');

    const binMap = {};
    for (const b of bins) binMap[b.name] = b;

    function tooltipHtml(d) {
      const b = binMap[d.data.name];
      if (!b) return `<strong>${d.data.name}</strong>`;
      let html = `<strong>${b.name}</strong>`;
      if (b.classification) html += `<br>${b.classification}`;
      if (b.completeness != null) html += `<br>Completeness: ${b.completeness}%`;
      if (b.contamination != null) html += `<br>Contamination: ${b.contamination}%`;
      if (b.ani != null) html += `<br>ANI: ${b.ani}`;
      if (b.method) html += `<br>Method: ${b.method}`;
      return html;
    }

    if (layout === 'radial') {
      const radius = width / 2 - 80;

      const cluster = d3.cluster()
        .size([2 * Math.PI, radius])
        .separation((a, b) => 1);

      cluster(root);

      svg.attr('viewBox', `${-width / 2} ${-width / 2} ${width} ${width}`);

      const g = svg.append('g');

      // Links
      g.selectAll('.link')
        .data(root.links())
        .join('path')
        .attr('class', 'link')
        .attr('fill', 'none')
        .attr('stroke', '#475569')
        .attr('stroke-width', 1)
        .attr('d', d => {
          return `M${radialPoint(d.target.x, d.target.y)}
                  L${radialPoint(d.target.x, d.source.y)}
                  L${radialPoint(d.source.x, d.source.y)}`;
        });

      // Nodes
      const nodes = g.selectAll('.node')
        .data(root.descendants())
        .join('g')
        .attr('class', 'node')
        .attr('transform', d => `translate(${radialPoint(d.x, d.y)})`);

      nodes.filter(d => !d.children)
        .append('circle')
        .attr('r', 4)
        .attr('fill', d => colorFn(d))
        .attr('stroke', d => {
          const b = binMap[d.data.name];
          return b?.binner === 'dastool' ? '#fff' : 'none';
        })
        .attr('stroke-width', d => {
          const b = binMap[d.data.name];
          return b?.binner === 'dastool' ? 1.5 : 0;
        })
        .on('mouseenter', (event, d) => {
          tooltip.html(tooltipHtml(d)).style('opacity', 1);
        })
        .on('mousemove', (event) => {
          const rect = container.getBoundingClientRect();
          tooltip
            .style('left', (event.clientX - rect.left + 10) + 'px')
            .style('top', (event.clientY - rect.top - 10) + 'px');
        })
        .on('mouseleave', () => tooltip.style('opacity', 0));

      // Labels for leaves
      nodes.filter(d => !d.children)
        .append('text')
        .attr('dy', '0.31em')
        .attr('x', d => d.x < Math.PI ? 6 : -6)
        .attr('text-anchor', d => d.x < Math.PI ? 'start' : 'end')
        .attr('transform', d => d.x >= Math.PI ? 'rotate(180)' : null)
        .attr('fill', '#94a3b8')
        .attr('font-size', '8px')
        .text(d => d.data.name);

    } else {
      // Rectangular layout
      const innerWidth = width - margin.left - margin.right;
      const innerHeight = height - margin.top - margin.bottom;

      const tree = d3.tree()
        .size([innerHeight, innerWidth]);

      tree(root);

      svg.attr('viewBox', `0 0 ${width} ${height}`);

      const g = svg.append('g')
        .attr('transform', `translate(${margin.left},${margin.top})`);

      // Links (elbow connectors)
      g.selectAll('.link')
        .data(root.links())
        .join('path')
        .attr('class', 'link')
        .attr('fill', 'none')
        .attr('stroke', '#475569')
        .attr('stroke-width', 1)
        .attr('d', d => {
          return `M${d.target.y},${d.target.x}
                  H${d.source.y}
                  V${d.source.x}`;
        });

      // Nodes
      const nodes = g.selectAll('.node')
        .data(root.descendants())
        .join('g')
        .attr('class', 'node')
        .attr('transform', d => `translate(${d.y},${d.x})`);

      nodes.filter(d => !d.children)
        .append('circle')
        .attr('r', 4)
        .attr('fill', d => colorFn(d))
        .attr('stroke', d => {
          const b = binMap[d.data.name];
          return b?.binner === 'dastool' ? '#fff' : 'none';
        })
        .attr('stroke-width', d => {
          const b = binMap[d.data.name];
          return b?.binner === 'dastool' ? 1.5 : 0;
        })
        .on('mouseenter', (event, d) => {
          tooltip.html(tooltipHtml(d)).style('opacity', 1);
        })
        .on('mousemove', (event) => {
          const rect = container.getBoundingClientRect();
          tooltip
            .style('left', (event.clientX - rect.left + 10) + 'px')
            .style('top', (event.clientY - rect.top - 10) + 'px');
        })
        .on('mouseleave', () => tooltip.style('opacity', 0));

      // Labels
      nodes.filter(d => !d.children)
        .append('text')
        .attr('dy', '0.31em')
        .attr('x', 8)
        .attr('text-anchor', 'start')
        .attr('fill', '#94a3b8')
        .attr('font-size', '9px')
        .text(d => d.data.name);
    }
  }

  function radialPoint(angle, radius) {
    return `${radius * Math.cos(angle - Math.PI / 2)},${radius * Math.sin(angle - Math.PI / 2)}`;
  }
</script>

<div bind:this={container} class="relative w-full" style="min-height: 300px;"></div>
