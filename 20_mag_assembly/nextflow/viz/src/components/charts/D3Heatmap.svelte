<script>
  import { onMount } from 'svelte';
  import * as d3 from 'd3';

  let { data = null, onCellClick = null, onRowClick = null,
        tooltipFormat = null, legendLabel = null, selectedRow = null } = $props();
  let container;

  onMount(() => {
    if (data) render(data);
  });

  $effect(() => {
    if (container && data) {
      container.innerHTML = '';
      render(data);
    }
  });

  function render(data) {
    const { mag_ids, module_ids, module_names, matrix, row_order, col_order } = data;
    if (!mag_ids.length || !module_ids.length) return;

    // Apply clustering order
    const rOrder = row_order && row_order.length === mag_ids.length
      ? row_order : mag_ids.map((_, i) => i);
    const cOrder = col_order && col_order.length === module_ids.length
      ? col_order : module_ids.map((_, i) => i);

    const orderedMags = rOrder.map(i => mag_ids[i]);
    const orderedMods = cOrder.map(i => module_ids[i]);
    const orderedNames = cOrder.map(i => module_names[i]);
    const orderedMatrix = rOrder.map(ri => cOrder.map(ci => matrix[ri][ci]));

    const cellW = 14;
    const cellH = 16;
    const marginLeft = 180;
    const marginTop = 260;
    const marginBottom = 20;
    const marginRight = 100;
    const width = marginLeft + orderedMods.length * cellW + marginRight;
    const height = marginTop + orderedMags.length * cellH + marginBottom;

    const colorScale = d3.scaleSequential(d3.interpolateViridis).domain([0, 1]);

    const defaultTooltip = (val, magId, modId, modName) =>
      `<strong>${modId}</strong><br>${modName}<br>MAG: ${magId}<br>Completeness: ${(val * 100).toFixed(0)}%`;
    const formatTooltip = tooltipFormat || defaultTooltip;

    const svg = d3.select(container)
      .append('svg')
      .attr('width', width)
      .attr('height', height);

    // Tooltip
    const tooltip = d3.select(container)
      .append('div')
      .attr('class', 'absolute bg-slate-800 text-slate-200 text-xs px-2 py-1 rounded border border-slate-600 pointer-events-none opacity-0 z-10')
      .style('transition', 'opacity 0.15s');

    // Cells
    for (let i = 0; i < orderedMags.length; i++) {
      for (let j = 0; j < orderedMods.length; j++) {
        const val = orderedMatrix[i][j];
        svg.append('rect')
          .attr('x', marginLeft + j * cellW)
          .attr('y', marginTop + i * cellH)
          .attr('width', cellW - 1)
          .attr('height', cellH - 1)
          .attr('fill', val > 0 ? colorScale(val) : '#0f172a')
          .attr('rx', 1)
          .style('cursor', 'pointer')
          .on('mouseover', function (event) {
            d3.select(this).attr('stroke', '#22d3ee').attr('stroke-width', 1.5);
            tooltip.style('opacity', 1)
              .html(formatTooltip(val, orderedMags[i], orderedMods[j], orderedNames[j]));
          })
          .on('mousemove', function (event) {
            const rect = container.getBoundingClientRect();
            tooltip
              .style('left', (event.clientX - rect.left + 12) + 'px')
              .style('top', (event.clientY - rect.top - 10) + 'px');
          })
          .on('mouseout', function () {
            d3.select(this).attr('stroke', 'none');
            tooltip.style('opacity', 0);
          })
          .on('click', function () {
            if (onCellClick) onCellClick(orderedMags[i], orderedMods[j], val);
            if (onRowClick) onRowClick(orderedMags[i]);
          });
      }

      // Row highlight for selected MAG
      if (selectedRow && orderedMags[i] === selectedRow) {
        svg.append('rect')
          .attr('x', marginLeft - 2)
          .attr('y', marginTop + i * cellH - 1)
          .attr('width', orderedMods.length * cellW + 4)
          .attr('height', cellH + 1)
          .attr('fill', 'none')
          .attr('stroke', '#22d3ee')
          .attr('stroke-width', 1.5)
          .attr('rx', 2)
          .attr('pointer-events', 'none');
      }
    }

    // Row labels (MAG names)
    for (let i = 0; i < orderedMags.length; i++) {
      const isSelected = selectedRow && orderedMags[i] === selectedRow;
      svg.append('text')
        .attr('x', marginLeft - 4)
        .attr('y', marginTop + i * cellH + cellH / 2)
        .attr('text-anchor', 'end')
        .attr('dominant-baseline', 'central')
        .attr('fill', isSelected ? '#22d3ee' : '#94a3b8')
        .attr('font-size', '9px')
        .attr('font-weight', isSelected ? 'bold' : 'normal')
        .attr('font-family', 'JetBrains Mono, monospace')
        .text(orderedMags[i])
        .style('cursor', 'pointer')
        .on('click', () => onRowClick && onRowClick(orderedMags[i]));
    }

    // Column labels (module names, rotated)
    for (let j = 0; j < orderedMods.length; j++) {
      svg.append('text')
        .attr('x', 0)
        .attr('y', 0)
        .attr('transform', `translate(${marginLeft + j * cellW + cellW / 2}, ${marginTop - 6}) rotate(-70)`)
        .attr('text-anchor', 'start')
        .attr('dominant-baseline', 'hanging')
        .attr('fill', '#64748b')
        .attr('font-size', '8px')
        .attr('font-family', 'Inter, sans-serif')
        .text(orderedNames[j]);
    }

    // Color legend
    const legendWidth = 120;
    const legendHeight = 10;
    const legendX = marginLeft;
    const legendY = height - 4;

    const defs = svg.append('defs');
    const gradient = defs.append('linearGradient').attr('id', 'heatmap-grad');
    for (let i = 0; i <= 10; i++) {
      const t = i / 10;
      gradient.append('stop')
        .attr('offset', `${i * 10}%`)
        .attr('stop-color', colorScale(t));
    }

    svg.append('rect')
      .attr('x', legendX).attr('y', legendY)
      .attr('width', legendWidth).attr('height', legendHeight)
      .attr('fill', 'url(#heatmap-grad)').attr('rx', 2);

    const legendLow = legendLabel ? legendLabel[0] : '0%';
    const legendHigh = legendLabel ? legendLabel[1] : '100%';
    svg.append('text').attr('x', legendX).attr('y', legendY - 2)
      .attr('fill', '#64748b').attr('font-size', '9px').text(legendLow);
    svg.append('text').attr('x', legendX + legendWidth).attr('y', legendY - 2)
      .attr('text-anchor', 'end').attr('fill', '#64748b').attr('font-size', '9px').text(legendHigh);
  }
</script>

<div bind:this={container} class="relative w-full overflow-x-auto"></div>
