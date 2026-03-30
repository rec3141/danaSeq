<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';

  let { data = null, colorBy = 'sample', sizeBy = 'fixed', sizeScale = 1.0, mode = 'tsne', colorMap = {}, coordExtents = null, onselect = null, onclick = null, searchMatchIds = null, idField = 'id', exportName = 'danaseq_scatter' } = $props();
  let container;
  let listenersAttached = false;

  const darkLayout = {
    paper_bgcolor: '#1e293b', plot_bgcolor: '#0f172a',
    font: { color: '#94a3b8', family: 'Inter, sans-serif', size: 12 },
    xaxis: { showgrid: false, zeroline: false, showticklabels: false, title: '' },
    yaxis: { showgrid: false, zeroline: false, showticklabels: false, title: '' },
    margin: { t: 10, r: 10, b: 10, l: 10 },
    modebar: { orientation: 'v', bgcolor: 'transparent', color: '#64748b', activecolor: '#22d3ee' },
    legend: { bgcolor: 'rgba(15,23,42,0.7)', font: { color: '#94a3b8', size: 10 },
              x: 1, y: 1, xanchor: 'right', yanchor: 'top' },
    dragmode: 'pan',
    hovermode: 'closest',
  };

  const palette = ['#22d3ee','#34d399','#fbbf24','#f87171','#a78bfa','#fb923c',
                   '#2dd4bf','#818cf8','#f472b6','#4ade80','#e879f9','#38bdf8',
                   '#94a3b8','#d4d4d8','#78716c'];

  // Precomputed sizing normalization (set in getTraces)
  let sizeNorm = null;

  function markerSize(c, sizeBy, scale) {
    if (sizeBy === 'fixed') return 8 * scale;
    const val = typeof c[sizeBy] === 'number' ? c[sizeBy] : 0;
    if (val <= 0 || !sizeNorm) return 4 * scale;
    const t = (Math.log1p(val) - sizeNorm.logMin) / sizeNorm.range;
    return (4 + t * 20) * scale; // 4px to 24px full range
  }

  function hoverText(c) {
    const parts = [`<b>${c[idField] || c.id}</b>`];
    if (c.flowcell) parts.push(`Flowcell: ${c.flowcell}`);
    if (c.barcode) parts.push(`Barcode: ${c.barcode}`);
    if (c.read_count != null) parts.push(`Reads: ${c.read_count.toLocaleString()}`);
    if (c.total_bases != null) parts.push(`Bases: ${(c.total_bases / 1e6).toFixed(1)} Mbp`);
    if (c.avg_length != null) parts.push(`Avg length: ${Math.round(c.avg_length).toLocaleString()} bp`);
    // Show the value of the active colorBy field
    const colorVal = c[colorBy];
    if (colorVal != null && colorVal !== '') {
      const label = colorBy.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
      parts.push(`${label}: ${colorVal}`);
    }
    return parts.join('<br>');
  }

  function isContinuousField(points, field) {
    let numCount = 0, total = 0;
    for (const c of points) {
      if (c[field] != null) { total++; if (typeof c[field] === 'number') numCount++; }
    }
    return total > 0 && numCount > total * 0.8;
  }

  function getTraces() {
    if (!data?.points?.length) return [];
    const points = data.points;
    const xKey = `${mode}_x`, yKey = `${mode}_y`;

    // Precompute size normalization from data range
    sizeNorm = null;
    if (sizeBy !== 'fixed') {
      const vals = points.map(c => c[sizeBy]).filter(v => typeof v === 'number' && v > 0);
      if (vals.length > 1) {
        const logMin = Math.log1p(Math.min(...vals));
        const logMax = Math.log1p(Math.max(...vals));
        sizeNorm = { logMin, logMax, range: logMax - logMin || 1 };
      }
    }

    if (isContinuousField(points, colorBy)) {
      const fn = c => typeof c[colorBy] === 'number' ? c[colorBy] : 0;
      const sorted = [...points].sort((a, b) => fn(a) - fn(b));
      const label = colorBy.replace(/_/g, ' ').replace(/\b\w/g, l => l.toUpperCase());
      return [{
        type: 'scattergl', mode: 'markers',
        x: sorted.map(c => c[xKey]), y: sorted.map(c => c[yKey]),
        text: sorted.map(c => hoverText(c)),
        customdata: sorted.map(c => c[idField] || c.id),
        marker: {
          size: sorted.map(c => markerSize(c, sizeBy, sizeScale)),
          color: sorted.map(fn), opacity: 0.75,
          colorscale: 'Viridis', showscale: true,
          colorbar: {
            title: { text: label, font: { color: '#94a3b8' } },
            tickfont: { color: '#94a3b8' },
            x: 0.99, xanchor: 'right', y: 0.98, yanchor: 'top',
            len: 0.4, thickness: 12, bgcolor: 'rgba(15,23,42,0.7)', borderwidth: 0,
          },
          line: { width: 0 },
        },
        hoverinfo: 'text',
      }];
    }

    // Categorical
    const groups = {};
    for (const c of points) {
      const key = c[colorBy] || 'Unknown';
      if (!groups[key]) groups[key] = { x: [], y: [], text: [], sizes: [], ids: [] };
      groups[key].x.push(c[xKey]);
      groups[key].y.push(c[yKey]);
      groups[key].text.push(hoverText(c));
      groups[key].sizes.push(markerSize(c, sizeBy, sizeScale));
      groups[key].ids.push(c[idField] || c.id);
    }

    const useParentMap = colorMap && Object.keys(colorMap).length > 0;
    const effectiveColorMap = useParentMap ? colorMap : (() => {
      const allNames = Object.keys(groups).filter(n => n !== 'Unknown').sort();
      const m = {};
      allNames.forEach((name, i) => { m[name] = palette[i % palette.length]; });
      return m;
    })();

    const entries = Object.entries(groups).sort((a, b) => {
      if (a[0] === 'Unknown') return -1;
      if (b[0] === 'Unknown') return 1;
      return a[0].localeCompare(b[0]);
    });

    return entries.map(([name, g]) => {
      const isBg = name === 'Unknown';
      const opacity = searchMatchIds
        ? g.ids.map(id => searchMatchIds.has(id) ? (isBg ? 0.25 : 0.75) : 0.005)
        : (isBg ? 0.25 : 0.75);
      return {
        type: 'scattergl', mode: 'markers',
        name: name.length > 25 ? name.slice(0, 23) + '..' : name,
        x: g.x, y: g.y, text: g.text, customdata: g.ids,
        marker: {
          size: g.sizes,
          color: isBg ? '#475569' : (effectiveColorMap[name] || '#94a3b8'),
          opacity, line: { width: 0 },
        },
        hoverinfo: 'text',
      };
    });
  }

  function plot() {
    const traces = getTraces();
    const ext = coordExtents;
    const pad = 0.05;
    const xAxisCfg = { ...darkLayout.xaxis, uirevision: mode };
    const yAxisCfg = { ...darkLayout.yaxis, uirevision: mode, scaleanchor: 'x', scaleratio: 1 };
    if (ext) {
      const xPad = (ext.xMax - ext.xMin) * pad;
      const yPad = (ext.yMax - ext.yMin) * pad;
      xAxisCfg.range = [ext.xMin - xPad, ext.xMax + xPad];
      yAxisCfg.range = [ext.yMin - yPad, ext.yMax + yPad];
      xAxisCfg.autorange = false;
      yAxisCfg.autorange = false;
    }
    const layout = { ...darkLayout, xaxis: xAxisCfg, yaxis: yAxisCfg, uirevision: mode + '-' + colorBy };
    Plotly.react(container, traces, layout, {
      responsive: true, displaylogo: false, scrollZoom: true,
      doubleClick: false, displayModeBar: 'hover',
      toImageButtonOptions: { format: 'svg', filename: exportName },
      modeBarButtonsToRemove: ['autoScale2d','toggleSpikelines','hoverCompareCartesian','hoverClosestCartesian','zoomIn2d','zoomOut2d'],
    });

    if (!listenersAttached) {
      // Shift held = lasso mode, shift released = back to pan
      document.addEventListener('keydown', (e) => {
        if (e.key === 'Shift' && container.layout?.dragmode !== 'lasso') {
          Plotly.relayout(container, { dragmode: 'lasso' });
        }
      });
      document.addEventListener('keyup', (e) => {
        if (e.key === 'Shift' && container.layout?.dragmode === 'lasso') {
          Plotly.relayout(container, { dragmode: 'pan' });
        }
      });

      container.on('plotly_doubleclick', () => {
        const lr = container.layout;
        const xr = lr.xaxis.range, yr = lr.yaxis.range;
        if (!xr || !yr) return;
        const cx = (xr[0]+xr[1])/2, cy = (yr[0]+yr[1])/2;
        const dx = (xr[1]-xr[0])/4, dy = (yr[1]-yr[0])/4;
        Plotly.relayout(container, { 'xaxis.range': [cx-dx,cx+dx], 'yaxis.range': [cy-dy,cy+dy] });
      });
      container.on('plotly_selected', (eventData) => {
        if (onselect && eventData?.points?.length) onselect(eventData.points.map(p => p.customdata).filter(Boolean));
      });
      container.on('plotly_deselect', () => { if (onselect) onselect(null); });
      container.on('plotly_click', (eventData) => {
        if (onclick && eventData?.points?.length) {
          const id = eventData.points[0].customdata;
          if (id) onclick(id);
        }
      });
      listenersAttached = true;
    }
  }

  onMount(() => {
    if (data) plot();
    requestAnimationFrame(() => { if (container) Plotly.Plots.resize(container); });
  });

  $effect(() => {
    const _deps = [data, colorBy, sizeBy, sizeScale, mode, colorMap, coordExtents, searchMatchIds];
    if (container && data) plot();
  });
</script>

<div bind:this={container} class="w-full plotly-modebar-left flex-1 min-h-0"></div>

<style>
  .plotly-modebar-left :global(.modebar) { left: 0 !important; right: auto !important; }
</style>
