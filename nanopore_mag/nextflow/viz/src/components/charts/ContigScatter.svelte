<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';

  let { data = null, colorBy = 'bin', sizeBy = 'length', sizeScale = 1.0, mode = 'pca', colorMap = {}, onselect = null, onclick = null } = $props();
  let container;
  let listenersAttached = false;

  const darkLayout = {
    paper_bgcolor: '#1e293b',
    plot_bgcolor: '#0f172a',
    font: { color: '#94a3b8', family: 'Inter, sans-serif', size: 12 },
    xaxis: { showgrid: false, zeroline: false, showticklabels: false, title: '' },
    yaxis: { showgrid: false, zeroline: false, showticklabels: false, title: '' },
    margin: { t: 10, r: 10, b: 10, l: 10 },
    modebar: { orientation: 'v', bgcolor: 'transparent', color: '#64748b', activecolor: '#22d3ee' },
    legend: { bgcolor: 'rgba(15,23,42,0.7)', font: { color: '#94a3b8', size: 10 },
              x: 1, y: 1, xanchor: 'right', yanchor: 'top' },
    dragmode: 'pan',
  };

  function markerSize(c, sizeBy, scale) {
    let s;
    if (sizeBy === 'depth') s = Math.max(1.5, Math.sqrt(c.depth) * 0.8);
    else if (sizeBy === 'fixed') s = 3;
    else s = Math.max(1.5, Math.pow(c.length, 0.3) * 0.25);
    // length: 1KB→2px, 10KB→4px, 100KB→8px, 500KB→14px, 1MB→20px
    return s * scale;
  }

  const BIN_LABELS = { bin: 'DAS Tool', semibin_bin: 'SemiBin2', metabat_bin: 'MetaBAT2', maxbin_bin: 'MaxBin2', lorbin_bin: 'LorBin', comebin_bin: 'COMEBin' };

  function hoverText(c, colorBy) {
    let tax = c.kaiju_phylum || c.kraken2_phylum || c.rrna_phylum || '?';
    // If coloring by a taxonomy field, show that field's value
    if (colorBy && !colorBy.endsWith('_bin') && colorBy !== 'bin' && c[colorBy] !== undefined) {
      tax = c[colorBy] || '?';
    }
    const isBin = colorBy === 'bin' || colorBy.endsWith('_bin');
    const binLine = isBin ? `<br>${BIN_LABELS[colorBy] || colorBy}: ${c[colorBy] || 'none'}` : '';
    return `${c.id}<br>Length: ${c.length.toLocaleString()} bp<br>Depth: ${c.depth}<br>GC: ${c.gc ?? '?'}%${binLine}<br>Tax: ${tax}`;
  }

  function getTraces(data, colorBy, sizeBy, sizeScale, mode) {
    if (!data || !data.contigs.length) return [];

    const contigs = data.contigs;
    let xKey = 'pca_x', yKey = 'pca_y';
    if (mode === 'tsne' && data.has_tsne) { xKey = 'tsne_x'; yKey = 'tsne_y'; }
    if (mode === 'umap' && data.has_umap) { xKey = 'umap_x'; yKey = 'umap_y'; }

    const palette = ['#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa', '#fb923c',
                     '#2dd4bf', '#818cf8', '#f472b6', '#4ade80', '#e879f9', '#38bdf8',
                     '#94a3b8', '#d4d4d8', '#78716c'];

    // Continuous color modes
    const continuousModes = {
      depth:  { fn: c => Math.log10(c.depth + 0.01), title: 'log\u2081\u2080(depth)', scale: 'Viridis' },
      length: { fn: c => Math.log10(c.length + 1),   title: 'log\u2081\u2080(length)', scale: 'Viridis' },
      gc:     { fn: c => c.gc,                        title: 'GC%',                     scale: 'Viridis' },
    };
    if (continuousModes[colorBy]) {
      const cm = continuousModes[colorBy];
      const x = contigs.map(c => c[xKey]);
      const y = contigs.map(c => c[yKey]);
      const text = contigs.map(c => hoverText(c, colorBy));
      const customdata = contigs.map(c => c.id);
      const sizes = contigs.map(c => markerSize(c, sizeBy, sizeScale));
      const colors = contigs.map(cm.fn);

      const markerOpts = {
        size: sizes,
        color: colors,
        colorscale: cm.scale,
        showscale: true,
        colorbar: {
          title: { text: cm.title, font: { color: '#94a3b8' } },
          tickfont: { color: '#94a3b8' },
          x: 0.99, xanchor: 'right', xpad: 0,
          y: 0.98, yanchor: 'top', ypad: 0,
          len: 0.4, thickness: 12,
          bgcolor: 'rgba(15,23,42,0.7)',
          borderwidth: 0,
        },
        line: { width: 0 },
      };
      return [{
        type: 'scattergl',
        mode: 'markers',
        x, y, text, customdata,
        marker: markerOpts,
        hoverinfo: 'text',
      }];
    }

    // Categorical color modes: colorBy matches the JSON field name directly
    const isBinField = colorBy === 'bin' || colorBy.endsWith('_bin');
    const groups = {};
    for (const c of contigs) {
      let key = c[colorBy];
      if (isBinField) key = key || 'unbinned';
      else key = key || 'Unknown';

      if (!groups[key]) groups[key] = { x: [], y: [], text: [], sizes: [], ids: [] };
      groups[key].x.push(c[xKey]);
      groups[key].y.push(c[yKey]);
      groups[key].text.push(hoverText(c, colorBy));
      groups[key].sizes.push(markerSize(c, sizeBy, sizeScale));
      groups[key].ids.push(c.id);
    }

    // Use stable colorMap from parent (derived from full unfiltered dataset) if available,
    // otherwise fall back to computing from visible groups
    const useParentMap = colorMap && Object.keys(colorMap).length > 0;
    const effectiveColorMap = useParentMap ? colorMap : (() => {
      const allNames = Object.keys(groups).filter(n => n !== 'unbinned' && n !== 'Unknown').sort();
      const m = {};
      allNames.forEach((name, i) => { m[name] = palette[i % palette.length]; });
      return m;
    })();

    // Render order: unbinned/Unknown first (bottom layer), then smallest groups on top
    const entries = Object.entries(groups).sort((a, b) => {
      if (a[0] === 'unbinned' || a[0] === 'Unknown') return -1;
      if (b[0] === 'unbinned' || b[0] === 'Unknown') return 1;
      return a[1].x.length - b[1].x.length;
    });

    return entries.map(([name, g]) => {
      const isBackground = name === 'unbinned' || name === 'Unknown';
      return {
        type: 'scattergl',
        mode: 'markers',
        name: name.length > 20 ? name.slice(0, 18) + '..' : name,
        x: g.x,
        y: g.y,
        text: g.text,
        customdata: g.ids,
        marker: {
          size: g.sizes,
          color: isBackground ? '#475569' : (effectiveColorMap[name] || '#94a3b8'),
          opacity: isBackground ? 0.25 : 0.75,
          line: { width: 0 },
        },
        hoverinfo: 'text',
      };
    });
  }

  onMount(() => {
    if (data) plot();
    // Plotly needs a resize after layout settles to fit its container
    requestAnimationFrame(() => {
      if (container) Plotly.Plots.resize(container);
    });
  });

  $effect(() => {
    // Read all reactive props so Svelte tracks them as dependencies
    const _deps = [data, colorBy, sizeBy, sizeScale, mode, colorMap];
    if (container && data) plot();
  });

  function plot() {
    const traces = getTraces(data, colorBy, sizeBy, sizeScale, mode);
    const layout = {
      ...darkLayout,
      xaxis: { ...darkLayout.xaxis, uirevision: mode },
      yaxis: { ...darkLayout.yaxis, uirevision: mode, scaleanchor: 'x', scaleratio: 1 },
      uirevision: mode + '-' + colorBy,
    };
    Plotly.react(container, traces, layout, {
      responsive: true,
      displaylogo: false,
      scrollZoom: true,
      toImageButtonOptions: { scale: 4 },
      doubleClick: false,
      displayModeBar: 'hover',
      modeBarButtonsToRemove: [
        'autoScale2d', 'toggleSpikelines',
        'hoverCompareCartesian', 'hoverClosestCartesian',
        'zoomIn2d', 'zoomOut2d',
        'zoom2d', 'select2d',
      ],
    });

    if (!listenersAttached) {
      // Double-click zooms in 2x at click point
      container.on('plotly_doubleclick', () => {
        // Get current axis ranges from the live layout
        const layout = container.layout;
        const xr = layout.xaxis.range;
        const yr = layout.yaxis.range;
        if (!xr || !yr) return;
        const cx = (xr[0] + xr[1]) / 2;
        const cy = (yr[0] + yr[1]) / 2;
        const dx = (xr[1] - xr[0]) / 4;
        const dy = (yr[1] - yr[0]) / 4;
        Plotly.relayout(container, {
          'xaxis.range': [cx - dx, cx + dx],
          'yaxis.range': [cy - dy, cy + dy],
        });
      });

    // Attach selection listeners once
      container.on('plotly_selected', (eventData) => {
        if (onselect && eventData?.points?.length) {
          const ids = eventData.points.map(p => p.customdata).filter(Boolean);
          onselect(ids);
        }
      });
      container.on('plotly_deselect', () => {
        if (onselect) onselect(null);
      });
      container.on('plotly_click', (eventData) => {
        if (onclick && eventData?.points?.length) {
          const id = eventData.points[0].customdata;
          if (id) onclick(id);
        }
      });
      listenersAttached = true;
    }
  }
</script>

<div bind:this={container} class="w-full plotly-modebar-left flex-1 min-h-0"></div>

<style>
  .plotly-modebar-left :global(.modebar) {
    left: 0 !important;
    right: auto !important;
  }
</style>
