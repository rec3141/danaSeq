<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';

  let { data = null, colorBy = 'bin', mode = 'pca' } = $props();
  let container;

  const darkLayout = {
    paper_bgcolor: 'transparent',
    plot_bgcolor: 'transparent',
    font: { color: '#94a3b8', family: 'Inter, sans-serif', size: 12 },
    xaxis: { gridcolor: '#1e293b', zerolinecolor: '#334155', title: 'Component 1' },
    yaxis: { gridcolor: '#1e293b', zerolinecolor: '#334155', title: 'Component 2' },
    margin: { t: 30, r: 20, b: 50, l: 60 },
    legend: { bgcolor: 'transparent', font: { color: '#94a3b8', size: 10 } },
    dragmode: 'pan',
  };

  function markerSize(length) {
    return Math.max(2, Math.log10(length + 1) * 2);
  }

  function hoverText(c) {
    return `${c.id}<br>Length: ${c.length.toLocaleString()} bp<br>Depth: ${c.depth}<br>Bin: ${c.bin || 'none'}<br>Phylum: ${c.phylum || '?'}`;
  }

  function getTraces(data, colorBy, mode) {
    if (!data || !data.contigs.length) return [];

    const contigs = data.contigs;
    let xKey = 'pca_x', yKey = 'pca_y';
    if (mode === 'tsne' && data.has_tsne) { xKey = 'tsne_x'; yKey = 'tsne_y'; }
    if (mode === 'umap' && data.has_umap) { xKey = 'umap_x'; yKey = 'umap_y'; }

    const palette = ['#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa', '#fb923c',
                     '#2dd4bf', '#818cf8', '#f472b6', '#4ade80', '#e879f9', '#38bdf8',
                     '#94a3b8', '#d4d4d8', '#78716c'];

    // Continuous color modes: depth and length
    if (colorBy === 'depth' || colorBy === 'length') {
      const x = contigs.map(c => c[xKey]);
      const y = contigs.map(c => c[yKey]);
      const text = contigs.map(c => hoverText(c));
      const sizes = contigs.map(c => markerSize(c.length));

      let colors, cbarTitle;
      if (colorBy === 'depth') {
        colors = contigs.map(c => Math.log10(c.depth + 0.01));
        cbarTitle = 'log\u2081\u2080(depth)';
      } else {
        colors = contigs.map(c => Math.log10(c.length + 1));
        cbarTitle = 'log\u2081\u2080(length)';
      }

      return [{
        type: 'scattergl',
        mode: 'markers',
        x, y, text,
        marker: {
          size: sizes,
          color: colors,
          colorscale: 'Viridis',
          showscale: true,
          colorbar: {
            title: { text: cbarTitle, font: { color: '#94a3b8' } },
            tickfont: { color: '#94a3b8' },
          },
          line: { width: 0 },
        },
        hoverinfo: 'text',
      }];
    }

    // Categorical color modes
    const categoryField = {
      'bin': 'bin', 'phylum': 'phylum', 'domain': 'domain',
      'kaiju_phylum': 'kaiju_phylum', 'kraken2_phylum': 'kraken2_phylum', 'rrna_phylum': 'rrna_phylum',
    };
    const groups = {};
    for (const c of contigs) {
      let key;
      const field = categoryField[colorBy];
      if (colorBy === 'bin') key = c.bin || 'unbinned';
      else if (field) key = c[field] || 'Unknown';
      else key = 'all';

      if (!groups[key]) groups[key] = { x: [], y: [], text: [], sizes: [] };
      groups[key].x.push(c[xKey]);
      groups[key].y.push(c[yKey]);
      groups[key].text.push(hoverText(c));
      groups[key].sizes.push(markerSize(c.length));
    }

    // Stable color map: assign colors alphabetically so each group always gets the same color
    const allNames = Object.keys(groups).filter(n => n !== 'unbinned' && n !== 'Unknown').sort();
    const colorMap = {};
    allNames.forEach((name, i) => { colorMap[name] = palette[i % palette.length]; });

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
        marker: {
          size: g.sizes,
          color: isBackground ? '#475569' : colorMap[name],
          opacity: isBackground ? 0.25 : 0.75,
          line: { width: 0 },
        },
        hoverinfo: 'text',
      };
    });
  }

  onMount(() => {
    if (data) plot();
  });

  $effect(() => {
    if (container && data) plot();
  });

  function plot() {
    const traces = getTraces(data, colorBy, mode);
    const layout = {
      ...darkLayout,
      xaxis: { ...darkLayout.xaxis, title: mode === 'tsne' ? 't-SNE 1' : mode === 'umap' ? 'UMAP 1' : 'PC 1', uirevision: mode, constrain: 'domain' },
      yaxis: { ...darkLayout.yaxis, title: mode === 'tsne' ? 't-SNE 2' : mode === 'umap' ? 'UMAP 2' : 'PC 2', uirevision: mode, scaleanchor: 'x', scaleratio: 1, constrain: 'domain' },
      uirevision: mode + '-' + colorBy,
    };
    Plotly.react(container, traces, layout, {
      responsive: true,
      displaylogo: false,
      scrollZoom: true,
      displayModeBar: 'hover',
      modeBarButtonsToRemove: [
        'autoScale2d', 'toggleSpikelines',
        'hoverCompareCartesian', 'hoverClosestCartesian',
      ],
    });
  }
</script>

<div bind:this={container} class="w-full" style="min-height: 500px;"></div>
