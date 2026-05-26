<script>
  /**
   * Shared Plotly wrapper matching the patterns used in mag_analysis and
   * nanopore_live: dark theme, scroll-zoom, drag-pan, mode-bar cleaned up.
   * Pass `data`, `layout`, `config` to override.
   */
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';

  let { data = [], layout = {}, config = {}, onclick = null, height = '500px' } = $props();
  let container;

  const darkLayout = {
    paper_bgcolor: '#0f172a',
    plot_bgcolor: '#0f172a',
    font: { color: '#cbd5e1', family: 'Inter, sans-serif', size: 11 },
    xaxis: { gridcolor: '#1e293b', zerolinecolor: '#334155' },
    yaxis: { gridcolor: '#1e293b', zerolinecolor: '#334155' },
    margin: { t: 30, r: 30, b: 50, l: 60 },
    legend: { bgcolor: 'transparent', font: { color: '#cbd5e1' } },
    dragmode: 'pan',
    hovermode: 'closest',
  };

  const defaultConfig = {
    responsive: true,
    displaylogo: false,
    scrollZoom: true,
    displayModeBar: 'hover',
    modeBarButtonsToRemove: [
      'select2d', 'lasso2d', 'autoScale2d',
      'toggleSpikelines', 'hoverCompareCartesian', 'hoverClosestCartesian',
    ],
  };

  function mergeDeep(target, source) {
    const result = { ...target };
    for (const key of Object.keys(source)) {
      const sv = source[key];
      if (sv && typeof sv === 'object' && !Array.isArray(sv)) {
        result[key] = mergeDeep(result[key] || {}, sv);
      } else {
        result[key] = sv;
      }
    }
    return result;
  }

  onMount(() => {
    const merged = mergeDeep(darkLayout, layout);
    const cfg = { ...defaultConfig, ...config };
    Plotly.newPlot(container, data, merged, cfg);
    if (onclick) container.on('plotly_click', onclick);
    return () => Plotly.purge(container);
  });

  $effect(() => {
    if (container && data) {
      const merged = mergeDeep(darkLayout, layout);
      const cfg = { ...defaultConfig, ...config };
      Plotly.react(container, data, merged, cfg);
    }
  });
</script>

<div bind:this={container} class="w-full" style="height: {height};"></div>
