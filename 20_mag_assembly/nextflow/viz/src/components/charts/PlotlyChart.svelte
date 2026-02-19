<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';

  let { data = [], layout = {}, config = {}, onclick = null } = $props();
  let container;

  const darkLayout = {
    paper_bgcolor: '#1e293b',
    plot_bgcolor: '#1e293b',
    font: { color: '#94a3b8', family: 'Inter, sans-serif', size: 12 },
    xaxis: { gridcolor: '#1e293b', zerolinecolor: '#334155' },
    yaxis: { gridcolor: '#1e293b', zerolinecolor: '#334155' },
    margin: { t: 40, r: 20, b: 50, l: 60 },
    legend: { bgcolor: 'transparent', font: { color: '#94a3b8' } },
    dragmode: 'pan',
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
    modeBarButtonsToAdd: [],
  };

  function mergeDeep(target, source) {
    const result = { ...target };
    for (const key of Object.keys(source)) {
      if (source[key] && typeof source[key] === 'object' && !Array.isArray(source[key])) {
        result[key] = mergeDeep(result[key] || {}, source[key]);
      } else {
        result[key] = source[key];
      }
    }
    return result;
  }

  onMount(() => {
    const merged = mergeDeep(darkLayout, layout);
    const cfg = { ...defaultConfig, ...config };
    Plotly.newPlot(container, data, merged, cfg);

    if (onclick) {
      container.on('plotly_click', onclick);
    }

    return () => {
      Plotly.purge(container);
    };
  });

  $effect(() => {
    if (container && data) {
      const merged = mergeDeep(darkLayout, layout);
      const cfg = { ...defaultConfig, ...config };
      Plotly.react(container, data, merged, cfg);
    }
  });
</script>

<div bind:this={container} class="w-full" style="min-height: 350px;"></div>
