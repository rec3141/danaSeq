<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';

  let { traces = [], layout = {}, config = {} } = $props();
  let container;

  const darkDefaults = {
    paper_bgcolor: '#1e293b', plot_bgcolor: '#0f172a',
    font: { color: '#94a3b8', family: 'Inter, sans-serif', size: 12 },
    margin: { t: 30, r: 20, b: 40, l: 50 },
    modebar: { orientation: 'v', bgcolor: 'transparent', color: '#64748b', activecolor: '#22d3ee' },
    xaxis: { gridcolor: '#1e293b', zerolinecolor: '#334155' },
    yaxis: { gridcolor: '#1e293b', zerolinecolor: '#334155' },
  };

  function plot() {
    if (!container || !traces.length) return;
    const mergedLayout = { ...darkDefaults, ...layout,
      xaxis: { ...darkDefaults.xaxis, ...(layout.xaxis || {}) },
      yaxis: { ...darkDefaults.yaxis, ...(layout.yaxis || {}) },
    };
    const mergedConfig = { responsive: true, displaylogo: false, displayModeBar: 'hover', ...config };
    Plotly.react(container, traces, mergedLayout, mergedConfig);
  }

  onMount(() => {
    plot();
    requestAnimationFrame(() => { if (container) Plotly.Plots.resize(container); });
  });

  $effect(() => {
    const _deps = [traces, layout, config];
    if (container) plot();
  });
</script>

<div bind:this={container} class="w-full min-h-0"></div>
