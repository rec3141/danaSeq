<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';

  let { traces = [], layout = {}, config = {}, onclick = null, onselect = null, exportName = 'danaseq_plot' } = $props();
  let container;
  let listenersAttached = false;

  const darkDefaults = {
    paper_bgcolor: '#1e293b', plot_bgcolor: '#0f172a',
    font: { color: '#94a3b8', family: 'Inter, sans-serif', size: 12 },
    margin: { t: 30, r: 20, b: 40, l: 50 },
    modebar: { orientation: 'v', bgcolor: 'transparent', color: '#64748b', activecolor: '#22d3ee' },
    xaxis: { gridcolor: '#1e293b', zerolinecolor: '#334155' },
    yaxis: { gridcolor: '#1e293b', zerolinecolor: '#334155' },
    hovermode: 'closest',
  };

  function plot() {
    if (!container || !traces.length) return;
    const mergedLayout = { ...darkDefaults, ...layout,
      xaxis: { ...darkDefaults.xaxis, ...(layout.xaxis || {}) },
      yaxis: { ...darkDefaults.yaxis, ...(layout.yaxis || {}) },
    };
    const mergedConfig = {
      responsive: true, displaylogo: false, displayModeBar: 'hover',
      toImageButtonOptions: { format: 'svg', filename: exportName },
      ...config,
    };
    Plotly.react(container, traces, mergedLayout, mergedConfig);
    // Force resize after react to prevent overflow in grid/flex containers
    requestAnimationFrame(() => { if (container) Plotly.Plots.resize(container); });

    if (!listenersAttached) {
      if (onclick) {
        container.on('plotly_click', (eventData) => {
          if (eventData?.points?.length) {
            const id = eventData.points[0].customdata;
            if (id) onclick(id);
          }
        });
      }
      if (onselect) {
        // Shift-drag for lasso select
        document.addEventListener('keydown', (e) => {
          if (e.key === 'Shift' && container._fullLayout?.dragmode !== 'lasso') {
            Plotly.relayout(container, { dragmode: 'lasso' });
          }
        });
        document.addEventListener('keyup', (e) => {
          if (e.key === 'Shift' && container._fullLayout?.dragmode === 'lasso') {
            Plotly.relayout(container, { dragmode: 'zoom' });
          }
        });
        container.on('plotly_selected', (eventData) => {
          if (eventData?.points?.length) onselect(eventData.points.map(p => p.customdata).filter(Boolean));
        });
        container.on('plotly_deselect', () => onselect(null));
      }
      listenersAttached = true;
    }
  }

  onMount(() => {
    plot();
    const ro = new ResizeObserver(() => {
      if (container) Plotly.Plots.resize(container);
    });
    ro.observe(container);
    return () => ro.disconnect();
  });

  $effect(() => {
    const _deps = [traces, layout, config];
    if (container) plot();
  });
</script>

<div bind:this={container} class="w-full min-h-0"></div>
