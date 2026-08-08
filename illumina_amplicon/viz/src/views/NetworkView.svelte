<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import {
    store, countsBySample,
    GROUP_HEX,
    buildTaxColorMap, getAsvColor, getEffectiveColorLevel, getClusterColor,
    makeTaxonMatcher,
  } from '../stores/data.svelte.js';

  let { filters = {} } = $props();

  let plotDiv;
  let hasPlot = false;

  let taxonMatches = $derived.by(() =>
    makeTaxonMatcher(filters.taxonFilter, filters.taxonFilterLevel)
  );

  let filteredAsvs = $derived.by(() => {
    const matches = taxonMatches;
    const gf = filters.groupFlags || {};
    return store.asvs.filter(a => {
      if ((a.n_samples ?? 0) < (filters.minPrevalence || 0)) return false;
      if ((a.total_reads ?? 0) < (filters.minReads || 0)) return false;
      const group = a.group ?? 'unknown';
      if (gf[group] === false) return false;
      if (matches && !matches(a.id ?? '', a.taxonomy ?? '')) return false;
      return true;
    });
  });

  let selectedAsvObj = $derived(
    store.selectedAsv != null ? store.asvs[store.selectedAsv] : null
  );

  $effect(() => {
    if (!plotDiv || filteredAsvs.length === 0) return;

    const colorLevel = filters.colorMode === 'group' ? 'group' : getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter, filters.taxonFilterLevel);
    const cmap = colorLevel !== 'group' ? buildTaxColorMap(colorLevel, filters.taxonFilter, filters.taxonFilterLevel).colorMap : null;

    const colors = filteredAsvs.map(a => {
      if (filters.colorMode === 'cluster') return getClusterColor(a.id, 'asvCluster', filters.asvClusterK);
      if (cmap) return getAsvColor(a.id, colorLevel, cmap);
      return GROUP_HEX[a.group ?? 'prokaryote'] ?? GROUP_HEX.unknown;
    });

    const trace = {
      x: filteredAsvs.map(a => a.x ?? 0),
      y: filteredAsvs.map(a => a.y ?? 0),
      mode: 'markers',
      type: 'scattergl',
      marker: {
        size: (() => {
          const s = filters.networkPointScale ?? 10;
          const nSamples = store.samples.length || 1;
          return filteredAsvs.map(a => {
            const meanRA = (a.total_reads ?? 0) / nSamples;
            return Math.min(60, Math.sqrt(meanRA) * s * 0.3);
          });
        })(),
        color: colors,
        opacity: 0.7,
        line: { width: 0 },
      },
      text: filteredAsvs.map(a =>
        `${a.id}<br>${a.taxonomy ?? ''}<br>${(a.total_reads ?? 0).toLocaleString()} reads<br>${a.n_samples ?? 0} samples`
      ),
      hoverinfo: 'text',
      showlegend: false,
    };

    const layout = {
      dragmode: 'pan',
      xaxis: { title: '', zeroline: false, showgrid: false, showticklabels: false },
      yaxis: { title: '', zeroline: false, showgrid: false, showticklabels: false, scaleanchor: 'x' },
      plot_bgcolor: 'rgba(2, 6, 15, 1)',
      paper_bgcolor: 'rgba(2, 6, 15, 1)',
      font: { color: '#94a3b8' },
      margin: { l: 20, r: 20, t: 10, b: 20 },
      title: { text: `${filteredAsvs.length} ASVs`, font: { size: 12, color: '#64748b' }, x: 0.01, y: 0.99 },
    };

    const config = { scrollZoom: true, displayModeBar: false, doubleClick: 'reset+autosize' };

    if (!hasPlot) {
      Plotly.newPlot(plotDiv, [trace], layout, config);
      hasPlot = true;

      plotDiv.on('plotly_click', (data) => {
        if (data.points?.[0]) {
          const idx = data.points[0].pointNumber;
          const oi = store.asvs.indexOf(filteredAsvs[idx]);
          store.selectedAsv = oi >= 0 ? oi : null;
        }
      });

    } else {
      const curLayout = plotDiv.layout;
      if (curLayout?.xaxis?.range) layout.xaxis.range = curLayout.xaxis.range;
      if (curLayout?.yaxis?.range) layout.yaxis.range = curLayout.yaxis.range;
      Plotly.react(plotDiv, [trace], layout, config);
    }
  });

  function handleKey(e) {
    if (!plotDiv || !hasPlot) return;
    if (e.key === 'Shift') {
      Plotly.relayout(plotDiv, { dragmode: e.type === 'keydown' ? 'lasso' : 'pan' });
    } else if (e.key === 'Escape' && e.type === 'keydown') {
      store.selectedAsv = null;
      Plotly.restyle(plotDiv, { selectedpoints: [null] });
    }
  }

  onMount(() => {
    document.addEventListener('keydown', handleKey);
    document.addEventListener('keyup', handleKey);
    return () => {
      document.removeEventListener('keydown', handleKey);
      document.removeEventListener('keyup', handleKey);
      if (plotDiv && hasPlot) Plotly.purge(plotDiv);
    };
  });
</script>

<div class="flex h-full flex-col">
  <div class="flex-1 relative">
    <div bind:this={plotDiv} class="absolute inset-0"></div>
  </div>

  {#if selectedAsvObj}
    <div class="border-t border-slate-800 bg-slate-900/80 p-4">
      <div class="flex items-center justify-between">
        <h3 class="text-sm font-semibold text-slate-200">
          {selectedAsvObj.id ?? 'ASV'} &mdash; {selectedAsvObj.group ?? ''}
        </h3>
        <button class="text-xs text-slate-500 hover:text-slate-300" onclick={() => store.selectedAsv = null}>Close</button>
      </div>
      <p class="mt-1 text-xs text-slate-400">{selectedAsvObj.taxonomy ?? 'No taxonomy'}</p>
      <p class="text-xs text-slate-500">
        {(selectedAsvObj.total_reads ?? 0).toLocaleString()} total reads |
        Prevalence: {selectedAsvObj.n_samples ?? 0}
      </p>
    </div>
  {/if}
</div>
