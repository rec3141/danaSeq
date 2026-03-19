<script>
  import DataTable from '../components/ui/DataTable.svelte';
  import { onMount } from 'svelte';
  import { ecosystemServices, loadEcosystemServices } from '../stores/data.js';

  let esData = $derived($ecosystemServices);

  onMount(() => { loadEcosystemServices(); });

  // === ES Summary Cards ===
  let esSummary = $derived.by(() => {
    if (!esData?.heatmap) return [];
    const { es_codes, es_names, matrix, mags } = esData.heatmap;
    return es_codes.map((code, i) => {
      const scores = matrix.map(row => row[i]).filter(v => v !== 0);
      const total = scores.reduce((a, b) => a + b, 0);
      const mean = scores.length ? total / scores.length : 0;
      const max = scores.length ? Math.max(...scores) : 0;
      const nMags = scores.length;
      return {
        code, name: es_names[code] || code,
        total: total.toFixed(1), mean: mean.toFixed(3),
        max: max.toFixed(3), nMags
      };
    }).sort((a, b) => b.total - a.total);
  });

  // === SDG Goal Scores ===
  let sdgGoals = $derived.by(() => {
    if (!esData?.sdg?.goals) return [];
    // Aggregate by goal across all entities
    const goals = {};
    for (const g of esData.sdg.goals) {
      const key = g.sdg_goal;
      if (!goals[key]) goals[key] = { goal: key, name: g.goal_name, score: 0 };
      goals[key].score += g.score;
    }
    return Object.values(goals).sort((a, b) => b.score - a.score);
  });

  let maxSdgScore = $derived(sdgGoals.length ? Math.max(...sdgGoals.map(g => g.score)) : 1);

  // SDG colors (official UN colors)
  const sdgColors = {
    1: '#E5243B', 2: '#DDA63A', 3: '#4C9F38', 6: '#26BDE2',
    7: '#FCC30B', 8: '#A21942', 9: '#FD6925', 11: '#FD9D24',
    12: '#BF8B2E', 13: '#3F7E44', 14: '#0A97D9', 15: '#56C02B'
  };

  // === Heatmap data ===
  let heatmapData = $derived.by(() => {
    if (!esData?.heatmap) return null;
    const { es_codes, es_names, matrix, mags } = esData.heatmap;
    // Filter to MAGs with any non-zero score
    const activeRows = [];
    const activeLabels = [];
    for (let i = 0; i < mags.length; i++) {
      if (matrix[i].some(v => v !== 0)) {
        activeRows.push(matrix[i]);
        activeLabels.push(mags[i]);
      }
    }
    return {
      rows: activeRows,
      rowLabels: activeLabels,
      colLabels: es_codes.map(c => es_names[c] || c),
      colCodes: es_codes
    };
  });

  // Color scale for heatmap
  function scoreColor(val, max) {
    if (val === 0) return 'rgb(15, 23, 42)'; // slate-900
    const t = Math.min(val / (max || 1), 1);
    // cyan gradient
    const r = Math.round(15 + t * (6 - 15));
    const g = Math.round(23 + t * (182 - 23));
    const b = Math.round(42 + t * (212 - 42));
    return `rgb(${r}, ${g}, ${b})`;
  }

  let heatmapMax = $derived.by(() => {
    if (!heatmapData) return 1;
    return Math.max(...heatmapData.rows.flat().filter(v => v > 0), 0.001);
  });

  // === Catalog table ===
  let catalogRows = $derived.by(() => {
    if (!esData?.catalog) return [];
    return esData.catalog.slice(0, 500).map(row => ({
      protein_id: row.protein_id,
      contig_id: row.contig_id,
      gene_id: row.gene_id,
      gene_type: row.gene_id_type,
      es_code: row.es_code,
      es_name: row.es_name,
      confidence: row.confidence,
      role: row.functional_role,
      method: row.detection_method
    }));
  });

  // ES category colors
  const esColors = {
    '2.1.1.1': '#f59e0b', // amber - waste decomposition
    '2.1.1.2': '#a855f7', // purple - pollutant sequestration
    '2.3.3.2': '#ef4444', // red - disease control
    '2.3.4.2': '#84cc16', // lime - soil OM
    '2.3.5.1': '#06b6d4', // cyan - freshwater quality
    '2.3.5.2': '#3b82f6', // blue - marine quality
    '2.3.6.1': '#10b981', // emerald - climate regulation
    '2.3.6.2': '#8b5cf6', // violet - air quality
  };
</script>

<div class="space-y-6">
  <h2 class="text-xl font-semibold text-slate-200">Ecosystem Services</h2>
  <p class="text-sm text-slate-400">
    Gene-to-ecosystem-service mapping via ECOSSDB ({esData ? `${esData.catalog?.length || 0} gene hits` : 'loading...'})
  </p>

  {#if !esData}
    <div class="text-slate-500 text-center py-12">
      <p>No ecosystem services data available.</p>
      <p class="text-xs mt-2">Run the pipeline with <code class="bg-slate-800 px-1 rounded">--run_ecossdb true</code></p>
    </div>
  {:else}
    <!-- ES Summary Cards -->
    <div class="grid grid-cols-2 md:grid-cols-4 gap-3">
      {#each esSummary as es}
        <div class="bg-slate-800/50 border border-slate-700 rounded-lg p-3"
             style="border-left: 3px solid {esColors[es.code] || '#64748b'}">
          <div class="text-xs text-slate-400 truncate" title={es.name}>{es.name}</div>
          <div class="text-lg font-semibold text-slate-200 mt-1">{es.total}</div>
          <div class="text-xs text-slate-500 mt-0.5">{es.nMags} MAGs · mean {es.mean}</div>
        </div>
      {/each}
    </div>

    <!-- SDG Goals -->
    {#if sdgGoals.length > 0}
      <div class="bg-slate-800/30 border border-slate-700 rounded-lg p-4">
        <h3 class="text-sm font-semibold text-slate-300 mb-3">UN Sustainable Development Goals</h3>
        <div class="space-y-1.5">
          {#each sdgGoals as goal}
            <div class="flex items-center gap-2">
              <div class="w-6 h-6 rounded flex items-center justify-center text-xs font-bold text-white shrink-0"
                   style="background-color: {sdgColors[goal.goal] || '#64748b'}">
                {goal.goal}
              </div>
              <div class="text-xs text-slate-400 w-40 truncate shrink-0" title={goal.name}>{goal.name}</div>
              <div class="flex-1 h-5 bg-slate-900 rounded overflow-hidden">
                <div class="h-full rounded transition-all"
                     style="width: {(goal.score / maxSdgScore * 100).toFixed(1)}%; background-color: {sdgColors[goal.goal] || '#64748b'}; opacity: 0.7">
                </div>
              </div>
              <div class="text-xs text-slate-400 w-16 text-right">{goal.score.toFixed(1)}</div>
            </div>
          {/each}
        </div>
      </div>
    {/if}

    <!-- MAG × ES Heatmap -->
    {#if heatmapData && heatmapData.rows.length > 0}
      <div class="bg-slate-800/30 border border-slate-700 rounded-lg p-4">
        <h3 class="text-sm font-semibold text-slate-300 mb-3">
          MAG × Ecosystem Service Heatmap
          <span class="text-slate-500 font-normal ml-2">({heatmapData.rowLabels.length} MAGs × {heatmapData.colLabels.length} ES categories)</span>
        </h3>
        <div class="overflow-x-auto">
          <table class="text-xs">
            <thead>
              <tr>
                <th class="text-left text-slate-500 px-1 py-0.5 sticky left-0 bg-slate-900/90">MAG</th>
                {#each heatmapData.colLabels as label, i}
                  <th class="text-center text-slate-500 px-1 py-0.5 max-w-20 truncate"
                      title={label}
                      style="border-bottom: 2px solid {esColors[heatmapData.colCodes[i]] || '#64748b'}">
                    {label.split(' ').slice(0, 3).join(' ')}
                  </th>
                {/each}
              </tr>
            </thead>
            <tbody>
              {#each heatmapData.rowLabels.slice(0, 50) as mag, i}
                <tr>
                  <td class="text-slate-400 px-1 py-0.5 sticky left-0 bg-slate-900/90 truncate max-w-32" title={mag}>
                    {mag}
                  </td>
                  {#each heatmapData.rows[i] as val}
                    <td class="px-1 py-0.5 text-center"
                        style="background-color: {scoreColor(val, heatmapMax)}"
                        title={val.toFixed(3)}>
                      {val > 0 ? val.toFixed(2) : ''}
                    </td>
                  {/each}
                </tr>
              {/each}
            </tbody>
          </table>
          {#if heatmapData.rowLabels.length > 50}
            <p class="text-xs text-slate-500 mt-2">Showing 50 of {heatmapData.rowLabels.length} MAGs</p>
          {/if}
        </div>
      </div>
    {/if}

    <!-- Gene Catalog Detail Table -->
    {#if catalogRows.length > 0}
      <div class="bg-slate-800/30 border border-slate-700 rounded-lg p-4">
        <h3 class="text-sm font-semibold text-slate-300 mb-3">Gene Catalog (top 500)</h3>
        <div class="overflow-x-auto max-h-96 overflow-y-auto">
          <table class="w-full text-xs">
            <thead class="sticky top-0 bg-slate-900">
              <tr class="text-slate-400 border-b border-slate-700">
                <th class="text-left px-2 py-1">Protein</th>
                <th class="text-left px-2 py-1">Gene ID</th>
                <th class="text-left px-2 py-1">Type</th>
                <th class="text-left px-2 py-1">ES Code</th>
                <th class="text-left px-2 py-1">ES Name</th>
                <th class="text-center px-2 py-1">Conf.</th>
                <th class="text-left px-2 py-1">Role</th>
              </tr>
            </thead>
            <tbody>
              {#each catalogRows as row}
                <tr class="border-b border-slate-800 hover:bg-slate-800/50">
                  <td class="px-2 py-0.5 text-slate-300 font-mono">{row.protein_id}</td>
                  <td class="px-2 py-0.5 text-cyan-400">{row.gene_id}</td>
                  <td class="px-2 py-0.5 text-slate-500">{row.gene_type}</td>
                  <td class="px-2 py-0.5" style="color: {esColors[row.es_code] || '#94a3b8'}">{row.es_code}</td>
                  <td class="px-2 py-0.5 text-slate-400 truncate max-w-48">{row.es_name}</td>
                  <td class="px-2 py-0.5 text-center text-slate-300">{row.confidence}</td>
                  <td class="px-2 py-0.5 text-slate-400">
                    {#if row.role === 'producer'}
                      <span class="text-emerald-400">▲ {row.role}</span>
                    {:else if row.role === 'inhibitor'}
                      <span class="text-rose-400">▼ {row.role}</span>
                    {:else if row.role === 'transformer'}
                      <span class="text-amber-400">◆ {row.role}</span>
                    {:else}
                      <span class="text-sky-400">○ {row.role}</span>
                    {/if}
                  </td>
                </tr>
              {/each}
            </tbody>
          </table>
        </div>
      </div>
    {/if}
  {/if}
</div>
