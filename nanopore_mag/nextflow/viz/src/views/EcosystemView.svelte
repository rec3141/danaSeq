<script>
  import DataTable from '../components/ui/DataTable.svelte';
  import D3Heatmap from '../components/charts/D3Heatmap.svelte';
  import { onMount } from 'svelte';
  import { ecosystemServices, loadEcosystemServices, mags as magsStore, binQuality, loadBinQuality } from '../stores/data.js';
  import { selectedMag } from '../stores/selection.js';
  import * as d3 from 'd3';

  let esData = $derived($ecosystemServices);
  let magsData = $derived($magsStore);
  let selected = $derived($selectedMag);

  onMount(() => { loadEcosystemServices(); loadBinQuality(); });

  // === Abundance mode toggle ===
  let abundanceMode = $state('even'); // 'even' | 'relative'

  // === Binner filter (MAG heatmap) ===
  const binnerDefs = [
    { key: 'dastool', label: 'DAS Tool' },
    { key: 'binette', label: 'Binette' },
    { key: 'magscot', label: 'MAGScoT' },
    { key: 'comebin', label: 'COMEBin' },
    { key: 'metabat', label: 'MetaBAT2' },
    { key: 'maxbin',  label: 'MaxBin2' },
    { key: 'lorbin',  label: 'LorBin' },
    { key: 'semibin', label: 'SemiBin2' },
    { key: 'vamb',    label: 'VAMB' },
    { key: 'vamb_tax', label: 'VAMB-tax' },
  ];
  let activeBinners = $state(new Set(['dastool']));

  function toggleBinner(key) {
    const next = new Set(activeBinners);
    if (next.has(key)) {
      if (next.size > 1) next.delete(key);
    } else {
      next.add(key);
    }
    activeBinners = next;
  }

  let binnerCounts = $derived.by(() => {
    const pw = esData?.pathway_heatmap;
    if (pw?.binners) {
      const counts = {};
      for (const b of pw.binners) counts[b] = (counts[b] || 0) + 1;
      return counts;
    }
    if (!magsData) return {};
    const counts = {};
    for (const m of Object.values(magsData)) {
      const bk = m.binner || 'unknown';
      counts[bk] = (counts[bk] || 0) + 1;
    }
    return counts;
  });

  // === Sort definitions ===
  const magSortDefs = [
    { key: 'clustered', label: 'Clustered' },
    { key: 'name', label: 'A-Z' },
    { key: 'total_score', label: 'Total Score' },
    { key: 'completeness', label: 'Completeness' },
    { key: 'contamination', label: 'Contamination' },
    { key: 'size', label: 'Genome Size' },
    { key: 'binner', label: 'Source Binner' },
  ];
  let magSortBy = $state('clustered');

  function cycleMagSort() {
    const keys = magSortDefs.map(d => d.key);
    magSortBy = keys[(keys.indexOf(magSortBy) + 1) % keys.length];
  }

  const sampleSortDefs = [
    { key: 'clustered', label: 'Clustered' },
    { key: 'name', label: 'A-Z' },
    { key: 'total_score', label: 'Total Score' },
    { key: 'top_es', label: 'Top ES' },
    { key: 'diversity', label: 'ES Diversity' },
  ];
  let sampleSortBy = $state('clustered');

  function cycleSampleSort() {
    const keys = sampleSortDefs.map(d => d.key);
    sampleSortBy = keys[(keys.indexOf(sampleSortBy) + 1) % keys.length];
  }

  // === Bin metadata lookup (all binners via checkm2_all) ===
  let allBins = $derived($binQuality);

  let magMeta = $derived.by(() => {
    const map = {};
    // Primary source: checkm2_all (covers all binners)
    if (allBins) {
      for (const b of allBins) map[b.name] = b;
    }
    // Fallback: mags.json (DAS Tool only, has coverage)
    if (magsData) {
      for (const m of Object.values(magsData)) {
        if (!map[m.name]) map[m.name] = m;
        else map[m.name].coverage = m.coverage; // merge coverage into checkm2 data
      }
    }
    return map;
  });

  // Quality tier helper
  function qualityTier(m) {
    if (m.completeness >= 90 && m.contamination < 5) return 'HQ';
    if (m.completeness >= 50 && m.contamination < 10) return 'MQ';
    return 'LQ';
  }

  // === Per-sample pathway matrix (relative abundance mode) ===
  // Uses pathway_heatmap if available, otherwise falls back to basic heatmap
  let perSampleES = $derived.by(() => {
    if (!magsData) return null;
    const pw = esData?.pathway_heatmap;
    const hm = esData?.heatmap;
    if (!pw && !hm) return null;

    // Source data
    const srcMags = pw ? pw.mags : hm.mags;
    const srcMatrix = pw ? pw.matrix : hm.matrix;
    const nCols = pw ? pw.columns.length : hm.es_codes.length;

    const magIdx = {};
    for (let i = 0; i < srcMags.length; i++) magIdx[srcMags[i]] = i;
    const magList = Object.values(magsData);
    if (!magList.length || !magList[0].coverage) return null;

    const sampleSet = new Set();
    for (const m of magList) {
      if (m.coverage) for (const s of Object.keys(m.coverage)) sampleSet.add(s);
    }
    const samples = [...sampleSet].sort();
    const rows = [];
    for (const sample of samples) {
      const row = new Array(nCols).fill(0);
      let totalCov = 0;
      for (const m of magList) {
        const cov = m.coverage?.[sample] || 0;
        if (cov === 0) continue;
        const idx = magIdx[m.name];
        if (idx === undefined) continue;
        totalCov += cov;
        for (let e = 0; e < nCols; e++) {
          row[e] += srcMatrix[idx][e] * cov;
        }
      }
      if (totalCov > 0) {
        for (let e = 0; e < nCols; e++) row[e] /= totalCov;
      }
      rows.push(row);
    }
    return { samples, rows, columns: pw ? pw.columns : null, es_codes: hm?.es_codes };
  });

  // === ES Summary Cards (always uses the 8-code heatmap for summary) ===
  let esSummary = $derived.by(() => {
    if (!esData?.heatmap) return [];
    const { es_codes, es_names, matrix } = esData.heatmap;
    if (abundanceMode === 'relative' && perSampleES) {
      // Aggregate per-sample pathway scores back to ES code level
      const pw = perSampleES.columns;
      return es_codes.map((code, codeIdx) => {
        let scores;
        if (pw) {
          // Sum pathway columns belonging to this ES code per sample
          const colIndices = pw.map((c, i) => c.es_code === code ? i : -1).filter(i => i >= 0);
          scores = perSampleES.rows.map(row => colIndices.reduce((s, i) => s + row[i], 0)).filter(v => v > 0);
        } else {
          scores = perSampleES.rows.map(row => row[codeIdx]).filter(v => v > 0);
        }
        const total = scores.reduce((a, b) => a + b, 0);
        const mean = scores.length ? total / scores.length : 0;
        const max = scores.length ? Math.max(...scores) : 0;
        return {
          code, name: es_names[code] || code,
          total: total.toFixed(1), mean: mean.toFixed(3),
          max: max.toFixed(3), count: scores.length,
          subtitle: `${scores.length} samples`
        };
      }).sort((a, b) => b.total - a.total);
    }
    return es_codes.map((code, i) => {
      const scores = matrix.map(row => row[i]).filter(v => v !== 0);
      const total = scores.reduce((a, b) => a + b, 0);
      const mean = scores.length ? total / scores.length : 0;
      const max = scores.length ? Math.max(...scores) : 0;
      return {
        code, name: es_names[code] || code,
        total: total.toFixed(1), mean: mean.toFixed(3),
        max: max.toFixed(3), count: scores.length,
        subtitle: `${scores.length} MAGs`
      };
    }).sort((a, b) => b.total - a.total);
  });

  // === SDG Goal Scores ===
  let sdgGoals = $derived.by(() => {
    if (!esData?.sdg?.goals) return [];
    const goals = {};
    for (const g of esData.sdg.goals) {
      const key = g.sdg_goal;
      if (!goals[key]) goals[key] = { goal: key, name: g.goal_name, score: 0 };
      goals[key].score += g.score;
    }
    return Object.values(goals).sort((a, b) => b.score - a.score);
  });

  let maxSdgScore = $derived(sdgGoals.length ? Math.max(...sdgGoals.map(g => g.score)) : 1);

  const sdgColors = {
    1: '#E5243B', 2: '#DDA63A', 3: '#4C9F38', 6: '#26BDE2',
    7: '#FCC30B', 8: '#A21942', 9: '#FD6925', 11: '#FD9D24',
    12: '#BF8B2E', 13: '#3F7E44', 14: '#0A97D9', 15: '#56C02B'
  };

  // === ES color scale (log-scaled cyan gradient) ===
  let esColorScale = $derived.by(() => {
    // Find max across the active heatmap data source
    let maxVal = 1;
    const src = abundanceMode === 'even'
      ? (esData?.pathway_heatmap?.matrix || esData?.heatmap?.matrix)
      : (perSampleES?.rows);
    if (src) {
      for (const row of src) {
        for (const v of row) if (v > maxVal) maxVal = v;
      }
    }
    // Log scale for hit counts (large dynamic range)
    const logMax = Math.log10(maxVal + 1);
    return (val) => {
      if (val <= 0) return '#0f172a';
      const t = Math.min(Math.log10(val + 1) / logMax, 1);
      const r = Math.round(15 + t * (6 - 15));
      const g = Math.round(23 + t * (182 - 23));
      const b = Math.round(42 + t * (212 - 42));
      return `rgb(${r}, ${g}, ${b})`;
    };
  });

  let esMaxVal = $derived.by(() => {
    let maxVal = 1;
    const src = abundanceMode === 'even'
      ? (esData?.pathway_heatmap?.matrix || esData?.heatmap?.matrix)
      : (perSampleES?.rows);
    if (src) {
      for (const row of src) {
        for (const v of row) if (v > maxVal) maxVal = v;
      }
    }
    return maxVal;
  });

  // === ES code filter for pathway heatmap ===
  let activeEsCodes = $state(null); // null = show all

  let availableEsCodes = $derived.by(() => {
    if (!esData?.pathway_heatmap) return [];
    const codes = [...new Set(esData.pathway_heatmap.columns.map(c => c.es_code))].sort();
    return codes;
  });

  function toggleEsCode(code) {
    if (!activeEsCodes) {
      // Currently showing all — switch to just this one
      activeEsCodes = new Set([code]);
    } else {
      const next = new Set(activeEsCodes);
      if (next.has(code)) {
        if (next.size > 1) next.delete(code);
        else activeEsCodes = null; // back to all
      } else {
        next.add(code);
      }
      if (next.size === availableEsCodes.length) {
        activeEsCodes = null; // all selected = show all
      } else {
        activeEsCodes = next;
      }
    }
  }

  // === MAG Heatmap D3 data (even mode) — uses pathway_heatmap if available ===
  let magD3Data = $derived.by(() => {
    // Prefer pathway heatmap (ES code + pathway columns)
    const pw = esData?.pathway_heatmap;
    const hm = esData?.heatmap;
    if (!pw && !hm) return null;

    let mags, matrix, moduleIds, moduleNames, colEsCodes, colIndices;
    if (pw) {
      // Filter columns by active ES codes
      colIndices = pw.columns.map((_, i) => i);
      if (activeEsCodes) {
        colIndices = colIndices.filter(i => activeEsCodes.has(pw.columns[i].es_code));
      }
      mags = pw.mags;
      matrix = pw.matrix.map(row => colIndices.map(i => row[i]));
      moduleIds = colIndices.map(i => `${pw.columns[i].es_code} ${pw.columns[i].pathway}`);
      moduleNames = colIndices.map(i => pw.columns[i].pathway);
      colEsCodes = colIndices.map(i => pw.columns[i].es_code);
    } else {
      const { es_codes, es_names, matrix: hMatrix, mags: hMags } = hm;
      mags = hMags;
      matrix = hMatrix;
      moduleIds = es_codes;
      moduleNames = es_codes.map(c => es_names[c] || c);
      colEsCodes = es_codes;
    }

    // Filter by binner
    const binners = pw?.binners;
    const filteredIndices = [];
    for (let i = 0; i < mags.length; i++) {
      if (!matrix[i].some(v => v !== 0)) continue;
      const bk = binners ? binners[i] : (magMeta[mags[i]]?.binner || 'unknown');
      if (!activeBinners.has(bk)) continue;
      filteredIndices.push(i);
    }

    const magIds = filteredIndices.map(i => mags[i]);
    const matrixFiltered = filteredIndices.map(i => matrix[i]);

    // Sort rows
    let order = filteredIndices.map((_, i) => i);
    if (magSortBy === 'clustered' && pw?.row_order) {
      // Map precomputed order through binner filter
      const filteredSet = new Set(filteredIndices);
      const clusterFiltered = pw.row_order.filter(i => filteredSet.has(i));
      order = clusterFiltered.map(origIdx => filteredIndices.indexOf(origIdx));
    } else if (magSortBy === 'name') {
      order.sort((a, b) => magIds[a].localeCompare(magIds[b]));
    } else if (magSortBy === 'total_score') {
      order.sort((a, b) => matrixFiltered[b].reduce((s, v) => s + v, 0) - matrixFiltered[a].reduce((s, v) => s + v, 0));
    } else if (magSortBy === 'completeness') {
      order.sort((a, b) => (magMeta[magIds[b]]?.completeness ?? 0) - (magMeta[magIds[a]]?.completeness ?? 0));
    } else if (magSortBy === 'contamination') {
      order.sort((a, b) => (magMeta[magIds[a]]?.contamination ?? 100) - (magMeta[magIds[b]]?.contamination ?? 100));
    } else if (magSortBy === 'size') {
      order.sort((a, b) => (magMeta[magIds[b]]?.size ?? 0) - (magMeta[magIds[a]]?.size ?? 0));
    } else if (magSortBy === 'binner') {
      const getBinner = (i) => binners ? binners[filteredIndices[i]] : (magMeta[magIds[i]]?.binner || '');
      order.sort((a, b) => getBinner(a).localeCompare(getBinner(b)));
    }

    // Column order: use Bray-Curtis precomputed order, filtered to active columns
    let colOrder;
    if (pw?.col_order && !activeEsCodes) {
      // Map precomputed full col_order through the colIndices filter
      colOrder = pw.col_order.filter(i => colIndices.includes(i)).map(origIdx => colIndices.indexOf(origIdx));
    } else if (pw?.col_order && activeEsCodes) {
      // Filtered columns — remap precomputed order
      const colSet = new Set(colIndices);
      const filtered = pw.col_order.filter(i => colSet.has(i));
      colOrder = filtered.map(origIdx => colIndices.indexOf(origIdx));
    } else {
      colOrder = moduleIds.map((_, i) => i);
    }

    return {
      mag_ids: magIds,
      module_ids: moduleIds,
      module_names: moduleNames,
      matrix: matrixFiltered,
      row_order: order,
      col_order: colOrder,
    };
  });

  // MAG heatmap row annotations
  let magRowAnnotations = $derived.by(() => {
    if (!magD3Data) return null;
    const ids = magD3Data.mag_ids;
    const qualityColor = (v) => v === 'HQ' ? '#34d399' : v === 'MQ' ? '#fbbf24' : '#64748b';
    const compColor = (v) => {
      const n = parseFloat(v);
      return n >= 90 ? '#34d399' : n >= 50 ? '#fbbf24' : '#64748b';
    };
    const contColor = (v) => {
      const n = parseFloat(v);
      return n < 5 ? '#34d399' : n < 10 ? '#fbbf24' : '#ef4444';
    };
    const scoreColor = (v) => {
      const n = parseFloat(v);
      return n > 0 ? '#22d3ee' : '#64748b';
    };
    return [
      {
        label: 'MAG', width: 130, align: 'start',
        values: ids.map(id => id),
      },
      {
        label: 'Qual', width: 30,
        values: ids.map(id => { const m = magMeta[id]; return m ? qualityTier(m) : '-'; }),
        colorFn: qualityColor,
      },
      {
        label: 'Comp%', width: 40,
        values: ids.map(id => { const m = magMeta[id]; return m ? m.completeness.toFixed(0) : '-'; }),
        colorFn: compColor,
      },
      {
        label: 'Cont%', width: 40,
        values: ids.map(id => { const m = magMeta[id]; return m ? m.contamination.toFixed(1) : '-'; }),
        colorFn: contColor,
      },
      {
        label: 'Score', width: 40,
        values: ids.map(id => {
          const idx = magD3Data.mag_ids.indexOf(id);
          if (idx < 0) return '-';
          return magD3Data.matrix[idx].reduce((s, v) => s + v, 0).toFixed(1);
        }),
        colorFn: scoreColor,
      },
    ];
  });

  // MAG heatmap tooltip
  function magTooltip(val, magId, colId, colName) {
    const m = magMeta[magId];
    const qual = m ? qualityTier(m) : '-';
    const comp = m ? m.completeness.toFixed(1) + '%' : '-';
    // colId is "es_code pathway" for pathway heatmap, or just es_code
    const isPathway = colId.includes(' ');
    const esCode = isPathway ? colId.split(' ')[0] : colId;
    const pathway = isPathway ? colName : '';
    const valLabel = Number.isInteger(val) || val > 10 ? `${val} hits` : `Score: ${val.toFixed(3)}`;
    return `<strong>${esCode}</strong>${pathway ? '<br>' + pathway : ''}<br>MAG: ${magId}<br>${valLabel}<br>Quality: ${qual} (${comp} comp)`;
  }

  // === Sample Heatmap D3 data (relative mode, pathway-expanded) ===
  let sampleD3Data = $derived.by(() => {
    if (!perSampleES || !esData?.heatmap) return null;
    const { es_codes, es_names } = esData.heatmap;
    const pw = esData?.pathway_heatmap;
    const hasPw = pw && perSampleES.columns;

    // Build module IDs/names from pathway columns or ES codes
    let allModuleIds, allModuleNames;
    if (hasPw) {
      // Filter by active ES codes
      let colIndices = pw.columns.map((_, i) => i);
      if (activeEsCodes) {
        colIndices = colIndices.filter(i => activeEsCodes.has(pw.columns[i].es_code));
      }
      allModuleIds = colIndices.map(i => `${pw.columns[i].es_code} ${pw.columns[i].pathway}`);
      allModuleNames = colIndices.map(i => pw.columns[i].pathway);
      // Filter sample matrix columns too
      var filteredRows = perSampleES.rows.map(row => colIndices.map(i => row[i]));
      var usedColIndices = colIndices;
    } else {
      allModuleIds = es_codes;
      allModuleNames = es_codes.map(c => es_names[c] || c);
      var filteredRows = perSampleES.rows;
      var usedColIndices = null;
    }

    const activeSamples = [];
    const activeRows = [];
    for (let i = 0; i < perSampleES.samples.length; i++) {
      const row = filteredRows[i];
      if (row.some(v => v > 0)) {
        activeSamples.push(perSampleES.samples[i]);
        activeRows.push(row);
      }
    }

    // Sort rows
    let order = activeSamples.map((_, i) => i);
    if (sampleSortBy === 'clustered' && activeRows.length > 2) {
      // Client-side greedy nearest-neighbour (Bray-Curtis)
      const n = activeRows.length;
      const bc = (a, b) => {
        let num = 0, den = 0;
        for (let k = 0; k < a.length; k++) { num += Math.abs(a[k] - b[k]); den += a[k] + b[k]; }
        return den > 0 ? num / den : 0;
      };
      const sums = activeRows.map(r => r.reduce((s, v) => s + v, 0));
      let start = 0;
      for (let i = 1; i < n; i++) if (sums[i] > sums[start]) start = i;
      const visited = new Array(n).fill(false);
      order = [start]; visited[start] = true;
      for (let step = 1; step < n; step++) {
        const curr = order[order.length - 1];
        let bestJ = -1, bestD = Infinity;
        for (let j = 0; j < n; j++) {
          if (!visited[j]) { const d = bc(activeRows[curr], activeRows[j]); if (d < bestD) { bestD = d; bestJ = j; } }
        }
        order.push(bestJ); visited[bestJ] = true;
      }
    } else if (sampleSortBy === 'name') {
      order.sort((a, b) => activeSamples[a].localeCompare(activeSamples[b]));
    } else if (sampleSortBy === 'total_score') {
      order.sort((a, b) => activeRows[b].reduce((s, v) => s + v, 0) - activeRows[a].reduce((s, v) => s + v, 0));
    } else if (sampleSortBy === 'top_es') {
      order.sort((a, b) => Math.max(...activeRows[b]) - Math.max(...activeRows[a]));
    } else if (sampleSortBy === 'diversity') {
      order.sort((a, b) => {
        const da = activeRows[a].filter(v => v > 0).length;
        const db = activeRows[b].filter(v => v > 0).length;
        return db - da || activeRows[b].reduce((s, v) => s + v, 0) - activeRows[a].reduce((s, v) => s + v, 0);
      });
    }

    // Column order: use precomputed Bray-Curtis from pathway_heatmap
    let colOrder;
    if (hasPw && pw.col_order) {
      const colSet = new Set(usedColIndices);
      colOrder = pw.col_order.filter(i => colSet.has(i)).map(origIdx => usedColIndices.indexOf(origIdx));
    } else {
      colOrder = allModuleIds.map((_, i) => i);
    }

    return {
      mag_ids: activeSamples,
      module_ids: allModuleIds,
      module_names: allModuleNames,
      matrix: activeRows,
      row_order: order,
      col_order: colOrder,
    };
  });

  // Sample heatmap row annotations
  let sampleRowAnnotations = $derived.by(() => {
    if (!sampleD3Data) return null;
    const ids = sampleD3Data.mag_ids;
    const scoreColor = (v) => {
      const n = parseFloat(v);
      return n > 0 ? '#22d3ee' : '#64748b';
    };
    const divColor = (v) => {
      const n = parseInt(v);
      return n >= 6 ? '#34d399' : n >= 3 ? '#fbbf24' : '#64748b';
    };
    return [
      {
        label: 'Sample', width: 150, align: 'start',
        values: ids.map(id => id),
      },
      {
        label: 'Score', width: 40,
        values: ids.map((id, i) => sampleD3Data.matrix[i].reduce((s, v) => s + v, 0).toFixed(2)),
        colorFn: scoreColor,
      },
      {
        label: 'ES#', width: 30,
        values: ids.map((id, i) => String(sampleD3Data.matrix[i].filter(v => v > 0).length)),
        colorFn: divColor,
      },
    ];
  });

  // Sample heatmap tooltip
  function sampleTooltip(val, sampleId, colId, colName) {
    const isPathway = colId.includes(' ');
    const esCode = isPathway ? colId.split(' ')[0] : colId;
    const pathway = isPathway ? colName : '';
    return `<strong>${esCode}</strong>${pathway ? '<br>' + pathway : ''}<br>Sample: ${sampleId}<br>Weighted score: ${val.toFixed(4)}`;
  }

  // Active D3 data
  let activeD3Data = $derived(abundanceMode === 'relative' ? sampleD3Data : magD3Data);
  let activeRowAnnotations = $derived(abundanceMode === 'relative' ? sampleRowAnnotations : magRowAnnotations);

  function handleRowClick(magId) {
    if (abundanceMode === 'even') {
      selectedMag.set(magId);
    }
  }

  // === Catalog table (unique gene × ES, sortable) ===
  const catalogColumns = [
    { key: 'gene_id',         label: 'Gene ID' },
    { key: 'gene_type',       label: 'Type' },
    { key: 'es_code',         label: 'ES Code' },
    { key: 'es_name',         label: 'ES Name' },
    { key: 'confidence',      label: 'Conf.' },
    { key: 'role',            label: 'Role' },
    { key: 'n_hits',          label: 'Hits' },
    { key: 'pathway_context', label: 'Pathway' },
  ];
  let catalogSortCol = $state('n_hits');
  let catalogSortAsc = $state(false);

  function handleCatalogSort(key) {
    if (catalogSortCol === key) {
      catalogSortAsc = !catalogSortAsc;
    } else {
      catalogSortCol = key;
      catalogSortAsc = (key === 'confidence' || key === 'n_hits') ? false : true;
    }
  }

  let catalogRows = $derived.by(() => {
    if (!esData?.catalog) return [];
    let rows = esData.catalog.map(row => ({
      gene_id: row.gene_id,
      gene_type: row.gene_id_type,
      es_code: row.es_code,
      es_name: row.es_name,
      confidence: row.confidence,
      role: row.functional_role,
      n_hits: row.n_hits,
      pathway_context: row.pathway_context || '',
      notes: row.notes || '',
    }));
    if (catalogSortCol) {
      const key = catalogSortCol;
      rows = [...rows].sort((a, b) => {
        let va = a[key], vb = b[key];
        if (key === 'confidence' || key === 'n_hits') {
          va = parseFloat(va) || 0;
          vb = parseFloat(vb) || 0;
          return catalogSortAsc ? va - vb : vb - va;
        }
        va = String(va || '');
        vb = String(vb || '');
        return catalogSortAsc ? va.localeCompare(vb) : vb.localeCompare(va);
      });
    }
    return rows;
  });

  // ES category colors
  const esColors = {
    '2.1.1.1': '#f59e0b',
    '2.1.1.2': '#a855f7',
    '2.3.3.2': '#ef4444',
    '2.3.4.2': '#84cc16',
    '2.3.5.1': '#06b6d4',
    '2.3.5.2': '#3b82f6',
    '2.3.6.1': '#10b981',
    '2.3.6.2': '#8b5cf6',
  };

  // === TSV Export ===
  function downloadTsv(filename, header, bodyRows) {
    const content = header + '\n' + bodyRows.join('\n');
    const blob = new Blob([content], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    a.click();
    URL.revokeObjectURL(url);
  }

  function exportHeatmapTsv() {
    if (!activeD3Data) return;
    const isRelative = abundanceMode === 'relative';
    const rowType = isRelative ? 'Sample' : 'MAG';
    const { mag_ids, module_ids, module_names, matrix, row_order } = activeD3Data;
    const header = [rowType, ...module_ids.map((id, i) => `${id} ${module_names[i]}`)].join('\t');
    const body = row_order.map(ri =>
      [mag_ids[ri], ...matrix[ri].map(v => v.toFixed(4))].join('\t')
    );
    const filename = isRelative
      ? 'ecossdb_sample_relative_abundance_es.tsv'
      : 'ecossdb_mag_even_abundance_es.tsv';
    downloadTsv(filename, header, body);
  }

  function exportCatalogTsv() {
    if (!catalogRows.length) return;
    const header = ['Gene_ID', 'Gene_Type', 'ES_Code', 'ES_Name', 'Confidence', 'Functional_Role', 'Protein_Hits', 'Pathway_Context', 'Notes'].join('\t');
    const body = catalogRows.map(r =>
      [r.gene_id, r.gene_type, r.es_code, r.es_name, r.confidence, r.role, r.n_hits, r.pathway_context, r.notes].join('\t')
    );
    downloadTsv('ecossdb_gene_es_mappings.tsv', header, body);
  }
</script>

<div class="space-y-6">
  <div class="flex items-center justify-between flex-wrap gap-2">
    <div>
      <h2 class="text-xl font-semibold text-slate-200">Ecosystem Services</h2>
      <p class="text-sm text-slate-400">
        Gene-to-ecosystem-service mapping via ECOSSDB ({esData ? `${esData.catalog?.length || 0} gene-ES mappings, ${(esData.total_protein_hits || 0).toLocaleString()} protein hits` : 'loading...'})
      </p>
    </div>
    {#if esData && magsData}
      <div class="flex items-center gap-1 bg-slate-800 rounded-lg p-0.5">
        <button class="px-3 py-1 text-xs rounded-md transition-colors {abundanceMode === 'even' ? 'bg-cyan-600 text-white' : 'text-slate-400 hover:text-slate-200'}"
                onclick={() => abundanceMode = 'even'}>
          Even Abundance
        </button>
        <button class="px-3 py-1 text-xs rounded-md transition-colors {abundanceMode === 'relative' ? 'bg-cyan-600 text-white' : 'text-slate-400 hover:text-slate-200'}"
                onclick={() => abundanceMode = 'relative'}>
          Relative Abundance
        </button>
      </div>
    {/if}
  </div>

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
          <div class="text-xs text-slate-400 truncate" title={es.name}>{es.code} {es.name}</div>
          <div class="text-lg font-semibold text-slate-200 mt-1">{es.total}</div>
          <div class="text-xs text-slate-500 mt-0.5">{es.subtitle} · mean {es.mean}</div>
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

    <!-- D3 Heatmap -->
    {#if activeD3Data && activeD3Data.mag_ids.length > 0}
      <div class="bg-slate-800/30 border border-slate-700 rounded-lg p-4">
        <div class="flex items-center justify-between flex-wrap gap-2 mb-2">
          <h3 class="text-sm font-medium text-slate-400">
            {abundanceMode === 'relative' ? 'Sample' : 'MAG'} × Ecosystem Service Heatmap
            ({activeD3Data.mag_ids.length} {abundanceMode === 'relative' ? 'samples' : 'MAGs'} × {activeD3Data.module_ids.length} {esData?.pathway_heatmap ? 'pathways' : 'ES categories'})
          </h3>
          <div class="flex items-center gap-2">
            {#if abundanceMode === 'even'}
              <button
                class="text-xs px-2 py-0.5 rounded border border-cyan-400 bg-cyan-400/10 text-cyan-400"
                onclick={cycleMagSort}
                title={`Sort rows: ${magSortDefs.map(d => d.label).join(' → ')}`}
              >Sort: {magSortDefs.find(d => d.key === magSortBy)?.label} &#x25BE;</button>
            {:else}
              <button
                class="text-xs px-2 py-0.5 rounded border border-cyan-400 bg-cyan-400/10 text-cyan-400"
                onclick={cycleSampleSort}
                title={`Sort rows: ${sampleSortDefs.map(d => d.label).join(' → ')}`}
              >Sort: {sampleSortDefs.find(d => d.key === sampleSortBy)?.label} &#x25BE;</button>
            {/if}
            <button
              class="text-xs px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:border-slate-500"
              onclick={exportHeatmapTsv}
              title="Export heatmap as TSV"
            >TSV</button>
          </div>
        </div>

        <!-- Binner filter (even mode) -->
        {#if abundanceMode === 'even'}
          <div class="flex items-center gap-2 flex-wrap text-xs mb-2">
            <span class="text-slate-400">Source binner:</span>
            {#each binnerDefs as { key, label }}
              {#if binnerCounts[key]}
                <button
                  class="px-2.5 py-0.5 rounded-md border transition-colors
                    {activeBinners.has(key) ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
                  onclick={() => toggleBinner(key)}
                >
                  {label} <span class="text-slate-500">({binnerCounts[key]})</span>
                </button>
              {/if}
            {/each}
            <span class="text-slate-500">{activeD3Data.mag_ids.length} MAGs</span>
          </div>
        {/if}

        <!-- ES code filter (pathway heatmap — both modes) -->
        {#if esData?.pathway_heatmap && availableEsCodes.length > 1}
          <div class="flex items-center gap-2 flex-wrap text-xs mb-2">
            <span class="text-slate-400">ES filter:</span>
            <button
              class="px-2 py-0.5 rounded-md border transition-colors
                {!activeEsCodes ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
              onclick={() => activeEsCodes = null}
            >All</button>
            {#each availableEsCodes as code}
              <button
                class="px-2 py-0.5 rounded-md border transition-colors"
                style="border-color: {!activeEsCodes || activeEsCodes.has(code) ? (esColors[code] || '#22d3ee') : '#475569'};
                       color: {!activeEsCodes || activeEsCodes.has(code) ? (esColors[code] || '#22d3ee') : '#64748b'};
                       background-color: {!activeEsCodes || activeEsCodes.has(code) ? (esColors[code] || '#22d3ee') + '18' : 'transparent'}"
                onclick={() => toggleEsCode(code)}
                title={esData.heatmap.es_names[code] || code}
              >{code}</button>
            {/each}
            <span class="text-slate-500">{activeD3Data.module_ids.length} pathways</span>
          </div>
        {/if}

        <div class="es-heatmap-scroll" style="max-height: 700px;">
          {#if abundanceMode === 'even'}
            <D3Heatmap
              data={magD3Data}
              onRowClick={handleRowClick}
              tooltipFormat={magTooltip}
              colorScale={esColorScale}
              legendLabel={['0', `${esMaxVal} hits`]}
              selectedRow={selected}
              rowAnnotations={magRowAnnotations}
            />
          {:else}
            <D3Heatmap
              data={sampleD3Data}
              tooltipFormat={sampleTooltip}
              colorScale={esColorScale}
              legendLabel={['0', esMaxVal.toFixed(2)]}
              rowAnnotations={sampleRowAnnotations}
            />
          {/if}
        </div>
      </div>
    {/if}

    <!-- ES Gene Mappings Table -->
    {#if catalogRows.length > 0}
      <div class="bg-slate-800/30 border border-slate-700 rounded-lg p-4">
        <div class="flex items-center justify-between mb-3">
          <h3 class="text-sm font-semibold text-slate-300">
            ES Gene Mappings
            <span class="text-slate-500 font-normal ml-2">({catalogRows.length} gene-ES pairs from {esData?.total_protein_hits?.toLocaleString() || 0} protein hits)</span>
          </h3>
          <button
            class="text-xs px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:border-slate-500"
            onclick={exportCatalogTsv}
            title="Export gene mappings as TSV"
          >TSV</button>
        </div>
        <div class="overflow-x-auto overflow-y-auto" style="max-height: 32rem">
          <table class="w-full text-xs">
            <thead class="sticky top-0 bg-slate-900 z-10">
              <tr class="border-b border-slate-700">
                {#each catalogColumns as col}
                  <th class="{col.key === 'confidence' || col.key === 'n_hits' ? 'text-center' : 'text-left'} px-2 py-1 text-slate-400 cursor-pointer hover:text-cyan-400 whitespace-nowrap select-none"
                      onclick={() => handleCatalogSort(col.key)}>
                    {col.label}{#if catalogSortCol === col.key}<span class="ml-1">{catalogSortAsc ? '↑' : '↓'}</span>{/if}
                  </th>
                {/each}
              </tr>
            </thead>
            <tbody>
              {#each catalogRows as row}
                <tr class="border-b border-slate-800 hover:bg-slate-800/50" title={row.notes}>
                  <td class="px-2 py-0.5 text-cyan-400 font-mono">{row.gene_id}</td>
                  <td class="px-2 py-0.5 text-slate-500">{row.gene_type}</td>
                  <td class="px-2 py-0.5" style="color: {esColors[row.es_code] || '#94a3b8'}">{row.es_code}</td>
                  <td class="px-2 py-0.5 text-slate-400 truncate max-w-40">{row.es_name}</td>
                  <td class="px-2 py-0.5 text-center text-slate-300">{row.confidence}</td>
                  <td class="px-2 py-0.5 text-slate-400">
                    {#if row.role === 'producer'}
                      <span class="text-emerald-400">▲ prod</span>
                    {:else if row.role === 'inhibitor'}
                      <span class="text-rose-400">▼ inh</span>
                    {:else if row.role === 'transformer'}
                      <span class="text-amber-400">◆ trans</span>
                    {:else}
                      <span class="text-sky-400">○ {row.role}</span>
                    {/if}
                  </td>
                  <td class="px-2 py-0.5 text-center text-slate-300 font-mono">{row.n_hits.toLocaleString()}</td>
                  <td class="px-2 py-0.5 text-slate-500 truncate max-w-48" title={row.pathway_context}>{row.pathway_context}</td>
                </tr>
              {/each}
            </tbody>
          </table>
        </div>
      </div>
    {/if}
  {/if}
</div>

<style>
  /* Force scrollbar to always show on heatmap container */
  .es-heatmap-scroll {
    overflow: scroll !important;
  }
  .es-heatmap-scroll :global(.relative.w-full.overflow-x-auto) {
    overflow: visible !important;
  }
</style>
