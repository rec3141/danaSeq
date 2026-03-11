<script>
  import StatCard from '../components/layout/StatCard.svelte';
  import PlotlyChart from '../components/charts/PlotlyChart.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { biosynthetic, mags, loadBiosynthetic } from '../stores/data.js';
  import { onMount } from 'svelte';

  onMount(() => { loadBiosynthetic(); });

  let data = $derived($biosynthetic);
  let magsData = $derived($mags);

  // BGC type bar chart
  let typeBarData = $derived.by(() => {
    if (!data?.type_counts) return [];
    const entries = Object.entries(data.type_counts).sort((a, b) => b[1] - a[1]);
    return [{
      type: 'bar',
      x: entries.map(([k]) => k),
      y: entries.map(([, v]) => v),
      marker: { color: '#34d399', opacity: 0.85 },
      hovertemplate: '%{x}: %{y}<extra></extra>',
    }];
  });

  // Per-MAG BGC count summary
  let magTableColumns = [
    { key: 'name', label: 'MAG' },
    { key: 'total', label: 'BGC Regions' },
    { key: 'types', label: 'Types' },
  ];

  let magTableRows = $derived.by(() => {
    if (!data?.per_bin) return [];
    return Object.entries(data.per_bin)
      .filter(([k]) => k !== '_unbinned')
      .map(([name, info]) => ({
        name,
        total: info.total,
        types: Object.entries(info.types).map(([t, c]) => `${t} (${c})`).join(', '),
      }))
      .sort((a, b) => b.total - a.total);
  });

  // Unbinned BGC count
  let unbinnedCount = $derived(data?.per_bin?.['_unbinned']?.total ?? 0);

  // All regions table
  let regionColumns = [
    { key: 'contig', label: 'Contig' },
    { key: 'products_str', label: 'Type' },
    { key: 'location', label: 'Location' },
    { key: 'length_str', label: 'Length' },
    { key: 'mag', label: 'MAG' },
    { key: 'best_match', label: 'Best Known Match' },
    { key: 'similarity', label: 'Similarity' },
  ];

  let regionRows = $derived.by(() => {
    if (!data?.regions) return [];
    return data.regions.map(r => {
      const best = r.known_matches?.[0];
      return {
        contig: r.contig,
        products_str: r.products?.join(', ') || '?',
        location: `${r.start.toLocaleString()}-${r.end.toLocaleString()}`,
        length_str: `${(r.length / 1000).toFixed(1)} kb`,
        mag: r.mag === '_unbinned' ? '-' : r.mag,
        best_match: best?.description || '-',
        similarity: best?.similarity != null ? `${best.similarity}%` : '-',
      };
    }).sort((a, b) => {
      // Sort by similarity descending (known matches first)
      const simA = parseFloat(a.similarity) || 0;
      const simB = parseFloat(b.similarity) || 0;
      return simB - simA;
    });
  });

  // Sorting state
  let sortKey = $state('similarity_val');
  let sortDir = $state(-1); // -1 = desc, 1 = asc

  function toggleSort(key) {
    if (sortKey === key) sortDir = -sortDir;
    else { sortKey = key; sortDir = key === 'contig' || key === 'products_str' || key === 'mag' ? 1 : -1; }
  }

  function sortIcon(key) {
    if (sortKey !== key) return '';
    return sortDir === 1 ? ' \u25B2' : ' \u25BC';
  }

  // Sorted regions for the table
  let sortedRegions = $derived.by(() => {
    if (!data?.regions) return [];
    const rows = data.regions.map((r, i) => {
      const best = r.known_matches?.[0];
      const badge = noveltyBadge(r);
      return {
        _idx: i,
        _region: r,
        contig: r.contig,
        products_str: r.products?.join(', ') || '?',
        length: r.length,
        mag: r.mag === '_unbinned' ? '' : r.mag,
        novelty: badge.label,
        best_match: best?.description || '',
        similarity_val: best?.similarity ?? -1,
      };
    });
    rows.sort((a, b) => {
      let va = a[sortKey], vb = b[sortKey];
      if (typeof va === 'string') return va.localeCompare(vb) * sortDir;
      return ((va ?? 0) - (vb ?? 0)) * sortDir;
    });
    return rows;
  });

  let expandedRegion = $state(null);

  let expandedDetail = $derived.by(() => {
    if (expandedRegion == null || !data?.regions) return null;
    return data.regions[expandedRegion] || null;
  });

  // TSV export
  function exportTSV() {
    if (!data?.regions) return;
    const header = ['contig', 'start', 'end', 'length_bp', 'type', 'mag', 'novelty', 'best_known_match', 'similarity_pct', 'match_accession', 'match_type', 'protoclusters', 'gene_biosynthetic', 'gene_transport', 'gene_regulatory', 'gene_other'].join('\t');
    const lines = sortedRegions.map(row => {
      const r = row._region;
      const best = r.known_matches?.[0];
      const badge = noveltyBadge(r);
      const gk = r.gene_kinds || {};
      return [
        r.contig, r.start, r.end, r.length,
        r.products?.join(';') || '',
        r.mag === '_unbinned' ? '' : r.mag,
        badge.label,
        best?.description || '', best?.similarity ?? '', best?.accession || '', best?.cluster_type || '',
        (r.protoclusters || []).map(p => p.product).join(';'),
        gk['biosynthetic'] || 0, gk['transport'] || 0, gk['regulatory'] || 0, gk['other'] || 0,
      ].join('\t');
    });
    const blob = new Blob([header + '\n' + lines.join('\n') + '\n'], { type: 'text/tab-separated-values' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = 'antismash_bgc_regions.tsv';
    a.click();
    URL.revokeObjectURL(url);
  }

  // Novelty classification
  function noveltyBadge(region) {
    if (!region.known_matches?.length) return { label: 'Novel', color: 'text-emerald-400 bg-emerald-400/10 border-emerald-400/30' };
    const sim = region.known_matches[0].similarity;
    if (sim >= 75) return { label: 'Known', color: 'text-slate-400 bg-slate-400/10 border-slate-400/30' };
    if (sim >= 40) return { label: 'Similar', color: 'text-amber-400 bg-amber-400/10 border-amber-400/30' };
    return { label: 'Distant', color: 'text-cyan-400 bg-cyan-400/10 border-cyan-400/30' };
  }

  // Novelty summary counts
  let noveltyCounts = $derived.by(() => {
    if (!data?.regions) return { novel: 0, distant: 0, similar: 0, known: 0 };
    const counts = { novel: 0, distant: 0, similar: 0, known: 0 };
    for (const r of data.regions) {
      const badge = noveltyBadge(r);
      if (badge.label === 'Novel') counts.novel++;
      else if (badge.label === 'Distant') counts.distant++;
      else if (badge.label === 'Similar') counts.similar++;
      else counts.known++;
    }
    return counts;
  });
</script>

{#if !data}
  <div class="flex items-center justify-center py-12 text-slate-500">
    <div class="w-6 h-6 border-2 border-slate-700 border-t-cyan-400 rounded-full animate-spin mr-3"></div>
    Loading biosynthetic gene cluster data...
  </div>
{:else if data.n_regions === 0}
  <div class="bg-slate-800 rounded-lg p-8 border border-slate-700 text-center">
    <p class="text-slate-400 text-lg">No biosynthetic gene clusters detected</p>
    <p class="text-slate-500 text-sm mt-2">antiSMASH did not identify any BGC regions in this assembly, or antiSMASH was not run.</p>
  </div>
{:else}
  <!-- Summary cards -->
  <div class="grid grid-cols-2 md:grid-cols-5 gap-4 mb-6">
    <StatCard label="BGC Regions" value={data.n_regions} color="emerald" />
    <StatCard label="BGC Types" value={Object.keys(data.type_counts).length} color="cyan" />
    <StatCard label="Novel" value={noveltyCounts.novel} color="emerald" />
    <StatCard label="Similar/Known" value={noveltyCounts.similar + noveltyCounts.known} color="amber" />
    <StatCard label="Unbinned" value={unbinnedCount} color="slate" />
  </div>

  <!-- Charts row -->
  <div class="grid grid-cols-1 lg:grid-cols-2 gap-6 mb-6">
    <!-- BGC type distribution -->
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <PlotlyChart
        data={typeBarData}
        layout={{
          title: { text: 'BGC Type Distribution', font: { size: 14 } },
          xaxis: { tickangle: -45 },
          yaxis: { title: 'Count' },
        }}
      />
    </div>

    <!-- Per-MAG BGC table -->
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <h3 class="text-sm font-medium text-slate-400 mb-2">BGCs per MAG</h3>
      {#if magTableRows.length}
        <DataTable
          columns={magTableColumns}
          rows={magTableRows}
          maxHeight="300px"
        />
      {:else}
        <p class="text-slate-500 text-sm">No MAG-assigned BGC regions found.</p>
      {/if}
    </div>
  </div>

  <!-- All regions table -->
  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
    <div class="flex items-center justify-between mb-2">
      <h3 class="text-sm font-medium text-slate-400">All BGC Regions ({data.n_regions})</h3>
      <button
        class="px-2 py-0.5 rounded border border-slate-600 text-slate-400 hover:text-slate-200 hover:border-slate-500 transition-colors text-xs"
        onclick={exportTSV}
        title="Export BGC regions as TSV"
      >TSV</button>
    </div>
    <div class="overflow-x-auto max-h-[500px] overflow-y-auto">
      <table class="w-full text-xs">
        <thead class="sticky top-0 bg-slate-800 z-10">
          <tr class="text-slate-500 border-b border-slate-700">
            <th class="text-left py-1.5 px-2 font-medium cursor-pointer select-none hover:text-slate-300" onclick={() => toggleSort('contig')}>Contig{sortIcon('contig')}</th>
            <th class="text-left py-1.5 px-2 font-medium cursor-pointer select-none hover:text-slate-300" onclick={() => toggleSort('products_str')}>Type{sortIcon('products_str')}</th>
            <th class="text-right py-1.5 px-2 font-medium cursor-pointer select-none hover:text-slate-300" onclick={() => toggleSort('length')}>Length{sortIcon('length')}</th>
            <th class="text-left py-1.5 px-2 font-medium cursor-pointer select-none hover:text-slate-300" onclick={() => toggleSort('mag')}>MAG{sortIcon('mag')}</th>
            <th class="text-center py-1.5 px-2 font-medium cursor-pointer select-none hover:text-slate-300" onclick={() => toggleSort('novelty')}>Novelty{sortIcon('novelty')}</th>
            <th class="text-left py-1.5 px-2 font-medium cursor-pointer select-none hover:text-slate-300" onclick={() => toggleSort('best_match')}>Best Known Match{sortIcon('best_match')}</th>
            <th class="text-right py-1.5 px-2 font-medium cursor-pointer select-none hover:text-slate-300" onclick={() => toggleSort('similarity_val')}>Similarity{sortIcon('similarity_val')}</th>
          </tr>
        </thead>
        <tbody>
          {#each sortedRegions as row}
            {@const region = row._region}
            {@const badge = noveltyBadge(region)}
            {@const best = region.known_matches?.[0]}
            <tr
              class="border-b border-slate-700/50 hover:bg-slate-700/30 cursor-pointer {expandedRegion === row._idx ? 'bg-slate-700/40' : ''}"
              onclick={() => { expandedRegion = expandedRegion === row._idx ? null : row._idx; }}
            >
              <td class="py-1 px-2 text-slate-300 font-mono text-[10px]">{region.contig}</td>
              <td class="py-1 px-2">
                <div class="flex flex-wrap gap-1">
                  {#each region.products || [] as prod}
                    <span class="px-1.5 py-0.5 rounded text-[10px] font-medium bg-emerald-400/10 text-emerald-400 border border-emerald-400/30">
                      {prod}
                    </span>
                  {/each}
                </div>
              </td>
              <td class="py-1 px-2 text-right text-slate-400 font-mono">{(region.length / 1000).toFixed(1)} kb</td>
              <td class="py-1 px-2 text-slate-400">{row.mag || '-'}</td>
              <td class="py-1 px-2 text-center">
                <span class="px-1.5 py-0.5 rounded text-[10px] font-medium border {badge.color}">{badge.label}</span>
              </td>
              <td class="py-1 px-2 text-slate-300 truncate max-w-[250px]" title={best?.description || ''}>
                {best?.description || '-'}
              </td>
              <td class="py-1 px-2 text-right font-mono {best?.similarity >= 75 ? 'text-slate-400' : best?.similarity >= 40 ? 'text-amber-400' : best?.similarity > 0 ? 'text-cyan-400' : 'text-slate-600'}">
                {best?.similarity != null ? `${best.similarity}%` : '-'}
              </td>
            </tr>
          {/each}
        </tbody>
      </table>
    </div>
  </div>

  <!-- Expanded region detail -->
  {#if expandedDetail}
    <div class="mt-4 bg-slate-800 rounded-lg p-4 border border-cyan-400/30">
      <h4 class="text-cyan-400 text-sm font-semibold mb-2">
        {expandedDetail.contig} : {expandedDetail.start.toLocaleString()}-{expandedDetail.end.toLocaleString()}
        ({(expandedDetail.length / 1000).toFixed(1)} kb)
      </h4>

      <div class="grid grid-cols-1 md:grid-cols-3 gap-4 text-xs">
        <!-- Products -->
        <div>
          <h5 class="text-slate-500 mb-1 font-medium">Products</h5>
          <div class="flex flex-wrap gap-1">
            {#each expandedDetail.products || [] as prod}
              <span class="px-2 py-0.5 rounded bg-emerald-400/10 text-emerald-400 border border-emerald-400/30">{prod}</span>
            {/each}
          </div>
        </div>

        <!-- Protoclusters -->
        {#if expandedDetail.protoclusters?.length}
          <div>
            <h5 class="text-slate-500 mb-1 font-medium">Protoclusters ({expandedDetail.protoclusters.length})</h5>
            {#each expandedDetail.protoclusters as pc}
              <div class="text-slate-300 mb-0.5">
                <span class="text-amber-400">{pc.product}</span>
                {#if pc.category}<span class="text-slate-500"> ({pc.category})</span>{/if}
                <span class="text-slate-500 font-mono"> {pc.core_start.toLocaleString()}-{pc.core_end.toLocaleString()}</span>
              </div>
            {/each}
          </div>
        {/if}

        <!-- Known matches -->
        {#if expandedDetail.known_matches?.length}
          <div>
            <h5 class="text-slate-500 mb-1 font-medium">Known Cluster Matches</h5>
            {#each expandedDetail.known_matches as km}
              <div class="text-slate-300 mb-1">
                <div class="truncate" title={km.description}>{km.description}</div>
                <div class="text-slate-500">
                  {#if km.accession}<span class="text-cyan-400">{km.accession}</span> — {/if}
                  Similarity: <span class="{km.similarity >= 75 ? 'text-slate-400' : km.similarity >= 40 ? 'text-amber-400' : 'text-cyan-400'}">{km.similarity}%</span>
                  {#if km.cluster_type} ({km.cluster_type}){/if}
                </div>
              </div>
            {/each}
          </div>
        {/if}

        <!-- Gene function summary -->
        {#if expandedDetail.gene_kinds && Object.keys(expandedDetail.gene_kinds).length}
          <div>
            <h5 class="text-slate-500 mb-1 font-medium">Gene Functions</h5>
            {#each Object.entries(expandedDetail.gene_kinds) as [kind, count]}
              <div class="text-slate-300">
                <span class="{kind === 'biosynthetic' ? 'text-emerald-400' : kind === 'transport' ? 'text-blue-400' : kind === 'regulatory' ? 'text-amber-400' : 'text-slate-400'}">{kind}</span>: {count}
              </div>
            {/each}
          </div>
        {/if}
      </div>
    </div>
  {/if}
{/if}
