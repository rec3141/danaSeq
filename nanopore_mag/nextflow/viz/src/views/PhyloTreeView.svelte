<script>
  import D3PhyloTree from '../components/charts/D3PhyloTree.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { onMount } from 'svelte';
  import { phyloTree, loadPhyloTree } from '../stores/data.js';

  let treeData = $derived($phyloTree);

  onMount(() => { loadPhyloTree(); });

  // Layout toggle (cycles)
  let layout = $state('radial');
  const layoutDefs = [
    { key: 'radial', label: 'Radial' },
    { key: 'rectangular', label: 'Rectangular' },
  ];

  function cycleLayout() {
    const keys = layoutDefs.map(l => l.key);
    layout = keys[(keys.indexOf(layout) + 1) % keys.length];
  }
  function layoutLabel() {
    return layoutDefs.find(l => l.key === layout)?.label || layout;
  }

  // Color-by toggle (cycles)
  let colorBy = $state('binner');
  const colorDefs = [
    { key: 'binner', label: 'Binner' },
    { key: 'domain', label: 'Domain' },
    { key: 'phylum', label: 'Phylum' },
    { key: 'quality', label: 'Quality' },
  ];

  function cycleColor() {
    const keys = colorDefs.map(c => c.key);
    colorBy = keys[(keys.indexOf(colorBy) + 1) % keys.length];
  }
  function colorLabel() {
    return colorDefs.find(c => c.key === colorBy)?.label || colorBy;
  }

  // Table columns
  const tableColumns = [
    { key: 'name', label: 'Bin' },
    { key: 'binner', label: 'Binner' },
    { key: 'domain', label: 'Domain' },
    { key: 'phylum', label: 'Phylum' },
    { key: 'class', label: 'Class' },
    { key: 'order', label: 'Order' },
    { key: 'family', label: 'Family' },
    { key: 'genus', label: 'Genus' },
    { key: 'species', label: 'Species' },
    { key: 'completeness', label: 'Comp%' },
    { key: 'contamination', label: 'Cont%' },
    { key: 'method', label: 'Method' },
  ];

  let tableRows = $derived.by(() => {
    if (!treeData?.bins?.length) return [];
    return treeData.bins.map(b => ({
      name: b.name,
      binner: b.binner,
      domain: b.lineage?.domain || '-',
      phylum: b.lineage?.phylum || '-',
      class: b.lineage?.class || '-',
      order: b.lineage?.order || '-',
      family: b.lineage?.family || '-',
      genus: b.lineage?.genus || '-',
      species: b.lineage?.species || '-',
      completeness: b.completeness != null ? b.completeness : '-',
      contamination: b.contamination != null ? b.contamination : '-',
      method: b.method || '-',
    }));
  });

  let hasData = $derived(treeData?.hierarchy != null && treeData?.bins?.length > 0);
</script>

<div class="sticky top-14 z-10 bg-slate-950 pb-2 -mx-4 px-4 pt-1">
<div class="flex items-center gap-3 flex-wrap text-xs">
  <span class="text-slate-400">Layout:</span>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
    style="min-width: 6rem"
    onclick={cycleLayout}
    title={`Click to cycle: ${layoutDefs.map(l => l.label).join(' → ')}`}
  >
    {layoutLabel()} &#x25BE;
  </button>
  <span class="text-slate-400">Color by:</span>
  <button
    class="px-3 py-1 rounded-md border transition-colors text-center border-cyan-400 bg-cyan-400/10 text-cyan-400"
    style="min-width: 5rem"
    onclick={cycleColor}
    title={`Click to cycle: ${colorDefs.map(c => c.label).join(' → ')}`}
  >
    {colorLabel()} &#x25BE;
  </button>
  {#if treeData?.bins}
    <span class="text-slate-500">{treeData.bins.length} bins classified</span>
  {/if}
</div>
</div>

{#if !treeData}
  <div class="text-slate-500 text-sm py-20 text-center">Loading phylogenetic data...</div>
{:else if !hasData}
  <div class="bg-slate-800 rounded-lg p-8 border border-slate-700 text-center">
    <h3 class="text-lg font-medium text-slate-300 mb-2">No GTDB-Tk Data</h3>
    <p class="text-slate-400 text-sm">
      Run the pipeline with <code class="bg-slate-700 px-1.5 py-0.5 rounded text-cyan-400">--gtdbtk_db /path/to/gtdbtk_data</code> to enable phylogenetic classification.
    </p>
  </div>
{:else}
  <div class="grid grid-cols-1 gap-6 mb-6">
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
      <h3 class="text-sm font-medium text-slate-400 mb-2">GTDB-Tk Phylogenetic Tree ({layoutLabel()})</h3>
      <D3PhyloTree
        data={treeData.hierarchy}
        {layout}
        {colorBy}
        bins={treeData.bins}
      />
    </div>
  </div>

  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
    <h3 class="text-sm font-medium text-slate-400 mb-2">Per-Bin Classification ({treeData.bins.length} bins)</h3>
    <DataTable
      columns={tableColumns}
      rows={tableRows}
      maxHeight="500px"
    />
  </div>
{/if}
