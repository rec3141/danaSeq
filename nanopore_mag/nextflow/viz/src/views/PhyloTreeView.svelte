<script>
  import PhylocanvasTree from '../components/charts/PhylocanvasTree.svelte';
  import DataTable from '../components/ui/DataTable.svelte';
  import { onMount } from 'svelte';
  import { phyloTree, loadPhyloTree, checkm2All, loadCheckm2All } from '../stores/data.js';
  import { selectedMag } from '../stores/selection.js';

  let treeData = $derived($phyloTree);

  let allBins = $derived($checkm2All);

  onMount(() => { loadPhyloTree(); loadCheckm2All(); });

  // ---- Tree selector ----
  let treeNames = $derived.by(() => {
    if (!treeData?.newick) return [];
    return Object.keys(treeData.newick).sort((a, b) => {
      const aBack = a.includes('backbone') ? 0 : 1;
      const bBack = b.includes('backbone') ? 0 : 1;
      if (aBack !== bBack) return aBack - bBack;
      return a.localeCompare(b);
    });
  });

  let selectedTree = $state('');

  $effect(() => {
    if (treeNames.length && !selectedTree) {
      selectedTree = treeNames.find(n => n.includes('backbone')) || treeNames[0];
    }
  });

  let currentNewick = $derived(treeData?.newick?.[selectedTree] || '');

  function treeLabel(name) {
    const meta = treeData?.tree_metadata?.[name];
    const binCount = meta?.user_bins ?? '?';
    const short = name.replace(/\.tree$/, '').replace(/\.nwk$/, '').replace(/^gtdbtk\./, '');
    return `${short} (${binCount} bins)`;
  }

  // ---- Layout cycle ----
  let treeType = $state('rc');
  const layoutDefs = [
    { key: 'rc', label: 'Rectangular' },
    { key: 'rd', label: 'Radial' },
    { key: 'cr', label: 'Circular' },
    { key: 'dg', label: 'Diagonal' },
    { key: 'hr', label: 'Hierarchical' },
  ];
  function cycleLayout() {
    const keys = layoutDefs.map(l => l.key);
    treeType = keys[(keys.indexOf(treeType) + 1) % keys.length];
  }
  function layoutLabel() {
    return layoutDefs.find(l => l.key === treeType)?.label || treeType;
  }

  // ---- Color-by cycle ----
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

  // ---- Filters ----
  let filterBinner = $state('all');
  let filterPhylum = $state('all');
  let filterQuality = $state('all');
  let collapseRefs = $state(false);

  // ---- Palette ----
  const palette = ['#22d3ee', '#34d399', '#fbbf24', '#f87171', '#a78bfa', '#fb923c',
                   '#2dd4bf', '#818cf8', '#f472b6', '#4ade80', '#e879f9', '#38bdf8',
                   '#94a3b8', '#d4d4d8', '#78716c', '#64748b', '#475569'];

  // ---- Derived: bin lookups ----
  let userBinSet = $derived(new Set(treeData?.bins?.map(b => b.name) || []));

  // Reverse lookup: name -> phylotree bin object
  let binByName = $derived.by(() => {
    const map = {};
    for (const b of (treeData?.bins || [])) {
      map[b.name] = b;
    }
    return map;
  });

  // Cross-reference: name -> checkm2All quality data (completeness, gc, genome_size, n50, etc.)
  let qualityByName = $derived.by(() => {
    const map = {};
    if (!allBins) return map;
    for (const b of allBins) {
      map[b.name] = b;
    }
    return map;
  });

  let uniqueBinners = $derived.by(() => {
    if (!treeData?.bins?.length) return [];
    return [...new Set(treeData.bins.map(b => b.binner))].sort();
  });

  let uniquePhyla = $derived.by(() => {
    if (!treeData?.bins?.length) return [];
    return [...new Set(treeData.bins.map(b => b.lineage?.phylum).filter(Boolean))].sort();
  });

  // ---- Filtered bins ----
  let filteredBins = $derived.by(() => {
    if (!treeData?.bins?.length) return [];
    return treeData.bins.filter(b => {
      if (filterBinner !== 'all' && b.binner !== filterBinner) return false;
      if (filterPhylum !== 'all' && b.lineage?.phylum !== filterPhylum) return false;
      if (filterQuality !== 'all') {
        const comp = b.completeness ?? 0;
        const cont = b.contamination ?? 100;
        if (filterQuality === 'high' && !(comp >= 90 && cont < 5)) return false;
        if (filterQuality === 'medium' && !(comp >= 50 && cont < 10)) return false;
        if (filterQuality === 'low' && (comp >= 50 && cont < 10)) return false;
      }
      return true;
    });
  });

  let filteredBinSet = $derived(new Set(filteredBins.map(b => b.name)));

  // ---- Color function ----
  function hexToRgba(hex) {
    const r = parseInt(hex.slice(1, 3), 16);
    const g = parseInt(hex.slice(3, 5), 16);
    const b = parseInt(hex.slice(5, 7), 16);
    return [r, g, b, 255];
  }

  function brighten(rgba, factor = 1.4) {
    return [
      Math.min(255, Math.round(rgba[0] * factor)),
      Math.min(255, Math.round(rgba[1] * factor)),
      Math.min(255, Math.round(rgba[2] * factor)),
      255,
    ];
  }

  function getColor(bin) {
    if (colorBy === 'binner') {
      const idx = uniqueBinners.indexOf(bin.binner);
      return palette[idx >= 0 ? idx % palette.length : palette.length - 1];
    }
    if (colorBy === 'domain') {
      const domains = [...new Set(treeData.bins.map(b => b.lineage?.domain).filter(Boolean))];
      const idx = domains.indexOf(bin.lineage?.domain);
      return idx >= 0 ? palette[idx % palette.length] : '#475569';
    }
    if (colorBy === 'phylum') {
      const idx = uniquePhyla.indexOf(bin.lineage?.phylum);
      return idx >= 0 ? palette[idx % palette.length] : '#475569';
    }
    if (colorBy === 'quality') {
      const comp = bin.completeness ?? 0;
      const cont = bin.contamination ?? 100;
      if (comp >= 90 && cont < 5) return '#34d399';
      if (comp >= 50 && cont < 10) return '#fbbf24';
      return '#f87171';
    }
    return '#22d3ee';
  }

  // Valid binomial species name: "Genus epithet" where epithet is not sp\d+ or all-caps
  const spPlaceholder = /^sp\d/;
  function isValidBinomial(name) {
    const parts = name.split(' ');
    if (parts.length < 2) return false;
    const epithet = parts[1];
    if (spPlaceholder.test(epithet)) return false;
    if (epithet === epithet.toUpperCase() && epithet.length > 1) return false;
    return true;
  }

  const GOLD = [251, 191, 36, 255]; // amber-400

  // ---- Node styles for Phylocanvas (keyed by display_label) ----
  let nodeStyles = $derived.by(() => {
    const styles = {};
    // Style user bins
    for (const bin of (treeData?.bins || [])) {
      const key = bin.name;
      if (filteredBinSet.has(key)) {
        styles[key] = {
          fillColour: brighten(hexToRgba(getColor(bin))),
          shape: bin.binner === 'dastool' ? 'diamond' : 'circle',
          nodeSize: 1,
        };
      } else {
        styles[key] = {
          fillColour: [71, 85, 105, 128],
          shape: 'circle',
          nodeSize: 0.5,
        };
      }
    }
    // Style reference taxa with valid binomial names in gold
    const leafInfo = treeData?.leaf_info;
    if (leafInfo) {
      for (const name in leafInfo) {
        if (leafInfo[name].user_bin) continue;
        if (isValidBinomial(name)) {
          styles[name] = { fillColour: GOLD };
        }
      }
    }
    return styles;
  });

  // ---- Collapse ref-only subtrees ----
  let collapsedIds = $derived.by(() => {
    if (!collapseRefs || !currentNewick || !userBinSet.size) return [];
    return computeRefOnlyNodes(currentNewick, userBinSet);
  });

  function computeRefOnlyNodes(nwk, binKeySet) {
    const ids = [];
    let pos = 0;
    let internalId = 0;

    function parseNode() {
      if (pos >= nwk.length) return { hasUser: false, id: null };

      if (nwk[pos] === '(') {
        pos++;
        const myId = internalId++;
        let hasUser = false;
        const children = [];

        while (pos < nwk.length && nwk[pos] !== ')') {
          const child = parseNode();
          children.push(child);
          hasUser = hasUser || child.hasUser;
          if (pos < nwk.length && nwk[pos] === ',') pos++;
        }
        if (pos < nwk.length) pos++;

        while (pos < nwk.length && nwk[pos] !== ',' && nwk[pos] !== ')' && nwk[pos] !== ';') pos++;

        if (!hasUser && children.length > 0) {
          ids.push(String(myId));
        }
        return { hasUser, id: String(myId) };
      } else {
        let name = '';
        while (pos < nwk.length && nwk[pos] !== ':' && nwk[pos] !== ',' && nwk[pos] !== ')' && nwk[pos] !== ';') {
          name += nwk[pos];
          pos++;
        }
        if (pos < nwk.length && nwk[pos] === ':') {
          pos++;
          while (pos < nwk.length && nwk[pos] !== ',' && nwk[pos] !== ')' && nwk[pos] !== ';') pos++;
        }
        const cleanName = name.replace(/'/g, '').trim();
        return { hasUser: binKeySet.has(cleanName), id: null };
      }
    }

    try { parseNode(); } catch (e) { return []; }
    return ids;
  }

  // ---- Click info panel ----
  let clickedNode = $state(null);

  function handleNodeClick(nodeId, node) {
    // Look up bin by display_label or name
    const bin = binByName[nodeId];
    if (bin) {
      const quality = qualityByName[nodeId];
      clickedNode = { type: 'bin', nodeId, bin, quality };
      selectedMag.set(bin.name);
    } else {
      // Reference genome — check leaf_info
      const info = treeData?.leaf_info?.[nodeId];
      if (info) {
        clickedNode = { type: 'ref', nodeId, info };
      } else if (node && !node.isLeaf) {
        // Internal node
        clickedNode = { type: 'internal', nodeId, node };
      } else {
        clickedNode = { type: 'unknown', nodeId };
      }
    }
  }

  function closeInfoPanel() { clickedNode = null; }

  // Parse a GTDB lineage string like "d__Bacteria;p__Proteo;..." into readable rows
  function parseLineage(lineageStr) {
    if (!lineageStr) return [];
    const ranks = { d: 'Domain', p: 'Phylum', c: 'Class', o: 'Order', f: 'Family', g: 'Genus', s: 'Species' };
    return lineageStr.split(';').map(seg => {
      const [prefix, ...rest] = seg.trim().split('__');
      const name = rest.join('__');
      return { rank: ranks[prefix] || prefix, name: name || '-' };
    }).filter(r => r.name && r.name !== '-');
  }

  // ---- Table ----
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
    if (!filteredBins.length) return [];
    return filteredBins.map(b => ({
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

  let hasData = $derived(
    treeData?.bins?.length > 0 && treeData?.newick && Object.keys(treeData.newick).length > 0
  );
</script>

<div class="sticky top-14 z-10 bg-slate-950 pb-2 -mx-4 px-4 pt-1">
<div class="flex items-center gap-3 flex-wrap text-xs">
  {#if treeNames.length > 1}
    <span class="text-slate-400">Tree:</span>
    <select
      class="bg-slate-800 border border-slate-600 text-slate-200 text-xs rounded px-2 py-1"
      bind:value={selectedTree}
    >
      {#each treeNames as name}
        <option value={name}>{treeLabel(name)}</option>
      {/each}
    </select>
  {/if}

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

  <span class="text-slate-400">Binner:</span>
  <select class="bg-slate-800 border border-slate-600 text-slate-200 text-xs rounded px-2 py-1" bind:value={filterBinner}>
    <option value="all">All</option>
    {#each uniqueBinners as b}
      <option value={b}>{b}</option>
    {/each}
  </select>

  <span class="text-slate-400">Phylum:</span>
  <select class="bg-slate-800 border border-slate-600 text-slate-200 text-xs rounded px-2 py-1" bind:value={filterPhylum}>
    <option value="all">All</option>
    {#each uniquePhyla as p}
      <option value={p}>{p}</option>
    {/each}
  </select>

  <span class="text-slate-400">Quality:</span>
  <select class="bg-slate-800 border border-slate-600 text-slate-200 text-xs rounded px-2 py-1" bind:value={filterQuality}>
    <option value="all">All</option>
    <option value="high">High (C≥90, X&lt;5)</option>
    <option value="medium">Medium (C≥50, X&lt;10)</option>
    <option value="low">Low</option>
  </select>

  <button
    class="px-3 py-1 rounded-md border transition-colors text-center {collapseRefs ? 'border-cyan-400 bg-cyan-400/20 text-cyan-400' : 'border-slate-600 bg-slate-800 text-slate-400'}"
    onclick={() => collapseRefs = !collapseRefs}
    title="Collapse reference-only subtrees"
  >
    Collapse refs
  </button>

  {#if treeData?.bins}
    <span class="text-slate-500">{filteredBins.length}/{treeData.bins.length} bins</span>
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
  <div class="flex gap-4 mb-6">
    <!-- Tree panel (2/3 or full width) -->
    <div class="bg-slate-800 rounded-lg p-4 border border-slate-700 {clickedNode ? 'w-2/3' : 'w-full'} transition-all">
      <h3 class="text-sm font-medium text-slate-400 mb-2">
        GTDB-Tk Placement Tree — {treeLabel(selectedTree)} ({layoutLabel()})
      </h3>
      <PhylocanvasTree
        newick={currentNewick}
        treeType={treeType}
        styles={nodeStyles}
        selectedIds={[]}
        {collapsedIds}
        onNodeClick={handleNodeClick}
      />
      <div class="text-[10px] text-slate-500 mt-1 flex gap-4">
        <span>Scroll: zoom</span>
        <span>Shift+scroll: y-axis zoom</span>
        <span>Ctrl+scroll: x-axis zoom</span>
        <span>Click+drag: pan</span>
        <span>Click node: info</span>
      </div>
    </div>

    <!-- Info panel (right 1/3) -->
    {#if clickedNode}
      <div class="w-1/3 bg-slate-800 rounded-lg border border-slate-700 flex flex-col overflow-hidden">
        <div class="flex items-center justify-between px-4 py-3 border-b border-slate-700 shrink-0">
          <button class="text-slate-500 hover:text-slate-300 ml-auto text-lg leading-none" onclick={closeInfoPanel}>&times;</button>
        </div>
        <div class="px-4 py-3 overflow-y-auto flex-1">
          {#if clickedNode.type === 'bin'}
            {@const b = clickedNode.bin}
            {@const q = clickedNode.quality}
            <h3 class="text-cyan-400 font-semibold mb-1 text-sm">{b.name}</h3>
            <p class="text-xs text-slate-500 mb-3">{b.binner}{q?.is_dastool ? ' (DAS Tool consensus)' : ''}</p>
            <dl class="space-y-2 text-sm">
              {#if q?.completeness != null}
                <div class="flex justify-between">
                  <dt class="text-slate-400">Completeness</dt>
                  <dd class="font-mono text-emerald-400">{q.completeness}%</dd>
                </div>
              {/if}
              {#if q?.contamination != null}
                <div class="flex justify-between">
                  <dt class="text-slate-400">Contamination</dt>
                  <dd class="font-mono text-amber-400">{q.contamination}%</dd>
                </div>
              {/if}
              {#if q?.gc != null}
                <div class="flex justify-between">
                  <dt class="text-slate-400">GC Content</dt>
                  <dd class="font-mono">{(q.gc * 100).toFixed(1)}%</dd>
                </div>
              {/if}
              {#if q?.genome_size != null}
                <div class="flex justify-between">
                  <dt class="text-slate-400">Genome Size</dt>
                  <dd class="font-mono">{(q.genome_size / 1e6).toFixed(2)} Mbp</dd>
                </div>
              {/if}
              {#if q?.n50 != null}
                <div class="flex justify-between">
                  <dt class="text-slate-400">N50</dt>
                  <dd class="font-mono">{(q.n50 / 1e3).toFixed(1)} Kb</dd>
                </div>
              {/if}
              {#if q?.n_contigs != null}
                <div class="flex justify-between">
                  <dt class="text-slate-400">Contigs</dt>
                  <dd class="font-mono">{q.n_contigs}</dd>
                </div>
              {/if}
              <div class="border-t border-slate-700 pt-2 mt-2">
                <dt class="text-slate-400 mb-1 font-medium">GTDB-Tk Closest Match</dt>
                {#if b.closest_ref || b.closest_ani != null}
                  {#if b.closest_ref}
                    <div class="flex justify-between">
                      <dt class="text-slate-400">Reference</dt>
                      <dd class="font-mono text-slate-200 text-xs">{b.closest_ref}</dd>
                    </div>
                  {/if}
                  {#if b.closest_ani != null}
                    <div class="flex justify-between">
                      <dt class="text-slate-400">ANI</dt>
                      <dd class="font-mono text-cyan-400">{b.closest_ani}%</dd>
                    </div>
                  {/if}
                  {#if b.closest_af != null}
                    <div class="flex justify-between">
                      <dt class="text-slate-400">Align. Fraction</dt>
                      <dd class="font-mono">{b.closest_af}</dd>
                    </div>
                  {/if}
                  {#if b.method}
                    <div class="flex justify-between">
                      <dt class="text-slate-400">Method</dt>
                      <dd class="font-mono">{b.method}</dd>
                    </div>
                  {/if}
                {:else}
                  <dd class="text-xs text-slate-500 italic">No ANI match (classified by {b.method || 'placement only'})</dd>
                {/if}
              </div>
              {#if b.classification}
                <div class="border-t border-slate-700 pt-2 mt-2">
                  <dt class="text-slate-400 mb-1 font-medium">GTDB Taxonomy</dt>
                  <dd class="font-mono text-xs">
                    {#each parseLineage(b.classification) as row}
                      <div><span class="text-slate-500">{row.rank}:</span> {row.name}</div>
                    {/each}
                  </dd>
                </div>
              {/if}
              {#if q}
                <div class="border-t border-slate-700 pt-2 mt-2">
                  <dt class="text-slate-400 mb-1 font-medium">MGE Summary</dt>
                  <dd class="font-mono text-xs grid grid-cols-2 gap-1">
                    <span>Virus: {q.n_virus ?? 0}</span>
                    <span>Plasmid: {q.n_plasmid ?? 0}</span>
                    <span>Defense: {q.n_defense ?? 0}</span>
                    <span>Integron: {q.n_integron ?? 0}</span>
                  </dd>
                </div>
              {/if}
            </dl>
          {:else if clickedNode.type === 'ref'}
            {@const info = clickedNode.info}
            <h3 class="text-cyan-400 font-semibold mb-1 text-sm">{clickedNode.nodeId}</h3>
            <p class="text-xs text-slate-500 mb-3">Reference genome</p>
            <dl class="space-y-2 text-sm">
              {#if info.id}
                <div class="flex justify-between">
                  <dt class="text-slate-400">Accession</dt>
                  <dd class="font-mono text-slate-200">{info.id}</dd>
                </div>
              {/if}
              {#if info.lineage}
                <div class="border-t border-slate-700 pt-2 mt-2">
                  <dt class="text-slate-400 mb-1 font-medium">GTDB Taxonomy</dt>
                  <dd class="font-mono text-xs">
                    {#each parseLineage(info.lineage) as row}
                      <div><span class="text-slate-500">{row.rank}:</span> {row.name}</div>
                    {/each}
                  </dd>
                </div>
              {/if}
            </dl>
          {:else if clickedNode.type === 'internal'}
            {@const n = clickedNode.node}
            <h3 class="text-cyan-400 font-semibold mb-1 text-sm">Internal Node</h3>
            <p class="text-xs text-slate-500 mb-3">Subtree</p>
            <dl class="space-y-2 text-sm">
              {#if n.totalLeaves != null}
                <div class="flex justify-between">
                  <dt class="text-slate-400">Total Leaves</dt>
                  <dd class="font-mono">{n.totalLeaves}</dd>
                </div>
              {/if}
              {#if n.totalNodes != null}
                <div class="flex justify-between">
                  <dt class="text-slate-400">Total Nodes</dt>
                  <dd class="font-mono">{n.totalNodes}</dd>
                </div>
              {/if}
              {#if n.label}
                <div class="flex justify-between">
                  <dt class="text-slate-400">Label</dt>
                  <dd class="font-mono text-slate-200">{n.label}</dd>
                </div>
              {/if}
            </dl>
          {:else}
            <h3 class="text-slate-300 font-semibold mb-1 text-sm">{clickedNode.nodeId}</h3>
            <p class="text-xs text-slate-500">Unknown node type</p>
          {/if}
        </div>
      </div>
    {/if}
  </div>

  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
    <h3 class="text-sm font-medium text-slate-400 mb-2">Per-Bin Classification ({filteredBins.length} bins)</h3>
    <DataTable
      columns={tableColumns}
      rows={tableRows}
      maxHeight="500px"
    />
  </div>
{/if}
