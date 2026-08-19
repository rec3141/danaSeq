<script>
  import PhylocanvasTree from '../components/PhylocanvasTree.svelte';
  import {
    store,
    GROUP_COLORS, GROUP_HEX,
    buildTaxColorMap, getAsvColor, getEffectiveColorLevel, getClusterColor, hexToRgba255,
    makeTaxonMatcher,
    lineageOf,
    activeDb,
  } from '../stores/data.svelte.js';
  import { chrome } from '../stores/theme.svelte.js';

  let { filters = {} } = $props();

  let treeType = $derived(filters.treeLayout || 'rc');

  // Parse newick into a tree, prune, and re-serialize
  function pruneNewick(nwk, keepIds) {
    if (!nwk || !keepIds.size) return nwk;

    // Parse newick into a tree structure
    function parse(s) {
      let pos = 0;
      function readNode() {
        const node = { children: [], label: '', branchLength: null };
        if (s[pos] === '(') {
          pos++; // skip (
          node.children.push(readNode());
          while (s[pos] === ',') {
            pos++; // skip ,
            node.children.push(readNode());
          }
          pos++; // skip )
        }
        // Read label
        let label = '';
        while (pos < s.length && s[pos] !== ':' && s[pos] !== ',' && s[pos] !== ')' && s[pos] !== ';') {
          label += s[pos++];
        }
        node.label = label;
        // Read branch length
        if (pos < s.length && s[pos] === ':') {
          pos++; // skip :
          let bl = '';
          while (pos < s.length && s[pos] !== ',' && s[pos] !== ')' && s[pos] !== ';') {
            bl += s[pos++];
          }
          node.branchLength = parseFloat(bl);
        }
        return node;
      }
      const tree = readNode();
      return tree;
    }

    // Prune: remove tips not in keepIds, collapse single-child internals
    function prune(node) {
      if (node.children.length === 0) {
        // Leaf: keep only if in keepIds
        return keepIds.has(node.label) ? node : null;
      }
      // Recurse
      node.children = node.children.map(prune).filter(Boolean);
      if (node.children.length === 0) return null;
      if (node.children.length === 1) {
        // Single child: merge branch lengths
        const child = node.children[0];
        if (node.branchLength != null && child.branchLength != null) {
          child.branchLength += node.branchLength;
        }
        return child;
      }
      return node;
    }

    // Serialize back to newick
    function toNewick(node) {
      if (node.children.length === 0) {
        return node.label + (node.branchLength != null ? ':' + node.branchLength : '');
      }
      const inner = node.children.map(toNewick).join(',');
      const bl = node.branchLength != null ? ':' + node.branchLength : '';
      return `(${inner})${node.label}${bl}`;
    }

    try {
      const tree = parse(nwk);
      const pruned = prune(tree);
      if (!pruned) return nwk;
      return toNewick(pruned) + ';';
    } catch {
      return nwk;
    }
  }

  let displayNewick = $derived.by(() => {
    if (!filters.treePrune || filteredAsvIds.size === 0) return store.treeNewick;
    return pruneNewick(store.treeNewick, filteredAsvIds);
  });

  // ---- Taxonomy data ----
  let taxDb = $derived(activeDb() || null);
  let taxByAsv = $derived.by(() => {
    if (!taxDb || !store.taxonomy[taxDb]?.assignments) return {};
    return store.taxonomy[taxDb].assignments;
  });
  let taxLevels = $derived(
    taxDb && store.taxonomy[taxDb]?.levels ? store.taxonomy[taxDb].levels : []
  );

  // ---- Taxonomy filter predicate ----
  let taxonMatches = $derived.by(() =>
    makeTaxonMatcher(filters.taxonFilter, filters.taxonFilterLevel)
  );

  // Bootstrap data
  let bootByAsv = $derived.by(() => {
    if (!taxDb || !store.taxonomy[taxDb]?.bootstraps) return {};
    return store.taxonomy[taxDb].bootstraps;
  });

  // ---- Filtered ASV set (by taxonomy, group, bootstrap) ----
  let filteredAsvIds = $derived.by(() => {
    const matches = taxonMatches;
    const gf = filters.groupFlags || {};
    const minBoot = filters.treeMinBootstrap || 0;
    const assayAsvs = store.assayAsvIds;
    const ids = new Set();
    for (const asv of store.asvs) {
      if (assayAsvs && !assayAsvs.has(asv.id)) continue;
      const group = asv.group ?? 'unknown';
      if (gf[group] === false) continue;
      if (matches && !matches(asv.id ?? '', lineageOf(asv.id, asv.taxonomy))) continue;
      // Check bootstrap at deepest assigned rank
      if (minBoot > 0) {
        const boot = bootByAsv[asv.id];
        const tax = taxByAsv[asv.id];
        if (boot && tax) {
          let deepestBoot = 0;
          for (let i = tax.length - 1; i >= 0; i--) {
            if (tax[i]) { deepestBoot = boot[i] || 0; break; }
          }
          if (deepestBoot < minBoot) continue;
        }
      }
      ids.add(asv.id);
    }
    return ids;
  });

  // ---- Node styles using shared color-by ----
  let effectiveColorLevel = $derived(
    filters.colorMode === 'group' ? 'group' : getEffectiveColorLevel(filters.colorByLevel, filters.taxonFilter, filters.taxonFilterLevel)
  );

  let taxColorMap = $derived.by(() => {
    if (effectiveColorLevel === 'group') return null;
    return buildTaxColorMap(effectiveColorLevel, filters.taxonFilter, filters.taxonFilterLevel);
  });

  // Build label for an ASV based on selected label levels.
  //
  // Bootstrap is a support value *for a rank*, not for the ASV, so it belongs
  // beside each rank shown rather than once at the end: a genus called at 62 and
  // the phylum above it called at 100 are two different statements, and a single
  // number cannot say which one it is about. Selecting Phylum, Genus and
  // bootstrap gives `Bacteroidota (100) | Polaribacter (62)`.
  function makeLabel(asvId) {
    const levels = filters.treeLabelLevels || ['id'];
    const tax = taxByAsv[asvId];
    const boot = bootByAsv[asvId];
    const withBoot = levels.includes('bootstrap');
    const ranks = levels.filter(l => l !== 'id' && l !== 'bootstrap');
    const parts = [];

    for (const level of levels) {
      if (level === 'id') {
        parts.push(asvId);
      } else if (level === 'bootstrap') {
        // Only meaningful on its own when no rank was chosen to attach it to;
        // then it still answers "how well is this classified at all?".
        if (boot && ranks.length === 0) {
          for (let i = boot.length - 1; i >= 0; i--) {
            if (tax && tax[i]) { parts.push(`(${boot[i]})`); break; }
          }
        }
      } else if (tax) {
        const idx = taxLevels.indexOf(level);
        if (idx >= 0 && tax[idx]) {
          const b = withBoot && boot ? boot[idx] : null;
          parts.push(b == null || b === '' ? tax[idx] : `${tax[idx]} (${b})`);
        }
      }
    }

    return parts.length > 0 ? parts.join(' | ') : asvId;
  }

  let nodeStyles = $derived.by(() => {
    const styles = {};
    const cmap = taxColorMap?.colorMap;

    for (const asv of store.asvs) {
      const id = asv.id;
      const isFiltered = filteredAsvIds.has(id);
      const label = makeLabel(id);

      if (!isFiltered) {
        styles[id] = {
          fillColour: chrome().collapsed,
          fontColour: chrome().collapsed,
          shape: 'circle',
          nodeSize: 0.3,
          label,
        };
        continue;
      }

      let hex;
      if (filters.colorMode === 'cluster') {
        hex = getClusterColor(id, 'asvCluster', filters.asvClusterK);
      } else if (cmap) {
        hex = getAsvColor(id, effectiveColorLevel, cmap);
      } else {
        hex = GROUP_HEX[asv.group] || GROUP_HEX.unknown;
      }

      styles[id] = {
        fillColour: hex,
        fontColour: hex,
        shape: 'circle',
        nodeSize: 1,
        label,
      };
    }
    return styles;
  });

  // ---- ASV lookup ----
  let asvById = $derived.by(() => {
    const map = {};
    for (const a of store.asvs) map[a.id] = a;
    return map;
  });

  // ---- Click handling ----
  let clickedNode = $state(null);

  function handleNodeClick(nodeId) {
    const asv = asvById[nodeId];
    if (asv) {
      const tax = taxByAsv[nodeId];
      clickedNode = { id: nodeId, asv, tax };
    } else {
      clickedNode = { id: nodeId, asv: null, tax: null };
    }
  }

  function closeInfoPanel() { clickedNode = null; }

  // Export what is on screen, not what was loaded: with pruning on, the tree you
  // are looking at is a subset, and that is the one worth taking away.
  function exportNewick() {
    const nwk = (displayNewick || '').trim();
    if (!nwk) return;
    const text = nwk.endsWith(';') ? nwk : nwk + ';';
    const blob = new Blob([text], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    const pruned = filters.treePrune && filteredAsvIds.size > 0;
    a.href = url;
    a.download = pruned ? 'tree.filtered.nwk' : 'tree.nwk';
    a.click();
    URL.revokeObjectURL(url);
  }
</script>

<div class="flex h-full flex-col">

  {#if store.treeNewick}
    <div class="flex items-center gap-2 px-3 pt-2">
      <button
        onclick={exportNewick}
        class="rounded border border-line px-2 py-1 text-xs text-fg2
               hover:border-line2 hover:text-strong"
        title="Download this tree in Newick format">
        Export Newick{filters.treePrune && filteredAsvIds.size > 0 ? ' (filtered)' : ''}
      </button>
      <span class="text-xs text-faint">
        {filteredAsvIds.size > 0 && filters.treePrune
          ? `${filteredAsvIds.size.toLocaleString()} tips`
          : `${store.asvs.length.toLocaleString()} tips`}
      </span>
    </div>
  {/if}

  <!-- Tree + info panel -->
  <div class="flex flex-1 overflow-hidden">
    {#if !displayNewick}
      <div class="flex-1 flex items-center justify-center">
        <div class="text-center">
          <p class="text-muted mb-2">No phylogenetic tree available</p>
          <p class="text-xs text-faint">Run with <code class="bg-raised px-1.5 py-0.5 rounded text-link">--run_phylogeny</code></p>
        </div>
      </div>
    {:else}
      <div class="{clickedNode ? 'w-2/3' : 'w-full'} transition-all">
        <PhylocanvasTree
          newick={displayNewick}
          {treeType}
          styles={nodeStyles}
          onNodeClick={handleNodeClick}
        />
      </div>

      {#if clickedNode}
        <div class="w-1/3 border-l border-line bg-surface/60 overflow-y-auto p-4">
          <div class="flex items-center justify-between mb-3">
            <h3 class="text-sm font-semibold text-link">{clickedNode.id}</h3>
            <button class="text-faint hover:text-fg2 text-lg" onclick={closeInfoPanel}>&times;</button>
          </div>

          {#if clickedNode.asv}
            <dl class="space-y-2 text-sm">
              <div class="flex justify-between">
                <dt class="text-muted">Group</dt>
                <dd>
                  <span class="inline-block h-2 w-2 rounded-full mr-1" style="background:{GROUP_HEX[clickedNode.asv.group] || GROUP_HEX.unknown}"></span>
                  {clickedNode.asv.group || 'unknown'}
                </dd>
              </div>
              <div class="flex justify-between">
                <dt class="text-muted">Total Reads</dt>
                <dd class="font-mono">{(clickedNode.asv.total_reads || 0).toLocaleString()}</dd>
              </div>
              <div class="flex justify-between">
                <dt class="text-muted">Prevalence</dt>
                <dd class="font-mono">{clickedNode.asv.n_samples || 0} samples</dd>
              </div>

              {#if clickedNode.tax}
                <div class="border-t border-line pt-2 mt-2">
                  <dt class="text-muted mb-1 font-medium">Taxonomy</dt>
                  <dd class="font-mono text-xs">
                    {#each taxLevels as level, i}
                      {#if clickedNode.tax[i]}
                        <div><span class="text-faint">{level}:</span> {clickedNode.tax[i]}</div>
                      {/if}
                    {/each}
                  </dd>
                </div>
              {/if}
            </dl>
          {:else}
            <p class="text-xs text-faint">Internal node</p>
          {/if}
        </div>
      {/if}
    {/if}
  </div>
</div>
