<script>
  import { onMount } from 'svelte';

  let {
    newick = '',
    treeType = 'rc',
    styles = {},
    selectedIds = [],
    collapsedIds = [],
    onNodeClick = null,
    onSelectionChange = null,
  } = $props();

  let container;
  let tree = null;
  let prevNewick = null;

  // Tooltip state
  let tipText = $state('');
  let tipX = $state(0);
  let tipY = $state(0);
  let tipShow = $state(false);

  // Dark theme colours (RGBA 0-255)
  const DARK_THEME = {
    fillColour: [148, 163, 184, 255],    // slate-400
    strokeColour: [71, 85, 105, 255],     // slate-600
    fontColour: [100, 116, 139, 255],      // slate-500 (dim default for unnamed refs)
    highlightColour: [220, 225, 235, 180], // soft white, visible on edges
  };

  function createTree() {
    if (!container || !newick) return;

    const PhylocanvasGL = window.phylocanvas?.PhylocanvasGL;
    if (!PhylocanvasGL) {
      console.error('PhylocanvasGL not found on window.phylocanvas');
      return;
    }

    if (tree) {
      tree.destroy();
      tree = null;
    }

    tree = new PhylocanvasGL(container, {
      source: newick,
      type: treeType,
      interactive: true,
      showLabels: true,
      showLeafLabels: true,
      showShapes: true,
      fontSize: 12,
      nodeSize: 6,
      lineWidth: 1,
      padding: 24,
      ...DARK_THEME,
      styles,
      selectedIds,
      collapsedIds,
    });

    // Collect all descendant leaf IDs of an internal node.
    // Phylocanvas nodes are pre-order indexed; descendants of a node span
    // a contiguous slice of the array.
    function getDescendantLeafIds(node) {
      const leaves = [];
      function walk(n) {
        if (!n.children || n.children.length === 0 || n.isLeaf) {
          leaves.push(n.id);
          return;
        }
        for (const c of n.children) walk(c);
      }
      if (node.children) { for (const c of node.children) walk(c); }
      return leaves;
    }

    // Cache hovered node — user always hovers before clicking, so this gives
    // us the node object without calling pickNodeFromLayer before origHandleClick
    // (which would consume the pick and prevent Phylocanvas edge highlighting).
    let lastHoveredNode = null;

    // Override handleHover for our tooltip only (skip Phylocanvas's built-in one)
    tree.handleHover = function(info, event) {
      const node = tree.pickNodeFromLayer(info);
      if (node) {
        lastHoveredNode = node;
        const label = node.label ?? node.id;
        tipText = node.isLeaf ? label : `${label} (${node.totalLeaves} leaves)`;
        const rect = container.getBoundingClientRect();
        tipX = (event?.srcEvent?.clientX ?? 0) - rect.left;
        tipY = (event?.srcEvent?.clientY ?? 0) - rect.top - 16;
        tipShow = true;
      } else {
        tipShow = false;
      }
    };

    // Override handleClick — let origHandleClick run first so Phylocanvas
    // highlights edges, then use cached hover node for our callbacks.
    const origHandleClick = tree.handleClick.bind(tree);
    tree.handleClick = function(info, event) {
      origHandleClick(info, event);
      // For leaves, pickNodeFromLayer still works after origHandleClick;
      // for internal nodes it returns null, so fall back to last hovered node.
      const node = tree.pickNodeFromLayer(info) || lastHoveredNode;
      if (node && onNodeClick) {
        onNodeClick(node.id, node);
      }
      if (onSelectionChange) {
        if (node && !node.isLeaf) {
          const leafIds = getDescendantLeafIds(node);
          onSelectionChange(leafIds);
        } else {
          const sel = tree.props?.selectedIds || [];
          onSelectionChange([...sel]);
        }
      }
    };

    prevNewick = newick;
  }

  onMount(() => {
    // Resize tree canvas when container changes size (e.g. info panel opens)
    const ro = new ResizeObserver(() => {
      if (tree && container) {
        tree.resize(container.offsetWidth, container.offsetHeight);
      }
    });
    if (container) ro.observe(container);

    return () => {
      ro.disconnect();
      if (tree) { tree.destroy(); tree = null; }
    };
  });

  $effect(() => {
    const _deps = [newick, treeType, styles, collapsedIds];
    if (!container || !newick) return;

    if (!tree || newick !== prevNewick) {
      createTree();
    } else {
      // Don't pass selectedIds here — Phylocanvas manages selection internally
      // via clicks; pushing [] would clear the highlight on every re-render.
      tree.setProps({
        type: treeType,
        styles,
        collapsedIds,
      });
    }
  });
</script>

<div class="w-full relative overflow-hidden" style="height: 75vh; background: #0f172a;">
  <div bind:this={container} class="absolute inset-0 w-full h-full"></div>
  {#if tipShow}
    <div
      class="absolute pointer-events-none bg-slate-900/95 text-slate-200 text-xs px-2 py-1 rounded shadow-lg border border-slate-600 whitespace-nowrap z-10"
      style="left: {tipX}px; top: {tipY}px; transform: translate(-50%, -100%)"
    >
      {tipText}
    </div>
  {/if}
</div>
