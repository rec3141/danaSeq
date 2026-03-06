<script>
  import { onMount } from 'svelte';

  let {
    newick = '',
    treeType = 'rc',
    styles = {},
    selectedIds = [],
    collapsedIds = [],
    onNodeClick = null,
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
    highlightColour: [34, 211, 238, 255], // cyan-400
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

    // Override handleClick to call our callback after default behavior
    const origHandleClick = tree.handleClick.bind(tree);
    tree.handleClick = function(info, event) {
      // Run default: pickNodeFromLayer + selectNode
      origHandleClick(info, event);
      // Then fire our callback with the picked node
      if (onNodeClick) {
        const node = tree.pickNodeFromLayer(info);
        if (node) {
          onNodeClick(node.id, node);
        }
      }
    };

    // Override handleHover for tooltip
    const origHandleHover = tree.handleHover.bind(tree);
    tree.handleHover = function(info, event) {
      origHandleHover(info, event);
      const node = tree.pickNodeFromLayer(info);
      if (node) {
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
    const _deps = [newick, treeType, styles, selectedIds, collapsedIds];
    if (!container || !newick) return;

    if (!tree || newick !== prevNewick) {
      createTree();
    } else {
      tree.setProps({
        type: treeType,
        styles,
        selectedIds,
        collapsedIds,
      });
    }
  });
</script>

<div class="w-full relative overflow-hidden" style="min-height: 600px; background: #0f172a;">
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
