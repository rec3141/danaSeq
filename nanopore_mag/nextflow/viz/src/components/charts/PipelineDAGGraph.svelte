<script>
  import { onMount } from 'svelte';

  let { processes = {} } = $props();

  // Real pipeline DAG: nodes and directed edges from main.nf
  // Each node: { id, label, x, y } — positioned by layout below
  // Each edge: { from, to }
  const nodes = [
    // Tier 0: Input
    { id: 'CONCAT_READS', label: 'Concat Reads', tier: 0 },
    // Tier 1: Assembly
    { id: 'ASSEMBLY_FLYE', label: 'Assembly', tier: 1 },
    // Tier 2: Post-assembly fan-out
    { id: 'MAP_READS', label: 'Mapping', tier: 2 },
    { id: 'CALCULATE_TNF', label: 'TNF', tier: 2 },
    { id: 'BAKTA_BASIC', label: 'Bakta', tier: 2 },
    { id: 'GENOMAD_CLASSIFY', label: 'geNomad', tier: 2 },
    { id: 'INTEGRONFINDER', label: 'Integron', tier: 2 },
    { id: 'KRAKEN2_CLASSIFY', label: 'Kraken2', tier: 2 },
    { id: 'SENDSKETCH_CLASSIFY', label: 'BBSketch', tier: 2 },
    { id: 'TIARA_CLASSIFY', label: 'Tiara', tier: 2 },
    { id: 'WHOKARYOTE_CLASSIFY', label: 'Whokaryote', tier: 2 },
    // Tier 3: Depends on mapping or annotation
    { id: 'CALCULATE_DEPTHS', label: 'Depths', tier: 3 },
    { id: 'CALCULATE_GENE_DEPTHS', label: 'Gene Depths', tier: 3 },
    { id: 'BAKTA_EXTRA', label: 'Bakta Extra', tier: 3 },
    { id: 'RNA_CLASSIFY', label: 'rRNA', tier: 3 },
    { id: 'CHECKV_QUALITY', label: 'CheckV', tier: 3 },
    { id: 'KAIJU_CLASSIFY', label: 'Kaiju', tier: 3 },
    { id: 'METAEUK_PREDICT', label: 'MetaEuk', tier: 3 },
    // Tier 4: Depends on annotation
    { id: 'KAIJU_CONTIG_CLASSIFY', label: 'Kaiju Ctg', tier: 4 },
    { id: 'ISLANDPATH_DIMOB', label: 'IslandPath', tier: 4 },
    { id: 'MACSYFINDER', label: 'MacSy', tier: 4 },
    { id: 'DEFENSEFINDER', label: 'Defense', tier: 4 },
    { id: 'KOFAMSCAN', label: 'KofamScan', tier: 4 },
    { id: 'EMAPPER', label: 'eggNOG', tier: 4 },
    { id: 'DBCAN', label: 'dbCAN', tier: 4 },
    { id: 'MARFERRET_CLASSIFY', label: 'MarFERReT', tier: 4 },
    // Tier 5: Binning (depends on depths + TNF)
    { id: 'BIN_SEMIBIN2', label: 'SemiBin2', tier: 5 },
    { id: 'BIN_METABAT2', label: 'MetaBAT2', tier: 5 },
    { id: 'BIN_MAXBIN2', label: 'MaxBin2', tier: 5 },
    { id: 'BIN_LORBIN', label: 'LorBin', tier: 5 },
    { id: 'BIN_COMEBIN', label: 'COMEBin', tier: 5 },
    { id: 'MERGE_ANNOTATIONS', label: 'Merge Annot', tier: 5 },
    // Tier 6: Consensus + downstream
    { id: 'DASTOOL_CONSENSUS', label: 'DAS Tool', tier: 6 },
    { id: 'MAP_TO_BINS', label: 'Map2Bins', tier: 6 },
    { id: 'KEGG_MODULES', label: 'KEGG Mod', tier: 6 },
    // Tier 7: Final QC + Viz
    { id: 'CHECKM2', label: 'CheckM2', tier: 7 },
    { id: 'VIZ_STAGE1', label: 'Viz Stage 1', tier: 7 },
    { id: 'VIZ_STAGE2', label: 'Viz Stage 2', tier: 7 },
  ];

  // Directed edges: [from, to]
  const edges = [
    // Input → Assembly
    ['CONCAT_READS', 'ASSEMBLY_FLYE'],
    // Assembly → fan-out
    ['ASSEMBLY_FLYE', 'MAP_READS'],
    ['ASSEMBLY_FLYE', 'CALCULATE_TNF'],
    ['ASSEMBLY_FLYE', 'BAKTA_BASIC'],
    ['ASSEMBLY_FLYE', 'GENOMAD_CLASSIFY'],
    ['ASSEMBLY_FLYE', 'INTEGRONFINDER'],
    ['ASSEMBLY_FLYE', 'KRAKEN2_CLASSIFY'],
    ['ASSEMBLY_FLYE', 'SENDSKETCH_CLASSIFY'],
    ['ASSEMBLY_FLYE', 'TIARA_CLASSIFY'],
    ['ASSEMBLY_FLYE', 'WHOKARYOTE_CLASSIFY'],
    // Mapping → depths
    ['MAP_READS', 'CALCULATE_DEPTHS'],
    ['MAP_READS', 'CALCULATE_GENE_DEPTHS'],
    // Assembly → RNA, MetaEuk (need assembly only)
    ['ASSEMBLY_FLYE', 'RNA_CLASSIFY'],
    ['ASSEMBLY_FLYE', 'METAEUK_PREDICT'],
    // Annotation → downstream
    ['BAKTA_BASIC', 'BAKTA_EXTRA'],
    ['BAKTA_BASIC', 'KAIJU_CLASSIFY'],
    ['BAKTA_BASIC', 'ISLANDPATH_DIMOB'],
    ['BAKTA_BASIC', 'MACSYFINDER'],
    ['BAKTA_BASIC', 'DEFENSEFINDER'],
    ['BAKTA_BASIC', 'KOFAMSCAN'],
    ['BAKTA_BASIC', 'EMAPPER'],
    ['BAKTA_BASIC', 'DBCAN'],
    // Kaiju gene → contig
    ['KAIJU_CLASSIFY', 'KAIJU_CONTIG_CLASSIFY'],
    // geNomad → CheckV
    ['GENOMAD_CLASSIFY', 'CHECKV_QUALITY'],
    // MetaEuk → MarFERReT
    ['METAEUK_PREDICT', 'MARFERRET_CLASSIFY'],
    // Depths + TNF → Binning
    ['CALCULATE_DEPTHS', 'BIN_SEMIBIN2'],
    ['CALCULATE_DEPTHS', 'BIN_METABAT2'],
    ['CALCULATE_DEPTHS', 'BIN_MAXBIN2'],
    ['CALCULATE_DEPTHS', 'BIN_LORBIN'],
    ['CALCULATE_DEPTHS', 'BIN_COMEBIN'],
    // Metabolism merge
    ['KOFAMSCAN', 'MERGE_ANNOTATIONS'],
    ['EMAPPER', 'MERGE_ANNOTATIONS'],
    ['DBCAN', 'MERGE_ANNOTATIONS'],
    // Binning → consensus
    ['BIN_SEMIBIN2', 'DASTOOL_CONSENSUS'],
    ['BIN_METABAT2', 'DASTOOL_CONSENSUS'],
    ['BIN_MAXBIN2', 'DASTOOL_CONSENSUS'],
    ['BIN_LORBIN', 'DASTOOL_CONSENSUS'],
    ['BIN_COMEBIN', 'DASTOOL_CONSENSUS'],
    // Consensus → QC
    ['DASTOOL_CONSENSUS', 'CHECKM2'],
    // Map annotations to bins
    ['MERGE_ANNOTATIONS', 'MAP_TO_BINS'],
    ['DASTOOL_CONSENSUS', 'MAP_TO_BINS'],
    // KEGG modules
    ['MAP_TO_BINS', 'KEGG_MODULES'],
    // Viz barriers
    ['DASTOOL_CONSENSUS', 'VIZ_STAGE1'],
    ['CHECKM2', 'VIZ_STAGE1'],
    ['VIZ_STAGE1', 'VIZ_STAGE2'],
  ];

  // Only include nodes that have a status in the processes map
  let activeNodeIds = $derived(new Set(
    Object.keys(processes).length > 0
      ? nodes.filter(n => processes[n.id]).map(n => n.id)
      : nodes.map(n => n.id)
  ));

  let activeNodes = $derived(nodes.filter(n => activeNodeIds.has(n.id)));
  let activeEdges = $derived(edges.filter(([from, to]) => activeNodeIds.has(from) && activeNodeIds.has(to)));

  // Layout: position nodes by tier
  const nodeW = 78;
  const nodeH = 24;
  const tierGap = 44;
  const nodeGapX = 8;
  const padX = 16;
  const padY = 16;

  let layout = $derived.by(() => {
    // Group active nodes by tier
    const tierMap = {};
    for (const n of activeNodes) {
      (tierMap[n.tier] ??= []).push(n);
    }
    const tierKeys = Object.keys(tierMap).map(Number).sort((a, b) => a - b);

    const positions = {};
    let maxTierWidth = 0;
    for (const tk of tierKeys) {
      const tierNodes = tierMap[tk];
      const tw = tierNodes.length * nodeW + (tierNodes.length - 1) * nodeGapX;
      if (tw > maxTierWidth) maxTierWidth = tw;
    }

    const svgWidth = maxTierWidth + padX * 2;
    let y = padY;
    for (const tk of tierKeys) {
      const tierNodes = tierMap[tk];
      const tw = tierNodes.length * nodeW + (tierNodes.length - 1) * nodeGapX;
      const startX = padX + (maxTierWidth - tw) / 2;
      tierNodes.forEach((n, i) => {
        positions[n.id] = {
          x: startX + i * (nodeW + nodeGapX),
          y,
          cx: startX + i * (nodeW + nodeGapX) + nodeW / 2,
          cy: y + nodeH / 2,
        };
      });
      y += nodeH + tierGap;
    }
    const svgHeight = y - tierGap + padY;
    return { positions, svgWidth, svgHeight };
  });

  // Color by status
  function statusColor(proc) {
    const s = processes[proc];
    if (s === 'completed' || s === 'stored') return { fill: 'rgba(34,211,238,0.15)', stroke: '#22d3ee', text: '#22d3ee' };
    if (s === 'running')   return { fill: 'rgba(251,191,36,0.15)', stroke: '#fbbf24', text: '#fbbf24' };
    if (s === 'failed')    return { fill: 'rgba(248,113,113,0.15)', stroke: '#f87171', text: '#f87171' };
    return { fill: 'rgba(100,116,139,0.1)', stroke: '#475569', text: '#64748b' };
  }

  function edgeColor(from, to) {
    const sf = processes[from];
    const st = processes[to];
    if (sf === 'completed' || sf === 'stored') {
      if (st === 'completed' || st === 'stored') return '#22d3ee';
      if (st === 'running') return '#fbbf24';
    }
    return '#334155';
  }

  // Curved path between two nodes (top-center → bottom-center)
  function edgePath(from, to) {
    const p1 = layout.positions[from];
    const p2 = layout.positions[to];
    if (!p1 || !p2) return '';
    const x1 = p1.cx, y1 = p1.y + nodeH;
    const x2 = p2.cx, y2 = p2.y;
    const dy = (y2 - y1) * 0.5;
    return `M${x1},${y1} C${x1},${y1 + dy} ${x2},${y2 - dy} ${x2},${y2}`;
  }
</script>

<div class="overflow-x-auto">
  <svg width={layout.svgWidth} height={layout.svgHeight} class="mx-auto">
    <defs>
      <marker id="arrow" viewBox="0 0 6 6" refX="6" refY="3"
              markerWidth="5" markerHeight="5" orient="auto-start-reverse">
        <path d="M0,0 L6,3 L0,6 Z" fill="#475569" />
      </marker>
    </defs>

    <!-- Edges -->
    {#each activeEdges as [from, to]}
      <path
        d={edgePath(from, to)}
        fill="none"
        stroke={edgeColor(from, to)}
        stroke-width="1"
        opacity="0.5"
        marker-end="url(#arrow)"
      />
    {/each}

    <!-- Nodes -->
    {#each activeNodes as node}
      {@const pos = layout.positions[node.id]}
      {@const colors = statusColor(node.id)}
      {#if pos}
        <rect
          x={pos.x} y={pos.y} width={nodeW} height={nodeH} rx="4"
          fill={colors.fill} stroke={colors.stroke} stroke-width="1"
          class={processes[node.id] === 'running' ? 'animate-pulse' : ''}
        />
        <text
          x={pos.cx} y={pos.cy + 1}
          text-anchor="middle" dominant-baseline="middle"
          fill={colors.text}
          class="text-[8px] font-medium select-none"
        >{node.label}</text>
      {/if}
    {/each}
  </svg>
</div>

<style>
  @keyframes pulse {
    0%, 100% { opacity: 1; }
    50% { opacity: 0.5; }
  }
  .animate-pulse {
    animation: pulse 2s ease-in-out infinite;
  }
</style>
