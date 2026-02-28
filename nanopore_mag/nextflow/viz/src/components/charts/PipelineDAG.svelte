<script>
  let { processes = {} } = $props();

  // DAG tiers: each tier has a label and list of [processName, displayLabel]
  const tiers = [
    { label: 'Input',       nodes: [['CONCAT_READS', 'Concat Reads']] },
    { label: 'Assembly',    nodes: [['ASSEMBLY_FLYE', 'Assembly'], ['CALCULATE_TNF', 'TNF']] },
    { label: 'Map+Annot',   nodes: [['MAP_READS', 'Mapping'], ['CALCULATE_DEPTHS', 'Depths'],
                                     ['PROKKA_ANNOTATE', 'Prokka'], ['BAKTA_BASIC', 'Bakta'],
                                     ['BAKTA_EXTRA', 'Bakta Extra']] },
    { label: 'Classify',    nodes: [['KAIJU_CLASSIFY', 'Kaiju'], ['KAIJU_CONTIG_CLASSIFY', 'Kaiju Ctg'],
                                     ['KRAKEN2_CLASSIFY', 'Kraken2'], ['SENDSKETCH_CLASSIFY', 'BBSketch'],
                                     ['RNA_CLASSIFY', 'rRNA'], ['GENOMAD_CLASSIFY', 'geNomad'],
                                     ['TIARA_CLASSIFY', 'Tiara'], ['WHOKARYOTE_CLASSIFY', 'Whokaryote']] },
    { label: 'MGE+Euk',     nodes: [['CHECKV_QUALITY', 'CheckV'], ['INTEGRONFINDER', 'Integron'],
                                     ['ISLANDPATH_DIMOB', 'IslandPath'], ['MACSYFINDER', 'MacSy'],
                                     ['DEFENSEFINDER', 'Defense'], ['METAEUK_PREDICT', 'MetaEuk'],
                                     ['MARFERRET_CLASSIFY', 'MarFERReT']] },
    { label: 'Binning',     nodes: [['BIN_SEMIBIN2', 'SemiBin2'], ['BIN_METABAT2', 'MetaBAT2'],
                                     ['BIN_MAXBIN2', 'MaxBin2'], ['BIN_LORBIN', 'LorBin'],
                                     ['BIN_COMEBIN', 'COMEBin']] },
    { label: 'Consensus',   nodes: [['DASTOOL_CONSENSUS', 'DAS Tool'], ['CHECKM2', 'CheckM2']] },
    { label: 'Metabolism',  nodes: [['KOFAMSCAN', 'KofamScan'], ['EMAPPER', 'eggNOG'], ['DBCAN', 'dbCAN'],
                                     ['MERGE_ANNOTATIONS', 'Merge'], ['MAP_TO_BINS', 'Map2Bins'],
                                     ['KEGG_MODULES', 'KEGG Mod'], ['KEGG_DECODER', 'KEGG Dec'],
                                     ['MINPATH', 'MinPath']] },
  ];

  // Layout constants
  const nodeW = 80;
  const nodeH = 26;
  const tierGap = 36;
  const nodeGapX = 6;
  const padX = 12;
  const padY = 8;
  const labelW = 72;

  // Compute max nodes per tier for width
  const maxNodes = Math.max(...tiers.map(t => t.length));

  // Compute SVG dimensions
  function getTierWidth(tier) {
    return tier.nodes.length * nodeW + (tier.nodes.length - 1) * nodeGapX;
  }
  const svgWidth = Math.max(...tiers.map(getTierWidth)) + labelW + padX * 3;
  const svgHeight = tiers.length * (nodeH + tierGap) - tierGap + padY * 2;

  // Color map
  function statusColor(proc) {
    const s = processes[proc];
    if (s === 'completed') return { fill: 'rgba(34,211,238,0.15)', stroke: '#22d3ee', text: '#22d3ee' };
    if (s === 'running')   return { fill: 'rgba(251,191,36,0.15)', stroke: '#fbbf24', text: '#fbbf24' };
    if (s === 'failed')    return { fill: 'rgba(248,113,113,0.15)', stroke: '#f87171', text: '#f87171' };
    if (s === 'warning')   return { fill: 'rgba(251,191,36,0.15)', stroke: '#f59e0b', text: '#f59e0b' };
    // pending or unknown
    return { fill: 'rgba(100,116,139,0.1)', stroke: '#475569', text: '#64748b' };
  }
</script>

<div class="overflow-x-auto">
  <svg width={svgWidth} height={svgHeight} class="mx-auto">
    {#each tiers as tier, tierIdx}
      {@const y = padY + tierIdx * (nodeH + tierGap)}
      {@const tierW = getTierWidth(tier)}
      {@const startX = labelW + padX * 2 + (svgWidth - labelW - padX * 3 - tierW) / 2}

      <!-- Tier label -->
      <text x={padX + labelW - 4} y={y + nodeH / 2 + 1}
            text-anchor="end" dominant-baseline="middle"
            class="fill-slate-500 text-[10px]">{tier.label}</text>

      <!-- Connector lines to next tier -->
      {#if tierIdx < tiers.length - 1}
        {@const nextTier = tiers[tierIdx + 1]}
        {@const nextY = padY + (tierIdx + 1) * (nodeH + tierGap)}
        {@const nextTierW = getTierWidth(nextTier)}
        {@const nextStartX = labelW + padX * 2 + (svgWidth - labelW - padX * 3 - nextTierW) / 2}
        <!-- Simple center-to-center line between tiers -->
        <line
          x1={startX + tierW / 2} y1={y + nodeH}
          x2={nextStartX + nextTierW / 2} y2={nextY}
          stroke="#334155" stroke-width="1" opacity="0.4" />
      {/if}

      <!-- Nodes -->
      {#each tier.nodes as [proc, label], nodeIdx}
        {@const x = startX + nodeIdx * (nodeW + nodeGapX)}
        {@const colors = statusColor(proc)}
        <rect
          x={x} y={y} width={nodeW} height={nodeH} rx="4"
          fill={colors.fill} stroke={colors.stroke} stroke-width="1"
          class={processes[proc] === 'running' ? 'animate-pulse' : ''} />
        <text
          x={x + nodeW / 2} y={y + nodeH / 2 + 1}
          text-anchor="middle" dominant-baseline="middle"
          fill={colors.text}
          class="text-[9px] font-medium select-none">{label}</text>
      {/each}
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
