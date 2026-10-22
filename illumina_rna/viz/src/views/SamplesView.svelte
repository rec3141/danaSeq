<script>
  import DataTable from '../components/ui/DataTable.svelte';
  import StatCard from '../components/layout/StatCard.svelte';
  import { samples, fmtInt, fmtPct, median } from '../stores/data.js';

  let rows = $derived(($samples?.rows ?? []).map((r, i) => ({ ...r, genes_detected: $samples.genes_detected?.[i] })));

  const columns = [
    { key: 'sample', label: 'Sample' },
    { key: 'reads_raw', label: 'Raw reads', render: v => fmtInt(v) },
    { key: 'reads_filtered', label: 'After QC', render: v => fmtInt(v) },
    { key: 'pct_rrna', label: 'rRNA filtered', render: v => fmtPct(v) },
    { key: 'pct_rrna_total', label: 'rRNA total', render: v => fmtPct(v) },
    { key: 'reads_norrna', label: 'Non-rRNA', render: v => fmtInt(v) },
    { key: 'reads_nohuman', label: 'Non-human', render: v => fmtInt(v) },
    { key: 'reads_mapped', label: 'Mapped', render: v => fmtInt(v) },
    { key: 'pct_mapped', label: 'Mapped %', render: v => fmtPct(v) },
    { key: 'pairs_assigned', label: 'Pairs on features', render: v => fmtInt(v) },
    { key: 'pairs_cds', label: 'Pairs on CDS', render: v => fmtInt(v) },
    { key: 'pairs_rrna_features', label: 'Pairs on rRNA loci', render: v => fmtInt(v) },
    { key: 'pairs_no_feature', label: 'No gene', render: v => fmtInt(v) },
    { key: 'pct_cds', label: 'On CDS %', render: v => fmtPct(v) },
    { key: 'genes_detected', label: 'Genes detected', render: v => fmtInt(v) },
  ];
  const sum = k => rows.reduce((s, r) => s + (r[k] ?? 0), 0);
</script>

{#if rows.length}
  <div class="grid grid-cols-2 md:grid-cols-4 gap-4 mb-6">
    <StatCard label="Raw Reads" value={fmtInt(sum('reads_raw'))} sub="all samples" color="cyan" />
    <StatCard label="Mapped Reads" value={fmtInt(sum('reads_mapped'))} color="cyan" />
    <StatCard label="Pairs on Genes" value={fmtInt(sum('pairs_assigned'))} color="emerald" />
    <StatCard label="Median Genes Detected" value={fmtInt(median(rows.map(r => r.genes_detected)))} color="slate" />
  </div>
  <div class="bg-slate-800 rounded-lg p-4 border border-slate-700">
    <h3 class="text-sm font-medium text-slate-400 mb-2">Read Accounting ({rows.length} samples)</h3>
    <DataTable {columns} {rows} idKey="sample" maxHeight="700px" />
    <p class="mt-2 text-xs text-slate-500">Reads count both mates; "Pairs on genes" counts fragments assigned by featureCounts. Mapped % is of non-rRNA, non-human reads.</p>
  </div>
{/if}
