<script>
  import MatrixPanel from '../components/MatrixPanel.svelte';
  import { bins } from '../stores/data.js';

  // Last GTDB rank that is filled in, e.g. "g__Pelagibacter"
  const shortTax = c => {
    if (!c) return null;
    const ranks = c.split(';').filter(x => x && !/^[a-z]__$/.test(x));
    return ranks.length ? ranks[ranks.length - 1] : null;
  };

  let rows = $derived(($bins?.bins ?? []).map(b => ({
    ...b, id: b.bin, taxon: shortTax(b.classification),
    size_mb: b.length ? b.length / 1e6 : null,
  })));

  const columns = [
    { key: 'taxon', label: 'Taxonomy (GTDB)' },
    { key: 'completeness', label: 'Compl. %', align: 'right', fmt: v => v == null ? '—' : v.toFixed(1) },
    { key: 'contamination', label: 'Contam. %', align: 'right', fmt: v => v == null ? '—' : v.toFixed(1) },
    { key: 'size_mb', label: 'Size Mb', align: 'right', fmt: v => v == null ? '—' : v.toFixed(2) },
    { key: 'n_genes', label: 'Genes', align: 'right', fmt: v => v == null ? '—' : v.toLocaleString() },
  ];
</script>

<MatrixPanel {rows} samples={$bins?.samples ?? []} {columns} idLabel="MAG"
  labelOf={r => r.taxon ? `${r.id} · ${r.taxon}` : r.id}
  title="Most active MAGs"
  subtitle="Sum of gene TPM per MAG (DAS Tool bins); 'unbinned' is expression on contigs in no MAG. Click a row or cell for its profile across samples." />
