<script>
  import MatrixPanel from '../components/MatrixPanel.svelte';
  import { topGenes } from '../stores/data.js';

  let rows = $derived($topGenes?.genes ?? []);
  let columns = $derived([
    { key: 'product', label: 'Product', wrap: true },
    { key: 'gene', label: 'Gene' },
    ...(rows.some(r => r.ko) ? [{ key: 'ko', label: 'KO' }] : []),
    { key: 'bin', label: 'MAG' },
    { key: 'type', label: 'Type' },
  ]);
</script>

<MatrixPanel {rows} samples={$topGenes?.samples ?? []} {columns} idLabel="Locus tag"
  labelOf={r => `${r.id}${r.gene ? ' ' + r.gene : ''} · ${(r.product || '').slice(0, 40)}`}
  title="Most expressed genes"
  subtitle="Top genes by mean TPM across samples (Bakta locus tags; full matrices are in expression/<ref>/)." />
