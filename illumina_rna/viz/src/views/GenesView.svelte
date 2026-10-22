<script>
  import ExpressionMatrix from '../components/ExpressionMatrix.svelte';
  import { topGenes, fmtInt } from '../stores/data.js';

  let rows = $derived($topGenes?.genes ?? []);
  let hasKo = $derived(rows.some(r => r.ko));
  const TYPES = ['CDS', 'rRNA', 'tRNA', 'tmRNA', 'ncRNA'].map(k => ({ key: k, label: k }));

  let columns = $derived([
    { key: 'id', label: 'Locus tag' },
    { key: 'gene', label: 'Gene' },
    { key: 'product', label: 'Product', render: v => v ? (v.length > 60 ? v.slice(0, 58) + '…' : v) : '-' },
    ...(hasKo ? [{ key: 'ko', label: 'KO' }] : []),
    { key: 'bin', label: 'MAG' },
    { key: 'type', label: 'Type' },
  ]);
</script>

{#if rows.length}
  <ExpressionMatrix {rows} samples={$topGenes.samples} noun="genes" title="Most Expressed Genes"
    groups={{ label: 'Type', field: 'type', values: TYPES }}
    searchFields={['id', 'gene', 'product', 'ko', 'cog', 'bin', 'contig']}
    searchPlaceholder="locus tag, gene, product, KO, MAG..."
    {columns} idLabel="Locus tag" idWidth={130}
    labelOf={r => `${r.id}${r.gene ? ' ' + r.gene : ''} · ${(r.product || '').slice(0, 40)}`}
    subtitleOf={r => r.product || 'hypothetical protein'}
    annotations={shown => [
      { label: 'Gene', values: shown.map(r => r.gene || '-') },
      ...(hasKo ? [{ label: 'KO', values: shown.map(r => r.ko || '-') }] : []),
      { label: 'MAG', values: shown.map(r => r.bin || '-'), width: 160 },
      { label: 'Product', values: shown.map(r => r.product || '-'), align: 'start', width: 420 },
    ]}
    facts={r => [
      { label: 'Gene', value: r.gene }, { label: 'Type', value: r.type },
      { label: 'Length', value: r.length ? `${fmtInt(r.length)} bp` : null },
      { label: 'Contig', value: r.contig }, { label: 'MAG', value: r.bin, cls: 'text-emerald-400' },
      { label: 'KO', value: r.ko }, { label: 'COG', value: r.cog },
    ]}
    links={r => r.ko ? [{ label: 'KEGG', href: `https://www.kegg.jp/entry/${r.ko}` }] : []}
    note="The genes with the highest mean TPM (Bakta locus tags); the full gene × sample matrix is in expression/<ref>/." />
{:else}
  <div class="flex flex-col items-center justify-center py-24 text-center">
    <p class="text-slate-400 text-lg mb-2">No gene-level counts for this reference</p>
    <p class="text-slate-500 text-sm">Gene counts need a GFF next to the reference (Bakta / Prokka).</p>
  </div>
{/if}
