<script>
  import MatrixPanel from '../components/MatrixPanel.svelte';
  import { funcs } from '../stores/data.js';

  const NAMES = { ko: 'KEGG orthologs', cog: 'COG categories', cog_id: 'COGs', cazy: 'CAZy families' };
  const COG = {
    J: 'Translation, ribosome', A: 'RNA processing', K: 'Transcription', L: 'Replication & repair',
    B: 'Chromatin', D: 'Cell cycle', Y: 'Nuclear structure', V: 'Defense', T: 'Signal transduction',
    M: 'Cell wall/membrane', N: 'Motility', Z: 'Cytoskeleton', W: 'Extracellular structures',
    U: 'Secretion & trafficking', O: 'Protein turnover, chaperones', X: 'Mobilome',
    C: 'Energy production', G: 'Carbohydrate metabolism', E: 'Amino acid metabolism',
    F: 'Nucleotide metabolism', H: 'Coenzyme metabolism', I: 'Lipid metabolism',
    P: 'Inorganic ion transport', Q: 'Secondary metabolites', R: 'General function', S: 'Function unknown',
  };

  let kinds = $derived(Object.keys($funcs?.kinds ?? {}));
  let kind = $state(null);
  $effect(() => { if (!kind || !kinds.includes(kind)) kind = kinds[0] ?? null; });

  let rows = $derived((kind ? $funcs.kinds[kind] : []).map(f => ({
    ...f, name: kind === 'cog' ? COG[f.id] ?? null : null,
  })));
  let columns = $derived(kind === 'cog' ? [{ key: 'name', label: 'Category' }] : []);
  const keggLink = r => `https://www.kegg.jp/entry/${r.id}`;
</script>

{#if kinds.length}
  <div class="flex gap-2 mb-4">
    {#each kinds as k}
      <button class="px-3 py-1 rounded text-xs border {kind === k ? 'border-cyan-500 text-cyan-300' : 'border-slate-700 text-slate-400 hover:text-slate-200'}"
        onclick={() => (kind = k)}>{NAMES[k] ?? k}</button>
    {/each}
  </div>
  {#key kind}
    <MatrixPanel {rows} samples={$funcs.samples} idLabel={NAMES[kind] ?? kind}
      columns={kind === 'ko' ? [{ key: 'id', label: 'KEGG', link: keggLink, fmt: () => 'entry ↗' }] : columns}
      labelOf={r => r.name ? `${r.id} · ${r.name}` : r.id}
      title="{NAMES[kind] ?? kind} by expression"
      subtitle="Sum of gene TPM per function; each gene counts toward its first assigned function." />
  {/key}
{:else}
  <div class="bg-slate-900 border border-slate-800 rounded-lg p-6 text-sm text-slate-400">
    No functional annotation for this reference. Run mag_analysis with its metabolism module (KofamScan / eggNOG) to get KO, COG and CAZy tables here.
  </div>
{/if}
