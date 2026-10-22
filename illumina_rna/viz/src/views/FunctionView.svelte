<script>
  import ExpressionMatrix from '../components/ExpressionMatrix.svelte';
  import { funcs } from '../stores/data.js';

  const KINDS = { ko: 'KEGG orthologs', cog: 'COG categories', cog_id: 'COGs', cazy: 'CAZy families' };
  const COG = {
    J: 'Translation, ribosomal structure', A: 'RNA processing', K: 'Transcription', L: 'Replication & repair',
    B: 'Chromatin structure', D: 'Cell cycle control', Y: 'Nuclear structure', V: 'Defense mechanisms',
    T: 'Signal transduction', M: 'Cell wall/membrane/envelope', N: 'Cell motility', Z: 'Cytoskeleton',
    W: 'Extracellular structures', U: 'Secretion & vesicular transport', O: 'Protein turnover, chaperones',
    X: 'Mobilome', C: 'Energy production & conversion', G: 'Carbohydrate metabolism', E: 'Amino acid metabolism',
    F: 'Nucleotide metabolism', H: 'Coenzyme metabolism', I: 'Lipid metabolism', P: 'Inorganic ion transport',
    Q: 'Secondary metabolites', R: 'General function prediction', S: 'Function unknown',
  };
  const CAZY_CLASS = { GH: 'Glycoside hydrolase', GT: 'Glycosyltransferase', PL: 'Polysaccharide lyase',
    CE: 'Carbohydrate esterase', AA: 'Auxiliary activity', CBM: 'Carbohydrate-binding module' };

  let kinds = $derived(Object.keys($funcs?.kinds ?? {}));
  let kind = $state(null);
  $effect(() => { if (!kind || !kinds.includes(kind)) kind = kinds[0] ?? null; });

  let rows = $derived((kind ? $funcs.kinds[kind] : []).map(f => ({
    ...f,
    name: kind === 'cog' ? COG[f.id] ?? null : f.name ?? null,
    cls: kind === 'cazy' ? CAZY_CLASS[f.id.match(/^[A-Z]+/)?.[0]] ?? null : null,
  })));

  let hasNames = $derived(kind === 'ko' && rows.some(r => r.name));
  let total = k => $funcs?.n_total?.[k] ?? $funcs.kinds[k].length;

  let columns = $derived([
    { key: 'id', label: KINDS[kind] ?? kind },
    ...(kind === 'cog' ? [{ key: 'name', label: 'Category' }] : []),
    ...(kind === 'ko' && hasNames ? [{ key: 'name', label: 'Definition', render: v => v ? (v.length > 70 ? v.slice(0, 68) + '…' : v) : '-' }] : []),
    ...(kind === 'cazy' ? [{ key: 'cls', label: 'Class' }] : []),
  ]);
  const link = r => kind === 'ko' ? [{ label: 'KEGG', href: `https://www.kegg.jp/entry/${r.id}` }]
    : kind === 'cazy' ? [{ label: 'CAZy', href: `http://www.cazy.org/${r.id}.html` }]
    : kind === 'cog_id' ? [{ label: 'COG', href: `https://www.ncbi.nlm.nih.gov/research/cog/cog/${r.id}/` }] : [];
</script>

{#if kinds.length}
  <div class="flex items-center gap-3 text-xs mb-2">
    <span class="text-slate-400">Annotation:</span>
    {#each kinds as k}
      <button
        class="px-3 py-1 rounded-md border transition-colors
          {kind === k ? 'border-cyan-400 bg-cyan-400/10 text-cyan-400' : 'border-slate-600 text-slate-400 hover:border-slate-500'}"
        onclick={() => kind = k}
      >{KINDS[k] ?? k} <span class="text-slate-500">({total(k).toLocaleString()})</span></button>
    {/each}
    {#if kind && total(kind) > $funcs.kinds[kind].length}
      <span class="text-slate-500">top {$funcs.kinds[kind].length.toLocaleString()} of {total(kind).toLocaleString()} by mean TPM; all in expression/&lt;ref&gt;/</span>
    {/if}
  </div>
  {#key kind}
    <ExpressionMatrix {rows} samples={$funcs.samples} noun={KINDS[kind] ?? kind} title={KINDS[kind] ?? kind}
      searchFields={['id', 'name', 'cls']} searchPlaceholder={kind === 'ko' ? (hasNames ? 'K number or definition...' : 'K number...') : kind === 'cog' ? 'letter or category...' : 'family or class...'}
      {columns} idLabel={kind === 'ko' ? 'KO' : kind === 'cazy' ? 'Family' : 'COG'} idWidth={70}
      labelOf={r => r.name ? `${r.id} · ${r.name}` : r.id}
      subtitleOf={r => r.name ?? r.cls ?? KINDS[kind]}
      annotations={shown => kind === 'cog' ? [{ label: 'Category', values: shown.map(r => r.name ?? '-'), align: 'start', width: 260 }]
        : kind === 'ko' && hasNames ? [{ label: 'Definition', values: shown.map(r => r.name ?? '-'), align: 'start', width: 520 }]
        : kind === 'cazy' ? [{ label: 'Class', values: shown.map(r => r.cls ?? '-'), align: 'start', width: 220 }] : null}
      facts={r => [{ label: kind === 'ko' ? 'Definition' : 'Category', value: r.name }, { label: 'Class', value: r.cls }]}
      links={link}
      note="Sum of gene TPM per function; each gene counts toward its first assigned function. Log-scaled colour." />
  {/key}
{:else}
  <div class="flex flex-col items-center justify-center py-24 text-center">
    <p class="text-slate-400 text-lg mb-2">No functional annotation for this reference</p>
    <p class="text-slate-500 text-sm">Run mag_analysis with its metabolism module (KofamScan / eggNOG / dbCAN) to get KO, COG and CAZy tables here.</p>
  </div>
{/if}
