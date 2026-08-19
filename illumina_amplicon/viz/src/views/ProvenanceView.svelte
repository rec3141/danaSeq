<script>
  import { onMount } from 'svelte';
  import Plotly from 'plotly.js-dist-min';
  import {
    store, listAssays, assayKey, assayHeading, assayMatchNote, sampleInAssays,
  } from '../stores/data.svelte.js';
  import { chrome } from '../stores/theme.svelte.js';

  let { filters = {} } = $props();

  let plotDiv = $state(null);
  let heatDiv = $state(null);
  let selected = $state('__all__');

  const prov = $derived(store.provenance);
  const stages = $derived(prov?.stages ?? []);

  const assayFilterOn = $derived((filters.assays?.size ?? 0) > 0);

  // Sample ids the assay selection admits, or null when it admits everything.
  // Null rather than "all of them" so every consumer below can skip the
  // restriction instead of testing membership once per sample per redraw.
  const allowedIds = $derived.by(() => {
    if (!assayFilterOn) return null;
    const ids = new Set();
    for (const s of store.samples) if (sampleInAssays(s, filters.assays)) ids.add(s.id);
    return ids;
  });

  const trackedIds = $derived(Object.keys(prov?.samples ?? {}));
  const sampleIds = $derived(
    trackedIds.filter(id => !allowedIds || allowedIds.has(id)).sort()
  );

  // A sample the filter no longer admits cannot stay selected — its funnel would
  // keep drawing while the sidebar says that assay is hidden.
  $effect(() => {
    if (selected !== '__all__' && !sampleIds.includes(selected)) selected = '__all__';
  });

  // Counts for whatever is selected: the combined total, or one sample.
  const counts = $derived.by(() => {
    if (!prov) return null;
    if (selected !== '__all__') return prov.samples?.[selected];
    // prov.total covers every sample the run tracked, so under an assay filter
    // the combined column has to be re-summed from the admitted samples.
    if (!allowedIds) return prov.total;
    const sum = {};
    for (const id of sampleIds) {
      const c = prov.samples[id] || {};
      for (const s of stages) sum[s.id] = (sum[s.id] ?? 0) + (c[s.id] ?? 0);
    }
    return sum;
  });

  // Retention relative to the raw FASTQ — the number that matters, since a
  // sample can look "shallow" when in fact primer removal discarded it.
  const rows = $derived.by(() => {
    if (!counts || !stages.length) return [];
    const raw = counts[stages[0].id] || 0;
    return stages.map(s => {
      const n = counts[s.id] ?? 0;
      return { id: s.id, label: s.label, n, pct: raw ? (100 * n) / raw : 0 };
    });
  });

  // Samples whose reads mostly died at primer removal — worth surfacing, it
  // usually means that run used a different primer pair than the one applied.
  const lostSamples = $derived.by(() => {
    if (!prov?.samples || stages.length < 2) return [];
    const [rawKey, primerKey] = [stages[0].id, stages[1].id];
    return sampleIds
      .map(id => ({ id, raw: prov.samples[id][rawKey] ?? 0, kept: prov.samples[id][primerKey] ?? 0 }))
      .filter(s => s.raw >= 1000 && s.kept / s.raw < 0.5)
      .map(s => ({ ...s, pct: (100 * s.kept) / s.raw }))
      .sort((a, b) => a.pct - b.pct);
  });

  // ── Primers ───────────────────────────────────────────────────────────────

  /**
   * The primer sequences cutadapt worked with, per assay.
   *
   * The sequences come from the sample records rather than the assay key: a
   * placed assay is identified by where it sits on the reference and keeps the
   * primers out of its key, precisely because samples of one assay differ by a
   * few bases. So one assay can hold more than one pair, and each is listed
   * with the number of samples that carried it.
   */
  const primerSets = $derived.by(() => {
    const listed = listAssays().filter(a => !assayFilterOn || filters.assays.has(a.key));
    const byKey = new Map(listed.map(a => [a.key, { ...a, pairs: new Map() }]));
    for (const s of store.samples) {
      const rec = byKey.get(assayKey(s));
      if (!rec) continue;
      const fwd = (s.assay_primer_fwd || '').trim();
      const rev = (s.assay_primer_rev || '').trim();
      const source = (s.assay_source || '').trim();
      const id = `${fwd}|${rev}|${source}`;
      const pair = rec.pairs.get(id) || { fwd, rev, source, n: 0 };
      pair.n++;
      rec.pairs.set(id, pair);
    }
    return [...byKey.values()].map(a => ({
      ...a,
      pairs: [...a.pairs.values()].sort((x, y) => y.n - x.n),
    }));
  });

  // ── Parameters ────────────────────────────────────────────────────────────

  const params = $derived(store.manifest?.parameters ?? null);

  function fmt(v) {
    if (v === true) return 'on';
    if (v === false) return 'off';
    if (v === null || v === undefined) return '—';
    return typeof v === 'number' ? v.toLocaleString() : String(v);
  }

  /**
   * The resolved parameters, as sentences in the order a read meets them.
   *
   * Each line declares the keys it speaks for, and whatever is left over is
   * listed verbatim under "Other" — a parameter the pipeline gains later has to
   * appear somewhere rather than be quietly dropped by a hard-coded list.
   */
  const paramLines = $derived.by(() => {
    const p = params;
    if (!p) return [];
    const has = k => p[k] !== undefined && p[k] !== null;
    const used = new Set();
    const lines = [];
    const line = (label, text, ...keys) => {
      keys.forEach(k => used.add(k));
      if (text) lines.push({ label, text });
    };

    let primer = '';
    if (p.skip_primer_removal) {
      primer = 'skipped — the reads keep whatever they arrived with';
    } else if (has('primer_error_rate') || has('primer_detect_reads')) {
      const bits = ['cutadapt'];
      if (has('primer_error_rate')) bits.push(`≤ ${(p.primer_error_rate * 100).toFixed(0)}% mismatch`);
      if (has('primer_detect_reads')) bits.push(`${fmt(p.primer_detect_reads)} reads read per file to detect the primer`);
      primer = bits.join(' · ');
    }
    line('Primer removal', primer, 'skip_primer_removal', 'primer_error_rate', 'primer_detect_reads');

    let trunc = '';
    if (p.auto_trim) {
      trunc = `automatic${has('auto_trim_min_quality') ? `, down to Q${fmt(p.auto_trim_min_quality)}` : ''}`
            + `${has('trunc_policy') ? ` · policy ${p.trunc_policy}` : ''}`;
    } else if (p.truncLen_fwd || p.truncLen_rev) {
      trunc = `fixed · forward ${fmt(p.truncLen_fwd)} bp, reverse ${fmt(p.truncLen_rev)} bp`;
    } else if (has('auto_trim')) {
      trunc = 'off — reads keep their full length';
    }
    line('Truncation', trunc, 'auto_trim', 'auto_trim_min_quality', 'trunc_policy',
         'truncLen_fwd', 'truncLen_rev');

    const qc = [];
    if (has('maxEE')) qc.push(`maxEE ${fmt(p.maxEE)}`);
    if (has('truncQ')) qc.push(`truncQ ${fmt(p.truncQ)}`);
    if (has('maxN')) qc.push(`maxN ${fmt(p.maxN)}`);
    line('Quality filter', qc.join(' · '), 'maxEE', 'truncQ', 'maxN');

    line('Merging', has('min_overlap') ? `≥ ${fmt(p.min_overlap)} bp overlap` : '', 'min_overlap');

    const asv = [];
    if (has('min_seq_length')) asv.push(`≥ ${fmt(p.min_seq_length)} bp`);
    if (has('min_seqs')) asv.push(`≥ ${fmt(p.min_seqs)} reads`);
    if (has('min_samples')) asv.push(`in ≥ ${fmt(p.min_samples)} sample${p.min_samples === 1 ? '' : 's'}`);
    line('ASVs kept', asv.join(' · '), 'min_seq_length', 'min_seqs', 'min_samples');

    line('Samples kept', has('min_reads') ? `≥ ${fmt(p.min_reads)} reads` : '', 'min_reads');
    line('Error-model group', has('min_group_samples')
      ? `≥ ${fmt(p.min_group_samples)} samples` : '', 'min_group_samples');

    const other = Object.entries(p).filter(([k]) => !used.has(k));
    if (other.length) {
      lines.push({ label: 'Other', text: other.map(([k, v]) => `${k} ${fmt(v)}`).join(' · ') });
    }
    return lines;
  });

  // Flags as they were typed, so options that never reach the parameters block
  // (--min_prevalence, --metadata, an explicit primer FASTA) are still readable.
  const cliFlags = $derived.by(() => {
    const cmd = store.manifest?.command_line;
    if (!cmd) return [];
    const toks = String(cmd).match(/'[^']*'|"[^"]*"|\S+/g) ?? [];
    const out = [];
    for (let i = 0; i < toks.length; i++) {
      if (!toks[i].startsWith('-')) continue;
      const next = toks[i + 1];
      const value = next && !next.startsWith('-') ? next.replace(/^['"]|['"]$/g, '') : '';
      out.push({ flag: toks[i], value });
      if (value) i++;
    }
    return out;
  });

  const cliMap = $derived(new Map(cliFlags.map(f => [f.flag.replace(/^-+/, ''), f.value])));

  const requestLines = $derived.by(() => {
    if (!cliMap.size) return [];
    const val = k => cliMap.get(k) ?? '';
    const lines = [];
    if (val('input')) lines.push({ label: 'Reads', text: val('input'), mono: true });
    if (val('primers_fwd') || val('primers_rev')) {
      lines.push({ label: 'Primers given', mono: true,
                   text: [val('primers_fwd'), val('primers_rev')].filter(Boolean).join(' · ') });
    } else if (!params?.skip_primer_removal) {
      lines.push({ label: 'Primers given', text: 'none — detected from the reads' });
    }
    if (val('metadata')) {
      lines.push({ label: 'Metadata', mono: true,
                   text: `${val('metadata')}${val('sample_id_column') ? ` · id column ${val('sample_id_column')}` : ''}` });
    }
    if (val('min_prevalence')) {
      lines.push({ label: 'Network', text: `ASVs in ≥ ${val('min_prevalence')} samples` });
    }
    if (cliMap.has('run_phylogeny')) lines.push({ label: 'Phylogeny', text: 'alignment and tree built' });
    return lines;
  });

  // ── Reference databases ───────────────────────────────────────────────────

  /**
   * The release a reference file states, or null when its name states none.
   *
   * A reference carries no version anywhere but its filename, and only some
   * naming conventions put one there: SILVA_138.2_SSURef_NR99.fasta,
   * pr2_version_5.1.1_SSU_dada2.fasta, silva_nr_v132_train_set.fa. A name that
   * matches none of these gets no version rather than a number lifted from
   * whichever digits it happens to contain.
   */
  function versionFromFilename(file) {
    const stem = String(file || '')
      .replace(/\.gz$/i, '')
      .replace(/\.(fasta|fa|fna|fas|faa|tsv|txt)$/i, '');
    const tagged = stem.match(/(?:^|[_.\-])v(?:ersion)?[_.\-]?(\d+(?:\.\d+)*)(?=[_.\-]|$)/i);
    if (tagged) return tagged[1];
    const dotted = stem.match(/(?:^|[_.\-])(\d+\.\d+(?:\.\d+)*)(?=[_.\-]|$)/);
    if (dotted) return dotted[1];
    const bare = stem.match(/(?:^|[_.\-])(\d{2,4})(?=[_.\-]|$)/);
    if (bare) return bare[1];
    return null;
  }

  /**
   * `name:path:Rank1,Rank2;name:path:...` as the run was given it.
   *
   * The path is delimited on both sides rather than split on every colon, since
   * a path may hold one of its own; a trailing field that looks like a path is
   * treated as one, which is what an entry naming no ranks looks like.
   */
  function parseReferences(spec) {
    if (!spec) return [];
    return String(spec).split(';').map(e => e.trim()).filter(Boolean).map((entry, i) => {
      const first = entry.indexOf(':');
      const last = entry.lastIndexOf(':');
      const name = first > 0 ? entry.slice(0, first) : entry;
      const rest = first > 0 ? entry.slice(first + 1) : '';
      const tail = last > first ? entry.slice(last + 1) : '';
      const hasRanks = tail !== '' && !tail.includes('/');
      const path = hasRanks ? entry.slice(first + 1, last) : rest;
      const ranks = hasRanks ? tail.split(',').map(r => r.trim()).filter(Boolean) : [];
      const file = path.split('/').pop();
      return { name, path, file, ranks, version: versionFromFilename(file), primary: i === 0 };
    });
  }

  const references = $derived(parseReferences(store.manifest?.reference_databases));

  // ── Plots ─────────────────────────────────────────────────────────────────

  function draw() {
    if (!plotDiv || !rows.length) return;
    // Reading the palette here is what redraws the plot when the theme flips.
    const c = chrome();
    Plotly.react(
      plotDiv,
      [{
        type: 'bar',
        x: rows.map(r => r.label),
        y: rows.map(r => r.n),
        marker: { color: rows.map((_, i) => ['#38bdf8', '#818cf8', '#34d399'][i % 3]) },
        text: rows.map(r => `${r.n.toLocaleString()}<br>${r.pct.toFixed(1)}%`),
        textposition: 'auto',
        hovertemplate: '%{x}<br>%{y:,} reads<extra></extra>',
      }],
      {
        margin: { l: 70, r: 20, t: 30, b: 60 },
        plot_bgcolor: c.plotBg,
        paper_bgcolor: c.plotBg,
        font: { color: c.axis, size: 12 },
        xaxis: { gridcolor: c.grid },
        yaxis: { title: 'reads', gridcolor: c.grid, rangemode: 'tozero' },
        showlegend: false,
      },
      { responsive: true, displaylogo: false },
    );
  }

  // All-samples heatmap: rows = samples, columns = stages, colour = log10(reads).
  // Gives a one-glance view of the whole dataset's funnel — a sample that
  // collapses at a stage (e.g. an outlier stripped by the prevalence filter)
  // shows up as a dark cell in that column. Rows are ordered by overall
  // retention so the worst-retained samples cluster together.
  function drawHeat() {
    if (!heatDiv || selected !== '__all__' || !stages.length || !sampleIds.length) return;
    // Reading the palette here is what redraws the plot when the theme flips.
    const c = chrome();
    const rawId = stages[0].id, finalId = stages[stages.length - 1].id;
    const retention = (id) => {
      const r = prov.samples[id][rawId] || 0;
      return r ? (prov.samples[id][finalId] || 0) / r : 0;
    };
    const order = [...sampleIds].sort((a, b) => retention(a) - retention(b));
    const x = stages.map(s => s.label);
    const z = order.map(id => stages.map(s => {
      const n = prov.samples[id][s.id] ?? 0;
      return n > 0 ? Math.log10(n) : null;   // null → blank cell for zero reads
    }));
    const raw = order.map(id => stages.map(s => prov.samples[id][s.id] ?? 0));
    Plotly.react(
      heatDiv,
      [{
        type: 'heatmap',
        x, y: order, z,
        customdata: raw,
        colorscale: 'Viridis',
        hoverongaps: false,
        hovertemplate: '%{y}<br>%{x}<br>%{customdata:,} reads<extra></extra>',
        colorbar: { title: { text: 'log₁₀(reads)', side: 'right' }, thickness: 12,
                    tickfont: { color: c.axis, size: 10 } },
      }],
      {
        margin: { l: 110, r: 10, t: 10, b: 70 },
        plot_bgcolor: c.plotBg,
        paper_bgcolor: c.plotBg,
        font: { color: c.axis, size: 11 },
        xaxis: { tickangle: -30, automargin: true },
        yaxis: { automargin: true, autorange: 'reversed' },
      },
      { responsive: true, displaylogo: false },
    );
  }

  onMount(() => { draw(); drawHeat(); });
  $effect(() => { rows; draw(); });
  $effect(() => { selected; stages; sampleIds; drawHeat(); });
</script>

<div class="flex h-full flex-col gap-3 overflow-y-auto p-3">
  <!-- What this run is and what produced it, above the read funnel. Both come
       from files the pipeline and the portal each write for their own reasons;
       neither is much use without the other when you are trying to work out
       why two runs of one accession disagree. -->
  {#if store.runInfo || store.manifest}
    <div class="rounded-lg border border-line bg-surface/40 p-3">
      {#if store.runInfo?.title}
        <h2 class="text-sm font-semibold text-strong">{store.runInfo.title}</h2>
      {/if}
      <dl class="mt-2 grid grid-cols-[auto_1fr] gap-x-4 gap-y-1 text-xs">
        {#if store.runInfo?.bioproject}
          <dt class="text-faint">BioProject</dt>
          <dd><a class="text-link hover:text-link2"
                 href="https://www.ncbi.nlm.nih.gov/bioproject/{store.runInfo.bioproject}"
                 target="_blank" rel="noopener">{store.runInfo.bioproject}</a></dd>
        {/if}
        {#if store.runInfo?.registered}
          <dt class="text-faint">Released</dt>
          <dd class="text-fg2">{store.runInfo.registered}
            {#if store.runInfo.updated && store.runInfo.updated !== store.runInfo.registered}
              <span class="text-faint">· updated {store.runInfo.updated}</span>
            {/if}
          </dd>
        {/if}
        {#if store.runInfo?.slug}
          <dt class="text-faint">Run</dt>
          <dd class="font-mono text-fg2">
            {#if store.runInfo.portal_url}
              <a class="text-link hover:text-link2" href={store.runInfo.portal_url}
                 target="_blank" rel="noopener">{store.runInfo.slug}</a>
            {:else}{store.runInfo.slug}{/if}
          </dd>
        {/if}
        {#if store.runInfo?.pipeline}
          <dt class="text-faint">Pipeline</dt><dd class="text-fg2">{store.runInfo.pipeline}</dd>
        {/if}
        {#if store.runInfo?.cluster}
          <dt class="text-faint">Cluster</dt><dd class="text-fg2">{store.runInfo.cluster}</dd>
        {/if}
        {#if store.manifest?.pipeline}
          <dt class="text-faint">Version</dt><dd class="text-fg2">{store.manifest.pipeline}</dd>
        {/if}
        {#if store.manifest?.commit_id}
          <dt class="text-faint">Build</dt>
          <dd class="font-mono text-fg2">
            <a class="text-link hover:text-link2"
               href="https://github.com/rec3141/danaSeq/commit/{store.manifest.commit_id}"
               target="_blank" rel="noopener">{store.manifest.commit_id.slice(0, 7)}</a>
            {#if store.manifest.container_built}
              <span class="text-faint">· built {store.manifest.container_built.slice(0, 10)}</span>
            {/if}
          </dd>
        {/if}
        {#if store.manifest?.nextflow_version || store.manifest?.denoise}
          <dt class="text-faint">Engine</dt>
          <dd class="text-fg2">
            {#if store.manifest.nextflow_version}Nextflow {store.manifest.nextflow_version}{/if}
            {#if store.manifest.denoise}<span class="text-faint">· {store.manifest.denoise}</span>{/if}
          </dd>
        {/if}
        {#if store.manifest?.completed}
          <dt class="text-faint">Finished</dt>
          <dd class="text-fg2">{store.manifest.completed.slice(0, 10)}
            {#if store.manifest.duration}<span class="text-faint">· took {store.manifest.duration}</span>{/if}
            {#if store.manifest.success === false}
              <span class="text-warn">· the run reported failures, so these outputs may be partial</span>
            {/if}
          </dd>
        {/if}
        <dt class="text-faint">Samples</dt>
        <dd class="text-fg2">
          {#if allowedIds}
            {allowedIds.size.toLocaleString()} of {store.samples.length.toLocaleString()}
            · {(store.assayAsvIds?.size ?? store.asvs.length).toLocaleString()} ASVs
            <span class="text-faint">· assay filter on</span>
          {:else}
            {store.samples.length.toLocaleString()} · {store.asvs.length.toLocaleString()} ASVs
          {/if}
        </dd>
      </dl>

      {#if store.manifest?.tool_versions}
        <details class="mt-2">
          <summary class="cursor-pointer text-xs text-muted hover:text-fg">Tool versions</summary>
          <pre class="mt-2 max-h-64 overflow-auto rounded bg-sunken/60 p-2 text-[11px] leading-relaxed text-muted">{JSON.stringify(store.manifest.tool_versions, null, 2)}</pre>
        </details>
      {/if}

      {#if references.length}
        <details class="mt-2">
          <summary class="cursor-pointer text-xs text-muted hover:text-fg">
            Database versions ({references.length})
          </summary>
          <div class="mt-2 space-y-2">
            {#if references.length > 1}
              <p class="text-[11px] text-faint">
                The first is the primary: the ASV-level labels and the renormalisation come from it,
                the rest are alternatives the taxonomy views can be switched to.
              </p>
            {/if}
            {#each references as r}
              <div class="rounded bg-sunken/60 p-2 text-[11px]">
                <p>
                  <span class="font-semibold text-fg2">{r.name}</span>
                  {#if r.version}
                    <span class="text-fg2">{r.version}</span>
                  {:else}
                    <span class="text-caution">version undeclared</span>
                  {/if}
                  {#if r.primary && references.length > 1}
                    <span class="text-faint">· primary</span>
                  {/if}
                </p>
                {#if r.ranks.length}
                  <p class="mt-1 text-muted">
                    {r.ranks.length} ranks: <span class="text-fg2">{r.ranks.join(' › ')}</span>
                  </p>
                {:else}
                  <p class="mt-1 text-caution">no ranks declared for this reference</p>
                {/if}
                <p class="mt-1 break-all font-mono text-faint">{r.path}</p>
              </div>
            {/each}
          </div>
        </details>
      {/if}
    </div>
  {/if}

  <div class="grid gap-3 lg:grid-cols-2">
    <!-- What was asked for. The parameters block is what the run resolved; the
         command line is what a person typed, and carries the options that never
         reach the parameters block. -->
    {#if paramLines.length || requestLines.length}
      <div class="rounded-lg border border-line bg-surface/40 p-3">
        <h3 class="text-xs font-semibold uppercase tracking-wider text-muted">How this run was configured</h3>
        <dl class="mt-2 grid grid-cols-[auto_1fr] gap-x-4 gap-y-1 text-xs">
          {#each paramLines as l}
            <dt class="whitespace-nowrap text-faint">{l.label}</dt>
            <dd class="text-fg2">{l.text}</dd>
          {/each}
          {#each requestLines as l}
            <dt class="whitespace-nowrap text-faint">{l.label}</dt>
            <dd class="break-all text-fg2 {l.mono ? 'font-mono' : ''}">{l.text}</dd>
          {/each}
        </dl>

        {#if params}
          <details class="mt-2">
            <summary class="cursor-pointer text-xs text-muted hover:text-fg">All parameters</summary>
            <pre class="mt-2 max-h-64 overflow-auto rounded bg-sunken/60 p-2 text-[11px] leading-relaxed text-muted">{JSON.stringify(params, null, 2)}</pre>
          </details>
        {/if}
        {#if store.manifest?.command_line}
          <details class="mt-2">
            <summary class="cursor-pointer text-xs text-muted hover:text-fg">Command line</summary>
            <pre class="mt-2 max-h-64 overflow-auto whitespace-pre-wrap break-all rounded bg-sunken/60 p-2 text-[11px] leading-relaxed text-muted">{store.manifest.command_line}</pre>
          </details>
        {/if}
      </div>
    {/if}

    <!-- The sequences themselves, in full: a truncated primer is not a primer,
         and this is the field people come to this tab to copy. -->
    <div class="rounded-lg border border-line bg-surface/40 p-3">
      <h3 class="text-xs font-semibold uppercase tracking-wider text-muted">Primer sequences</h3>
      {#if params?.skip_primer_removal}
        <p class="mt-1 text-[11px] text-caution">
          Primer removal was skipped for this run — nothing was trimmed off the reads.
        </p>
      {/if}
      {#if !primerSets.length}
        <p class="mt-2 text-xs text-faint">
          {assayFilterOn
            ? 'No assay in the current selection records a primer.'
            : 'This run records no assay per sample, so there are no primer sequences to show.'}
        </p>
      {:else}
        <div class="mt-2 space-y-2">
          {#each primerSets as a}
            {@const match = assayMatchNote(a.meanMatch)}
            <div class="border-t border-line pt-2 first:border-t-0 first:pt-0">
              <p class="text-xs text-fg2">
                {assayHeading(a.gene, a.region, a.fwd, a.rev, a.place)}
                {#if a.span}<span class="text-faint">· {a.span}</span>{/if}
                <span class="text-faint">· {a.n} sample{a.n === 1 ? '' : 's'}</span>
              </p>
              {#if match}
                <p class="text-[10px] text-faint" title={match.title}>{match.text}</p>
              {/if}
              {#each a.pairs as p}
                {#if p.fwd || p.rev}
                  <dl class="mt-1 grid grid-cols-[4rem_1fr] gap-x-2 text-[11px]">
                    <dt class="text-faint">forward</dt>
                    <dd class="break-all font-mono text-fg2">{p.fwd || 'none recorded'}</dd>
                    <dt class="text-faint">reverse</dt>
                    <dd class="break-all font-mono text-fg2">{p.rev || 'none recorded'}</dd>
                  </dl>
                  <p class="text-[10px] text-faint">
                    {p.n} sample{p.n === 1 ? '' : 's'}{p.source ? ` · via ${p.source}` : ''}
                  </p>
                {:else}
                  <p class="mt-1 text-[11px] text-faint">
                    {p.n} sample{p.n === 1 ? '' : 's'} carry no primer sequence — their reads
                    arrived with the primers already off{p.source ? ` (via ${p.source})` : ''}.
                  </p>
                {/if}
              {/each}
            </div>
          {/each}
        </div>
      {/if}
    </div>
  </div>

  {#if !prov}
    <div class="flex flex-1 items-center justify-center text-sm text-faint">
      No read-tracking data — <code class="mx-1">data/provenance.json</code> was not produced for this run.
    </div>
  {:else if !sampleIds.length}
    <div class="flex flex-1 items-center justify-center text-center text-sm text-faint">
      No sample in the selected assay{filters.assays?.size === 1 ? '' : 's'} was tracked through the read funnel.
    </div>
  {:else}
    <div class="flex items-center gap-3">
      <label for="prov-sample" class="text-xs text-muted">Sample</label>
      <select
        id="prov-sample"
        bind:value={selected}
        class="rounded border border-line bg-surface px-2 py-1 text-xs text-fg"
      >
        <option value="__all__">All samples combined ({sampleIds.length})</option>
        {#each sampleIds as id}
          <option value={id}>{id}</option>
        {/each}
      </select>
      {#if rows.length}
        <span class="text-xs text-faint">
          {rows[0].n.toLocaleString()} raw &rarr; {rows[rows.length - 1].n.toLocaleString()} retained
          ({rows[rows.length - 1].pct.toFixed(1)}%)
        </span>
      {/if}
      {#if allowedIds}
        <span class="text-xs text-faint">
          {sampleIds.length} of {trackedIds.length} tracked samples pass the assay filter
        </span>
      {/if}
    </div>

    <div bind:this={plotDiv} class="h-72 shrink-0"></div>

    {#if selected === '__all__' && sampleIds.length}
      <div class="shrink-0">
        <p class="mb-1 text-xs text-muted">
          All samples &middot; log<sub>10</sub>(reads) at each stage
          <span class="text-faint">— rows sorted by retention (lowest first)</span>
        </p>
        <div
          bind:this={heatDiv}
          style="height: {Math.min(700, 70 + sampleIds.length * 16)}px"
        ></div>
      </div>
    {/if}

    <div class="max-h-40 overflow-y-auto">
      <table class="w-full text-xs">
        <thead class="text-muted">
          <tr><th class="p-1 text-left">Step</th><th class="p-1 text-right">Reads</th><th class="p-1 text-right">% of raw</th></tr>
        </thead>
        <tbody class="text-fg2">
          {#each rows as r}
            <tr class="border-t border-line">
              <td class="p-1">{r.label}</td>
              <td class="p-1 text-right">{r.n.toLocaleString()}</td>
              <td class="p-1 text-right" class:text-warn={r.pct < 50}>{r.pct.toFixed(1)}%</td>
            </tr>
          {/each}
        </tbody>
      </table>

      {#if selected === '__all__' && lostSamples.length}
        <p class="mt-2 text-xs text-caution">
          {lostSamples.length} sample{lostSamples.length === 1 ? '' : 's'} lost over half their reads at
          primer removal — usually a different primer pair than the one applied:
        </p>
        <table class="w-full text-xs">
          <tbody class="text-fg2">
            {#each lostSamples.slice(0, 20) as s}
              <tr class="border-t border-line">
                <td class="p-1">{s.id}</td>
                <td class="p-1 text-right">{s.raw.toLocaleString()} &rarr; {s.kept.toLocaleString()}</td>
                <td class="p-1 text-right text-warn">{s.pct.toFixed(1)}%</td>
              </tr>
            {/each}
          </tbody>
        </table>
      {/if}
    </div>
  {/if}
</div>
