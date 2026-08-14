#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// dānaSeq Illumina Amplicon Pipeline - Nextflow DSL2
// ============================================================================
//
// Amplicon sequencing analysis pipeline: demultiplexing, primer removal,
// denoising, chimera removal, taxonomy assignment, phylogeny,
// normalization, clustering, network analysis, and visualization.
//
// Usage:
//   nextflow run main.nf --input /path/to/reads -resume
//
// ============================================================================

// ============================================================================
// Help message
// ============================================================================

def helpMessage() {
    log.info """
    =========================================
     dānaSeq Illumina Amplicon Pipeline
     https://github.com/rec3141/danaSeq
    =========================================

    Usage:
      nextflow run main.nf --input /path/to/reads [options] -resume

    Required:
      --input DIR          Input directory containing demultiplexed paired-end
                           *.fastq.gz files

    Primer removal:
      --primer_auto        Auto-select primer pair by filename prefix [default: true]
      --primers_fwd PATH   Forward primer FASTA (overrides auto-selection)
      --primers_rev PATH   Reverse primer FASTA (overrides auto-selection)
      --primer_error_rate  Cutadapt error rate [default: 0.12]
      --primer_detect_reads N  Reads sampled for primer detection [default: 500]
      --min_group_samples N    Min samples before a run/assay group is pooled [default: 3]

    Demultiplexing (optional):
      --run_demultiplex    Enable demultiplexing with Mr_Demuxy [default: false]
      --forward_bcs PATH   Forward barcode file
      --reverse_bcs PATH   Reverse barcode file

    Denoising parameters:
      --maxEE N            Max expected errors [default: 2]
      --truncQ N           Truncation quality [default: 11]
      --maxN N             Max Ns [default: 0]
      --truncLen_fwd N     Truncate forward reads at position N [default: 0 = off]
      --truncLen_rev N     Truncate reverse reads at position N [default: 0 = off]
      --min_overlap N      Min overlap for pair merging [default: 10]

    QC filtering:
      --min_seq_length N   Min sequence length [default: 50]
      --min_reads N        Min reads per sample [default: 100]
      --min_samples N      Min samples per ESV [default: 2]
      --min_seqs N         Min total reads per ESV [default: 3]

    Output:
      --outdir DIR         Output directory [default: results]
      --store_dir DIR      Persistent cache directory [default: off]

    Taxonomy:
      --ref_databases STR  Reference databases (semicolon-separated):
                           "name:path:Level1,Level2,..." [default: none]
      --run_phylogeny      Build phylogenetic tree [default: false]

    Metadata:
      --metadata PATH      MIMARKS-compliant TSV/CSV metadata file [default: none]
      --sample_id_column S Column name matching sample IDs [default: sample_name]

    Network:
      --min_prevalence N   Min samples per ASV for SparCC [default: 5]

    Resources:
      --threads N          General thread count [default: 8]
      --denoise_cpus N     CPUs for denoising processes [default: 8]
      --denoise_memory S   Memory for denoising processes [default: '16 GB']
    """.stripIndent()
}

// ============================================================================
// Parameter validation
// ============================================================================

def validateParams() {
    if (!params.input && !params.samplesheet) {
        log.error "ERROR: --input or --samplesheet is required. Run with --help for usage."
        System.exit(1)
    }
    if (params.run_demultiplex && (!params.forward_bcs || !params.reverse_bcs)) {
        log.error "ERROR: --forward_bcs and --reverse_bcs are required when --run_demultiplex is enabled."
        System.exit(1)
    }
}

// ============================================================================
// Module imports
// ============================================================================

// Stage A: Preprocessing and denoising
include { DEMULTIPLEX }       from './modules/demultiplex'
include { DETECT_PRIMERS }      from './modules/primers'
include { RESOLVE_PRIMER_SET }  from './modules/primers'
include { REMOVE_PRIMERS }      from './modules/primers'
include { PRIMER_ASSIGNMENT } from './modules/primers'
include { AUTO_TRIM }         from './modules/denoise'
include { TRUNC_POLICY }      from './modules/denoise'
include { FILTER_TRIM }       from './modules/denoise'
include { LEARN_ERRORS }      from './modules/denoise'
include { DENOISE }           from './modules/denoise'
include { MERGE_SEQTABS }     from './modules/merge'
include { REMOVE_CHIMERAS }   from './modules/merge'
include { FILTER_SEQTAB }     from './modules/merge'

// Stage B: Taxonomy, phylogeny, renormalization
include { ASSIGN_TAXONOMY }   from './modules/taxonomy'
include { BUILD_PHYLOGENY }   from './modules/phylogeny'
include { RENORMALIZE }       from './modules/renormalize'

// Stage C: Metadata, clustering, network, visualization
include { LOAD_METADATA }     from './modules/metadata'
include { CLUSTER_TSNE }      from './modules/cluster'
include { NETWORK_SPARCC }    from './modules/network'
include { BUILD_VIZ }         from './modules/shiny'
include { BUNDLE_VIZ_SITE }  from './modules/shiny'

// ============================================================================
// Main workflow
// ============================================================================


// ── Sequencing-run recovery ──────────────────────────────────────────────────
//
// `plate` groups samples for TRUNC_POLICY (which collapses the group's per-sample
// truncation lengths by an explicit policy) and LEARN_ERRORS
// (one DADA2 error model per group). The grouping has to be the sequencing run:
// one sample per group fits a parametric error model on as few as ~1.4k reads,
// and picks a truncation length from one sample's ragged length distribution.
//
// The grouping is recoverable. SRA rewrites the read identifier, but
// fasterq-dump keeps the submitter's original name in field 2:
//
//   @SRR38958117.1 M04437:232:000000000-BNWCH:1:1101:16513:1933 length=251
//                  └── instrument:run:flowcell:lane ──┘
//
// which is exactly the unit DADA2 wants an error model for. Present in 822/822
// runs across every dataset checked (NCBI and DDBJ), but it depends on the
// submitter having uploaded original names, so this always falls back.
def runIdFromReads(reads) {
    try {
        // Indexable unless it is a single path: asking "is it an array?" needs
        // Object[], which the strict parser rejects, so ask the other way round.
        def single = reads instanceof CharSequence || reads instanceof java.nio.file.Path
        def first = single ? reads : reads[0]
        def path = first as java.nio.file.Path
        def raw = java.nio.file.Files.newInputStream(path)
        def stream = path.toString().endsWith('.gz') ? new java.util.zip.GZIPInputStream(raw) : raw
        def line = new BufferedReader(new InputStreamReader(stream)).readLine()
        stream.close()
        if (!line) return null
        def fields = line.trim().split(/\s+/)
        if (fields.size() < 2) return null
        def parts = fields[1].split(':')
        // instrument:run:flowcell:lane — anything shorter isn't an Illumina header
        if (parts.size() < 4) return null
        // The group name becomes part of a filename (e.g. <plate>_errF.pkl), and
        // Nextflow escapes colons in paths, which then don't resolve. Keep it
        // filesystem-safe while still readable as instrument_run_flowcell_lane.
        return parts[0..3].join(':').replaceAll(/[^A-Za-z0-9._-]/, '_')
    }
    catch (Exception e) {
        return null
    }
}

// Grouping key for AUTO_TRIM and LEARN_ERRORS: the sequencing run when we can
// recover it, else the historical accession/plate split.
//
// The assay is part of the key because one flowcell can carry more than one
// amplicon, and those differ in both length and error structure — a 16S V3-V4
// product and an 18S V4 product on the same run want their own truncation and
// their own error model. Where a run carries a single assay this is a no-op.
def plateFor(meta, reads, assay = null) {
    def runId = runIdFromReads(reads)
    def base = runId
    if (!base) {
        def parts = meta.id.split('_')
        base = parts.size() > 2 ? parts[0..1].join('_') : parts[0]
    }
    if (assay && assay != 'none') {
        base = "${base}__${assay}".replaceAll(/[^A-Za-z0-9._-]/, '_')
    }
    return base
}

workflow {

    // --help lives here rather than at the top of the file: a statement outside
    // a process, workflow or function is not script-level syntax any more.
    if (params.help) {
        helpMessage()
        System.exit(0)
    }
    validateParams()

    // 1. Discover input reads
    if (params.samplesheet) {
        // ── Samplesheet mode: CSV with sample_id, run, plate_id, r1_path, r2_path, primer_pair ──
        ch_samplesheet = Channel.fromPath(params.samplesheet)
            .splitCsv(header: true)
            .filter { row ->
                // Only include selected rows
                (row.selected ?: '1') == '1'
            }
            .map { row ->
                def meta = [
                    id: "${row.run}_${row.sample_id}",
                    sample_id: row.sample_id,
                    run: row.run,
                    plate: "${row.run}_${row.plate_id}",
                    primer_pair: row.primer_pair ?: '',
                ]
                [meta, file(row.r1_path), file(row.r2_path)]
            }

        // Filter by primer pair if specified
        if (params.primers) {
            ch_samplesheet = ch_samplesheet.filter { meta, r1, r2 ->
                meta.primer_pair == params.primers
            }
        }

        ch_demuxed = ch_samplesheet

    } else if (params.input) {
        def input_dir = file(params.input)
        if (!input_dir.isDirectory()) {
            error "ERROR: --input directory does not exist: ${params.input}\nRun with --help for usage."
        }

        // 2. Optional demultiplexing
        if (params.run_demultiplex) {
            ch_raw = Channel.fromPath("${params.input}/*.fastq.gz")
                .map { fastq ->
                    def name = fastq.baseName.replace('.fastq', '')
                    def lane = name.split('_')[0]
                    [[id: lane], fastq]
                }
                .groupTuple()

            DEMULTIPLEX(ch_raw,
                        file(params.forward_bcs),
                        file(params.reverse_bcs))
            ch_demuxed = DEMULTIPLEX.out.reads.flatten()
                .filter { it.name.endsWith('.fastq.gz') }
                .map { fastq ->
                    def name = fastq.baseName.replace('.fastq', '').replace('.ctrimmed', '')
                    [[id: name], fastq]
                }
        } else {
            // Pair up R1/R2 files by sample prefix
            ch_demuxed = Channel.fromFilePairs("${params.input}/*_{R1,R2,1,2}*.fastq.gz", flat: true)
                .map { sample_id, r1, r2 ->
                    [[id: sample_id], r1, r2]
                }
        }
    }

    // 3. Remove primers with cutadapt (skip if input is already trimmed)
    // Cutadapt's own per-adapter counts are the record of which assay each sample
    // really carries, so collect the logs from whichever branch ran. Stays empty
    // when primer removal is skipped, and BUILD_VIZ handles that.
    ch_cutadapt_logs = Channel.empty()
    // Whether the assay is ribosomal at all. Only primer detection can tell —
    // when primers are supplied we are not in a position to second-guess them,
    // so assume they mean what an rRNA reference database expects.
    ch_ribosomal = Channel.value(true)
    if (params.skip_primer_removal || params.samplesheet) {
        // Samplesheet mode: primer removal happens per primer_pair, or data is pre-trimmed
        ch_trimmed = ch_demuxed
            .map { meta, r1, r2 ->
                if (!meta.plate) {
                    meta = meta + [plate: plateFor(meta, r1)]
                }
                [meta, r1, r2]
            }

        // If samplesheet has primer_pair and we have primer files, run cutadapt
        if (params.samplesheet && !params.skip_primer_removal) {
            ch_with_primers = ch_trimmed.filter { meta, r1, r2 -> meta.primer_pair }
                .map { meta, r1, r2 ->
                    def pf = file("${projectDir}/primers/primers-${meta.primer_pair}.fa")
                    [meta, r1, r2, pf]
                }

            REMOVE_PRIMERS(ch_with_primers)
            ch_cutadapt_logs = REMOVE_PRIMERS.out.log
            ch_trimmed = REMOVE_PRIMERS.out.reads
                .map { meta, r1, r2, assay ->
                    [meta + [plate: plateFor(meta, r1, assay)], r1, r2]
                }
        }
    } else if (params.primers_fwd && params.primers_rev) {
        // REMOVE_PRIMERS applies one primer file to both reads (cutadapt -g/-G
        // selects the matching primer per read), exactly like the samplesheet and
        // auto-detect paths whose primers-*.fa hold both fwd and rev. Passing only
        // primers_fwd here trimmed R2 with the forward primer, so cutadapt's
        // --discard-untrimmed dropped ~every pair. Combine fwd + rev into one file.
        ch_primers_combined = Channel
            .fromPath([params.primers_fwd, params.primers_rev])
            .collectFile(name: 'primers_combined.fa')
        ch_with_primers = ch_demuxed.combine(ch_primers_combined)
        REMOVE_PRIMERS(ch_with_primers)
        ch_cutadapt_logs = REMOVE_PRIMERS.out.log
        ch_trimmed = REMOVE_PRIMERS.out.reads
            .map { meta, r1, r2, assay ->
                [meta + [plate: plateFor(meta, r1, assay)], r1, r2]
            }
    } else {
        // Detect per sample, then trim every sample with the same resolved set,
        // so which assay a sample belongs to is decided by what cutadapt matches
        // rather than by which sample's own detection it happened to be given.
        DETECT_PRIMERS(ch_demuxed)
        // report is optional, so an all-samples-failed detection would leave this
        // empty, RESOLVE_PRIMER_SET would never run, and the pipeline would end
        // with no output and no error. Say so instead.
        ch_detections = DETECT_PRIMERS.out.report
            .collect(sort: true)
            .ifEmpty { error "no primers were detected in any sample — nothing to trim with" }
        RESOLVE_PRIMER_SET(ch_detections)
        ch_ribosomal = RESOLVE_PRIMER_SET.out.report
            .map { rep -> new groovy.json.JsonSlurper().parse(rep.toFile()).ribosomal != false }
            .first()
        // Trim from ch_demuxed rather than from DETECT_PRIMERS' own output: every
        // sample is trimmed with the resolved set, so one whose detection found
        // nothing must still reach cutadapt.
        //
        // .first() makes the resolved set a value channel, so every sample is
        // paired with it. Without it the set is a queue channel of one item, the
        // pairing stops as soon as that item is consumed, and every sample after
        // the first is silently dropped.
        ch_with_primers = ch_demuxed
            .combine(RESOLVE_PRIMER_SET.out.primers.first())
        REMOVE_PRIMERS(ch_with_primers)
        ch_cutadapt_logs = REMOVE_PRIMERS.out.log
        ch_trimmed = REMOVE_PRIMERS.out.reads
            .map { meta, r1, r2, assay ->
                [meta + [plate: plateFor(meta, r1, assay)], r1, r2]
            }
    }

    // A group too small to fit an error model is pooled rather than fitted alone.
    // DADA2 fits a parametric model over quality bins; below a few samples most
    // bins are noise, and a noisy model is worse than a slightly mismatched one
    // shared with neighbours. Pool and say so, rather than fit in silence (#7).
    ch_trimmed = ch_trimmed
        .map { meta, r1, r2 -> [meta.plate, meta, r1, r2] }
        .groupTuple(by: 0)
        .flatMap { plate, metas, r1s, r2s ->
            def pooled = metas.size() < params.min_group_samples
            if (pooled) {
                log.warn "Group '${plate}' has only ${metas.size()} sample(s) — " +
                         "pooling it for truncation and error learning " +
                         "(--min_group_samples ${params.min_group_samples})"
            }
            def key = pooled ? 'pooled' : plate
            (0..<metas.size()).collect { i -> [metas[i] + [plate: key], r1s[i], r2s[i]] }
        }

    // 4a. Auto-detect truncation lengths (or use explicit params)
    if (params.auto_trim && params.truncLen_fwd == 0 && params.truncLen_rev == 0) {
        AUTO_TRIM(ch_trimmed)

        // Collapse each group's per-sample values by --trunc_policy.
        // Read lengths ride along so TRUNC_POLICY can cap the floor at the reads
        // that exist: a floor above the read length zeroes every sample (issue #4).
        ch_trunc_policy_in = AUTO_TRIM.out.trim_params
            .map { meta, trunc_fwd, trunc_rev, read_fwd, read_rev ->
                [meta.plate, meta.id, trunc_fwd as Integer, trunc_rev as Integer,
                 read_fwd as Integer, read_rev as Integer]
            }
            .groupTuple(by: 0)

        // Primer sequences let TRUNC_POLICY resolve the expected fragment length
        // from the primer database. Only the explicit --primers_* path can supply
        // them here; without them the overlap check falls back to the structural
        // floor and says it is unchecked.
        ch_trunc_primers = (params.primers_fwd && params.primers_rev)
            ? Channel.fromPath([params.primers_fwd, params.primers_rev]).collect()
            : Channel.value([])

        TRUNC_POLICY(ch_trunc_policy_in, ch_trunc_primers)

        // Join trim params back to individual samples by plate
        ch_filter_input = ch_trimmed
            .map { meta, r1, r2 -> [meta.plate, meta, r1, r2] }
            .combine(TRUNC_POLICY.out.trim_params, by: 0)
            .map { plate, meta, r1, r2, trunc_fwd, trunc_rev ->
                [meta, r1, r2, trunc_fwd, trunc_rev]
            }

        FILTER_TRIM(ch_filter_input)
    } else {
        // Use explicit params — add them to each sample tuple
        ch_filter_input = ch_trimmed.map { meta, r1, r2 ->
            [meta, r1, r2, params.truncLen_fwd, params.truncLen_rev]
        }
        FILTER_TRIM(ch_filter_input)
    }

    // 4b. Group filtered reads by plate
    ch_by_plate = FILTER_TRIM.out.reads
        .filter { meta, r1, r2 -> r1.size() > 0 && r2.size() > 0 }
        .map { meta, r1, r2 -> [meta.plate, r1, r2] }
        .groupTuple(by: 0)
        .map { plate, r1s, r2s ->
            [[id: plate, plate: plate], r1s, r2s]
        }

    // 4c. Learn error rates per-plate
    LEARN_ERRORS(ch_by_plate)

    // 4d. Denoise per-plate using learned error rates
    //     Join filtered reads with error models by plate name
    ch_plate_reads = ch_by_plate
        .map { meta, r1s, r2s -> [meta.plate, meta, r1s, r2s] }

    ch_plate_errors = LEARN_ERRORS.out.error_models
        .map { plate, meta, errF, errR -> [plate, errF, errR] }

    ch_denoise_input = ch_plate_reads
        .join(ch_plate_errors, by: 0)
        .map { plate, meta, r1s, r2s, errF, errR ->
            [meta, r1s, r2s, errF, errR]
        }

    DENOISE(ch_denoise_input)


    // One row per sample: which primer actually matched, and how many reads.
    PRIMER_ASSIGNMENT(ch_cutadapt_logs.collect())
    ch_primer_assignment = PRIMER_ASSIGNMENT.out.tsv

    // 5. Merge all sequence tables and remove chimeras
    ch_all_seqtabs = DENOISE.out.seqtab.map { meta, rds -> rds }.collect()

    MERGE_SEQTABS(ch_all_seqtabs)

    // Per-plate chimera removal already happens in DENOISE; the post-merge pass
    // catches cross-plate chimeras.
    REMOVE_CHIMERAS(MERGE_SEQTABS.out.seqtab)
    FILTER_SEQTAB(REMOVE_CHIMERAS.out.seqtab)

    // ======================================================================
    // Stage B: Taxonomy, phylogeny, renormalization
    // ======================================================================

    // 6. Assign taxonomy — one task per reference database (parallel)
    //    Build a channel of (db_name, db_path, tax_levels) tuples from params
    if (params.ref_databases) {
        def db_entries = params.ref_databases.tokenize(';')
        ch_databases = Channel.from(db_entries).map { entry ->
            // Format: "name:path:Level1,Level2,..." or "name:path" (no levels)
            def parts = entry.tokenize(':')
            def db_name = parts[0].trim()
            def db_path = file(parts[1].trim())
            def tax_levels = parts.size() > 2 ? parts[2].trim() : null
            [db_name, db_path, tax_levels]
        }
        // An rRNA database has nothing to say about a marker that is not rRNA,
        // but it answers anyway: PRJNA1287210 amplifies a denitrification gene,
        // and SILVA assigned 1047 of its 1053 ASVs to Eukaryota in a study of
        // bacteria. Wrong labels that carry confidence are worse than none, so
        // where the primers match neither SSU reference, skip taxonomy and leave
        // the ASV table to stand on its own.
        .combine(ch_ribosomal)
        .branch { _n, _p, _l, ribosomal ->
            assign: ribosomal
            skip:   true
        }

        ch_databases.skip.view { n, _p, _l, _r ->
            "[WARN] skipping taxonomy against ${n}: the detected primers match " +
            "neither the 16S nor the 18S reference, so this is not a ribosomal " +
            "assay and an rRNA database would label it confidently and wrongly"
        }

        ASSIGN_TAXONOMY(ch_databases.assign.map { n, p, l, _r -> [n, p, l] },
                        FILTER_SEQTAB.out.seqtab)

        // 9. Renormalize using the primary taxonomy database: the first one named
        //    in --ref_databases, selected by that name. Taking whichever task
        //    emitted first would work only while the process emits in input
        //    order, which is a scheduler directive set in nextflow.config — so
        //    which database classified the groups depended on a setting three
        //    files away, and would silently change to whichever of SILVA or PR2
        //    finished first without it.
        def primary_db = db_entries[0].tokenize(':')[0].trim()
        ch_primary_tax = ASSIGN_TAXONOMY.out.taxonomy
            .filter { db_name, _tax, _boot -> db_name == primary_db }
            .first()
        RENORMALIZE(FILTER_SEQTAB.out.seqtab, ch_primary_tax)
    }

    // 7. Build phylogeny (optional)
    if (params.run_phylogeny) {
        BUILD_PHYLOGENY(FILTER_SEQTAB.out.seqtab)
    }

    // ======================================================================
    // Stage C: Metadata, clustering, network, visualization
    // ======================================================================

    // 8. Load metadata (optional — if metadata file provided)
    if (params.metadata) {
        LOAD_METADATA(FILTER_SEQTAB.out.seqtab,
                      file(params.metadata),
                      params.sample_id_column)
        ch_metadata = LOAD_METADATA.out.metadata
    } else {
        ch_metadata = Channel.empty()
    }

    // 10. t-SNE clustering (samples + ASVs)
    CLUSTER_TSNE(FILTER_SEQTAB.out.seqtab)

    // 11. SparCC network analysis (requires renormalized data)
    if (params.ref_databases) {
        NETWORK_SPARCC(RENORMALIZE.out.merged, params.min_prevalence)

        // 12. Build visualization (bundle all outputs)
        ch_tax_files = ASSIGN_TAXONOMY.out.taxonomy
            .map { db_name, tax, boot -> [tax, boot] }
            .flatten()
            .collect()

        // t-SNE and network are optional enhancements running under
        // errorStrategy 'ignore'. Feed sentinels when their channels are empty
        // (the process failed/was skipped) so BUILD_VIZ still runs — the viz
        // just loses scatter coords / network edges instead of vanishing.
        BUILD_VIZ(
            FILTER_SEQTAB.out.seqtab,
            RENORMALIZE.out.merged,
            ch_tax_files,
            ch_metadata.ifEmpty(file('NO_METADATA')),
            CLUSTER_TSNE.out.sample_tsne.ifEmpty(file('NO_SAMPLE_TSNE')),
            CLUSTER_TSNE.out.seq_tsne.ifEmpty(file('NO_SEQ_TSNE')),
            NETWORK_SPARCC.out.correlations.ifEmpty(file('NO_NETWORK')),
            ch_primer_assignment.ifEmpty(file('NO_PRIMER_ASSIGNMENT'))
        )

        // 13. Optional: bundle static viz site (requires Node.js)
        if (params.build_viz_site) {
            BUNDLE_VIZ_SITE(
                BUILD_VIZ.out.json.collect(),
                BUILD_VIZ.out.json_gz.collect()
            )
        }
    }

    // ============================================================================
    // Pipeline completion handler
    // ============================================================================

    workflow.onComplete = {
        def msg = """\
            Pipeline completed at : ${workflow.complete}
            Duration              : ${workflow.duration}
            Success               : ${workflow.success}
            Exit status           : ${workflow.exitStatus}
            Output directory      : ${params.outdir}
            """.stripIndent()

        println msg

        // Run manifest — record what was actually run (pipeline version, resolved parameters,
        // per-process tool images, reference databases) so downstream reporting (omc-platform
        // Methods drafting + agents) can state it instead of "not specified in the outputs".
        // Published into viz/ alongside provenance.json.
        try {
            def vizDir = new File("${params.outdir}/viz")
            if (vizDir.exists()) {
                // Provenance that survives the container. workflow.commitId is null for
                // every run launched from the .sif (there is no git repo inside), which is
                // how this pipeline actually runs — so fall back to the SHA baked in at
                // build time. scriptId is the MD5 of main.nf and pins the source even when
                // both of those are missing.
                def sysEnv = System.getenv()
                def bakedSha = sysEnv['DANASEQ_GIT_SHA']
                def toolVersions = null
                try {
                    def probe = new File("${workflow.projectDir}/bin/tool_versions.sh")
                    if (probe.exists()) {
                        def proc = ["bash", probe.absolutePath].execute()
                        def out = new StringBuffer()
                        def err = new StringBuffer()
                        proc.consumeProcessOutput(out, err)
                        proc.waitForOrKill(120000)
                        toolVersions = new groovy.json.JsonSlurper().parseText(out.toString())
                    }
                } catch (Exception e) {
                    println "[WARN] tool version probe failed: ${e.message}"
                }

                def manifest = [
                    pipeline            : "${workflow.manifest.name} v${workflow.manifest.version}",
                    // Null when run from the .sif; bakedSha is then the only real answer.
                    revision            : (workflow.revision ?: workflow.commitId ?: null),
                    commit_id           : (workflow.commitId ?: bakedSha),
                    commit_source       : (workflow.commitId ? 'git checkout'
                                           : (bakedSha ? 'baked into container at build' : 'UNKNOWN')),
                    container_git_ref   : sysEnv['DANASEQ_GIT_REF'],
                    container_built     : sysEnv['DANASEQ_BUILD_DATE'],
                    // MD5 of main.nf — identifies the source even with no git and no build arg.
                    script_id           : workflow.scriptId,
                    // What each tool actually reported, not what the env spec asked for.
                    tool_versions       : toolVersions,
                    nextflow_version    : "${nextflow.version}",
                    command_line        : workflow.commandLine,
                    started             : "${workflow.start}",
                    completed           : "${workflow.complete}",
                    duration            : "${workflow.duration}",
                    success             : workflow.success,
                    // process -> container image; image tags give tool versions
                    containers          : workflow.container,
                    reference_databases : params.ref_databases,
                    denoise             : 'papa2 (DADA2-compatible)',
                    parameters          : [
                        skip_primer_removal  : params.skip_primer_removal,
                        primer_error_rate    : params.primer_error_rate,
                        // added on main while this branch was open — a manifest that omits
                        // a parameter is the same failure as one that reports a stale version
                        primer_detect_reads  : params.primer_detect_reads,
                        min_group_samples    : params.min_group_samples,
                        auto_trim            : params.auto_trim,
                        auto_trim_min_quality: params.auto_trim_min_quality,
                        maxEE                : params.maxEE,
                        truncQ               : params.truncQ,
                        maxN                 : params.maxN,
                        truncLen_fwd         : params.truncLen_fwd,
                        truncLen_rev         : params.truncLen_rev,
                        trunc_policy         : params.trunc_policy,
                        min_overlap          : params.min_overlap,
                        min_seq_length       : params.min_seq_length,
                        min_reads            : params.min_reads,
                        min_samples          : params.min_samples,
                        min_seqs             : params.min_seqs,
                    ],
                ]
                new File(vizDir, 'run_manifest.json').text =
                    groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(manifest))
                println "[INFO] Wrote run manifest: ${params.outdir}/viz/run_manifest.json"
            }
        } catch (Exception e) {
            println "[WARN] Could not write run_manifest.json: ${e.message}"
        }

        if (!workflow.success) {
            println "[WARNING] Pipeline completed with errors. Check .nextflow.log for details."
        }
    }

    workflow.onError = {
        println "[ERROR] Pipeline failed: ${workflow.errorMessage}"
    }
}
