#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Microscape Amplicon Pipeline - Nextflow DSL2
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
     Microscape Amplicon Pipeline
     https://github.com/rec3141/microscape
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

    Language:
      --lang STR           'R' or 'python' [default: python]

    Resources:
      --threads N          General thread count [default: 8]
      --denoise_cpus N     CPUs for denoising processes [default: 8]
      --denoise_memory S   Memory for denoising processes [default: '16 GB']
    """.stripIndent()
}

if (params.help) {
    helpMessage()
    System.exit(0)
}

// ============================================================================
// Parameter validation
// ============================================================================

if (!params.input && !params.samplesheet) {
    log.error "ERROR: --input or --samplesheet is required. Run with --help for usage."
    System.exit(1)
}

if (params.run_demultiplex && (!params.forward_bcs || !params.reverse_bcs)) {
    log.error "ERROR: --forward_bcs and --reverse_bcs are required when --run_demultiplex is enabled."
    System.exit(1)
}

// ============================================================================
// Module imports
// ============================================================================

// Stage A: Preprocessing and denoising
include { DEMULTIPLEX }       from './modules/demultiplex'
include { DETECT_PRIMERS }    from './modules/primers'
include { REMOVE_PRIMERS }    from './modules/primers'
include { PRIMER_ASSIGNMENT } from './modules/primers'
include { AUTO_TRIM }         from './modules/denoise'
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

workflow {

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
    if (params.skip_primer_removal || params.samplesheet) {
        // Samplesheet mode: primer removal happens per primer_pair, or data is pre-trimmed
        ch_trimmed = ch_demuxed
            .map { meta, r1, r2 ->
                if (!meta.plate) {
                    def parts = meta.id.split('_'); def plate = parts.size() > 2 ? parts[0..1].join('_') : parts[0]
                    meta = meta + [plate: plate]
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
                .map { meta, r1, r2 -> [meta, r1, r2] }
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
            .map { meta, r1, r2 ->
                def parts = meta.id.split('_'); def plate = parts.size() > 2 ? parts[0..1].join('_') : parts[0]
                def new_meta = meta + [plate: plate]
                [new_meta, r1, r2]
            }
    } else {
        // Auto-detect best primer pair per sample
        ch_primer_files = Channel.fromPath("${projectDir}/primers/primers-*.fa").collect()
        DETECT_PRIMERS(ch_demuxed, ch_primer_files)
        REMOVE_PRIMERS(DETECT_PRIMERS.out.detected)
        ch_cutadapt_logs = REMOVE_PRIMERS.out.log
        ch_trimmed = REMOVE_PRIMERS.out.reads
            .map { meta, r1, r2 ->
                def parts = meta.id.split('_'); def plate = parts.size() > 2 ? parts[0..1].join('_') : parts[0]
                def new_meta = meta + [plate: plate]
                [new_meta, r1, r2]
            }
    }

    // 4a. Auto-detect truncation lengths per plate (or use explicit params)
    if (params.auto_trim && params.truncLen_fwd == 0 && params.truncLen_rev == 0) {
        // Group trimmed reads by plate for per-plate quality profiling
        ch_plate_reads_for_trim = ch_trimmed
            .map { meta, r1, r2 -> [meta.plate, r1, r2] }
            .groupTuple(by: 0)
            .map { plate, r1s, r2s -> [plate, r1s.flatten() + r2s.flatten()] }

        AUTO_TRIM(ch_plate_reads_for_trim)

        // Join trim params back to individual samples by plate
        ch_filter_input = ch_trimmed
            .map { meta, r1, r2 -> [meta.plate, meta, r1, r2] }
            .combine(AUTO_TRIM.out.trim_params, by: 0)
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

    // Per-plate chimera removal already happens in DENOISE (both R and Python).
    // Post-merge chimera pass catches cross-plate chimeras.
    def effectiveEngine = params.denoise_engine ?: params.lang
    if (effectiveEngine == 'python') {
        REMOVE_CHIMERAS(MERGE_SEQTABS.out.seqtab)
        FILTER_SEQTAB(REMOVE_CHIMERAS.out.seqtab)
    } else {
        FILTER_SEQTAB(MERGE_SEQTABS.out.seqtab)
    }

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

        ASSIGN_TAXONOMY(ch_databases, FILTER_SEQTAB.out.seqtab)

        // 9. Renormalize using the primary taxonomy database
        //    Use the first database as the primary for group classification
        ch_primary_tax = ASSIGN_TAXONOMY.out.taxonomy.first()
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
}

// ============================================================================
// Pipeline completion handler
// ============================================================================

workflow.onComplete {
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
    // Published into viz/ alongside provenance.json. (microscape-nf #5)
    try {
        def vizDir = new File("${params.outdir}/viz")
        if (vizDir.exists()) {
            def manifest = [
                pipeline            : "${workflow.manifest.name} v${workflow.manifest.version}",
                revision            : (workflow.revision ?: workflow.commitId ?: null),
                commit_id           : workflow.commitId,
                nextflow_version    : "${nextflow.version}",
                command_line        : workflow.commandLine,
                started             : "${workflow.start}",
                completed           : "${workflow.complete}",
                duration            : "${workflow.duration}",
                success             : workflow.success,
                // process -> container image; image tags give tool versions
                containers          : workflow.container,
                reference_databases : params.ref_databases,
                denoise             : (params.denoise_engine ?: params.lang) == 'python' ? 'papa2 (DADA2-compatible)' : 'DADA2 (R)',
                parameters          : [
                    lang                 : params.lang,
                    denoise_engine       : params.denoise_engine,
                    skip_primer_removal  : params.skip_primer_removal,
                    primer_error_rate    : params.primer_error_rate,
                    auto_trim            : params.auto_trim,
                    auto_trim_min_quality: params.auto_trim_min_quality,
                    maxEE                : params.maxEE,
                    truncQ               : params.truncQ,
                    maxN                 : params.maxN,
                    truncLen_fwd         : params.truncLen_fwd,
                    truncLen_rev         : params.truncLen_rev,
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

workflow.onError {
    println "[ERROR] Pipeline failed: ${workflow.errorMessage}"
}
