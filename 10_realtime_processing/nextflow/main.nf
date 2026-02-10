#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Dana Pipeline - Real-Time Nanopore Processing (Nextflow DSL2)
// ============================================================================
//
// Nextflow implementation of the 24_process_reads_optimized.sh pipeline.
// Provides DAG-based parallelism, native resume (-resume), and resource-aware
// scheduling without manual checkpoint logic or semaphores.
//
// Usage:
//   nextflow run main.nf --input /path/to/data --run_kraken --run_prokka -resume
//
// Watch mode for live sequencing:
//   nextflow run main.nf --input /path/to/data --watch --run_db_integration -resume
//
// ============================================================================

// Import modules
include { VALIDATE_FASTQ }   from './modules/validate'
include { QC_BBDUK }         from './modules/qc'
include { QC_FILTLONG }      from './modules/qc'
include { CONVERT_TO_FASTA } from './modules/qc'
include { KRAKEN2_CLASSIFY }  from './modules/kraken'
include { PROKKA_ANNOTATE }   from './modules/prokka'
include { SENDSKETCH }        from './modules/sketch'
include { HMM_SEARCH }        from './modules/hmm'
include { TETRAMER_FREQ }     from './modules/tetramer'
include { DB_INTEGRATION }    from './modules/db_integration'
include { DB_SYNC }           from './modules/db_integration'

// ============================================================================
// Input channel: discover FASTQ files from Nanopore output structure
// ============================================================================
// Expected filename format: FLOWCELL_pass_barcodeNN_*.fastq.gz
// Extracts metadata (flowcell, barcode) from filename, carried through entire DAG

def create_fastq_channel() {
    def pattern = "${params.input}/**/fastq_pass/barcode*/*.fastq.gz"

    // watchPath monitors for new files during live sequencing runs;
    // 'create' event fires when MinKNOW atomically moves completed files into barcode dirs
    def ch_raw = params.watch
        ? Channel.watchPath(pattern, 'create')
        : Channel.fromPath(pattern, checkIfExists: true)

    ch_raw
        .filter { it.size() >= params.min_file_size }
        .map { fastq ->
            def name = fastq.baseName.replace('.fastq', '')
            def parts = name.split('_')
            def flowcell = parts[0]
            def barcode = parts[2]
            def meta = [
                id:       name,
                flowcell: flowcell,
                barcode:  barcode,
                sample:   "${flowcell}_${barcode}"
            ]
            [ meta, fastq ]
        }
}

// ============================================================================
// Main workflow
// ============================================================================

workflow {

    // Discover and validate input files
    ch_input = create_fastq_channel()

    // Stage 1: Validate FASTQ integrity (repair corrupted gzip)
    VALIDATE_FASTQ(ch_input)

    // Stage 2: Quality control pipeline
    QC_BBDUK(VALIDATE_FASTQ.out.validated)
    QC_FILTLONG(QC_BBDUK.out.trimmed)
    CONVERT_TO_FASTA(QC_FILTLONG.out.filtered)

    // FASTA channel feeds all downstream analyses
    ch_fasta = CONVERT_TO_FASTA.out.fasta

    // Stage 3: Optional downstream analyses (controlled by params)

    // Kraken2 taxonomic classification (maxForks=1 in process definition)
    if (params.run_kraken) {
        KRAKEN2_CLASSIFY(ch_fasta)
    }

    // Sendsketch profiling
    if (params.run_sketch) {
        SENDSKETCH(ch_fasta)
    }

    // Tetranucleotide frequency
    if (params.run_tetra) {
        TETRAMER_FREQ(ch_fasta)
    }

    // Prokka annotation
    if (params.run_prokka) {
        PROKKA_ANNOTATE(ch_fasta)

        // HMM search on Prokka proteins (requires Prokka output)
        if (params.hmm_databases) {
            // Build channel of HMM database files
            ch_hmm_dbs = Channel
                .of(params.hmm_databases.split(','))
                .map { it.trim() }
                .map { dbpath ->
                    def dbfile = file(dbpath)
                    def dbname = dbfile.baseName
                    [ dbname, dbfile ]
                }

            // Cartesian product: each protein file x each HMM database
            ch_hmm_input = PROKKA_ANNOTATE.out.proteins
                .combine(ch_hmm_dbs)
                .map { meta, faa, dbname, hmm_db ->
                    [ meta, faa, dbname, hmm_db ]
                }

            HMM_SEARCH(ch_hmm_input)
        }
    }

    // Stage 4: DuckDB integration (post-pipeline)
    // Collects all published output directories and loads into DuckDB
    // Uses string paths since R scripts operate on the host filesystem
    if (params.run_db_integration) {
        def abs_outdir = file(params.outdir).toAbsolutePath().toString()

        if (params.watch) {
            // Watch mode: channel never closes so .collect() would block forever.
            // Instead, use a periodic timer to scan output dirs and sync DB.
            // R scripts are idempotent (import_log tracks what's loaded).
            ch_timer = Channel.interval(params.db_sync_minutes * 60 * 1000)
            DB_SYNC(ch_timer, abs_outdir, params.danadir)
        } else {
            // Batch mode: barrier approach — wait for all processes to finish
            // .collect() blocks until ALL mixed channels are drained
            ch_done = ch_fasta.map { meta, f -> "${abs_outdir}/${meta.flowcell}/${meta.barcode}" }
            if (params.run_kraken)  { ch_done = ch_done.mix(KRAKEN2_CLASSIFY.out.parsed.map  { meta, f -> "${abs_outdir}/${meta.flowcell}/${meta.barcode}" }) }
            if (params.run_sketch)  { ch_done = ch_done.mix(SENDSKETCH.out.sketch.map        { meta, f -> "${abs_outdir}/${meta.flowcell}/${meta.barcode}" }) }
            if (params.run_tetra)   { ch_done = ch_done.mix(TETRAMER_FREQ.out.lrn.map        { meta, f -> "${abs_outdir}/${meta.flowcell}/${meta.barcode}" }) }
            if (params.run_prokka)  { ch_done = ch_done.mix(PROKKA_ANNOTATE.out.tsv.map      { meta, f -> "${abs_outdir}/${meta.flowcell}/${meta.barcode}" }) }
            if (params.run_prokka && params.hmm_databases) {
                ch_done = ch_done.mix(HMM_SEARCH.out.tsv.map { meta, f -> "${abs_outdir}/${meta.flowcell}/${meta.barcode}" })
            }

            // Wait for everything, then deduplicate to unique barcode dirs
            ch_barcode_dirs = ch_done
                .collect()
                .flatMap { it.unique() }

            DB_INTEGRATION(ch_barcode_dirs)
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

    if (!workflow.success) {
        println "[WARNING] Pipeline completed with errors. Check .nextflow.log for details."
    }
}

workflow.onError {
    println "[ERROR] Pipeline failed: ${workflow.errorMessage}"
}
