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
// Watch mode for live sequencing (--input = parent of run directory):
//   nextflow run main.nf --input /path/to/runs --watch --run_db_integration
//   # Or point --input at the run dir with a shallower glob:
//   nextflow run main.nf --input /path/to/run_dir --watch --watch_glob 'fastq_pass/barcode*/*.fastq.gz'
//
// ============================================================================

// ============================================================================
// Help message
// ============================================================================

def helpMessage() {
    log.info """
    =========================================
     Dana Pipeline - Real-Time Nanopore Processing
     https://github.com/rec3141/danaSeq
    =========================================

    Usage:
      nextflow run main.nf --input /path/to/data [options] -resume

    Required:
      --input DIR        Path to directory CONTAINING fastq_pass/
                         (e.g. /data/run1, where /data/run1/fastq_pass/barcode01/ exists)
      --outdir DIR       Output directory [default: results]

    Analysis flags:
      --run_kraken       Kraken2 taxonomic classification (requires --kraken_db)
      --run_prokka       Prokka gene annotation
      --run_sketch       Sendsketch taxonomic profiling
      --run_tetra        Tetranucleotide frequency analysis
      --hmm_databases    Comma-separated HMM file paths (requires --run_prokka)

    Database paths:
      --kraken_db DIR    Path to Kraken2 database directory
      --danadir DIR      Path to R scripts for DuckDB integration

    Integration:
      --run_db_integration  Load results into DuckDB after processing
      --cleanup             Compress/delete source files after DuckDB import

    Watch mode (live sequencing):
      --watch               Monitor for new FASTQ files continuously
      --db_sync_minutes N   DB sync interval in minutes [default: 10]

    QC parameters:
      --min_readlen N    Minimum read length in bp [default: 1500]
      --keep_percent N   Filtlong keep percent [default: 80]
      --min_file_size N  Minimum FASTQ size in bytes [default: 1000000]

    Examples:
      # Basic QC only
      nextflow run main.nf --input /data/run1 -resume

      # Full pipeline with Kraken2
      nextflow run main.nf --input /data/run1 \\
          --run_kraken --kraken_db /path/to/krakendb \\
          --run_prokka --run_sketch --run_tetra -resume

      # Kitchen sink — all modules, all options with defaults
      nextflow run main.nf --input /data/run1 --outdir results \\
          --run_kraken --kraken_db /path/to/krakendb \\
          --run_prokka \\
          --hmm_databases /path/to/CANT-HYD.hmm,/path/to/FOAM.hmm \\
          --run_sketch \\
          --run_tetra \\
          --run_db_integration \\
          --cleanup \\
          --min_readlen 1500 \\
          --keep_percent 80 \\
          --min_file_size 1000000 \\
          -resume

      # Kitchen sink — watch mode for live sequencing
      nextflow run main.nf --input /data/runs --outdir results \\
          --watch --db_sync_minutes 10 \\
          --run_kraken --kraken_db /path/to/krakendb \\
          --run_prokka \\
          --hmm_databases /path/to/CANT-HYD.hmm \\
          --run_sketch \\
          --run_tetra \\
          --run_db_integration

      # Using launcher script (local conda or Docker)
      ./run-realtime.sh --input /data/run1 --outdir /data/output \\
          --run_kraken --kraken_db /path/to/krakendb --run_prokka
      ./run-realtime.sh --docker --input /data/run1 --outdir /data/output \\
          --run_kraken --kraken_db /path/to/krakendb --run_prokka

      # Quick test with bundled test data
      nextflow run main.nf --input test-data -profile test -resume

    Input structure:
      --input must point to a directory CONTAINING fastq_pass/:
        /data/run1/                  <-- use this as --input
        └── fastq_pass/
            ├── barcode01/
            │   └── FLOWCELL_pass_barcode01_*.fastq.gz
            └── barcode02/
                └── ...

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    System.exit(0)
}

// ============================================================================
// Parameter validation
// ============================================================================

if (!params.input) {
    log.error "ERROR: --input is required. Provide path to directory containing fastq_pass/. Run with --help for usage."
    System.exit(1)
}
if (params.run_kraken && !params.kraken_db) {
    log.error "ERROR: --kraken_db is required when using --run_kraken. Provide path to Kraken2 database directory."
    System.exit(1)
}
if (params.run_db_integration && !params.danadir) {
    log.error "ERROR: --danadir is required when using --run_db_integration. Provide path to R scripts directory."
    System.exit(1)
}

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
include { CLEANUP }           from './modules/db_integration'

// ============================================================================
// Input channel: discover FASTQ files from Nanopore output structure
// ============================================================================
// Expected filename format: FLOWCELL_pass_barcodeNN_*.fastq.gz
// Extracts metadata (flowcell, barcode) from filename, carried through entire DAG

def create_fastq_channel() {
    // Batch mode uses ** recursive glob (fromPath handles this natively).
    // Watch mode uses a flat glob (watch_glob param) because Java's WatchService
    // cannot recursively monitor subdirectories via **.  Default watch_glob is
    // '*/fastq_pass/barcode*/*.fastq.gz' — set --input to the parent of the run
    // directory, or override --watch_glob to match your directory depth.

    def glob_pattern = "${params.input}/**/fastq_pass/barcode*/*.fastq.gz"
    def ch_raw
    if (params.watch) {
        ch_raw = Channel.watchPath("${params.input}/${params.watch_glob}", 'create')
    } else {
        // Validate input before creating channel
        def input_dir = file(params.input)
        if (!input_dir.isDirectory()) {
            error "ERROR: --input directory does not exist: ${params.input}\nRun with --help for usage."
        }

        // Try to find files and give targeted advice on failure
        def found = file("${params.input}/fastq_pass/barcode*/*.fastq.gz")
        if (!found) {
            found = file("${params.input}/*/fastq_pass/barcode*/*.fastq.gz")
        }
        if (!found) {
            def has_fastq_pass = file("${params.input}/fastq_pass").isDirectory()
            def has_barcode_dirs = file("${params.input}/barcode*/*.fastq.gz")

            def hint = ""
            if (has_fastq_pass) {
                hint = "    ${params.input} contains fastq_pass/ but no barcode*/*.fastq.gz files inside it.\n    Check: ls ${params.input}/fastq_pass/barcode*/*.fastq.gz"
            } else if (has_barcode_dirs) {
                hint = "    It looks like --input points at fastq_pass/ itself.\n    Point --input one level UP:\n      --input ${input_dir.getParent()}"
            } else {
                hint = "    Expected structure:\n      <input>/fastq_pass/barcode01/*.fastq.gz\n    Point --input at the directory CONTAINING fastq_pass/."
            }

            error "ERROR: No FASTQ files found.\n\n    Searched: ${glob_pattern}\n${hint}\n\n    Run with --help for full usage information."
        }

        ch_raw = Channel.fromPath(glob_pattern)
    }

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
            // Instead, DB_SYNC runs as a long-lived process with an internal
            // sleep loop that periodically scans output dirs and syncs DB.
            // R scripts are idempotent (import_log tracks what's loaded).
            // Cleanup runs inline after each sync cycle if enabled.
            def sync_secs = params.db_sync_minutes * 60
            def cleanup_flag = params.cleanup ? "true" : "false"
            DB_SYNC(abs_outdir, params.danadir, sync_secs, cleanup_flag)
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

            // Post-DB cleanup: gzip fa/, delete kraken/sketch/tetra, compress prokka
            if (params.cleanup) {
                CLEANUP(DB_INTEGRATION.out)
            }
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
