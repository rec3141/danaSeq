#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Dana MAG Assembly Pipeline - Nextflow DSL2
// ============================================================================
//
// Nextflow implementation of the MAG assembly and binning workflow.
// Co-assembles all reads with Flye, maps reads back, runs five binners
// (SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin) in parallel, integrates with DAS_Tool,
// detects mobile genetic elements (viruses, plasmids, proviruses) with geNomad + CheckV,
// and identifies anti-phage defense systems with DefenseFinder.
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
     Dana MAG Assembly Pipeline
     https://github.com/rec3141/danaSeq
    =========================================

    Usage:
      nextflow run main.nf --input /path/to/reads [options] -resume

    Required:
      --input DIR        Input directory (see Input section below)
      --outdir DIR       Output directory [default: results]

    Assembly:
      --min_overlap N    Flye --min-overlap [default: 1000]
      --polish           Enable Flye polishing iterations [default: true]
      --dedupe           Enable BBDuk deduplication before assembly
      --filtlong_size N  Filtlong target bases (e.g. 40000000000); skip if not set

    Binning:
      --run_maxbin       Include MaxBin2 in consensus binning [default: true]
      --run_lorbin       Include LorBin in consensus binning [default: true]
      --run_comebin      Include COMEBin in consensus binning [default: true]
      --metabat_min_cls N  MetaBAT2 minimum cluster size [default: 50000]
      --lorbin_min_length N LorBin minimum bin size in bp [default: 80000]

    Annotation:
      --run_prokka       Run Prokka gene annotation on co-assembly [default: true]

    Taxonomy:
      --run_kaiju        Run Kaiju protein-level taxonomy on Prokka proteins [default: true]
      --kaiju_db PATH    Path to Kaiju database (*.fmi + nodes.dmp + names.dmp); null = skip

    Mobile Genetic Elements:
      --run_genomad      Run geNomad virus + plasmid detection [default: true]
      --genomad_db PATH  Path to geNomad database; null = skip geNomad
      --run_checkv       Run CheckV viral quality assessment [default: true]
      --checkv_db PATH   Path to CheckV database; null = skip CheckV
      --run_integronfinder  Run IntegronFinder integron detection [default: true]
      --run_islandpath    Run IslandPath-DIMOB genomic island detection [default: true]
      --run_macsyfinder   Run MacSyFinder secretion/conjugation detection [default: true]
      --macsyfinder_models PATH  Path to MacSyFinder models dir; null = skip
      --run_defensefinder Run DefenseFinder anti-phage defense detection [default: true]
      --defensefinder_models PATH  Path to DefenseFinder models dir; null = auto-download

    Quality:
      --checkm2_db PATH  Path to CheckM2 DIAMOND database; null = skip CheckM2

    Bin Refinement (NCLB):
      --run_nclb          Run NCLB bin refinement after DAS_Tool [default: false]
      --nclb_dir PATH     Path to NCLB repository (contains bin/, lib/, envs/)
      --nclb_base_url URL LLM server URL for conversations [default: http://localhost:1234/v1]
      --nclb_model NAME   LLM model name [default: auto-detect from server]
      --nclb_with_ani     Run minimap2 ANI during Elder investigations [default: false]

    Resources:
      --assembly_cpus N     CPUs for assembly [default: 24]
      --assembly_memory STR Memory for assembly [default: '64 GB']

    Examples:
      # Basic run with all defaults
      nextflow run main.nf --input /data/reads -resume

      # With Filtlong pre-filtering and no MaxBin2
      nextflow run main.nf --input /data/reads \\
          --filtlong_size 40000000000 --run_maxbin false -resume

      # Kitchen sink — all options with defaults
      nextflow run main.nf --input /data/reads --outdir results \\
          --dedupe \\
          --filtlong_size 40000000000 \\
          --min_overlap 1000 \\
          --run_maxbin true \\
          --metabat_min_cls 50000 \\
          --checkm2_db /path/to/checkm2_db \\
          --assembly_cpus 24 \\
          --assembly_memory '64 GB' \\
          -resume

      # Using launcher script (local conda or Docker)
      ./run-mag.sh --input /data/reads --outdir /data/output
      ./run-mag.sh --docker --input /data/reads --outdir /data/output

      # Quick test with bundled test data
      nextflow run main.nf --input test-data -profile test -resume

    Input:
      --input can be either:
        (a) Directory with *.fastq.gz files (one per sample)
        (b) Nanopore run or project directory containing fastq_pass/barcode*/
            at any depth. Files are concatenated per barcode automatically.
            Example: --input /data/Ebb_Flow  (finds all runs underneath)

    Output:
      results/
      ├── assembly/
      │   └── assembly.fasta
      ├── mapping/
      │   ├── *.sorted.bam
      │   └── depths.txt
      ├── binning/
      │   ├── semibin/contig_bins.tsv
      │   ├── metabat/contig_bins.tsv
      │   ├── maxbin/contig_bins.tsv
      │   ├── lorbin/contig_bins.tsv
      │   ├── comebin/contig_bins.tsv
      │   ├── dastool/
      │   │   ├── bins/*.fa
      │   │   ├── contig2bin.tsv
      │   │   ├── allbins.fa
      │   │   ├── bin_quality.tsv
      │   │   └── summary.tsv
      │   └── checkm2/
      │       └── quality_report.tsv  (if --checkm2_db set)
      ├── mge/                         (if --genomad_db set)
      │   ├── genomad/
      │   │   ├── virus_summary.tsv
      │   │   ├── plasmid_summary.tsv
      │   │   ├── virus.fna
      │   │   ├── plasmid.fna
      │   │   ├── virus_proteins.faa
      │   │   ├── plasmid_proteins.faa
      │   │   ├── virus_genes.tsv
      │   │   ├── plasmid_genes.tsv
      │   │   ├── provirus.tsv
      │   │   ├── provirus.fna
      │   │   ├── taxonomy.tsv
      │   │   └── genomad_summary.tsv
      │   ├── checkv/                  (if --checkv_db set)
      │   │   ├── quality_summary.tsv
      │   │   ├── viruses.fna
      │   │   └── proviruses.fna
      │   ├── integrons/               (if --run_integronfinder)
      │   │   ├── integrons.tsv
      │   │   └── summary.tsv
      │   ├── genomic_islands/        (if --run_islandpath)
      │   │   └── genomic_islands.tsv
      │   ├── macsyfinder/            (if --macsyfinder_models set)
      │   │   ├── all_systems.tsv
      │   │   └── all_systems.txt
      │   └── defensefinder/         (if --run_defensefinder)
      │       ├── systems.tsv
      │       ├── genes.tsv
      │       └── hmmer.tsv
      ├── taxonomy/                    (if --kaiju_db set)
      │   └── kaiju/
      │       ├── kaiju_genes.tsv      Per-gene Kaiju classifications
      │       └── kaiju_contigs.tsv    Per-contig taxonomy (majority vote)
      ├── binning/nclb/               (if --run_nclb + --nclb_dir set)
      │   ├── communities/*.fa         Refined community FASTAs
      │   ├── gathering.json           Identity cards + resonance data
      │   ├── proposals.json           LLM conversation proposals
      │   ├── elder_reports.json       SCG redundancy investigations
      │   ├── chronicle.json           Machine-readable decision log
      │   ├── chronicle.md            Human-readable narrative
      │   ├── contig2community.tsv     Contig membership assignments
      │   ├── quality_report.tsv       Community quality metrics
      │   └── valence_report.tsv       Per-contig valence scores
      └── pipeline_info/

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
    log.error "ERROR: --input is required. Provide path to directory containing *.fastq.gz files. Run with --help for usage."
    System.exit(1)
}

// Import modules
include { CONCAT_READS }        from './modules/preprocess'
include { ASSEMBLY_FLYE }       from './modules/assembly'
include { CALCULATE_TNF }       from './modules/assembly'
include { MAP_READS }           from './modules/mapping'
include { CALCULATE_DEPTHS }    from './modules/mapping'
include { BIN_SEMIBIN2 }        from './modules/binning'
include { BIN_METABAT2 }        from './modules/binning'
include { BIN_MAXBIN2 }         from './modules/binning'
include { BIN_LORBIN }          from './modules/binning'
include { BIN_COMEBIN }         from './modules/binning'
include { DASTOOL_CONSENSUS }   from './modules/binning'
include { CHECKM2 }             from './modules/binning'
include { PROKKA_ANNOTATE }     from './modules/annotation'
include { KAIJU_CLASSIFY }      from './modules/taxonomy'
include { GENOMAD_CLASSIFY }    from './modules/mge'
include { CHECKV_QUALITY }      from './modules/mge'
include { INTEGRONFINDER }      from './modules/mge'
include { ISLANDPATH_DIMOB }    from './modules/mge'
include { MACSYFINDER }         from './modules/mge'
include { DEFENSEFINDER }       from './modules/mge'
include { NCLB_GATHER }         from './modules/refinement'
include { NCLB_CONVERSE }       from './modules/refinement'
include { NCLB_ELDERS }         from './modules/refinement'
include { NCLB_INTEGRATE }      from './modules/refinement'

// ============================================================================
// Main workflow
// ============================================================================

workflow {

    // 1. Discover input reads — auto-detect nanopore barcode vs flat directory
    def input_dir = file(params.input)
    if (!input_dir.isDirectory()) {
        error "ERROR: --input directory does not exist: ${params.input}\nRun with --help for usage."
    }

    // Check for nanopore barcode structure at any depth
    def barcode_dirs = file("${params.input}/**/fastq_pass/barcode*", type: 'dir')
    def flat_fastqs  = file("${params.input}/*.fastq.gz")

    if (barcode_dirs) {
        // Nanopore barcode structure detected — concat per barcode
        log.info "Detected nanopore barcode structure: ${barcode_dirs.size()} barcode directories"
        ch_barcode_raw = Channel.fromPath("${params.input}/**/fastq_pass/barcode*", type: 'dir')
            .flatMap { dir ->
                def barcode = dir.name
                // Extract flowcell ID from MinKNOW run dir name
                // e.g. 20251021_1541_X1_FAZ84459_9acd0eaa -> FAZ84459
                def run_name = dir.parent.parent.name
                def parts = run_name.tokenize('_')
                def flowcell = parts.size() >= 4 ? parts[3] : run_name
                def sample_id = "${flowcell}_${barcode}"
                file("${dir}/*.fastq.gz").collect { fq -> [sample_id, fq] }
            }
            .groupTuple()
            .map { sample_id, fastqs -> [[id: sample_id], fastqs] }
            .filter { meta, fastqs -> fastqs.size() > 0 }

        CONCAT_READS(ch_barcode_raw)
        ch_reads = CONCAT_READS.out.reads
            .filter { meta, fastq -> fastq.size() > 1024 }  // skip empty/tiny outputs
    } else if (flat_fastqs) {
        // Flat directory of FASTQ files — existing behavior
        log.info "Detected flat FASTQ directory: ${flat_fastqs.size()} files"
        ch_reads = Channel.fromPath("${params.input}/*.fastq.gz")
            .map { fastq ->
                def name = fastq.baseName.replace('.fastq', '')
                [[id: name], fastq]
            }
    } else {
        error "ERROR: No FASTQ files found in ${params.input}. Expected either *.fastq.gz or fastq_pass/barcode*/ structure.\nRun with --help for usage."
    }

    // 2. Co-assembly: fan-in all reads into one Flye assembly
    ch_all_reads = ch_reads.map { meta, fastq -> fastq }.collect()
    ASSEMBLY_FLYE(ch_all_reads)

    // 2b. Tetranucleotide frequencies from assembly
    CALCULATE_TNF(ASSEMBLY_FLYE.out.assembly)

    // 2c. Prokka gene annotation on co-assembly (optional)
    if (params.run_prokka) {
        PROKKA_ANNOTATE(ASSEMBLY_FLYE.out.assembly)
    }

    // 2c2. Kaiju protein-level taxonomy (requires Prokka .faa + .gff)
    if (params.run_kaiju && params.kaiju_db && params.run_prokka) {
        KAIJU_CLASSIFY(PROKKA_ANNOTATE.out.proteins, PROKKA_ANNOTATE.out.gff)
    }

    // 2d. Mobile genetic element detection (geNomad + CheckV)
    if (params.run_genomad && params.genomad_db) {
        GENOMAD_CLASSIFY(ASSEMBLY_FLYE.out.assembly)

        if (params.run_checkv && params.checkv_db) {
            CHECKV_QUALITY(GENOMAD_CLASSIFY.out.virus_fasta)
        }
    }

    // 2e. Integron detection (IntegronFinder)
    if (params.run_integronfinder) {
        INTEGRONFINDER(ASSEMBLY_FLYE.out.assembly)
    }

    // 2f. Genomic island detection (IslandPath-DIMOB, requires Prokka .gbk)
    if (params.run_islandpath && params.run_prokka) {
        ISLANDPATH_DIMOB(PROKKA_ANNOTATE.out.gbk)
    }

    // 2g. Secretion system + conjugation detection (MacSyFinder, requires Prokka .faa)
    if (params.run_macsyfinder && params.run_prokka && params.macsyfinder_models) {
        MACSYFINDER(PROKKA_ANNOTATE.out.proteins)
    }

    // 2h. Anti-phage defense system detection (DefenseFinder, requires Prokka .faa)
    if (params.run_defensefinder && params.run_prokka) {
        DEFENSEFINDER(PROKKA_ANNOTATE.out.proteins)
    }

    // 3. Map each sample back to assembly: fan-out
    ch_map_input = ch_reads.combine(ASSEMBLY_FLYE.out.assembly)
    MAP_READS(ch_map_input)

    // 4. Calculate depths from all BAMs: fan-in
    // Collect BAMs and BAIs together so they're staged in the same directory
    ch_bam_files = MAP_READS.out.bam
        .flatMap { meta, bam, bai -> [bam, bai] }
        .collect()
    CALCULATE_DEPTHS(ch_bam_files, ASSEMBLY_FLYE.out.assembly)

    // 5. Binners run in parallel; each emits [label, contig_bins.tsv]
    // New binners can be added by appending to ch_binner_results
    ch_binner_results = Channel.empty()

    // SemiBin2 (optional, BAM-based)
    if (params.run_semibin) {
        BIN_SEMIBIN2(ASSEMBLY_FLYE.out.assembly, ch_bam_files)
        ch_binner_results = ch_binner_results.mix(
            BIN_SEMIBIN2.out.bins.map { ['semibin', it] }
        )
    }

    // MetaBAT2 (depth-based)
    BIN_METABAT2(ASSEMBLY_FLYE.out.assembly, CALCULATE_DEPTHS.out.jgi_depth)
    ch_binner_results = ch_binner_results.mix(
        BIN_METABAT2.out.bins.map { ['metabat', it] }
    )

    // MaxBin2 (optional, depth-based)
    if (params.run_maxbin) {
        BIN_MAXBIN2(ASSEMBLY_FLYE.out.assembly, CALCULATE_DEPTHS.out.jgi_depth)
        ch_binner_results = ch_binner_results.mix(
            BIN_MAXBIN2.out.bins.map { ['maxbin', it] }
        )
    }

    // LorBin (optional, BAM-based — deep learning binner for long reads)
    if (params.run_lorbin) {
        BIN_LORBIN(ASSEMBLY_FLYE.out.assembly, ch_bam_files)
        ch_binner_results = ch_binner_results.mix(
            BIN_LORBIN.out.bins.map { ['lorbin', it] }
        )
    }

    // COMEBin (optional, BAM-based — contrastive multi-view deep learning binner)
    if (params.run_comebin) {
        BIN_COMEBIN(ASSEMBLY_FLYE.out.assembly, ch_bam_files)
        ch_binner_results = ch_binner_results.mix(
            BIN_COMEBIN.out.bins.map { ['comebin', it] }
        )
    }

    // 6. DAS Tool consensus -- collects all binner outputs dynamically
    ch_bin_labels = ch_binner_results.collect { it[0] }
    ch_bin_files  = ch_binner_results.collect { it[1] }

    DASTOOL_CONSENSUS(
        ASSEMBLY_FLYE.out.assembly,
        ch_bin_files,
        ch_bin_labels
    )

    // 7. Quality assessment with CheckM2 (optional — requires database path)
    if (params.checkm2_db) {
        ch_all_bins = BIN_METABAT2.out.fastas
        if (params.run_semibin) {
            ch_all_bins = ch_all_bins.mix(BIN_SEMIBIN2.out.fastas)
        }
        if (params.run_maxbin) {
            ch_all_bins = ch_all_bins.mix(BIN_MAXBIN2.out.fastas)
        }
        if (params.run_lorbin) {
            ch_all_bins = ch_all_bins.mix(BIN_LORBIN.out.fastas)
        }
        if (params.run_comebin) {
            ch_all_bins = ch_all_bins.mix(BIN_COMEBIN.out.fastas)
        }
        ch_all_bins = ch_all_bins
            .mix(DASTOOL_CONSENSUS.out.bins)
            .collect()

        CHECKM2(ch_all_bins)
    }

    // 8. NCLB bin refinement (optional — requires --nclb_dir + LLM server)
    if (params.run_nclb && params.nclb_dir) {
        // Build a "ready" signal: collect DAS Tool output (+ CheckM2 if available)
        // This ensures all upstream results are published before NCLB reads them
        ch_nclb_ready = DASTOOL_CONSENSUS.out.contig2bin
            .mix(CALCULATE_TNF.out.tnf)
            .mix(CALCULATE_DEPTHS.out.jgi_depth)
            .collect()

        NCLB_GATHER(ch_nclb_ready)
        NCLB_CONVERSE(NCLB_GATHER.out.gathering)
        NCLB_ELDERS(NCLB_GATHER.out.gathering)
        NCLB_INTEGRATE(NCLB_CONVERSE.out.proposals)
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
