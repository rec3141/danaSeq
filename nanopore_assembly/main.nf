#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Dana Nanopore Assembly Pipeline - Nextflow DSL2
// ============================================================================
//
// Nanopore-specific assembly pipeline: preprocessing, co-assembly (Flye,
// metaMDBG, or myloasm), read mapping, and depth calculation.
//
// Produces assembly.fasta + depths.txt + BAMs that can be fed into
// mag_analysis for downstream binning, annotation, taxonomy, etc.
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
     Dana Nanopore Assembly Pipeline
     https://github.com/rec3141/danaSeq
    =========================================

    Usage:
      nextflow run main.nf --input /path/to/reads [options] -resume

    Required:
      --input DIR        Input directory (*.fastq.gz or fastq_pass/barcode* structure)
      --outdir DIR       Output directory [default: results]

    Caching:
      --store_dir DIR    Persistent cache directory (storeDir); completed processes are
                         skipped across runs even after work/ cleanup. Off by default.

    Preprocessing:
      --run_remove_human Remove human reads via minimap2 [default: true]
      --human_ref PATH   Path to human reference FASTA or .mmi

    Assembly:
      --assembler STR    Assembler to use: 'flye', 'metamdbg', or 'myloasm' [default: flye]
      --min_overlap N    Flye --min-overlap [default: 1000]
      --polish           Enable Flye polishing iterations [default: true for flye]
      --dedupe           Enable BBDuk deduplication before assembly
      --filtlong_size N  Filtlong target bases (e.g. 40000000000); skip if not set

    Resources:
      --assembly_cpus N    CPUs for assembly [default: 16]
      --assembly_memory S  Memory for assembly [default: '60 GB']

    Output:
      results/assembly/assembly.fasta     Co-assembly
      results/assembly/tnf.tsv            Tetranucleotide frequencies
      results/mapping/depths.txt          CoverM depth table (MetaBAT2 format)
      results/mapping/*.sorted.bam        Per-sample alignments

    These outputs can be passed to mag_analysis:
      mag_analysis/run-mag-analysis.sh \\
          --assembly results/assembly/assembly.fasta \\
          --depths results/mapping/depths.txt \\
          --bam_dir results/mapping/
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

if (!(params.assembler in ['flye', 'metamdbg', 'myloasm'])) {
    log.error "ERROR: Invalid --assembler '${params.assembler}'. Choose from: flye, metamdbg, myloasm"
    System.exit(1)
}
log.info "Assembler: ${params.assembler}"

// Import modules
include { CONCAT_READS }        from './modules/preprocess'
include { PREPARE_READS }        from './modules/preprocess'
include { REMOVE_HUMAN }        from './modules/preprocess'
include { FLYE_ASSEMBLE }       from './modules/assembly'
include { FLYE_POLISH }         from './modules/assembly'
include { ASSEMBLY_METAMDBG }   from './modules/assembly'
include { ASSEMBLY_MYLOASM }    from './modules/assembly'
include { CALCULATE_TNF }       from './modules/assembly'
include { MAP_READS }           from './modules/mapping'
include { CALCULATE_DEPTHS }    from './modules/mapping'

// ============================================================================
// Main workflow
// ============================================================================

process WRITE_PROVENANCE {
    tag "provenance"
    label 'process_low'
    publishDir "${params.store_dir ?: params.outdir}/pipeline_info", mode: 'copy'

    output:
    path("versions.yml"), emit: versions

    script:
    """
    cat > versions.yml <<YAML
pipeline:
  name: ${workflow.manifest.name}
  version: ${workflow.manifest.version}
  revision: ${workflow.revision ?: 'unknown'}
  commitId: ${workflow.commitId ?: 'unknown'}
  sessionId: ${workflow.sessionId}

nextflow:
  version: ${nextflow.version}

params:
  assembler: ${params.assembler}
  polish: ${params.polish ?: 'auto'}
  dedupe: ${params.dedupe}
  filtlong_size: ${params.filtlong_size ?: 'none'}
  min_overlap: ${params.min_overlap}
  read_type: ${params.read_type}
  assembly_cpus: ${params.assembly_cpus}
  assembly_memory: ${params.assembly_memory}

tools:
YAML

    {
      echo "  flye: \$(flye --version 2>&1 || echo 'not found')"
      echo "  minimap2: \$(minimap2 --version 2>&1 || echo 'not found')"
      echo "  samtools: \$(samtools --version 2>&1 | head -1 | awk '{print \$2}' || echo 'not found')"
    } >> versions.yml 2>/dev/null || true
    """
}

workflow {

    WRITE_PROVENANCE()

    // 1. Discover input reads — auto-detect nanopore barcode vs flat directory
    def input_dir = file(params.input)
    if (!input_dir.isDirectory()) {
        error "ERROR: --input directory does not exist: ${params.input}\nRun with --help for usage."
    }

    def barcode_dirs = file("${params.input}/**/fastq_pass/barcode*", type: 'dir') +
                       file("${params.input}/fastq_pass/barcode*", type: 'dir')
    def flat_fastqs  = file("${params.input}/*.fastq.gz")

    if (barcode_dirs) {
        log.info "Detected nanopore barcode structure: ${barcode_dirs.size()} barcode directories"
        def all_pairs = []
        for (dir in barcode_dirs) {
            def barcode = dir.name
            def run_name = dir.parent.parent.name
            def parts = run_name.tokenize('_')
            def flowcell = parts.size() >= 4 ? parts[3] : run_name
            def sample_id = "${flowcell}_${barcode}"
            def fqs = file("${dir}/*.fastq.gz")
            if (fqs instanceof List) {
                for (fq in fqs) { all_pairs.add([sample_id, fq]) }
            } else if (fqs) {
                all_pairs.add([sample_id, fqs])
            }
        }
        log.info "Found ${all_pairs.size()} FASTQ files across ${barcode_dirs.size()} barcodes"
        ch_barcode_raw = Channel.from(all_pairs)
            .groupTuple()
            .map { sample_id, fastqs -> [[id: sample_id], fastqs] }
            .filter { meta, fastqs -> fastqs.size() > 0 }

        CONCAT_READS(ch_barcode_raw)
        ch_reads = CONCAT_READS.out.reads
            .filter { meta, fastq -> fastq.size() > 1024 }
    } else if (flat_fastqs) {
        log.info "Detected flat FASTQ directory: ${flat_fastqs.size()} files"
        ch_reads = Channel.fromPath("${params.input}/*.fastq.gz")
            .map { fastq ->
                def name = fastq.baseName.replace('.fastq', '')
                [[id: name], fastq]
            }
    } else {
        error "ERROR: No FASTQ files found in ${params.input}. Expected either *.fastq.gz or fastq_pass/barcode*/ structure.\nRun with --help for usage."
    }

    // 2. Concatenate + dedupe + optional filtlong -> single all_reads.fastq.gz
    ch_per_barcode = ch_reads.map { meta, fastq -> fastq }.collect()
    PREPARE_READS(ch_per_barcode)

    if (params.run_remove_human) {
        REMOVE_HUMAN(PREPARE_READS.out.reads)
        ch_asm_input = REMOVE_HUMAN.out.reads
    } else {
        ch_asm_input = PREPARE_READS.out.reads
    }

    if (params.assembler == 'flye') {
        FLYE_ASSEMBLE(ch_asm_input)
        ch_raw_assembly = FLYE_ASSEMBLE.out.assembly
        ch_asm_info     = FLYE_ASSEMBLE.out.info
        ch_asm_graph    = FLYE_ASSEMBLE.out.graph
    } else if (params.assembler == 'metamdbg') {
        ASSEMBLY_METAMDBG(ch_asm_input)
        ch_raw_assembly = ASSEMBLY_METAMDBG.out.assembly
        ch_asm_info     = ASSEMBLY_METAMDBG.out.info
        ch_asm_graph    = ASSEMBLY_METAMDBG.out.graph
    } else if (params.assembler == 'myloasm') {
        ASSEMBLY_MYLOASM(ch_asm_input)
        ch_raw_assembly = ASSEMBLY_MYLOASM.out.assembly
        ch_asm_info     = ASSEMBLY_MYLOASM.out.info
        ch_asm_graph    = ASSEMBLY_MYLOASM.out.graph
    }

    // Optional polishing (default: true for flye, false for others)
    def do_polish = (params.polish != null) ? params.polish : (params.assembler == 'flye')
    if (do_polish) {
        FLYE_POLISH(ch_raw_assembly, ch_asm_info, ch_asm_graph, ch_asm_input)
        ch_assembly  = FLYE_POLISH.out.assembly
        ch_asm_info  = FLYE_POLISH.out.info
        ch_asm_graph = FLYE_POLISH.out.graph
    } else {
        ch_assembly = ch_raw_assembly
    }

    // Tetranucleotide frequencies from assembly
    CALCULATE_TNF(ch_assembly)

    // 3. Map each sample back to assembly: fan-out
    ch_map_input = ch_reads.combine(ch_assembly)
    MAP_READS(ch_map_input)

    // 4. Calculate depths from all BAMs: fan-in
    ch_bam_files = MAP_READS.out.bam
        .flatMap { meta, bam, bai -> [bam, bai] }
        .collect()
    CALCULATE_DEPTHS(ch_bam_files, ch_assembly)
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

        Next step: run mag_analysis on these outputs:
          --assembly ${params.outdir}/assembly/assembly.fasta
          --depths   ${params.outdir}/mapping/depths.txt
          --bam_dir  ${params.outdir}/mapping/
        """.stripIndent()

    println msg

    if (!workflow.success) {
        println "[WARNING] Pipeline completed with errors. Check .nextflow.log for details."
    }
}

workflow.onError {
    println "[ERROR] Pipeline failed: ${workflow.errorMessage}"
}
