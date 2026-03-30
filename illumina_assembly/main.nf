#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Illumina Assembly Pipeline - Nextflow DSL2
// ============================================================================
//
// Illumina paired-end metagenomic assembly pipeline.
// Processes reads through BBTools QC (optical dedup, tile filtering, adapter
// trimming, artifact/PhiX removal, human decontamination), three-phase error
// correction, normalization, read merging, four parallel assemblers (Tadpole,
// Megahit, SPAdes, metaSPAdes), cascade deduplication, mapping, and depth
// calculation.
//
// Produces assembly FASTA + depths.txt + BAMs that can be fed into
// mag_analysis for downstream binning, annotation, taxonomy, etc.
//
// Supports per-sample (default) and co-assembly (--coassembly) modes.
//
// Usage:
//   nextflow run main.nf --input /path/to/reads --human_ref /path/to/ref -resume
//
// ============================================================================

// ============================================================================
// Help message
// ============================================================================

def helpMessage() {
    log.info """
    =========================================
     Illumina Assembly Pipeline
     https://github.com/rec3141/danaSeq
    =========================================

    Usage:
      nextflow run main.nf --input /path/to/reads --human_ref /path/to/ref [options] -resume

    Required:
      --input DIR        Input directory containing paired-end *_R1_*.fastq.gz reads
      --outdir DIR       Output directory [default: results]

    Caching:
      --store_dir DIR    Persistent cache directory (storeDir); completed processes are
                         skipped across runs even after work/ cleanup. Off by default.

    Mode:
      --coassembly       Co-assemble all samples together (default: per-sample assembly)

    Preprocessing:
      --min_readlen N    Minimum read length after trimming [default: 70]
      --run_remove_human Remove human reads (removehuman.sh) [default: true]
      --human_ref PATH   Path to BBTools human reference index
      --run_fastqc       Run FastQC on final preprocessed reads [default: true]

    Normalization:
      --run_normalize    Enable bbnorm coverage normalization [default: true]

    Assembly:
      --run_tadpole      Run Tadpole assembler [default: true]
      --run_megahit      Run Megahit assembler [default: true]
      --run_spades       Run SPAdes assembler [default: true]
      --run_metaspades   Run metaSPAdes assembler [default: true]

    Deduplication:
      --dedupe_identity N  Final deduplication identity threshold [default: 98]
      --min_contig_len N   Minimum contig length after deduplication [default: 500]

    Resources:
      --assembly_cpus N     CPUs for assembly/preprocessing [default: 24]
      --assembly_memory STR Memory for assembly [default: '250 GB']

    Output:
      results/assembly/<sample>/<sample>.dedupe.fasta   Final assembly
      results/mapping/<sample>/<sample>.depths.txt      Depth table (MetaBAT2 format)
      results/mapping/<sample>/<sample>.sorted.bam      Alignments

    These outputs can be passed to mag_analysis:
      mag_analysis/run-mag-analysis.sh \\
          --assembly results/assembly/<sample>/<sample>.dedupe.fasta \\
          --depths results/mapping/<sample>/<sample>.depths.txt \\
          --bam_dir results/mapping/<sample>/
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
    log.error "ERROR: --input is required. Provide path to directory containing *_R1_*.fastq.gz files. Run with --help for usage."
    System.exit(1)
}

if (params.run_remove_human && !params.human_ref) {
    log.error "ERROR: --run_remove_human requires --human_ref (path to BBTools human reference index directory)."
    System.exit(1)
}

// ============================================================================
// Import modules
// ============================================================================

include { CLUMPIFY }               from './modules/preprocess'
include { FILTER_BY_TILE }         from './modules/preprocess'
include { BBDUK_TRIM }             from './modules/preprocess'
include { BBDUK_FILTER }           from './modules/preprocess'
include { REMOVE_HUMAN }           from './modules/preprocess'
include { FASTQC }                 from './modules/preprocess'
include { ERROR_CORRECT_ECCO }     from './modules/error_correct'
include { ERROR_CORRECT_ECC }      from './modules/error_correct'
include { ERROR_CORRECT_TADPOLE }  from './modules/error_correct'
include { NORMALIZE_READS }        from './modules/normalize'
include { MERGE_READS }            from './modules/merge_reads'
include { QUALITY_TRIM }           from './modules/merge_reads'
include { ASSEMBLE_TADPOLE }       from './modules/assembly'
include { ASSEMBLE_MEGAHIT }       from './modules/assembly'
include { ASSEMBLE_SPADES }        from './modules/assembly'
include { ASSEMBLE_METASPADES }    from './modules/assembly'
include { DEDUPE_ASSEMBLIES }      from './modules/dedupe'
include { MAP_READS_BBMAP }        from './modules/mapping'
include { CALCULATE_DEPTHS }       from './modules/mapping'

// ============================================================================
// Main workflow
// ============================================================================

workflow {

    // 1. Discover input reads — scan for *_R1_*.fastq.gz and auto-derive R2
    def input_dir = file(params.input)
    if (!input_dir.isDirectory()) {
        error "ERROR: --input directory does not exist: ${params.input}\nRun with --help for usage."
    }

    def r1_files = file("${params.input}/*_R1_*.fastq.gz")
    if (!r1_files) {
        error "ERROR: No *_R1_*.fastq.gz files found in ${params.input}.\nRun with --help for usage."
    }

    log.info "Detected ${r1_files.size()} paired-end samples"
    log.info "Mode: ${params.coassembly ? 'co-assembly' : 'per-sample assembly'}"

    // Create channel of [meta, r1, r2] tuples
    ch_raw_reads = Channel.fromFilePairs("${params.input}/*_R{1,2}_*.fastq.gz", flat: true)
        .map { sample_id, r1, r2 ->
            def prefix = r1.name.toString().split('_')[0]
            [[id: prefix], r1, r2]
        }

    // ======================================================================
    // 2. Preprocessing (per-sample)
    // ======================================================================

    CLUMPIFY(ch_raw_reads)
    FILTER_BY_TILE(CLUMPIFY.out.reads)
    BBDUK_TRIM(FILTER_BY_TILE.out.reads)
    BBDUK_FILTER(BBDUK_TRIM.out.reads)

    if (params.run_remove_human) {
        REMOVE_HUMAN(BBDUK_FILTER.out.reads)
        ch_filtered = REMOVE_HUMAN.out.reads
    } else {
        ch_filtered = BBDUK_FILTER.out.reads
    }

    if (params.run_fastqc) {
        FASTQC(ch_filtered)
    }

    // ======================================================================
    // 3. Error correction (per-sample)
    // ======================================================================

    ERROR_CORRECT_ECCO(ch_filtered)
    ERROR_CORRECT_ECC(ERROR_CORRECT_ECCO.out.reads)
    ERROR_CORRECT_TADPOLE(ERROR_CORRECT_ECC.out.reads)

    ch_corrected = ERROR_CORRECT_TADPOLE.out.reads

    // ======================================================================
    // 4. Normalization (optional, per-sample)
    // ======================================================================

    if (params.run_normalize) {
        NORMALIZE_READS(ch_corrected)
        ch_for_merge = NORMALIZE_READS.out.reads
        ch_normalized = NORMALIZE_READS.out.reads
    } else {
        ch_for_merge = ch_corrected
        ch_normalized = ch_corrected
    }

    // ======================================================================
    // 5. Read merging (per-sample)
    // ======================================================================

    MERGE_READS(ch_for_merge)
    QUALITY_TRIM(MERGE_READS.out.unmerged)

    ch_merged   = MERGE_READS.out.merged
    ch_qtrimmed = QUALITY_TRIM.out.reads

    // ======================================================================
    // 6. Assembly
    // ======================================================================

    if (params.coassembly) {
        // Co-assembly mode: combine all samples' reads into one assembly
        ch_all_merged   = ch_merged.map { meta, reads -> reads }.collect()
        ch_all_qtrimmed = ch_qtrimmed.map { meta, reads -> reads }.collect()
        ch_all_normalized = ch_normalized.map { meta, reads -> reads }.collect()

        def coasm_meta = [id: 'coassembly']

        ch_coasm_merged = ch_all_merged.map { files -> [coasm_meta, files] }
        ch_coasm_qtrimmed = ch_all_qtrimmed.map { files -> [coasm_meta, files] }
        ch_coasm_normalized = ch_all_normalized.map { files -> [coasm_meta, files] }

        ch_assembly_contigs = Channel.empty()

        if (params.run_tadpole) {
            ASSEMBLE_TADPOLE(ch_coasm_merged.combine(ch_coasm_qtrimmed.map { m, f -> f }))
            ch_assembly_contigs = ch_assembly_contigs.mix(ASSEMBLE_TADPOLE.out.contigs)
        }
        if (params.run_megahit) {
            ASSEMBLE_MEGAHIT(ch_coasm_merged.combine(ch_coasm_qtrimmed.map { m, f -> f }))
            ch_assembly_contigs = ch_assembly_contigs.mix(ASSEMBLE_MEGAHIT.out.contigs)
        }
        if (params.run_spades) {
            ASSEMBLE_SPADES(ch_coasm_merged.combine(ch_coasm_qtrimmed.map { m, f -> f }))
            ch_assembly_contigs = ch_assembly_contigs.mix(ASSEMBLE_SPADES.out.contigs)
        }
        if (params.run_metaspades) {
            ASSEMBLE_METASPADES(ch_coasm_normalized)
            ch_assembly_contigs = ch_assembly_contigs.mix(ASSEMBLE_METASPADES.out.contigs)
        }

        ch_for_dedupe = ch_assembly_contigs
            .map { meta, contigs -> contigs }
            .collect()
            .map { contigs_list -> [coasm_meta, contigs_list] }

        DEDUPE_ASSEMBLIES(ch_for_dedupe)

        // Map each sample back to the co-assembly
        ch_map_input = ch_filtered.combine(DEDUPE_ASSEMBLIES.out.assembly.map { meta, asm -> asm })
            .map { meta, reads, asm -> [meta, reads, asm] }

        MAP_READS_BBMAP(ch_map_input)

        ch_all_bams = MAP_READS_BBMAP.out.bam
            .map { meta, bam, bai -> bam }
            .collect()

        CALCULATE_DEPTHS(
            DEDUPE_ASSEMBLIES.out.assembly
                .combine(ch_all_bams)
                .map { meta, asm, bams -> [meta, bams, asm] }
        )

    } else {
        // Per-sample assembly mode (default)

        ch_assembly_input = ch_merged.join(ch_qtrimmed)
            .map { meta, merged, qtrimmed -> [meta, merged, qtrimmed] }

        ch_assembly_contigs = Channel.empty()

        if (params.run_tadpole) {
            ASSEMBLE_TADPOLE(ch_assembly_input)
            ch_assembly_contigs = ch_assembly_contigs.mix(ASSEMBLE_TADPOLE.out.contigs)
        }
        if (params.run_megahit) {
            ASSEMBLE_MEGAHIT(ch_assembly_input)
            ch_assembly_contigs = ch_assembly_contigs.mix(ASSEMBLE_MEGAHIT.out.contigs)
        }
        if (params.run_spades) {
            ASSEMBLE_SPADES(ch_assembly_input)
            ch_assembly_contigs = ch_assembly_contigs.mix(ASSEMBLE_SPADES.out.contigs)
        }
        if (params.run_metaspades) {
            ASSEMBLE_METASPADES(ch_normalized)
            ch_assembly_contigs = ch_assembly_contigs.mix(ASSEMBLE_METASPADES.out.contigs)
        }

        ch_for_dedupe = ch_assembly_contigs
            .groupTuple()
            .map { meta, contigs_list -> [meta, contigs_list.flatten()] }

        DEDUPE_ASSEMBLIES(ch_for_dedupe)

        ch_map_input = ch_filtered.join(DEDUPE_ASSEMBLIES.out.assembly)
            .map { meta, reads, asm -> [meta, reads, asm] }

        MAP_READS_BBMAP(ch_map_input)

        ch_depth_input = MAP_READS_BBMAP.out.bam
            .join(DEDUPE_ASSEMBLIES.out.assembly)
            .map { meta, bam, bai, asm -> [meta, [bam], asm] }

        CALCULATE_DEPTHS(ch_depth_input)
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

        Next step: run mag_analysis on these outputs:
          --assembly <outdir>/assembly/<sample>/<sample>.dedupe.fasta
          --depths   <outdir>/mapping/<sample>/<sample>.depths.txt
          --bam_dir  <outdir>/mapping/<sample>/
        """.stripIndent()

    println msg

    if (!workflow.success) {
        println "[WARNING] Pipeline completed with errors. Check .nextflow.log for details."
    }
}

workflow.onError {
    println "[ERROR] Pipeline failed: ${workflow.errorMessage}"
}
