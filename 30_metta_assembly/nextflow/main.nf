#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// METTA Assembly Pipeline - Nextflow DSL2
// ============================================================================
//
// Illumina paired-end metagenomic assembly pipeline.
// Processes reads through BBTools QC, three-phase error correction,
// normalization, read merging, four parallel assemblers (Tadpole, Megahit,
// SPAdes, metaSPAdes), cascade deduplication, mapping, depth calculation,
// and MetaBAT2 binning.
//
// Supports per-sample (default) and co-assembly (--coassembly) modes.
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
     METTA Assembly Pipeline
     https://github.com/rec3141/danaSeq
    =========================================

    Usage:
      nextflow run main.nf --input /path/to/reads [options] -resume

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

    Normalization:
      --run_normalize    Enable bbnorm coverage normalization [default: true]

    Assembly:
      --run_tadpole      Run Tadpole assembler [default: true]
      --run_megahit      Run Megahit assembler [default: true]
      --run_spades       Run SPAdes assembler [default: true]
      --run_metaspades   Run metaSPAdes assembler [default: true]

    Deduplication:
      --dedupe_identity N  Final deduplication identity threshold [default: 98]

    Binning:
      --metabat_min_cls N  MetaBAT2 minimum cluster size [default: 2000]

    Resources:
      --assembly_cpus N     CPUs for assembly/preprocessing [default: 24]
      --assembly_memory STR Memory for assembly [default: '250 GB']

    SLURM:
      --slurm_account STR   SLURM --account [default: def-rec3141]
      --conda_path PATH     Path to conda/mamba bin/ dir for SLURM jobs (e.g. ~/scratch/miniforge3/bin)

    Examples:
      # Per-sample assembly (default)
      nextflow run main.nf --input /path/to/reads -resume

      # Co-assembly of all samples
      nextflow run main.nf --input /path/to/reads --coassembly -resume

      # Skip Megahit and SPAdes, just Tadpole + metaSPAdes
      nextflow run main.nf --input /path/to/reads \\
          --run_megahit false --run_spades false -resume

      # Using launcher script
      ./run-metta.sh --input /path/to/reads --outdir /path/to/output

      # SLURM profile
      nextflow run main.nf --input /path/to/reads -profile slurm -resume

    Input:
      --input must point to a directory containing paired-end Illumina reads
      with the naming convention *_R1_*.fastq.gz. R2 files are auto-detected
      by replacing _R1_ with _R2_. Sample ID is the first field before '_'.

    Output:
      results/
      \u251c\u2500\u2500 preprocess/<sample>/     QC'd reads
      \u251c\u2500\u2500 error_correct/<sample>/  Error-corrected reads
      \u251c\u2500\u2500 normalize/<sample>/      Normalized reads
      \u251c\u2500\u2500 merge/<sample>/          Merged + quality-trimmed reads
      \u251c\u2500\u2500 assembly/<sample>/       Per-assembler + deduplicated contigs
      \u251c\u2500\u2500 mapping/<sample>/        BAM + coverage stats
      \u251c\u2500\u2500 binning/<sample>/        MetaBAT2 bins
      \u2514\u2500\u2500 pipeline_info/           Nextflow reports

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

// ============================================================================
// Import modules
// ============================================================================

include { CLUMPIFY }               from './modules/preprocess'
include { FILTER_BY_TILE }         from './modules/preprocess'
include { BBDUK_TRIM }             from './modules/preprocess'
include { BBDUK_FILTER }           from './modules/preprocess'
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
include { BIN_METABAT2 }           from './modules/binning'

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
            // Extract prefix (first field before '_') as sample ID
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

    // Save filtered reads for mapping later
    ch_filtered = BBDUK_FILTER.out.reads

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

        // Concatenate merged reads
        ch_coasm_merged = ch_all_merged.map { files ->
            [coasm_meta, files]
        }
        ch_coasm_qtrimmed = ch_all_qtrimmed.map { files ->
            [coasm_meta, files]
        }
        ch_coasm_normalized = ch_all_normalized.map { files ->
            [coasm_meta, files]
        }

        // Run assemblers on combined data
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

        // Collect all assemblies and deduplicate
        ch_for_dedupe = ch_assembly_contigs
            .map { meta, contigs -> contigs }
            .collect()
            .map { contigs_list -> [coasm_meta, contigs_list] }

        DEDUPE_ASSEMBLIES(ch_for_dedupe)

        // Map each sample back to the co-assembly
        ch_map_input = ch_filtered.combine(DEDUPE_ASSEMBLIES.out.assembly.map { meta, asm -> asm })
            .map { meta, reads, asm -> [meta, reads, asm] }

        MAP_READS_BBMAP(ch_map_input)

        // Collect all BAMs for depth calculation
        ch_all_bams = MAP_READS_BBMAP.out.bam
            .map { meta, bam, bai -> bam }
            .collect()

        CALCULATE_DEPTHS(
            DEDUPE_ASSEMBLIES.out.assembly
                .combine(ch_all_bams)
                .map { meta, asm, bams -> [meta, bams, asm] }
        )

        // Binning on co-assembly
        BIN_METABAT2(
            DEDUPE_ASSEMBLIES.out.assembly
                .combine(CALCULATE_DEPTHS.out.depths.map { meta, depths -> depths })
                .map { meta, asm, depths -> [meta, asm, depths] }
        )

    } else {
        // Per-sample assembly mode (default)

        // Prepare assembly inputs: join merged + qtrimmed per sample
        ch_assembly_input = ch_merged.join(ch_qtrimmed)
            .map { meta, merged, qtrimmed -> [meta, merged, qtrimmed] }

        // Run assemblers in parallel (per-sample)
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
            // metaSPAdes uses normalized reads (interleaved)
            ASSEMBLE_METASPADES(ch_normalized)
            ch_assembly_contigs = ch_assembly_contigs.mix(ASSEMBLE_METASPADES.out.contigs)
        }

        // Collect all assemblies per sample and deduplicate
        ch_for_dedupe = ch_assembly_contigs
            .groupTuple()
            .map { meta, contigs_list -> [meta, contigs_list.flatten()] }

        DEDUPE_ASSEMBLIES(ch_for_dedupe)

        // Map filtered reads back to each sample's assembly
        ch_map_input = ch_filtered.join(DEDUPE_ASSEMBLIES.out.assembly)
            .map { meta, reads, asm -> [meta, reads, asm] }

        MAP_READS_BBMAP(ch_map_input)

        // Depth calculation: per-sample (single BAM per assembly)
        ch_depth_input = MAP_READS_BBMAP.out.bam
            .join(DEDUPE_ASSEMBLIES.out.assembly)
            .map { meta, bam, bai, asm -> [meta, [bam], asm] }

        CALCULATE_DEPTHS(ch_depth_input)

        // Binning: per-sample
        ch_bin_input = DEDUPE_ASSEMBLIES.out.assembly
            .join(CALCULATE_DEPTHS.out.depths)
            .map { meta, asm, depths -> [meta, asm, depths] }

        BIN_METABAT2(ch_bin_input)
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
