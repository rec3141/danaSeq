#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Illumina RNA Pipeline - Nextflow DSL2
// ============================================================================
//
// Paired-end Illumina RNA-seq pipeline tailored for metatranscriptomic data
// mapped against references from other danaSeq pipelines (assemblies from
// illumina_assembly/nanopore_assembly, or annotated MAGs from mag_analysis).
//
// Stages:
//   1. QA/QC      — BBTools (clumpify, filterbytile, bbduk_trim/filter, removehuman) + FastQC
//   2. rRNA       — SortMeRNA (optional, on by default)
//   3. Mapping    — BBmap to each reference FASTA in --references
//   4. Quantify   — featureCounts when a <ref>.gff is present, samtools idxstats always
//   5. Summarize  — merged gene_counts.tsv per reference + viz JSONs
//
// Usage:
//   nextflow run main.nf \\
//       --input /path/to/reads \\
//       --references /path/to/refs_dir \\
//       --human_ref /path/to/bbtools_human_index \\
//       --sortmerna_refs /path/to/sortmerna_fastas -resume
//
// ============================================================================

def helpMessage() {
    log.info """
    =========================================
     Illumina RNA Pipeline
     https://github.com/rec3141/danaSeq
    =========================================

    Usage:
      nextflow run main.nf --input DIR --references DIR [options] -resume

    Required:
      --input DIR        Directory of paired-end *_R1_*.fastq.gz
      --references DIR   Directory of reference <name>.fasta (+ optional <name>.gff)
      --outdir DIR       Output directory [default: results]

    Caching:
      --store_dir DIR    Persistent cache (storeDir); skips completed processes across runs

    Preprocessing:
      --min_readlen N        Minimum read length after trimming [default: 50]
      --run_remove_human     Remove human reads via removehuman.sh [default: true]
      --human_ref PATH       BBTools human reference index dir
      --run_fastqc           Run FastQC on preprocessed reads [default: true]

    rRNA removal:
      --run_remove_rrna      Strip rRNA via SortMeRNA [default: true]
      --sortmerna_refs DIR   Directory of rRNA FASTAs for SortMeRNA

    Mapping:
      --min_identity N       bbmap minid [default: 90]

    Quantification:
      --strandedness STR     'unstranded' | 'forward' | 'reverse' [default: unstranded]
      --feature_type STR     GFF feature (column 3) for featureCounts [default: CDS]
      --attr_type STR        GFF attribute for grouping [default: locus_tag]

    Resources:
      --rna_cpus N           CPUs per high-resource task [default: 16]
      --rna_memory STR       Memory per high-resource task [default: '64 GB']

    Output:
      results/preprocess/<sample>/         Cleaned reads + FastQC + SortMeRNA log
      results/mapping/<sample>/<ref>/      sorted.bam, idxstats, flagstat, covstats
      results/quantify/<sample>/<ref>/     featureCounts per-sample counts
      results/expression/<ref>/            merged <ref>.gene_counts.tsv (genes × samples)
      results/viz/                         JSONs for the Svelte viz
    """.stripIndent()
}

def validateParams() {
    if (params.help) {
        helpMessage()
        System.exit(0)
    }
    if (!params.input) {
        log.error "ERROR: --input is required. Run with --help for usage."
        System.exit(1)
    }
    if (!params.references) {
        log.error "ERROR: --references is required. Run with --help for usage."
        System.exit(1)
    }
    if (params.run_remove_human && !params.human_ref) {
        log.error "ERROR: --run_remove_human requires --human_ref."
        System.exit(1)
    }
}

// ============================================================================
// Import modules
// ============================================================================

include { CLUMPIFY }         from './modules/preprocess'
include { FILTER_BY_TILE }   from './modules/preprocess'
include { BBDUK_TRIM }       from './modules/preprocess'
include { BBDUK_FILTER }     from './modules/preprocess'
include { REMOVE_HUMAN }     from './modules/preprocess'
include { FASTQC }           from './modules/preprocess'
include { COUNT_READS as COUNT_READS_CLUMPED  } from './modules/preprocess'
include { COUNT_READS as COUNT_READS_TRIMMED  } from './modules/preprocess'
include { COUNT_READS as COUNT_READS_FILTERED } from './modules/preprocess'
include { COUNT_READS as COUNT_READS_NOHUMAN  } from './modules/preprocess'
include { COUNT_READS as COUNT_READS_NORRNA   } from './modules/preprocess'
include { REMOVE_RRNA }      from './modules/rrna'
include { BBMAP_INDEX }      from './modules/alignment'
include { MAP_READS_BBMAP }  from './modules/alignment'
include { FEATURECOUNTS }    from './modules/quantify'
include { MERGE_GENE_COUNTS } from './modules/summarize'
include { VIZ_PREPROCESS }    from './modules/summarize'

// ============================================================================
// Main workflow
// ============================================================================

workflow {

    main:

    validateParams()

    // ----- 1. Discover input reads -----
    // Supports both Illumina ("*_R1_*.fastq.gz") and BGI/MGI ("*_1.fq.gz") naming.
    def input_dir = file(params.input)
    if (!input_dir.isDirectory()) {
        error "ERROR: --input directory does not exist: ${params.input}"
    }

    def r1_files = file("${params.input}/*_R1_*.fastq.gz") +
                   file("${params.input}/*_1.fq.gz") +
                   file("${params.input}/*_1.fastq.gz")
    if (!r1_files) {
        error "ERROR: No paired-end reads found in ${params.input}.\n" +
              "  Looked for *_R1_*.fastq.gz (Illumina) and *_{1,2}.fq.gz/*.fastq.gz (BGI)."
    }
    log.info "Detected ${r1_files.size()} R1 file(s)"

    ch_raw_reads = Channel.fromPath([
            "${params.input}/*_R1_*.fastq.gz",
            "${params.input}/*_1.fq.gz",
            "${params.input}/*_1.fastq.gz",
        ])
        .map { r1 ->
            def name = r1.name.toString()
            def sample_id
            def r2_name
            if (name =~ /_R1_/) {
                sample_id = name.split('_')[0]
                r2_name   = name.replaceFirst('_R1_', '_R2_')
            } else if (name =~ /_1\.fq\.gz$/) {
                sample_id = name.replaceFirst(/_1\.fq\.gz$/, '')
                r2_name   = name.replaceFirst(/_1\.fq\.gz$/, '_2.fq.gz')
            } else {
                sample_id = name.replaceFirst(/_1\.fastq\.gz$/, '')
                r2_name   = name.replaceFirst(/_1\.fastq\.gz$/, '_2.fastq.gz')
            }
            def r2 = file("${params.input}/${r2_name}")
            [[id: sample_id], r1, r2]
        }
        .filter { meta, r1, r2 ->
            if (!r2.exists()) {
                log.warn "[WARNING] Missing R2 for ${meta.id} (expected ${r2.name}); skipping."
                return false
            }
            return true
        }

    // ----- 2. Discover references -----
    def ref_dir = file(params.references)
    if (!ref_dir.isDirectory()) {
        error "ERROR: --references directory does not exist: ${params.references}"
    }
    def ref_fastas = file("${params.references}/*.{fasta,fa,fna}")
    if (!ref_fastas) {
        error "ERROR: No *.fasta/*.fa/*.fna files found in ${params.references}."
    }
    log.info "Detected ${ref_fastas.size()} reference FASTAs"

    ch_references = Channel.fromPath("${params.references}/*.{fasta,fa,fna}")
        .map { fa ->
            def name = fa.baseName
            def gff  = file("${params.references}/${name}.gff")
            def has_gff = gff.exists()
            [[name: name, has_gff: has_gff], fa]
        }

    ch_ref_gffs = Channel.fromPath("${params.references}/*.{fasta,fa,fna}")
        .map { fa ->
            def name = fa.baseName
            def gff  = file("${params.references}/${name}.gff")
            [[name: name], gff.exists() ? gff : file("${projectDir}/assets/empty.gff", checkIfExists: false)]
        }

    // ----- 3. Preprocessing -----
    // Accumulate per-stage read counts into a single channel — much easier
    // than chasing conditionally-invoked process outputs in section 8.
    ch_readcounts = Channel.empty()

    CLUMPIFY(ch_raw_reads)
    COUNT_READS_CLUMPED(CLUMPIFY.out.reads.map { meta, r -> [meta, 'clumped', r] })
    ch_readcounts = ch_readcounts.mix(COUNT_READS_CLUMPED.out.count)

    FILTER_BY_TILE(CLUMPIFY.out.reads)
    BBDUK_TRIM(FILTER_BY_TILE.out.reads)
    COUNT_READS_TRIMMED(BBDUK_TRIM.out.reads.map { meta, r -> [meta, 'trimmed', r] })
    ch_readcounts = ch_readcounts.mix(COUNT_READS_TRIMMED.out.count)

    BBDUK_FILTER(BBDUK_TRIM.out.reads)
    COUNT_READS_FILTERED(BBDUK_FILTER.out.reads.map { meta, r -> [meta, 'filtered', r] })
    ch_readcounts = ch_readcounts.mix(COUNT_READS_FILTERED.out.count)

    if (params.run_remove_human) {
        REMOVE_HUMAN(BBDUK_FILTER.out.reads)
        ch_filtered = REMOVE_HUMAN.out.reads
        COUNT_READS_NOHUMAN(ch_filtered.map { meta, r -> [meta, 'nohuman', r] })
        ch_readcounts = ch_readcounts.mix(COUNT_READS_NOHUMAN.out.count)
    } else {
        ch_filtered = BBDUK_FILTER.out.reads
    }

    if (params.run_fastqc) {
        FASTQC(ch_filtered)
    }

    // ----- 4. rRNA removal -----
    if (params.run_remove_rrna) {
        REMOVE_RRNA(ch_filtered)
        ch_clean = REMOVE_RRNA.out.reads
        COUNT_READS_NORRNA(ch_clean.map { meta, r -> [meta, 'norrna', r] })
        ch_readcounts = ch_readcounts.mix(COUNT_READS_NORRNA.out.count)
    } else {
        ch_clean = ch_filtered
    }

    // ----- 5. Reference indexing + mapping -----
    BBMAP_INDEX(ch_references)

    // Cartesian product: every sample × every reference
    ch_map_input = ch_clean.combine(BBMAP_INDEX.out.index)
        .map { meta, reads, ref, ref_fa, ref_idx ->
            [meta, reads, ref, ref_fa, ref_idx]
        }

    MAP_READS_BBMAP(ch_map_input)

    // ----- 6. Quantification (featureCounts only when a GFF exists) -----
    ch_quant_input = MAP_READS_BBMAP.out.bam
        .map { meta, ref, bam, bai -> [[meta.id, ref.name], meta, ref, bam, bai] }
        .join(
            ch_ref_gffs.map { ref, gff -> [ref.name, gff] }
                .combine(MAP_READS_BBMAP.out.bam.map { meta, ref, bam, bai -> [meta.id, ref.name] })
                .map { ref_name, gff, sample_id, rn -> [[sample_id, ref_name], gff] }
        )
        .map { key, meta, ref, bam, bai, gff -> [meta, ref, bam, bai, gff] }
        .filter { meta, ref, bam, bai, gff -> ref.has_gff }

    FEATURECOUNTS(ch_quant_input)

    // ----- 7. Per-reference gene-count matrix -----
    ch_counts_by_ref = FEATURECOUNTS.out.counts
        .map { meta, ref, tsv -> [ref, tsv] }
        .groupTuple()

    MERGE_GENE_COUNTS(ch_counts_by_ref)

    // ----- 8. Viz JSONs -----
    ch_idxstats   = MAP_READS_BBMAP.out.idxstats.map { meta, ref, f -> f }.collect()
    ch_flagstat   = MAP_READS_BBMAP.out.flagstat.map { meta, ref, f -> f }.collect()
    ch_covstats   = MAP_READS_BBMAP.out.covstats.map { meta, ref, f -> f }.collect()
    ch_readcounts_collected = ch_readcounts
        .map { meta, stage, f -> f }
        .collect()
    ch_genecounts = MERGE_GENE_COUNTS.out.counts.map { ref, tsv -> tsv }.collect()

    VIZ_PREPROCESS(
        ch_idxstats,
        ch_flagstat,
        ch_covstats,
        ch_readcounts_collected,
        ch_genecounts
    )

    workflow.onComplete = {
        def msg = """\
            Pipeline completed at : ${workflow.complete}
            Duration              : ${workflow.duration}
            Success               : ${workflow.success}
            Exit status           : ${workflow.exitStatus}
            Output directory      : ${params.outdir}

            Viz JSONs              : ${params.outdir}/viz/
            Per-reference counts   : ${params.outdir}/expression/<ref>/<ref>.gene_counts.tsv
            """.stripIndent()
        println msg
        if (!workflow.success) {
            println "[WARNING] Pipeline completed with errors. Check .nextflow.log for details."
        }
    }

    workflow.onError = {
        println "[ERROR] Pipeline failed: ${workflow.errorMessage}"
    }
}
