#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Illumina RNA Pipeline - Nextflow DSL2
// ============================================================================
//
// Paired-end Illumina/BGI RNA-seq pipeline for metatranscriptomes, mapped
// against references from other danaSeq pipelines or supplied directly.
//
// Stages:
//   1. QC         — one streaming pass per sample: bbduk adapter trim | bbduk
//                   artifact/phiX/entropy filter | bbduk rRNA split (SILVA)
//   2. Human      — strobealign vs a masked human FASTA, removehuman.sh stringency
//   3. Mapping    — strobealign (default) or BBmap to each reference
//   4. Quantify   — featureCounts on genes (from the reference GFF) and contigs
//   5. Summarize  — genes/contigs/bins/functions × samples matrices + viz JSONs
//
// References (pick one):
//   --from_mag DIR       mag_analysis outdir: assembly + Bakta GFF + bins,
//                        CheckM2, GTDB-Tk and KO/COG annotations
//   --from_assembly DIR  nanopore_assembly / illumina_assembly outdir (contigs)
//   --references DIR     <name>.fasta (+ <name>.gff, <name>.bins.tsv)
//
// ============================================================================

def helpMessage() {
    log.info """
    =========================================
     Illumina RNA Pipeline
     https://github.com/rec3141/danaSeq
    =========================================

    Usage:
      nextflow run main.nf --input DIR --from_mag DIR [options] -resume

    Reads (one of):
      --input DIR        Directory searched recursively for paired reads
                         (*_R1_*.fastq.gz Illumina, *_1.fq.gz / *_1.fastq.gz BGI)
      --samplesheet CSV  sample,r1,r2 (header required)

    References (one of):
      --from_mag DIR       mag_analysis output dir
      --from_assembly DIR  nanopore_assembly / illumina_assembly output dir
      --references DIR     <name>.fasta (+ optional <name>.gff, <name>.bins.tsv)
      --ref_name STR       Name for a --from_mag/--from_assembly reference
                           [default: the directory name]
      --bin_set STR        mag_analysis bin set for per-MAG tables [default: dastool]
      --contig_explorer DIR  Contig t-SNE in mag_analysis viz format (contig_explorer_*)
                           [default: <from_mag>/viz/data when present]

    Output:
      --outdir DIR       Output directory [default: results]
      --store_dir DIR    Persistent cache (storeDir); per-sample results are
                         written here instead of published to --outdir

    Batching (large cohorts):
      --max_read_pairs N Use only the first N read pairs of each sample [default: off]
      --chunk I/N        Process only the I-th of N deterministic sample chunks
                         (1-based); pair with --summarize false and a shared --store_dir
      --summarize BOOL   Build the merged matrices + viz JSONs [default: true]
      --prepare_only     Build the reference index and features, then stop

    QC:
      --min_readlen N        Minimum read length after trimming [default: 50]
      --run_filterbytile     filterbytile.sh (Illumina tiles only) [default: false]
      --run_remove_human     Remove human reads via removehuman.sh [default: true]
      --human_ref PATH       Human reference FASTA (masked hg19 etc.; indexed once)
      --run_fastqc           FastQC on the final reads [default: false]

    rRNA removal:
      --rrna_method STR      bbduk | sortmerna | none [default: bbduk]
      --rrna_refs PATH       rRNA FASTA file or dir of FASTAs (SILVA SSU/LSU; U→T handled)

    Mapping:
      --mapper STR           strobealign | bbmap [default: strobealign]
      --read_len N           Read length for the strobealign index [default: 150]
      --min_identity N       bbmap minid [default: 90]

    Quantification:
      --strandedness STR     unstranded | forward | reverse [default: unstranded]
      --feature_types STR    GFF features counted [default: CDS,rRNA,tRNA,tmRNA,ncRNA]
      --min_mapq N           featureCounts -Q [default: 0]

    Resources:
      --rna_cpus N           CPUs per high-resource task [default: 16]
      --rna_memory STR       Memory per high-resource task [default: '64 GB']
      --map_cpus N           CPUs per mapping task [default: --rna_cpus]
      --map_memory STR       Memory per mapping/index task [default: --rna_memory]
    """.stripIndent()
}

def validateParams() {
    if (params.help) {
        helpMessage()
        System.exit(0)
    }
    if (workflow.profile.tokenize(',').contains('slurm') && !params.slurm_account) {
        log.error "ERROR: -profile slurm requires --slurm_account <account>. There is no default; SLURM rejects jobs submitted without a valid account."
        System.exit(1)
    }
    if (!params.prepare_only && !params.input && !params.samplesheet) {
        log.error "ERROR: --input or --samplesheet is required. Run with --help for usage."
        System.exit(1)
    }
    def n_ref_sources = [params.references, params.from_mag, params.from_assembly].count { it }
    if (n_ref_sources != 1) {
        log.error "ERROR: give exactly one of --references, --from_mag, --from_assembly."
        System.exit(1)
    }
    if (params.run_remove_human && !params.human_ref) {
        log.error "ERROR: --run_remove_human requires --human_ref."
        System.exit(1)
    }
    if (!(params.rrna_method in ['bbduk', 'sortmerna', 'none'])) {
        log.error "ERROR: --rrna_method must be bbduk, sortmerna or none."
        System.exit(1)
    }
    if (params.rrna_method != 'none' && !params.rrna_refs) {
        log.error "ERROR: --rrna_method ${params.rrna_method} requires --rrna_refs (or use --rrna_method none)."
        System.exit(1)
    }
    if (!(params.mapper in ['strobealign', 'bbmap'])) {
        log.error "ERROR: --mapper must be strobealign or bbmap."
        System.exit(1)
    }
    if (params.chunk && !(params.chunk.toString() ==~ /\d+\/\d+/)) {
        log.error "ERROR: --chunk must look like I/N (e.g. 2/6)."
        System.exit(1)
    }
}

// ----------------------------------------------------------------------------
// Read discovery
// ----------------------------------------------------------------------------

// Returns a sorted list of [sample_id, r1, r2]. Sorting makes --chunk stable
// across separate runs over the same input.
def discoverReads() {
    def samples = []
    if (params.samplesheet) {
        def sheet = file(params.samplesheet, checkIfExists: true)
        def lines = sheet.readLines().findAll { it.trim() && !it.startsWith('#') }
        def header = lines[0].split(',')*.trim()
        if (header[0..2] != ['sample', 'r1', 'r2']) {
            error "ERROR: --samplesheet header must start with sample,r1,r2 (got: ${lines[0]})"
        }
        lines.drop(1).each { line ->
            def f = line.split(',')*.trim()
            samples << [f[0], file(f[1], checkIfExists: true), file(f[2], checkIfExists: true)]
        }
    } else {
        def input_dir = file(params.input)
        if (!input_dir.isDirectory()) {
            error "ERROR: --input directory does not exist: ${params.input}"
        }
        def r1_files = file("${params.input}/**_R1_*.fastq.gz") +
                       file("${params.input}/**_1.fq.gz") +
                       file("${params.input}/**_1.fastq.gz")
        r1_files.each { r1 ->
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
            def r2 = r1.resolveSibling(r2_name)
            if (!r2.exists()) {
                log.warn "[WARNING] Missing R2 for ${sample_id} (expected ${r2}); skipping."
                return
            }
            samples << [sample_id, r1, r2]
        }
    }
    if (!samples) {
        error "ERROR: No paired-end reads found.\n" +
              "  Looked for *_R1_*.fastq.gz (Illumina) and *_{1,2}.fq.gz/*.fastq.gz (BGI)."
    }
    def dups = samples.countBy { it[0] }.findAll { k, v -> v > 1 }
    if (dups) {
        error "ERROR: Duplicate sample IDs: ${dups.keySet().join(', ')}"
    }
    return samples.sort { it[0] }
}

def selectChunk(samples) {
    if (!params.chunk) return samples
    def (i, n) = params.chunk.toString().tokenize('/')*.toInteger()
    if (i < 1 || i > n) error "ERROR: --chunk ${params.chunk}: I must be in 1..N"
    // Contiguous blocks keep each chunk's samples together in the store dir.
    def per = Math.ceil(samples.size() / n) as int
    return samples.drop((i - 1) * per).take(per)
}

// ----------------------------------------------------------------------------
// Reference discovery
// ----------------------------------------------------------------------------

// Each reference is a map:
//   name, fasta, gff, bins, checkm2, gtdbtk, annotations
// where everything after fasta may be null.
def firstExisting(paths) {
    paths.collect { file(it) }.find { it.exists() }
}

def referencesFromMag(dir_path) {
    def dir = file(dir_path)
    if (!dir.isDirectory()) error "ERROR: --from_mag directory does not exist: ${dir_path}"
    def fasta = firstExisting(["${dir}/assembly/assembly.fasta"])
    if (!fasta) error "ERROR: no assembly/assembly.fasta under --from_mag ${dir_path}"
    def gff = firstExisting([
        "${dir}/annotation/bakta/extra/annotation.gff3",
        "${dir}/annotation/bakta/basic/annotation.gff3",
    ]) ?: file("${dir}/annotation/prokka").with { it.isDirectory() ? it.listFiles().find { f -> f.name.endsWith('.gff') } : null }
    def bs = params.bin_set
    def bins = firstExisting([
        "${dir}/binning/${bs}/contig2bin.tsv",
        "${dir}/binning/${bs}/${bs}_bins.tsv",
    ])
    if (!bins) log.warn "[WARNING] No ${bs} contig→bin table under ${dir}/binning/${bs}; per-MAG tables will be skipped."
    if (!gff)  log.warn "[WARNING] No Bakta/Prokka GFF under ${dir}/annotation; quantification will be contig-level only."
    return [[
        name:        params.ref_name ?: dir.name,
        fasta:       fasta,
        gff:         gff,
        bins:        bins,
        checkm2:     firstExisting(["${dir}/binning/checkm2/quality_report.tsv"]),
        gtdbtk:      firstExisting(["${dir}/taxonomy/gtdbtk/gtdbtk_taxonomy.tsv"]),
        annotations: firstExisting(["${dir}/metabolism/merged/merged_annotations.tsv"]),
        kofam:       firstExisting(["${dir}/metabolism/kofamscan/kofamscan_results.tsv"]),
        rrna_genes:  firstExisting(["${dir}/taxonomy/rrna/rrna_genes.tsv"]),
        // contig t-SNE computed by the mag_analysis viz, reused for the RNA map
        contig_explorer: firstExisting([
            "${dir}/viz/data/contig_explorer_manifest.json.gz",
            "${dir}/viz/data/contig_explorer_manifest.json",
        ])?.parent,
    ]]
}

def referencesFromAssembly(dir_path) {
    def dir = file(dir_path)
    if (!dir.isDirectory()) error "ERROR: --from_assembly directory does not exist: ${dir_path}"
    // nanopore_assembly: assembly/assembly.fasta
    def single = file("${dir}/assembly/assembly.fasta")
    if (single.exists()) {
        return [[name: params.ref_name ?: dir.name, fasta: single]]
    }
    // illumina_assembly: assembly/<name>/<name>.dedupe.fasta (else .megahit.fasta)
    def refs = []
    file("${dir}/assembly").listFiles().findAll { it.isDirectory() }.sort { it.name }.each { sub ->
        def fa = firstExisting(["${sub}/${sub.name}.dedupe.fasta", "${sub}/${sub.name}.megahit.fasta"])
        if (fa) refs << [name: sub.name, fasta: fa]
    }
    if (!refs) error "ERROR: no assembly FASTA found under --from_assembly ${dir_path}/assembly"
    return refs
}

def referencesFromDir(dir_path) {
    def dir = file(dir_path)
    if (!dir.isDirectory()) error "ERROR: --references directory does not exist: ${dir_path}"
    def fastas = file("${dir}/*.{fasta,fa,fna}")
    if (!fastas) error "ERROR: No *.fasta/*.fa/*.fna files found in ${dir_path}."
    return fastas.sort { it.name }.collect { fa ->
        def name = fa.baseName
        [
            name: name,
            fasta: fa,
            gff:  firstExisting(["${dir}/${name}.gff", "${dir}/${name}.gff3"]),
            bins: firstExisting(["${dir}/${name}.bins.tsv"]),
        ]
    }
}

def discoverReferences() {
    def refs = params.from_mag      ? referencesFromMag(params.from_mag) :
               params.from_assembly ? referencesFromAssembly(params.from_assembly) :
                                      referencesFromDir(params.references)
    if (params.contig_explorer) {
        def ce = file(params.contig_explorer, checkIfExists: true)
        refs.each { r -> r.contig_explorer = ce }
    }
    refs.each { r ->
        if (!(r.name ==~ /[A-Za-z0-9._-]+/)) error "ERROR: reference name '${r.name}' must be alnum/._- only (set --ref_name)"
    }
    return refs
}

// ============================================================================
// Import modules
// ============================================================================

include { QC_READS }              from './modules/preprocess'
include { HUMAN_INDEX }           from './modules/preprocess'
include { REMOVE_HUMAN }          from './modules/preprocess'
include { FASTQC }                from './modules/preprocess'
include { PREPARE_RRNA_REF }      from './modules/rrna'
include { REMOVE_RRNA_SORTMERNA } from './modules/rrna'
include { PREPARE_RRNA_TAXDB }    from './modules/rrna'
include { RRNA_TAXONOMY }         from './modules/rrna'
include { REF_FEATURES }          from './modules/quantify'
include { STROBEALIGN_INDEX }     from './modules/alignment'
include { MAP_READS_STROBEALIGN } from './modules/alignment'
include { BBMAP_INDEX }           from './modules/alignment'
include { MAP_READS_BBMAP }       from './modules/alignment'
include { FEATURECOUNTS }         from './modules/quantify'
include { SUMMARIZE_EXPRESSION }  from './modules/summarize'
include { VIZ_PREPROCESS }        from './modules/summarize'

// ============================================================================
// Main workflow
// ============================================================================

workflow {

    main:

    validateParams()

    // ----- 1. References -----
    def refs = discoverReferences()
    refs.each { r ->
        log.info "Reference ${r.name}: ${r.fasta}" +
                 (r.gff  ? "\n  genes: ${r.gff}"  : "\n  genes: none (contig-level only)") +
                 (r.bins ? "\n  bins:  ${r.bins}" : "")
    }
    def empty = file("${projectDir}/assets/NO_FILE")
    // Only the name travels as the ref meta; files go through path() inputs so
    // they are staged (and bound into containers) properly.
    ch_ref_fasta = Channel.fromList(refs.collect { r -> [[name: r.name], r.fasta] })
    ch_ref_gff   = Channel.fromList(refs.collect { r -> [[name: r.name], r.fasta, r.gff ?: empty, r.rrna_genes ?: empty] })
    ch_ref_extra = Channel.fromList(refs.collect { r ->
        [r.name, r.bins ?: empty, r.checkm2 ?: empty, r.gtdbtk ?: empty, r.annotations ?: empty,
         r.contig_explorer ?: empty, r.kofam ?: empty]
    })

    REF_FEATURES(ch_ref_gff)

    if (params.mapper == 'strobealign') {
        STROBEALIGN_INDEX(ch_ref_fasta)
        ch_index = STROBEALIGN_INDEX.out.index
    } else {
        BBMAP_INDEX(ch_ref_fasta)
        ch_index = BBMAP_INDEX.out.index
    }

    // Shared per-run references: rRNA k-mers and the human index, built once.
    ch_rrna_ref = params.rrna_method == 'bbduk'
        ? PREPARE_RRNA_REF(Channel.value(file(params.rrna_refs, checkIfExists: true))).ref
        : Channel.value(empty)

    def do_rrna_tax = params.rrna_method == 'bbduk' && params.run_rrna_taxonomy
    ch_rrna_taxdb = do_rrna_tax
        ? PREPARE_RRNA_TAXDB(Channel.value(file(params.rrna_refs, checkIfExists: true))).db
        : Channel.empty()

    if (params.run_remove_human) {
        def human = file(params.human_ref, checkIfExists: true)
        ch_human_index = HUMAN_INDEX(Channel.value(human)).index
    }

    if (!params.prepare_only) {

    // ----- 2. Reads -----
    def all_samples = discoverReads()
    def samples = selectChunk(all_samples)
    log.info "Detected ${all_samples.size()} samples" +
             (params.chunk ? "; chunk ${params.chunk} has ${samples.size()} (${samples ? samples[0][0] + ' … ' + samples[-1][0] : 'none'})" : "")
    ch_raw_reads = Channel.fromList(samples.collect { id, r1, r2 -> [[id: id], r1, r2] })

    // ----- 3. QC (+ bbduk rRNA split) -----
    QC_READS(ch_raw_reads, ch_rrna_ref)
    ch_readcounts = QC_READS.out.counts
    ch_reads = QC_READS.out.reads

    // Taxonomy of the rRNA reads the QC split off
    if (do_rrna_tax) {
        RRNA_TAXONOMY(QC_READS.out.rrna.combine(ch_rrna_taxdb))
        ch_rrna_tax = RRNA_TAXONOMY.out.taxonomy.map { meta, f -> f }.collect()
    } else {
        ch_rrna_tax = Channel.value([])
    }

    if (params.rrna_method == 'sortmerna') {
        REMOVE_RRNA_SORTMERNA(ch_reads)
        ch_reads = REMOVE_RRNA_SORTMERNA.out.reads
        ch_readcounts = ch_readcounts.mix(REMOVE_RRNA_SORTMERNA.out.counts)
    }

    if (params.run_remove_human) {
        REMOVE_HUMAN(ch_reads, ch_human_index)
        ch_reads = REMOVE_HUMAN.out.reads
        ch_readcounts = ch_readcounts.mix(REMOVE_HUMAN.out.counts)
    }

    if (params.run_fastqc) {
        FASTQC(ch_reads)
    }

    // ----- 4. Mapping: every sample × every reference -----
    ch_map_input = ch_reads.combine(ch_index)
    if (params.mapper == 'strobealign') {
        MAP_READS_STROBEALIGN(ch_map_input)
        ch_bam = MAP_READS_STROBEALIGN.out.bam
        ch_flagstat = MAP_READS_STROBEALIGN.out.flagstat
    } else {
        MAP_READS_BBMAP(ch_map_input)
        ch_bam = MAP_READS_BBMAP.out.bam.map { meta, ref, bam, bai -> [meta, ref, bam] }
        ch_flagstat = MAP_READS_BBMAP.out.flagstat
    }

    // ----- 5. Quantification -----
    ch_quant_input = ch_bam
        .map { meta, ref, bam -> [ref.name, meta, ref, bam] }
        .combine(REF_FEATURES.out.saf.map { ref, genes, contigs -> [ref.name, genes, contigs] }, by: 0)
        .map { name, meta, ref, bam, genes, contigs -> [meta, ref, bam, genes, contigs] }

    FEATURECOUNTS(ch_quant_input)

    // ----- 6. Merged matrices + viz JSONs -----
    if (params.summarize) {
        ch_by_ref = FEATURECOUNTS.out.counts
            .map { meta, ref, genes, contigs, summary -> [ref.name, genes, contigs, summary] }
            .groupTuple()
        ch_flagstat_by_ref = ch_flagstat
            .map { meta, ref, f -> [ref.name, f] }
            .groupTuple()
        ch_all_counts = ch_readcounts.map { meta, f -> f }.collect()

        ch_summary_input = ch_by_ref
            .join(ch_flagstat_by_ref)
            .join(REF_FEATURES.out.annotation.map { ref, annot, lengths -> [ref.name, annot, lengths] })
            .join(ch_ref_extra)
            .combine(ch_all_counts.map { files -> [files] })
            .combine(ch_rrna_tax.ifEmpty([]).map { files -> [files] })

        SUMMARIZE_EXPRESSION(ch_summary_input)

        VIZ_PREPROCESS(
            SUMMARIZE_EXPRESSION.out.viz_tables.map { name, files -> files }.flatten().collect()
        )
    }

    } // !prepare_only

    workflow.onComplete = {
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

    workflow.onError = {
        println "[ERROR] Pipeline failed: ${workflow.errorMessage}"
    }
}
