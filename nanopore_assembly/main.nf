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
      --metamdbg_gfa     Also build metaMDBG's base-space assembly graph (realigns all
                         reads; about as slow as the assembly) [default: false]
      --read_type STR    Flye read mode, REQUIRED: nano-raw | nano-hq | nano-corr
                         (Flye: nano-hq for SUP reads under ~5% error, nano-raw
                         for older reads up to ~20%). 'auto' guesses from the
                         first 10,000 reads and is not recommended.
      --min_overlap N    Flye --min-overlap [default: 1000]
      --polish           Enable Flye polishing iterations [default: true for flye]
      --dedupe           Enable BBDuk deduplication before assembly
      --filtlong_size N  Filtlong target bases (e.g. 40000000000); skip if not set
      --onepass_filter   Select those bases with the legacy single-pass streaming
                         threshold. Cheaper on disk, but the kept subset depends on
                         input order and is biased toward whatever streams first.
                         Default is an order-independent two-pass bucket sort.

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

def validateParams() {
    if (params.help) {
        helpMessage()
        System.exit(0)
    }
    if (!params.input) {
        log.error "ERROR: --input is required. Provide path to directory containing *.fastq.gz files. Run with --help for usage."
        System.exit(1)
    }
    // Flye's read mode must be chosen, not inferred. Refusing here costs nothing:
    // no node has been allocated to any task yet.
    def readTypes = ['nano-raw', 'nano-hq', 'nano-corr', 'auto']
    if (params.assembler == 'flye' && !params.read_type) {
        log.error "ERROR: --read_type is required: nano-raw, nano-hq or nano-corr. " +
                  "It sets Flye's index, overlap settings and error model, so it is not inferred " +
                  "by default. Flye's criterion is error rate: nano-hq for Guppy5+/Dorado SUP reads " +
                  "under ~5% error, nano-raw for older reads up to ~20%. Judge by error-based " +
                  "read quality across the whole set, not the first reads of one file."
        System.exit(1)
    }
    if (params.read_type && !(params.read_type in readTypes)) {
        log.error "ERROR: Invalid --read_type '${params.read_type}'. Choose from: ${readTypes.join(', ')}"
        System.exit(1)
    }
    if (params.read_type == 'auto') {
        log.warn "--read_type auto: Flye's mode will be guessed from the first 10,000 reads " +
                 "of the stream, which can be one barcode and unrepresentative of the rest"
    }
    if (!(params.assembler in ['flye', 'metamdbg', 'myloasm'])) {
        log.error "ERROR: Invalid --assembler '${params.assembler}'. Choose from: flye, metamdbg, myloasm"
        System.exit(1)
    }
    log.info "Assembler: ${params.assembler}"
}

// Import modules
include { CONCAT_READS }        from './modules/preprocess'
include { PREPARE_READS }        from './modules/preprocess'
include { REMOVE_HUMAN }        from './modules/preprocess'
include { FLYE_ASSEMBLE }       from './modules/assembly'
include { FLYE_POLISH }         from './modules/assembly'
include { PUBLISH_UNPOLISHED } from './modules/assembly'
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

    main:

    validateParams()

    WRITE_PROVENANCE()

    // 1. Discover input reads — auto-detect nanopore barcode vs flat directory
    def input_dir = file(params.input)
    if (!input_dir.isDirectory()) {
        error "ERROR: --input directory does not exist: ${params.input}\nRun with --help for usage."
    }

    def barcode_dirs = file("${params.input}/**/fastq_pass/barcode*", type: 'dir') +
                       file("${params.input}/fastq_pass/barcode*", type: 'dir')
    // Flat mode accepts fastq(.gz) and fasta(.fa/.fasta[.gz]). FASTA input (e.g.
    // pre-QC'd reads exported from nanopore_live's fa/ store) is converted to
    // fastq with placeholder quality inside PREPARE_READS; pair with
    // --read_type to avoid quality-based Flye mode auto-detection.
    def flat_fastqs  = file("${params.input}/*.{fastq.gz,fq.gz,fa,fasta,fa.gz,fasta.gz}")

    if (barcode_dirs) {
        log.info "Detected nanopore barcode structure: ${barcode_dirs.size()} barcode directories"
        def all_pairs = []
        def empty_barcodes = []
        for (dir in barcode_dirs) {
            def barcode = dir.name
            def run_name = dir.parent.parent.name
            def parts = run_name.tokenize('_')
            def flowcell = parts.size() >= 4 ? parts[3] : run_name
            def sample_id = "${flowcell}_${barcode}"
            def fqs = file("${dir}/*.fastq.gz")
            if (fqs instanceof List) {
                if (fqs.isEmpty()) { empty_barcodes.add(dir) }
                for (fq in fqs) { all_pairs.add([sample_id, fq]) }
            } else if (fqs) {
                all_pairs.add([sample_id, fqs])
            } else {
                empty_barcodes.add(dir)
            }
        }
        log.info "Found ${all_pairs.size()} FASTQ files across ${barcode_dirs.size()} barcodes"
        // A barcode that contributes no file is dropped by groupTuple with no
        // error, so it is missing from both the assembly and the depth matrix.
        // The barcode-directory count above stays correct either way, so warn
        // explicitly; stale symlinks to renamed or concatenated files are the
        // usual cause.
        if (empty_barcodes) {
            log.warn "${empty_barcodes.size()} barcode director(ies) contain no *.fastq.gz and are ABSENT from the assembly and the depth matrix."
            log.warn "  A dangling symlink farm is the usual cause (targets renamed, moved or concatenated)."
            empty_barcodes.take(20).each { log.warn "    ${it}" }
            if (empty_barcodes.size() > 20) {
                log.warn "    ... and ${empty_barcodes.size() - 20} more"
            }
        }
        ch_barcode_raw = Channel.from(all_pairs)
            .groupTuple()
            .map { sample_id, fastqs -> [[id: sample_id], fastqs] }
            .filter { meta, fastqs -> fastqs.size() > 0 }

        CONCAT_READS(ch_barcode_raw)
        ch_reads = CONCAT_READS.out.reads
            .filter { meta, fastq -> fastq.size() > 1024 }
    } else if (flat_fastqs) {
        log.info "Detected flat read directory: ${flat_fastqs.size()} files (fastq/fasta)"
        ch_reads = Channel.fromPath("${params.input}/*.{fastq.gz,fq.gz,fa,fasta,fa.gz,fasta.gz}")
            .map { reads ->
                def name = reads.baseName.replaceAll(/\.(fastq|fq|fasta|fa)$/, '')
                [[id: name], reads]
            }
    } else {
        error "ERROR: No reads found in ${params.input}. Expected *.fastq.gz / *.fasta[.gz] or a fastq_pass/barcode*/ structure.\nRun with --help for usage."
    }

    // 2. Concatenate + dedupe + optional filtlong -> single all_reads.fastq.gz
    //
    // Order matters here and must not be left to chance. A plain collect()
    // emits in CONCAT_READS *completion* order, which varies run to run, and
    // fastq_filter's single-pass --target_bases mode decides accept/reject on
    // arrival against a threshold built only from the reads seen so far. It is
    // therefore lenient early and strict late, and reads with identical length
    // and quality are kept or dropped by their position in the stream.
    //
    // Largest file first, with the name as tiebreak so equal sizes cannot
    // reintroduce the nondeterminism. This pins the selection; it does not make
    // it unbiased; fastq_filter's default two-pass selection is what makes the
    // choice order-independent. Only --onepass still depends on this ordering.
    ch_per_barcode = ch_reads
        .map { meta, fastq -> fastq }
        .collect()
        .map { files -> files.toSorted { a, b -> (b.size() <=> a.size()) ?: (a.name <=> b.name) } }
    PREPARE_READS(ch_per_barcode)
    // Surface PREPARE_READS' warnings in the Nextflow log, not only in the
    // task's .command.err, where two truncated read sets went unread.
    PREPARE_READS.out.warnings.subscribe { f -> f.readLines().each { log.warn "PREPARE_READS: ${it}" } }

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
        // No polishing: the draft IS the final assembly, so republish it under the
        // plain name rather than leaving only draft_assembly.fasta behind.
        PUBLISH_UNPOLISHED(ch_raw_assembly, ch_asm_info, ch_asm_graph)
        ch_assembly  = PUBLISH_UNPOLISHED.out.assembly
        ch_asm_info  = PUBLISH_UNPOLISHED.out.info
        ch_asm_graph = PUBLISH_UNPOLISHED.out.graph
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
    // Every sample we asked MAP_READS to map. Passed to CALCULATE_DEPTHS so it
    // can name the ones that never produced a BAM -- errorStrategy 'ignore'
    // makes those failures invisible otherwise.
    ch_expected_ids = ch_reads.map { meta, fastq -> meta.id }.collect()
    CALCULATE_DEPTHS(ch_bam_files, ch_assembly, ch_expected_ids)

    workflow.onComplete = {
        // The draft is a rescue/-resume artifact, not a deliverable. Once the
        // polished assembly is in place, drop the draft copies from the outdir.
        if (workflow.success && do_polish) {
            def asmdir = file("${params.outdir}/assembly")
            def polished = asmdir.resolve('assembly.fasta')
            if (polished.exists() && polished.size() > 0) {
                ['draft_assembly.fasta', 'draft_assembly_info.txt',
                 'draft_assembly_graph.gfa'].each { n ->
                    def f = asmdir.resolve(n)
                    if (f.exists()) {
                        log.info "Removing draft artifact ${n} (polished assembly present)"
                        f.delete()
                    }
                }
            } else {
                log.warn "Polished assembly missing or empty - keeping draft_* for rescue"
            }
        }

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

        // Run manifest — record pipeline version, resolved parameters, and per-process tool
        // images so downstream reporting (omc-platform Methods drafting + agents) can state
        // what actually ran instead of "not specified in the outputs". (danaSeq #24)
        try {
            def outDir = new File("${params.outdir}")
            if (outDir.exists()) {
                def safeParams = params.collectEntries { k, v ->
                    [k, (v == null || v instanceof Number || v instanceof Boolean || v instanceof String || v instanceof List || v instanceof Map) ? v : "${v}"]
                }
                def sysEnv = System.getenv()
                def bakedSha = sysEnv['DANASEQ_GIT_SHA']
                def manifest = [
                    pipeline        : "${workflow.manifest.name ?: 'danaSeq'} v${workflow.manifest.version ?: 'dev'}",
                    revision        : (workflow.revision ?: workflow.commitId ?: null),
                    // Null for every run launched from the .sif — there is no git
                    // repo inside — so fall back to the SHA baked in at image build.
                    commit_id       : (workflow.commitId ?: bakedSha),
                    // Where that id came from, as a token a reader can switch on.
                    commit_source   : (workflow.commitId ? 'git-checkout'
                                       : (bakedSha ? 'container-build' : 'unknown')),
                    container_git_ref: sysEnv['DANASEQ_GIT_REF'],
                    container_built : sysEnv['DANASEQ_BUILD_DATE'],
                    // MD5 of main.nf — pins the source even with no git and no build arg.
                    script_id       : workflow.scriptId,
                    nextflow_version: "${nextflow.version}",
                    command_line    : workflow.commandLine,
                    started         : "${workflow.start}",
                    completed       : "${workflow.complete}",
                    duration        : "${workflow.duration}",
                    success         : workflow.success,
                    containers      : workflow.container,
                    parameters      : safeParams,
                ]
                new File(outDir, 'run_manifest.json').text =
                    groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(manifest))
                println "[INFO] Wrote run manifest: ${params.outdir}/run_manifest.json"
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
