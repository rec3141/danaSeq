#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// ============================================================================
// Dana MAG Analysis Pipeline - Nextflow DSL2
// ============================================================================
//
// Technology-agnostic downstream analysis pipeline for metagenome-assembled
// genomes (MAGs). Accepts a pre-computed assembly + depth table (+ optional
// BAMs) from any assembler (nanopore_assembly, illumina_assembly, or external).
//
// Runs binning (7 binners + 3 consensus), annotation, taxonomy, metabolic
// profiling, mobile genetic element detection, eukaryotic analysis, ecosystem
// services mapping, phylogenetics, and interactive visualization.
//
// Usage:
//   nextflow run main.nf --assembly /path/to/assembly.fasta \
//       --depths /path/to/depths.txt [options] -resume
//
// ============================================================================

// ============================================================================
// Help message
// ============================================================================

def helpMessage() {
    log.info """
    =========================================
     Dana MAG Analysis Pipeline
     https://github.com/rec3141/danaSeq
    =========================================

    Usage:
      nextflow run main.nf --assembly FILE --depths FILE [options] -resume

    Required:
      --assembly FILE    Path to assembly FASTA (e.g. assembly.fasta)
      --depths FILE      Path to MetaBAT2-format depth table (e.g. depths.txt)
      --outdir DIR       Output directory [default: results]

    Optional inputs:
      --bam_dir DIR      Directory containing *.sorted.bam + *.sorted.bam.bai files
                         (needed for BAM-based binners: SemiBin2, LorBin, COMEBin)
      --tnf FILE         Path to tetranucleotide frequency table (for viz t-SNE/UMAP)
      --assembly_info FILE  Assembler metadata (for viz)

    Caching:
      --store_dir DIR    Persistent cache directory (storeDir)

    Binning:
      --run_semibin      Include SemiBin2 (needs --bam_dir) [default: false]
      --run_maxbin       Include MaxBin2 [default: false]
      --run_lorbin       Include LorBin (needs --bam_dir) [default: false]
      --run_comebin      Include COMEBin (needs --bam_dir) [default: false]
      --run_vamb         Include VAMB [default: false]
      --run_vamb_tax     Include taxonomy-guided VAMB [default: false]
      --run_binette      Run Binette consensus (needs --checkm2_db) [default: false]
      --run_magscot      Run MAGScoT consensus [default: false]

    Annotation:
      --annotator STR    'prokka', 'bakta', or 'none' [default: bakta]

    Taxonomy:
      --run_kaiju, --run_kraken2, --run_sendsketch, --run_rrna

    MGE detection:
      --run_genomad, --run_checkv, --run_integronfinder, --run_islandpath,
      --run_macsyfinder, --run_defensefinder

    Metabolic profiling:
      --run_metabolism, --run_antismash, --run_ecossdb

    Eukaryotic:
      --run_eukaryotic, --run_metaeuk, --run_marferret

    Quality & Phylogenetics:
      --checkm2_db, --run_gtdbtk, --gtdbtk_db

    Visualization:
      --run_viz          Build interactive dashboard [default: false]

    Shortcut:
      --all              Enable all optional modules (still requires DB paths)

    """.stripIndent()
}

if (params.help) {
    helpMessage()
    System.exit(0)
}

// ============================================================================
// Parameter validation
// ============================================================================

if (!params.assembly) {
    log.error "ERROR: --assembly is required. Provide path to assembly FASTA."
    System.exit(1)
}

if (!params.depths) {
    log.error "ERROR: --depths is required. Provide path to MetaBAT2-format depth table."
    System.exit(1)
}

// Resolve annotator
def effective_annotator = params.annotator ?:
    (params.run_bakta ? 'bakta' : (params.run_prokka ? 'prokka' : 'none'))

if (effective_annotator == 'bakta' && !params.bakta_light_db) {
    log.error "ERROR: --bakta_light_db is required when using Bakta annotation."
    System.exit(1)
}
if (effective_annotator == 'bakta' && params.bakta_extra && !params.bakta_db) {
    log.error "ERROR: --bakta_db is required when using --bakta_extra."
    System.exit(1)
}

if (effective_annotator != 'none') {
    log.info "Annotator: ${effective_annotator}"
}

// ============================================================================
// Import modules
// ============================================================================

include { BIN_SEMIBIN2 }        from './modules/binning'
include { BIN_METABAT2 }        from './modules/binning'
include { BIN_MAXBIN2 }         from './modules/binning'
include { BIN_LORBIN }          from './modules/binning'
include { BIN_COMEBIN }         from './modules/binning'
include { BIN_VAMB }            from './modules/binning'
include { BIN_VAMB_TAX }        from './modules/binning'
include { DASTOOL_CONSENSUS }   from './modules/binning'
include { BINETTE_CONSENSUS }   from './modules/binning'
include { MAGSCOT_CONSENSUS }   from './modules/binning'
include { CHECKM2 }             from './modules/binning'
include { GTDBTK_CLASSIFY }     from './modules/phylogeny'
include { PROKKA_ANNOTATE }     from './modules/annotation'
include { BAKTA_BASIC }          from './modules/annotation'
include { BAKTA_EXTRA }          from './modules/annotation'
include { KAIJU_CONTIG_CLASSIFY } from './modules/taxonomy'
include { KAIJU_CLASSIFY }      from './modules/taxonomy'
include { KRAKEN2_CLASSIFY }   from './modules/taxonomy'
include { SENDSKETCH_CLASSIFY } from './modules/taxonomy'
include { RNA_CLASSIFY }       from './modules/rrna'
include { GENOMAD_CLASSIFY }    from './modules/mge'
include { CHECKV_QUALITY }      from './modules/mge'
include { INTEGRONFINDER }      from './modules/mge'
include { ISLANDPATH_DIMOB }    from './modules/mge'
include { MACSYFINDER }         from './modules/mge'
include { DEFENSEFINDER }       from './modules/mge'
include { TIARA_CLASSIFY }      from './modules/eukaryotic'
include { WHOKARYOTE_CLASSIFY } from './modules/eukaryotic'
include { METAEUK_PREDICT }     from './modules/eukaryotic'
include { MARFERRET_CLASSIFY }  from './modules/eukaryotic'
include { KOFAMSCAN }           from './modules/metabolism'
include { EMAPPER }             from './modules/metabolism'
include { DBCAN }               from './modules/metabolism'
include { MERGE_ANNOTATIONS }   from './modules/metabolism'
include { MAP_TO_BINS }         from './modules/metabolism'
include { KEGG_MODULES }        from './modules/metabolism'
include { MINPATH }             from './modules/metabolism'
include { KEGG_DECODER }        from './modules/metabolism'
include { ANTISMASH }           from './modules/metabolism'
include { ECOSSDB_MAP }         from './modules/metabolism'
include { ECOSSDB_SCORE }       from './modules/metabolism'
include { ECOSSDB_SDG }         from './modules/metabolism'
include { ECOSSDB_VIZ }         from './modules/metabolism'
include { CALCULATE_GENE_DEPTHS } from './modules/gene_depths'
include { VIZ_PREPROCESS as VIZ_STAGE1 } from './modules/viz'
include { VIZ_PREPROCESS as VIZ_STAGE2 } from './modules/viz'
include { VIZ_PREPROCESS as VIZ_STAGE3 } from './modules/viz'
include { VIZ_PREPROCESS as VIZ_STAGE4 } from './modules/viz'

// ============================================================================
// Helper: stage assembly inputs into outdir so viz preprocess can find them
// ============================================================================

process STAGE_INPUTS {
    tag "stage_inputs"
    label 'process_low'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path(assembly)
    path(depths)

    output:
    path("assembly/assembly.fasta")
    path("assembly/assembly_info.txt")
    path("assembly/gc.tsv")
    path("mapping/depths.txt")

    script:
    """
    mkdir -p assembly mapping
    cp ${assembly} assembly/assembly.fasta
    # Generate assembly_info.txt (contig lengths) from FASTA
    awk '/^>/{if(name) print name "\\t" len; name=substr(\$1,2); len=0; next} {len+=length(\$0)} END{if(name) print name "\\t" len}' ${assembly} \
        | sort -k2,2nr \
        | awk 'BEGIN{print "#seq_name\\tlength\\tcov.\\tcirc.\\trepeat\\tmult.\\talt_group\\tgraph_path"} {print \$1 "\\t" \$2 "\\tNA\\tN\\tN\\t1\\t*\\t*"}' \
        > assembly/assembly_info.txt
    # Compute per-contig GC%
    awk '/^>/{if(name) printf "%s\\t%.2f\\n", name, (gc/(gc+at))*100; name=substr(\$1,2); gc=0; at=0; next} {for(i=1;i<=length(\$0);i++){c=substr(\$0,i,1); if(c~/[GCgc]/)gc++; else at++}} END{if(name) printf "%s\\t%.2f\\n", name, (gc/(gc+at))*100}' ${assembly} \
        | awk 'BEGIN{print "contig_id\\tgc_pct"} {print}' \
        > assembly/gc.tsv
    cp ${depths} mapping/depths.txt
    """
}

// ============================================================================
// Main workflow
// ============================================================================

workflow {

    // 1. Create channels from file-based inputs
    ch_assembly = Channel.fromPath(params.assembly, checkIfExists: true)
    ch_depths   = Channel.fromPath(params.depths, checkIfExists: true)

    // Optional BAM files for BAM-based binners
    if (params.bam_dir) {
        ch_bam_files = Channel.fromPath("${params.bam_dir}/*.sorted.bam")
            .mix(Channel.fromPath("${params.bam_dir}/*.sorted.bam.bai"))
            .collect()
    } else {
        ch_bam_files = Channel.empty()
    }

    // Optional TNF for viz
    if (params.tnf) {
        ch_tnf = Channel.fromPath(params.tnf, checkIfExists: true)
    } else {
        ch_tnf = Channel.empty()
    }

    // 1b. Publish assembly inputs into outdir so viz preprocess can find them
    STAGE_INPUTS(ch_assembly, ch_depths)

    // 2. Gene annotation on assembly (optional)
    ch_proteins       = Channel.empty()
    ch_gff            = Channel.empty()
    ch_annotation_tsv = Channel.empty()

    if (effective_annotator == 'prokka') {
        PROKKA_ANNOTATE(ch_assembly)
        ch_proteins       = PROKKA_ANNOTATE.out.proteins
        ch_gff            = PROKKA_ANNOTATE.out.gff
        ch_annotation_tsv = PROKKA_ANNOTATE.out.tsv
    } else if (effective_annotator == 'bakta') {
        BAKTA_BASIC(ch_assembly)
        ch_proteins       = BAKTA_BASIC.out.proteins
        ch_gff            = BAKTA_BASIC.out.gff
        ch_annotation_tsv = BAKTA_BASIC.out.tsv

        if (params.bakta_extra) {
            ch_assembly_after_basic = BAKTA_BASIC.out.proteins
                .combine(ch_assembly)
                .map { prot, asm -> asm }
            BAKTA_EXTRA(ch_assembly_after_basic)
            ch_annotation_tsv = BAKTA_EXTRA.out.tsv
        }
    }

    // 2b. Kaiju taxonomy — six-frame translation on contigs
    if (params.run_kaiju && params.kaiju_db) {
        KAIJU_CONTIG_CLASSIFY(ch_assembly)
    }

    // 2b2. Kaiju per-gene (requires annotation)
    if (params.run_kaiju && params.kaiju_db && effective_annotator != 'none') {
        KAIJU_CLASSIFY(ch_proteins, ch_gff)
    }

    // 2b3. Kraken2 k-mer taxonomy
    if (params.run_kraken2 && params.kraken2_db) {
        KRAKEN2_CLASSIFY(ch_assembly)
    }

    // 2b4. BBSketch/sendsketch GTDB taxonomy
    if (params.run_sendsketch && params.sendsketch_address) {
        SENDSKETCH_CLASSIFY(ch_assembly)
    }

    // 2b5. RNA gene classification
    if (params.run_rrna && params.silva_ssu_db) {
        RNA_CLASSIFY(ch_assembly)
    }

    // 2c. Eukaryotic classification
    if (params.run_eukaryotic) {
        TIARA_CLASSIFY(ch_assembly)
        WHOKARYOTE_CLASSIFY(ch_assembly, ch_gff)

        if (params.run_metaeuk && params.metaeuk_db) {
            METAEUK_PREDICT(
                ch_assembly,
                TIARA_CLASSIFY.out.classifications,
                WHOKARYOTE_CLASSIFY.out.classifications
            )

            if (params.run_marferret && params.marferret_db) {
                MARFERRET_CLASSIFY(
                    METAEUK_PREDICT.out.proteins,
                    file(params.marferret_db)
                )
            }
        }
    }

    // 2d. Mobile genetic element detection
    if (params.run_genomad && params.genomad_db) {
        GENOMAD_CLASSIFY(ch_assembly)

        if (params.run_checkv && params.checkv_db) {
            CHECKV_QUALITY(GENOMAD_CLASSIFY.out.virus_fasta)
        }
    }

    if (params.run_integronfinder) {
        INTEGRONFINDER(ch_assembly)
    }

    if (params.run_islandpath && effective_annotator != 'none') {
        ISLANDPATH_DIMOB(ch_assembly, ch_gff, ch_proteins)
    }

    if (params.run_macsyfinder && effective_annotator != 'none' && params.macsyfinder_models) {
        MACSYFINDER(ch_proteins)
    }

    if (params.run_defensefinder && effective_annotator != 'none') {
        DEFENSEFINDER(ch_proteins, ch_gff)
    }

    // 2e. Metabolic profiling
    ch_annot_for_merge = Channel.empty()

    if (params.run_metabolism && effective_annotator != 'none') {
        if (params.kofam_db) {
            KOFAMSCAN(ch_proteins, file("${params.kofam_db}/profiles"), file("${params.kofam_db}/ko_list"))
            ch_annot_for_merge = ch_annot_for_merge.mix(
                KOFAMSCAN.out.ko_assignments.map { ['kofamscan', it] }
            )
        }

        if (params.eggnog_db) {
            EMAPPER(ch_proteins, file(params.eggnog_db))
            ch_annot_for_merge = ch_annot_for_merge.mix(
                EMAPPER.out.annotations.map { ['emapper', it] }
            )
        }

        if (params.dbcan_db) {
            DBCAN(ch_proteins, file(params.dbcan_db))
            ch_annot_for_merge = ch_annot_for_merge.mix(
                DBCAN.out.overview.map { ['dbcan', it] }
            )
        }
    }

    // 2f. Per-gene depths (requires annotation GFF + BAMs)
    if (effective_annotator != 'none' && params.bam_dir) {
        CALCULATE_GENE_DEPTHS(ch_bam_files, ch_gff)
    }

    // 3. Binners run in parallel
    ch_binner_results = Channel.empty()

    if (params.run_semibin && params.bam_dir) {
        BIN_SEMIBIN2(ch_assembly, ch_bam_files)
        ch_binner_results = ch_binner_results.mix(
            BIN_SEMIBIN2.out.bins.map { ['semibin', it] }
        )
    }

    // MetaBAT2 (depth-based, always runs)
    BIN_METABAT2(ch_assembly, ch_depths)
    ch_binner_results = ch_binner_results.mix(
        BIN_METABAT2.out.bins.map { ['metabat', it] }
    )

    if (params.run_maxbin) {
        BIN_MAXBIN2(ch_assembly, ch_depths)
        ch_binner_results = ch_binner_results.mix(
            BIN_MAXBIN2.out.bins.map { ['maxbin', it] }
        )
    }

    if (params.run_lorbin && params.bam_dir) {
        BIN_LORBIN(ch_assembly, ch_bam_files)
        ch_binner_results = ch_binner_results.mix(
            BIN_LORBIN.out.bins.map { ['lorbin', it] }
        )
    }

    if (params.run_comebin && params.bam_dir) {
        BIN_COMEBIN(ch_assembly, ch_bam_files)
        ch_binner_results = ch_binner_results.mix(
            BIN_COMEBIN.out.bins.map { ['comebin', it] }
        )
    }

    if (params.run_vamb) {
        BIN_VAMB(ch_assembly, ch_depths)
        ch_binner_results = ch_binner_results.mix(
            BIN_VAMB.out.bins.map { ['vamb', it] }
        )
    }

    if (params.run_vamb_tax && params.run_sendsketch && params.sendsketch_address) {
        BIN_VAMB_TAX(
            ch_assembly,
            ch_depths,
            SENDSKETCH_CLASSIFY.out.contig_taxonomy
        )
        ch_binner_results = ch_binner_results.mix(
            BIN_VAMB_TAX.out.bins.map { ['vamb_tax', it] }
        )
    }

    // 4. DAS Tool consensus
    ch_bin_labels = ch_binner_results.collect { it[0] }
    ch_bin_files  = ch_binner_results.collect { it[1] }

    DASTOOL_CONSENSUS(
        ch_assembly,
        ch_bin_files,
        ch_bin_labels
    )

    // 4a. Binette consensus refinement
    if (params.run_binette && params.checkm2_db) {
        ch_binner_fastas = BIN_METABAT2.out.fastas
        if (params.run_semibin && params.bam_dir) {
            ch_binner_fastas = ch_binner_fastas.mix(BIN_SEMIBIN2.out.fastas)
        }
        if (params.run_maxbin) {
            ch_binner_fastas = ch_binner_fastas.mix(BIN_MAXBIN2.out.fastas)
        }
        if (params.run_lorbin && params.bam_dir) {
            ch_binner_fastas = ch_binner_fastas.mix(BIN_LORBIN.out.fastas)
        }
        if (params.run_comebin && params.bam_dir) {
            ch_binner_fastas = ch_binner_fastas.mix(BIN_COMEBIN.out.fastas)
        }
        if (params.run_vamb) {
            ch_binner_fastas = ch_binner_fastas.mix(BIN_VAMB.out.fastas)
        }
        if (params.run_vamb_tax && params.run_sendsketch && params.sendsketch_address) {
            ch_binner_fastas = ch_binner_fastas.mix(BIN_VAMB_TAX.out.fastas)
        }
        BINETTE_CONSENSUS(ch_assembly, ch_binner_fastas.collect())
    }

    // 4b. MAGScoT consensus refinement
    if (params.run_magscot) {
        MAGSCOT_CONSENSUS(ch_assembly, ch_bin_files, ch_bin_labels)
    }

    // 4c. Metabolic merge + map to bins
    if (params.run_metabolism && effective_annotator != 'none') {
        ch_annot_labels = ch_annot_for_merge.collect { it[0] }
        ch_annot_files  = ch_annot_for_merge.collect { it[1] }

        MERGE_ANNOTATIONS(ch_annot_labels, ch_annot_files)

        ch_all_binner_bins = ch_binner_results.collect { it[1] }
        MAP_TO_BINS(MERGE_ANNOTATIONS.out.merged, DASTOOL_CONSENSUS.out.contig2bin, ch_all_binner_bins, ch_gff)

        KEGG_MODULES(MAP_TO_BINS.out.per_mag)
        MINPATH(MAP_TO_BINS.out.per_mag)
        KEGG_DECODER(MAP_TO_BINS.out.per_mag)

        // ECOSSDB ecosystem services
        if (params.run_ecossdb) {
            def ecossdb_db = "${projectDir}/ecossdb/db"
            ch_es_mapping  = Channel.fromPath("${ecossdb_db}/mappings/es_gene_mapping.tsv")
            ch_es_ontology = Channel.fromPath("${ecossdb_db}/ontology/cices_v5.2.tsv")

            ECOSSDB_MAP(
                MAP_TO_BINS.out.community,
                DASTOOL_CONSENSUS.out.contig2bin,
                ch_es_mapping,
                ch_es_ontology
            )

            ECOSSDB_SCORE(
                ECOSSDB_MAP.out.catalog,
                ch_es_mapping,
                params.ecossdb_role_weights
            )

            if (params.ecossdb_sdg) {
                ch_sdg_crosswalk = Channel.fromPath("${ecossdb_db}/ontology/sdg/cices_to_sdg.tsv")
                ch_sdg_targets   = Channel.fromPath("${ecossdb_db}/ontology/sdg/sdg_targets.tsv")
                ECOSSDB_SDG(
                    ECOSSDB_SCORE.out.scores,
                    ch_sdg_crosswalk,
                    ch_sdg_targets
                )
            }

            ch_es_hierarchy = Channel.fromPath("${ecossdb_db}/ontology/es_hierarchy.json")
            ECOSSDB_VIZ(
                ECOSSDB_SCORE.out.scores,
                ECOSSDB_MAP.out.catalog,
                ECOSSDB_MAP.out.mag_profiles,
                ch_es_hierarchy
            )
        }
    }

    // 4d. antiSMASH
    if (params.run_antismash) {
        def antismash_gff = ch_gff ?: file('NO_GFF')
        ANTISMASH(ch_assembly, antismash_gff)
    }

    // 5. Quality assessment with CheckM2
    if (params.checkm2_db) {
        ch_all_bins = BIN_METABAT2.out.fastas
        if (params.run_semibin && params.bam_dir) {
            ch_all_bins = ch_all_bins.mix(BIN_SEMIBIN2.out.fastas)
        }
        if (params.run_maxbin) {
            ch_all_bins = ch_all_bins.mix(BIN_MAXBIN2.out.fastas)
        }
        if (params.run_lorbin && params.bam_dir) {
            ch_all_bins = ch_all_bins.mix(BIN_LORBIN.out.fastas)
        }
        if (params.run_comebin && params.bam_dir) {
            ch_all_bins = ch_all_bins.mix(BIN_COMEBIN.out.fastas)
        }
        if (params.run_vamb) {
            ch_all_bins = ch_all_bins.mix(BIN_VAMB.out.fastas)
        }
        if (params.run_vamb_tax && params.run_sendsketch && params.sendsketch_address) {
            ch_all_bins = ch_all_bins.mix(BIN_VAMB_TAX.out.fastas)
        }
        ch_all_bins = ch_all_bins.mix(DASTOOL_CONSENSUS.out.bins)
        if (params.run_binette && params.checkm2_db) {
            ch_all_bins = ch_all_bins.mix(BINETTE_CONSENSUS.out.fastas)
        }
        if (params.run_magscot) {
            ch_all_bins = ch_all_bins.mix(MAGSCOT_CONSENSUS.out.fastas)
        }
        ch_all_bins = ch_all_bins.collect()

        CHECKM2(ch_all_bins)
    }

    // 6. GTDB-Tk phylogenetic classification
    if (params.run_gtdbtk && params.gtdbtk_db) {
        ch_gtdbtk_bins = BIN_METABAT2.out.fastas
        if (params.run_semibin && params.bam_dir) {
            ch_gtdbtk_bins = ch_gtdbtk_bins.mix(BIN_SEMIBIN2.out.fastas)
        }
        if (params.run_maxbin) {
            ch_gtdbtk_bins = ch_gtdbtk_bins.mix(BIN_MAXBIN2.out.fastas)
        }
        if (params.run_lorbin && params.bam_dir) {
            ch_gtdbtk_bins = ch_gtdbtk_bins.mix(BIN_LORBIN.out.fastas)
        }
        if (params.run_comebin && params.bam_dir) {
            ch_gtdbtk_bins = ch_gtdbtk_bins.mix(BIN_COMEBIN.out.fastas)
        }
        if (params.run_vamb) {
            ch_gtdbtk_bins = ch_gtdbtk_bins.mix(BIN_VAMB.out.fastas)
        }
        if (params.run_vamb_tax && params.run_sendsketch && params.sendsketch_address) {
            ch_gtdbtk_bins = ch_gtdbtk_bins.mix(BIN_VAMB_TAX.out.fastas)
        }
        ch_gtdbtk_bins = ch_gtdbtk_bins.mix(DASTOOL_CONSENSUS.out.bins)
        if (params.run_binette && params.checkm2_db) {
            ch_gtdbtk_bins = ch_gtdbtk_bins.mix(BINETTE_CONSENSUS.out.fastas)
        }
        if (params.run_magscot) {
            ch_gtdbtk_bins = ch_gtdbtk_bins.mix(MAGSCOT_CONSENSUS.out.fastas)
        }
        ch_gtdbtk_bins = ch_gtdbtk_bins.collect()

        GTDBTK_CLASSIFY(ch_gtdbtk_bins)
    }

    // 7. Viz dashboard
    if (params.run_viz) {
        // Stage 1: TNF done (or assembly-only snapshot)
        ch_viz1_input = ch_tnf.ifEmpty(ch_assembly).collect()
        VIZ_STAGE1(ch_viz1_input, false, false)

        // Stage 2: annotation done
        if (effective_annotator != 'none') {
            VIZ_STAGE2(ch_proteins.collect(), true, false)
        }

        // Stage 3: binning + quality done
        ch_viz_stage3 = DASTOOL_CONSENSUS.out.summary.collect()
        if (params.checkm2_db) {
            ch_viz_stage3 = ch_viz_stage3.mix(CHECKM2.out.report.collect())
        }
        VIZ_STAGE3(ch_viz_stage3.collect(), true, false)

        // Stage 4: final barrier
        ch_viz_stage4 = DASTOOL_CONSENSUS.out.summary.collect()
        if (params.checkm2_db) {
            ch_viz_stage4 = ch_viz_stage4.mix(CHECKM2.out.report.collect())
        }
        if (effective_annotator == 'bakta' && params.bakta_extra) {
            ch_viz_stage4 = ch_viz_stage4.mix(BAKTA_EXTRA.out.tsv.collect())
        }
        if (params.run_genomad && params.genomad_db) {
            ch_viz_stage4 = ch_viz_stage4.mix(GENOMAD_CLASSIFY.out.summary.collect())
        }
        if (params.run_integronfinder) {
            ch_viz_stage4 = ch_viz_stage4.mix(INTEGRONFINDER.out.integrons.collect())
        }
        if (params.run_islandpath && effective_annotator != 'none') {
            ch_viz_stage4 = ch_viz_stage4.mix(ISLANDPATH_DIMOB.out.islands.collect())
        }
        if (params.run_macsyfinder && effective_annotator != 'none' && params.macsyfinder_models) {
            ch_viz_stage4 = ch_viz_stage4.mix(MACSYFINDER.out.systems.collect())
        }
        if (params.run_defensefinder && effective_annotator != 'none') {
            ch_viz_stage4 = ch_viz_stage4.mix(DEFENSEFINDER.out.systems.collect())
        }
        if (params.run_eukaryotic) {
            ch_viz_stage4 = ch_viz_stage4.mix(WHOKARYOTE_CLASSIFY.out.classifications.collect())
            if (params.run_metaeuk && params.metaeuk_db && params.run_marferret && params.marferret_db) {
                ch_viz_stage4 = ch_viz_stage4.mix(MARFERRET_CLASSIFY.out.contigs.collect())
            }
        }
        if (params.run_kraken2 && params.kraken2_db) {
            ch_viz_stage4 = ch_viz_stage4.mix(KRAKEN2_CLASSIFY.out.contig_taxonomy.collect())
        }
        if (params.run_rrna && params.silva_ssu_db) {
            ch_viz_stage4 = ch_viz_stage4.mix(RNA_CLASSIFY.out.contig_classifications.collect())
        }
        if (params.run_sendsketch && params.sendsketch_address) {
            ch_viz_stage4 = ch_viz_stage4.mix(SENDSKETCH_CLASSIFY.out.contig_taxonomy.collect())
        }
        if (params.run_kaiju && params.kaiju_db) {
            ch_viz_stage4 = ch_viz_stage4.mix(KAIJU_CONTIG_CLASSIFY.out.contig_taxonomy.collect())
        }
        if (params.run_metabolism && effective_annotator != 'none') {
            ch_viz_stage4 = ch_viz_stage4.mix(KEGG_MODULES.out.modules.collect())
            ch_viz_stage4 = ch_viz_stage4.mix(MINPATH.out.pathways.collect())
            ch_viz_stage4 = ch_viz_stage4.mix(KEGG_DECODER.out.output.collect())
            if (params.run_ecossdb) {
                ch_viz_stage4 = ch_viz_stage4.mix(ECOSSDB_VIZ.out.json.collect())
            }
        }
        if (params.run_antismash) {
            ch_viz_stage4 = ch_viz_stage4.mix(ANTISMASH.out.summary.collect())
            ch_viz_stage4 = ch_viz_stage4.mix(ANTISMASH.out.json.collect())
        }
        if (params.run_gtdbtk && params.gtdbtk_db) {
            ch_viz_stage4 = ch_viz_stage4.mix(GTDBTK_CLASSIFY.out.taxonomy.collect())
        }
        VIZ_STAGE4(ch_viz_stage4.collect(), true, true)
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
