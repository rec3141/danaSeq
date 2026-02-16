// Metabolic profiling: KofamScan (KEGG Orthology), eggNOG-mapper (COG/GO/EC/Pfam),
// dbCAN3 (CAZyme detection), annotation merging, bin mapping, and KEGG module scoring.
//
// Architecture: "annotate once, map to bins" — all tools run on the full .faa,
// then annotations are partitioned by DAS_Tool contig2bin assignments.

process KOFAMSCAN {
    tag "kofamscan"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-kofamscan"
    publishDir "${params.outdir}/metabolism/kofamscan", mode: 'copy'

    input:
    path(proteins)
    path(kofam_profiles)
    path(ko_list)

    output:
    path("kofamscan_results.tsv"), emit: ko_assignments

    script:
    """
    # KofamScan: assign KEGG Orthology (KO) numbers using adaptive score thresholds
    # --format detail-tsv: output with significance mark, gene, KO, threshold, score, E-value
    # Profiles are HMMs from KEGG, ko_list provides per-KO adaptive thresholds

    if [ ! -s "${proteins}" ]; then
        echo "[WARNING] No protein sequences — skipping KofamScan" >&2
        printf '# gene name\\tKO\\tthrshld\\tscore\\tE-value\\tKO definition\\n' > kofamscan_results.tsv
        exit 0
    fi

    set +e
    exec_annotation \\
        --profile "${kofam_profiles}" \\
        --ko-list "${ko_list}" \\
        --cpu ${task.cpus} \\
        --format detail-tsv \\
        -o kofamscan_results.tsv \\
        "${proteins}"
    kofam_exit=\$?
    set -e

    if [ \$kofam_exit -ne 0 ] || [ ! -f kofamscan_results.tsv ]; then
        echo "[WARNING] KofamScan exited with code \$kofam_exit" >&2
        printf '# gene name\\tKO\\tthrshld\\tscore\\tE-value\\tKO definition\\n' > kofamscan_results.tsv
        exit 0
    fi
    """
}

process EMAPPER {
    tag "emapper"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-emapper"
    publishDir "${params.outdir}/metabolism/emapper", mode: 'copy'

    input:
    path(proteins)
    path(eggnog_db)

    output:
    path("emapper_results.emapper.annotations"), emit: annotations

    script:
    """
    # eggNOG-mapper: functional annotation via DIAMOND search against eggNOG database
    # Assigns COG categories, GO terms, EC numbers, KEGG KOs, Pfam domains, descriptions
    # -m diamond: fast search mode (default for large datasets)
    # --override: overwrite existing output files

    if [ ! -s "${proteins}" ]; then
        echo "[WARNING] No protein sequences — skipping eggNOG-mapper" >&2
        printf '## eggNOG-mapper\\n#query\\tseed_ortholog\\tevalue\\tscore\\teggNOG_OGs\\tmax_annot_lvl\\tCOG_category\\tDescription\\tPreferred_name\\tGOs\\tEC\\tKEGG_ko\\tKEGG_Pathway\\tKEGG_Module\\tKEGG_Reaction\\tKEGG_rclass\\tBRITE\\tKEGG_TC\\tCAZy\\tBiGG_Reaction\\tPFAMs\\n' > emapper_results.emapper.annotations
        exit 0
    fi

    set +e
    emapper.py \\
        -i "${proteins}" \\
        --data_dir "${eggnog_db}" \\
        -m diamond \\
        --cpu ${task.cpus} \\
        --output emapper_results \\
        --override
    emapper_exit=\$?
    set -e

    if [ \$emapper_exit -ne 0 ] || [ ! -f emapper_results.emapper.annotations ]; then
        echo "[WARNING] eggNOG-mapper exited with code \$emapper_exit" >&2
        printf '## eggNOG-mapper\\n#query\\tseed_ortholog\\tevalue\\tscore\\teggNOG_OGs\\tmax_annot_lvl\\tCOG_category\\tDescription\\tPreferred_name\\tGOs\\tEC\\tKEGG_ko\\tKEGG_Pathway\\tKEGG_Module\\tKEGG_Reaction\\tKEGG_rclass\\tBRITE\\tKEGG_TC\\tCAZy\\tBiGG_Reaction\\tPFAMs\\n' > emapper_results.emapper.annotations
        exit 0
    fi
    """
}

process DBCAN {
    tag "dbcan"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-dbcan"
    publishDir "${params.outdir}/metabolism/dbcan", mode: 'copy'

    input:
    path(proteins)
    path(dbcan_db)

    output:
    path("overview.tsv"),    emit: overview
    path("hmmer.out"),       emit: hmmer, optional: true
    path("diamond.out"),     emit: diamond, optional: true
    path("substrate.out"),   emit: substrate, optional: true

    script:
    """
    # dbCAN3 (run_dbcan): CAZyme annotation using three methods:
    #   1. HMMER against dbCAN HMM database
    #   2. DIAMOND against CAZy sequence database
    #   3. dbCAN-sub (HMM of CAZyme subfamilies)
    # Consensus: CAZyme families agreed upon by >= 2/3 methods

    if [ ! -s "${proteins}" ]; then
        echo "[WARNING] No protein sequences — skipping dbCAN" >&2
        printf 'Gene ID\\tEC#\\tdbCAN_hmm\\tdbCAN_sub\\tDIAMOND\\t#ofTools\\tRecommend Results\\tSubstrate\\n' > overview.tsv
        exit 0
    fi

    set +e
    run_dbcan CAZyme_annotation \\
        --input_raw_data "${proteins}" \\
        --mode protein \\
        --db_dir "${dbcan_db}" \\
        --output_dir dbcan_out \\
        --threads ${task.cpus}
    dbcan_exit=\$?
    set -e

    if [ \$dbcan_exit -ne 0 ] || [ ! -d dbcan_out ]; then
        echo "[WARNING] dbCAN exited with code \$dbcan_exit" >&2
        printf 'Gene ID\\tEC#\\tdbCAN_hmm\\tdbCAN_sub\\tDIAMOND\\t#ofTools\\tRecommend Results\\tSubstrate\\n' > overview.tsv
        exit 0
    fi

    # Copy overview (run_dbcan v4+ outputs .tsv, older versions output .txt)
    if [ -f dbcan_out/overview.tsv ]; then
        cp dbcan_out/overview.tsv .
    elif [ -f dbcan_out/overview.txt ]; then
        cp dbcan_out/overview.txt overview.tsv
    else
        printf 'Gene ID\\tEC#\\tdbCAN_hmm\\tdbCAN_sub\\tDIAMOND\\t#ofTools\\tRecommend Results\\tSubstrate\\n' > overview.tsv
    fi

    # Optional detailed outputs
    [ -f dbcan_out/hmmer.out ] && cp dbcan_out/hmmer.out . || true
    [ -f dbcan_out/diamond.out ] && cp dbcan_out/diamond.out . || true
    [ -f dbcan_out/substrate.out ] && cp dbcan_out/substrate.out . || true
    """
}

process MERGE_ANNOTATIONS {
    tag "merge"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-kofamscan"
    publishDir "${params.outdir}/metabolism/merged", mode: 'copy'

    input:
    val(labels)
    path(annotation_files)

    output:
    path("merged_annotations.tsv"), emit: merged

    script:
    // Build flags dynamically based on which annotation tools ran
    def kofam_file   = labels.contains('kofamscan') ? annotation_files.find { it.name.contains('kofamscan') } : null
    def emapper_file = labels.contains('emapper')   ? annotation_files.find { it.name.contains('emapper') }   : null
    def dbcan_file   = labels.contains('dbcan')      ? annotation_files.find { it.name.contains('overview') }  : null
    def flags = []
    if (kofam_file)   flags << "--kofamscan ${kofam_file}"
    if (emapper_file) flags << "--emapper ${emapper_file}"
    if (dbcan_file)   flags << "--dbcan ${dbcan_file}"
    """
    # Merge annotation sources into unified per-protein TSV
    # Tools that ran: ${labels.join(', ')}

    merge_annotations.py ${flags.join(' ')} -o merged_annotations.tsv
    """
}

process MAP_TO_BINS {
    tag "map_to_bins"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-kofamscan"
    publishDir "${params.outdir}/metabolism", mode: 'copy', pattern: 'per_mag'
    publishDir "${params.outdir}/metabolism/community", mode: 'copy', pattern: 'community_annotations.tsv'

    input:
    path(merged)
    path(contig2bin)
    path(gff)

    output:
    path("per_mag/"),                     emit: per_mag
    path("community_annotations.tsv"),    emit: community

    script:
    """
    # Map merged annotations to individual MAGs via DAS_Tool contig2bin.tsv
    # Uses GFF for accurate gene->contig mapping, then contig->bin from DAS_Tool
    # Outputs: per-MAG TSVs + community-wide table with bin_id column

    map_annotations_to_bins.py \\
        --annotations "${merged}" \\
        --contig2bin "${contig2bin}" \\
        --gff "${gff}" \\
        --outdir per_mag \\
        --community community_annotations.tsv
    """
}

process KEGG_MODULES {
    tag "kegg_modules"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-kofamscan"
    publishDir "${params.outdir}/metabolism/modules", mode: 'copy'

    input:
    path(per_mag_dir)

    output:
    path("module_completeness.tsv"), emit: modules
    path("module_heatmap.svg"),      emit: heatmap, optional: true

    script:
    """
    # Evaluate ~150 KEGG module definitions against KO presence per MAG
    # Outputs completeness matrix (MAG x module) + clustered SVG heatmap
    # Covers: carbon fixation, nitrogen, sulfur, methane, ETC, vitamins, central carbon

    kegg_module_completeness.py \\
        --input "${per_mag_dir}" \\
        --output module_completeness.tsv \\
        --heatmap module_heatmap.svg
    """
}
