// Metabolic profiling: KofamScan (KEGG Orthology), eggNOG-mapper (COG/GO/EC/Pfam),
// dbCAN3 (CAZyme detection), annotation merging, bin mapping, and KEGG module scoring.
//
// Architecture: "annotate once, map to bins" — all tools run on the full .faa,
// then annotations are partitioned by DAS_Tool contig2bin assignments.

process KOFAMSCAN {
    tag "kofamscan"
    label 'process_high'
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
    label 'process_high'
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
    path("overview.txt"),    emit: overview
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
        printf 'Gene ID\\tHMMER\\tdbCAN_sub\\tDIAMOND\\t#ofTools\\n' > overview.txt
        exit 0
    fi

    set +e
    run_dbcan \\
        "${proteins}" \\
        protein \\
        --db_dir "${dbcan_db}" \\
        --out_dir dbcan_out \\
        --dia_cpu ${task.cpus} \\
        --hmm_cpu ${task.cpus}
    dbcan_exit=\$?
    set -e

    if [ \$dbcan_exit -ne 0 ] || [ ! -d dbcan_out ]; then
        echo "[WARNING] dbCAN exited with code \$dbcan_exit" >&2
        printf 'Gene ID\\tHMMER\\tdbCAN_sub\\tDIAMOND\\t#ofTools\\n' > overview.txt
        exit 0
    fi

    # Copy outputs to working directory
    if [ -f dbcan_out/overview.txt ]; then
        cp dbcan_out/overview.txt .
    else
        printf 'Gene ID\\tHMMER\\tdbCAN_sub\\tDIAMOND\\t#ofTools\\n' > overview.txt
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
    path(kofamscan_tsv)
    path(emapper_tsv)
    path(dbcan_overview)

    output:
    path("merged_annotations.tsv"), emit: merged

    script:
    """
    # Merge KofamScan + eggNOG-mapper + dbCAN annotations into unified per-protein TSV
    # Joins on protein_id; outputs: protein_id, contig_id, KO, COG, GO, EC, Pfam, CAZy, desc

    merge_annotations.py \\
        --kofamscan "${kofamscan_tsv}" \\
        --emapper "${emapper_tsv}" \\
        --dbcan "${dbcan_overview}" \\
        -o merged_annotations.tsv
    """
}

process MAP_TO_BINS {
    tag "map_to_bins"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-kofamscan"
    publishDir "${params.outdir}/metabolism/per_mag", mode: 'copy', pattern: 'per_mag/*.tsv'
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
