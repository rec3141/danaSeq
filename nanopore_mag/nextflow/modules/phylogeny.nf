// Phylogenetic classification: GTDB-Tk classify_wf
//
// Assigns GTDB taxonomy to MAG bins using the GTDB-Tk classify workflow.
// Produces a merged taxonomy TSV and placement Newick trees.
// Follows the CHECKM2 pattern: collects all binner + consensus bin directories
// into a flat staging dir before running.

process GTDBTK_CLASSIFY {
    tag "gtdbtk"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-gtdbtk"
    publishDir "${params.outdir}/taxonomy/gtdbtk", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/taxonomy/gtdbtk" : null

    input:
    path(bins_dirs, stageAs: 'bins_?')   // collected list of bins/ directories

    output:
    path("gtdbtk_taxonomy.tsv"), emit: taxonomy
    path("gtdbtk_trees/"),       emit: trees

    script:
    """
    # Merge all bin directories into one flat directory
    mkdir -p all_bins gtdbtk_trees
    for d in bins_*; do
        [ -d "\$d" ] || continue
        cp "\$d"/*.fa all_bins/ 2>/dev/null || true
    done

    if [ -z "\$(ls all_bins/*.fa 2>/dev/null)" ]; then
        echo "[WARNING] No bin FASTAs found -- skipping GTDB-Tk" >&2
        printf 'user_genome\\tclassification\\tfastani_reference\\tfastani_reference_radius\\tfastani_taxonomy\\tfastani_ani\\tfastani_af\\tclosest_placement_reference\\tclosest_placement_radius\\tclosest_placement_taxonomy\\tclosest_placement_ani\\tclosest_placement_af\\tpplacer_taxonomy\\tclassification_method\\tnote\\tother_related_references\\tmsa_percent\\ttranslation_table\\tred_value\\twarnings\\n' > gtdbtk_taxonomy.tsv
        exit 0
    fi

    export GTDBTK_DATA_PATH="${params.gtdbtk_db}"

    set +e
    gtdbtk classify_wf \\
        --genome_dir all_bins \\
        --out_dir gtdbtk_out \\
        --extension fa \\
        --cpus ${task.cpus} \\
        --pplacer_cpus 1 \\
        --mash_db gtdbtk_mash_db
    gtdbtk_exit=\$?
    set -e

    if [ \$gtdbtk_exit -ne 0 ]; then
        echo "[WARNING] GTDB-Tk exited with code \$gtdbtk_exit" >&2
        printf 'user_genome\\tclassification\\tfastani_reference\\tfastani_reference_radius\\tfastani_taxonomy\\tfastani_ani\\tfastani_af\\tclosest_placement_reference\\tclosest_placement_radius\\tclosest_placement_taxonomy\\tclosest_placement_ani\\tclosest_placement_af\\tpplacer_taxonomy\\tclassification_method\\tnote\\tother_related_references\\tmsa_percent\\ttranslation_table\\tred_value\\twarnings\\n' > gtdbtk_taxonomy.tsv
        exit 0
    fi

    # Merge bacterial and archaeal summary TSVs
    BAC_SUMMARY=\$(find gtdbtk_out -name 'gtdbtk.bac120.summary.tsv' 2>/dev/null | head -1)
    AR_SUMMARY=\$(find gtdbtk_out -name 'gtdbtk.ar53.summary.tsv' 2>/dev/null | head -1)

    if [ -n "\$BAC_SUMMARY" ] && [ -s "\$BAC_SUMMARY" ]; then
        cp "\$BAC_SUMMARY" gtdbtk_taxonomy.tsv
        if [ -n "\$AR_SUMMARY" ] && [ -s "\$AR_SUMMARY" ]; then
            tail -n +2 "\$AR_SUMMARY" >> gtdbtk_taxonomy.tsv
        fi
    elif [ -n "\$AR_SUMMARY" ] && [ -s "\$AR_SUMMARY" ]; then
        cp "\$AR_SUMMARY" gtdbtk_taxonomy.tsv
    else
        echo "[WARNING] No GTDB-Tk summary files found" >&2
        printf 'user_genome\\tclassification\\tfastani_reference\\tfastani_reference_radius\\tfastani_taxonomy\\tfastani_ani\\tfastani_af\\tclosest_placement_reference\\tclosest_placement_radius\\tclosest_placement_taxonomy\\tclosest_placement_ani\\tclosest_placement_af\\tpplacer_taxonomy\\tclassification_method\\tnote\\tother_related_references\\tmsa_percent\\ttranslation_table\\tred_value\\twarnings\\n' > gtdbtk_taxonomy.tsv
    fi

    # Copy placement tree files
    find gtdbtk_out -name '*.classify.tree*' -exec cp {} gtdbtk_trees/ \\; 2>/dev/null || true
    # Ensure directory is not empty (Nextflow needs at least one file)
    if [ -z "\$(ls gtdbtk_trees/ 2>/dev/null)" ]; then
        touch gtdbtk_trees/.empty
    fi
    """
}
