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
    path("gtdbtk_markers/"),     emit: markers

    script:
    """
    # Merge all bin directories into one flat directory
    mkdir -p all_bins gtdbtk_trees gtdbtk_markers
    for d in bins_*; do
        [ -d "\$d" ] || continue
        cp "\$d"/*.fa all_bins/ 2>/dev/null || true
    done

    # Bins far larger than any prokaryotic genome are catch-all clusters; their
    # classification means nothing and they are a large share of the sequence
    # (file size approximates length).
    for f in all_bins/*.fa; do
        [ -e "\$f" ] || continue
        if [ "\$(stat -Lc %s "\$f")" -gt ${params.gtdbtk_max_bin_bp} ]; then
            echo "[WARNING] GTDB-Tk: skipping \$(basename "\$f") (\$(( \$(stat -Lc %s "\$f") / 1000000 )) MB)" >&2
            rm -f "\$f"
        fi
    done

    # Binners often return the same contig set; classify each set once and copy
    # the result to the other names afterwards (duplicates.tsv: copy <tab> kept).
    : > duplicates.tsv
    declare -A first_with
    for f in all_bins/*.fa; do
        [ -e "\$f" ] || continue
        key=\$(grep '>' "\$f" | cut -f1 -d' ' | sort | md5sum | cut -c1-32)
        if [ -n "\${first_with[\$key]:-}" ]; then
            printf '%s\t%s\n' "\$(basename "\$f" .fa)" "\${first_with[\$key]}" >> duplicates.tsv
            rm -f "\$f"
        else
            first_with[\$key]=\$(basename "\$f" .fa)
        fi
    done
    [ -s duplicates.tsv ] && echo "[INFO] GTDB-Tk: \$(wc -l < duplicates.tsv) bins duplicate another bin's contigs; classified once" >&2

    if [ -z "\$(ls all_bins/*.fa 2>/dev/null)" ]; then
        echo "[WARNING] No bin FASTAs found -- skipping GTDB-Tk" >&2
        printf 'user_genome\\tclassification\\tclosest_genome_reference\\tclosest_genome_reference_radius\\tclosest_genome_taxonomy\\tclosest_genome_ani\\tclosest_genome_af\\tclosest_placement_reference\\tclosest_placement_radius\\tclosest_placement_taxonomy\\tclosest_placement_ani\\tclosest_placement_af\\tpplacer_taxonomy\\tclassification_method\\tnote\\tother_related_references\\tmsa_percent\\ttranslation_table\\tred_value\\twarnings\\n' > gtdbtk_taxonomy.tsv
        touch gtdbtk_markers/.empty
        exit 0
    fi

    export GTDBTK_DATA_PATH="${params.gtdbtk_db}"

    set +e
    gtdbtk classify_wf \\
        --genome_dir all_bins \\
        --out_dir gtdbtk_out \\
        --extension fa \\
        --cpus ${task.cpus} \\
        --pplacer_cpus ${params.gtdbtk_pplacer_cpus}
    gtdbtk_exit=\$?
    set -e

    # A failed run must not leave outputs: storeDir would keep an empty result
    # and skip GTDB-Tk on every resume. The task is recorded as failed and the
    # pipeline carries on (errorStrategy).
    if [ \$gtdbtk_exit -ne 0 ]; then
        echo "[WARNING] GTDB-Tk exited with code \$gtdbtk_exit; no results stored, it reruns on resume" >&2
        exit \$gtdbtk_exit
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
        printf 'user_genome\\tclassification\\tclosest_genome_reference\\tclosest_genome_reference_radius\\tclosest_genome_taxonomy\\tclosest_genome_ani\\tclosest_genome_af\\tclosest_placement_reference\\tclosest_placement_radius\\tclosest_placement_taxonomy\\tclosest_placement_ani\\tclosest_placement_af\\tpplacer_taxonomy\\tclassification_method\\tnote\\tother_related_references\\tmsa_percent\\ttranslation_table\\tred_value\\twarnings\\n' > gtdbtk_taxonomy.tsv
    fi

    # Rows for the bins that were classified as their duplicate: each kept
    # bin's row is repeated under the names of its copies.
    add_duplicate_rows() {
        [ -s duplicates.tsv ] && [ -s "\$1" ] || return 0
        awk -F'\\t' -v OFS='\\t' 'NR==FNR { copies[\$2] = copies[\$2] " " \$1; next }
            { print } FNR > 1 && (\$1 in copies) { n = split(copies[\$1], c, " "); for (i = 1; i <= n; i++) { \$1 = c[i]; print } }' \\
            duplicates.tsv "\$1" > "\$1.tmp" && mv "\$1.tmp" "\$1"
    }
    add_duplicate_rows gtdbtk_taxonomy.tsv

    # Copy placement tree files
    find gtdbtk_out -name '*.classify.tree*' -exec cp {} gtdbtk_trees/ \\; 2>/dev/null || true
    # Ensure directory is not empty (Nextflow needs at least one file)
    if [ -z "\$(ls gtdbtk_trees/ 2>/dev/null)" ]; then
        touch gtdbtk_trees/.empty
    fi

    # Copy markers_summary TSVs for GTDB-Tk quality metrics
    find gtdbtk_out -name 'gtdbtk.bac120.markers_summary.tsv' -exec cp {} gtdbtk_markers/ \\; 2>/dev/null || true
    find gtdbtk_out -name 'gtdbtk.ar53.markers_summary.tsv' -exec cp {} gtdbtk_markers/ \\; 2>/dev/null || true
    find gtdbtk_out -name 'gtdbtk.bac120.tree.mapping.tsv' -exec cp {} gtdbtk_markers/ \\; 2>/dev/null || true
    for f in gtdbtk_markers/*.markers_summary.tsv; do
        [ -e "\$f" ] && add_duplicate_rows "\$f"
    done
    if [ -z "\$(ls gtdbtk_markers/*.tsv 2>/dev/null)" ]; then
        touch gtdbtk_markers/.empty
    fi
    """
}
