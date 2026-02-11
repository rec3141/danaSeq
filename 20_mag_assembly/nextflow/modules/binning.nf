// Binning: SemiBin2, MetaBAT2, MaxBin2, DAS_Tool consensus

process BIN_SEMIBIN2 {
    tag "semibin2"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-semibin"
    publishDir "${params.outdir}/binning/semibin", mode: 'copy', saveAs: { "contig_bins.tsv" }

    input:
    path(assembly)
    path(bams)

    output:
    path("semibin_bins.tsv"), emit: bins

    script:
    """
    # SemiBin2 can crash on very small datasets (0 bins → empty ORFs → hmmsearch fail)
    # Catch failures and produce an empty output so the pipeline can continue
    set +e
    SemiBin2 single_easy_bin \\
        -i "${assembly}" \\
        -b *.sorted.bam \\
        -o semibin_out
    semibin_exit=\$?
    set -e

    if [ \$semibin_exit -ne 0 ]; then
        echo "[WARNING] SemiBin2 exited with code \$semibin_exit (dataset may be too small)" >&2
        touch semibin_bins.tsv
    elif [ -f semibin_out/contig_bins.tsv ]; then
        # SemiBin2 outputs a header line -- remove it for DAS_Tool compatibility
        tail -n +2 semibin_out/contig_bins.tsv > semibin_bins.tsv
    else
        echo "[WARNING] SemiBin2 produced no contig_bins.tsv" >&2
        touch semibin_bins.tsv
    fi

    if [ ! -s semibin_bins.tsv ]; then
        echo "[WARNING] SemiBin2 produced no bins" >&2
    fi
    """
}

process BIN_METABAT2 {
    tag "metabat2"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-binning"
    publishDir "${params.outdir}/binning/metabat", mode: 'copy', saveAs: { "contig_bins.tsv" }

    input:
    path(assembly)
    path(jgi_depth)

    output:
    path("metabat_bins.tsv"), emit: bins

    script:
    """
    metabat2 \\
        -i "${assembly}" \\
        -o metabat_out/bin \\
        --saveCls \\
        --minClsSize ${params.metabat_min_cls} \\
        -a "${jgi_depth}"

    # Convert MetaBAT2 FASTA bins to contig_bins.tsv format
    > metabat_bins.tsv
    for bin_file in metabat_out/bin*.fa; do
        [ -e "\$bin_file" ] || continue
        bin_name=\$(basename "\$bin_file")
        grep '>' "\$bin_file" | tr -d '>' | cut -f1 -d' ' | while read -r contig; do
            printf '%s\\t%s\\n' "\$contig" "\$bin_name"
        done >> metabat_bins.tsv
    done

    if [ ! -s metabat_bins.tsv ]; then
        echo "[WARNING] MetaBAT2 produced no bins" >&2
    fi
    """
}

process BIN_MAXBIN2 {
    tag "maxbin2"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-binning"
    publishDir "${params.outdir}/binning/maxbin", mode: 'copy', saveAs: { "contig_bins.tsv" }

    input:
    path(assembly)
    path(jgi_depth)

    output:
    path("maxbin_bins.tsv"), emit: bins

    script:
    """
    # Extract coverage column from JGI depth table for MaxBin2 format
    # JGI format: contigName contigLen totalAvgDepth sample1.var sample2 ...
    # MaxBin2 wants: contigName avgDepth
    cut -f1,3 "${jgi_depth}" | tail -n +2 > coverage.txt

    mkdir -p maxbin_out
    run_MaxBin.pl \\
        -contig "${assembly}" \\
        -abund coverage.txt \\
        -out maxbin_out/bin \\
        -thread ${task.cpus}

    # Convert MaxBin2 FASTA bins to contig_bins.tsv format
    > maxbin_bins.tsv
    for bin_file in maxbin_out/bin*.fasta; do
        [ -e "\$bin_file" ] || continue
        bin_name=\$(basename "\$bin_file")
        grep '>' "\$bin_file" | tr -d '>' | cut -f1 -d' ' | while read -r contig; do
            printf '%s\\t%s\\n' "\$contig" "\$bin_name"
        done >> maxbin_bins.tsv
    done

    if [ ! -s maxbin_bins.tsv ]; then
        echo "[WARNING] MaxBin2 produced no bins" >&2
    fi
    """
}

process DASTOOL_CONSENSUS {
    tag "dastool"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-binning"
    publishDir "${params.outdir}/binning/consensus", mode: 'copy'

    input:
    path(assembly)
    path(bin_files)
    val(bin_labels)

    output:
    path("bins/"),           emit: bins
    path("contig2bin.tsv"),  emit: contig2bin
    path("scores.tsv"),     emit: scores
    path("allbins.fa"),     emit: allbins

    script:
    // Build comma-separated file and label lists from collected inputs
    // bin_files is a list of paths; bin_labels is a list of strings
    def files_csv = bin_files instanceof List ? bin_files.join(',') : bin_files
    def labels_csv = bin_labels instanceof List ? bin_labels.join(',') : bin_labels
    """
    # Filter out empty bin files (binners that produced 0 bins)
    FILTERED_FILES=""
    FILTERED_LABELS=""
    IFS=',' read -ra F_ARR <<< "${files_csv}"
    IFS=',' read -ra L_ARR <<< "${labels_csv}"
    for i in "\${!F_ARR[@]}"; do
        if [ -s "\${F_ARR[\$i]}" ]; then
            [ -n "\$FILTERED_FILES" ] && FILTERED_FILES+=","
            FILTERED_FILES+="\${F_ARR[\$i]}"
            [ -n "\$FILTERED_LABELS" ] && FILTERED_LABELS+=","
            FILTERED_LABELS+="\${L_ARR[\$i]}"
        else
            echo "[WARNING] Skipping empty binner output: \${L_ARR[\$i]}" >&2
        fi
    done

    if [ -z "\$FILTERED_FILES" ]; then
        echo "[WARNING] All binners produced empty output -- skipping DAS_Tool" >&2
        mkdir -p bins
        touch contig2bin.tsv scores.tsv allbins.fa
        exit 0
    fi

    # DAS_Tool exits 1 when no bins pass score_threshold (0.5 default)
    # This is expected with small datasets -- catch and produce empty output
    set +e
    DAS_Tool \\
        -i "\$FILTERED_FILES" \\
        -l "\$FILTERED_LABELS" \\
        -c "${assembly}" \\
        -o dastool_out/dastool \\
        --threads ${task.cpus} \\
        --write_bin_evals \\
        --write_bins
    dastool_exit=\$?
    set -e

    # Rename DAS_Tool verbose output to clean names
    mkdir -p bins

    if [ \$dastool_exit -ne 0 ]; then
        echo "[WARNING] DAS_Tool exited with code \$dastool_exit (no bins above score threshold)" >&2
    fi

    if [ -d dastool_out/dastool_DASTool_bins ]; then
        cp dastool_out/dastool_DASTool_bins/*.fa bins/ 2>/dev/null || true
    fi

    if [ -f dastool_out/dastool_DASTool_contig2bin.tsv ]; then
        cp dastool_out/dastool_DASTool_contig2bin.tsv contig2bin.tsv
    else
        touch contig2bin.tsv
    fi

    if [ -f dastool_out/dastool_DASTool_scores.tsv ]; then
        cp dastool_out/dastool_DASTool_scores.tsv scores.tsv
    else
        touch scores.tsv
    fi

    # Combine all bin FASTAs
    if ls bins/*.fa 1>/dev/null 2>&1; then
        cat bins/*.fa > allbins.fa
    else
        echo "[WARNING] DAS_Tool produced no consensus bins" >&2
        touch allbins.fa
    fi
    """
}
