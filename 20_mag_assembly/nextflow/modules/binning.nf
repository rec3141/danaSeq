// Binning: SemiBin2, MetaBAT2, MaxBin2, DAS_Tool consensus

process BIN_SEMIBIN2 {
    tag "semibin2"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-semibin"
    publishDir "${params.outdir}/binning/semibin", mode: 'copy', saveAs: { fn -> fn == 'semibin_bins.tsv' ? 'contig_bins.tsv' : fn }

    input:
    path(assembly)
    path(bams)

    output:
    path("semibin_bins.tsv"), emit: bins
    path("bins/"),            emit: fastas

    script:
    """
    mkdir -p bins

    # SemiBin2 can crash on very small datasets (0 bins → empty ORFs → hmmsearch fail)
    # Catch failures and produce an empty output so the pipeline can continue
    set +e
    SemiBin2 single_easy_bin \\
        -i "${assembly}" \\
        -b *.sorted.bam \\
        -o semibin_out \\
        --sequencing-type long_read
    semibin_exit=\$?
    set -e

    if [ \$semibin_exit -ne 0 ]; then
        echo "[WARNING] SemiBin2 exited with code \$semibin_exit (dataset may be too small)" >&2
        touch semibin_bins.tsv
    elif [ -f semibin_out/contig_bins.tsv ]; then
        # SemiBin2 outputs a header line -- remove it for DAS_Tool compatibility
        tail -n +2 semibin_out/contig_bins.tsv | \
            awk -F'\\t' '{printf "%s\\tsemibin_%03d\\n", \$1, \$2+1}' > semibin_bins.tsv

        # Rename intermediate bin FASTAs to standardized names
        # SemiBin2 outputs gzipped bins (.fa.gz)
        bin_num=0
        for bin_file in semibin_out/output_bins/*.fa.gz; do
            [ -e "\$bin_file" ] || continue
            bin_num=\$((bin_num + 1))
            zcat "\$bin_file" > "bins/\$(printf 'semibin_%03d.fa' \$bin_num)"
        done
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
    publishDir "${params.outdir}/binning/metabat", mode: 'copy', saveAs: { fn -> fn == 'metabat_bins.tsv' ? 'contig_bins.tsv' : fn }

    input:
    path(assembly)
    path(jgi_depth)

    output:
    path("metabat_bins.tsv"), emit: bins
    path("bins/"),            emit: fastas

    script:
    """
    mkdir -p bins

    metabat2 \\
        -i "${assembly}" \\
        -o metabat_out/bin \\
        --saveCls \\
        --minClsSize ${params.metabat_min_cls} \\
        -a "${jgi_depth}"

    # Convert MetaBAT2 FASTA bins to contig_bins.tsv format
    > metabat_bins.tsv
    bin_num=0
    for bin_file in metabat_out/bin*.fa; do
        [ -e "\$bin_file" ] || continue
        bin_num=\$((bin_num + 1))
        bin_name=\$(printf 'metabat_%03d' \$bin_num)
        cp "\$bin_file" "bins/\${bin_name}.fa"
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
    publishDir "${params.outdir}/binning/maxbin", mode: 'copy', saveAs: { fn -> fn == 'maxbin_bins.tsv' ? 'contig_bins.tsv' : fn }

    input:
    path(assembly)
    path(jgi_depth)

    output:
    path("maxbin_bins.tsv"), emit: bins
    path("bins/"),           emit: fastas

    script:
    """
    mkdir -p bins

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
    bin_num=0
    for bin_file in maxbin_out/bin*.fasta; do
        [ -e "\$bin_file" ] || continue
        bin_num=\$((bin_num + 1))
        bin_name=\$(printf 'maxbin_%03d' \$bin_num)
        cp "\$bin_file" "bins/\${bin_name}.fa"
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
    publishDir "${params.outdir}/binning/dastool", mode: 'copy'

    input:
    path(assembly)
    path(bin_files)
    val(bin_labels)

    output:
    path("bins/"),             emit: bins
    path("contig2bin.tsv"),    emit: contig2bin
    path("allbins.fa"),       emit: allbins
    path("bin_quality.tsv"),  emit: bin_quality
    path("summary.tsv"),     emit: summary

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
        touch contig2bin.tsv allbins.fa bin_quality.tsv summary.tsv
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

    # ---- Collect DAS_Tool evaluation and summary files ----
    if [ -f dastool_out/dastool_allBins.eval ]; then
        cp dastool_out/dastool_allBins.eval bin_quality.tsv
    else
        touch bin_quality.tsv
    fi

    if [ -f dastool_out/dastool_DASTool_summary.tsv ]; then
        cp dastool_out/dastool_DASTool_summary.tsv summary.tsv
    else
        touch summary.tsv
    fi

    if [ -f dastool_out/dastool_DASTool_contig2bin.tsv ]; then
        cp dastool_out/dastool_DASTool_contig2bin.tsv contig2bin.tsv
    else
        touch contig2bin.tsv
    fi

    if [ -d dastool_out/dastool_DASTool_bins ]; then
        for f in dastool_out/dastool_DASTool_bins/*.fa; do
            [ -e "\$f" ] || continue
            cp "\$f" "bins/dastool-\$(basename "\$f")"
        done
    fi

    # Prefix bin names in contig2bin.tsv and summary.tsv so they match the
    # dastool- prefixed FASTA filenames
    if [ -s contig2bin.tsv ]; then
        sed -i 's/\\t/\\tdastool-/' contig2bin.tsv
    fi
    if [ -s summary.tsv ]; then
        # Skip header line, prefix bin name in first column
        sed -i '2,\$s/^/dastool-/' summary.tsv
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

process CHECKM2 {
    tag "checkm2"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-checkm2"
    publishDir "${params.outdir}/binning/checkm2", mode: 'copy'

    input:
    path(bins_dirs)   // collected list of bins/ directories

    output:
    path("quality_report.tsv"), emit: report

    script:
    """
    # Merge all bin directories into one flat directory
    # Nextflow stages collected path("bins/") as bins, bins_2, bins_3, ...
    mkdir -p all_bins
    for d in bins*; do
        [ -d "\$d" ] || continue
        cp "\$d"/*.fa all_bins/ 2>/dev/null || true
    done

    if [ -z "\$(ls all_bins/*.fa 2>/dev/null)" ]; then
        echo "[WARNING] No bin FASTAs found -- skipping CheckM2" >&2
        printf 'Name\\tCompleteness\\tContamination\\n' > quality_report.tsv
        exit 0
    fi

    checkm2 predict \\
        --threads ${task.cpus} \\
        --input all_bins \\
        --output-directory checkm2_out \\
        -x fa \\
        --database_path ${params.checkm2_db}

    cp checkm2_out/quality_report.tsv .
    """
}
