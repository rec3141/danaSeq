// Binning: MetaBAT2 binning from assembly + depth matrix

process BIN_METABAT2 {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-metta-binning"
    publishDir "${params.outdir}/binning/${meta.id}", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/binning/${meta.id}" : null

    input:
    tuple val(meta), path(assembly), path(depths)

    output:
    tuple val(meta), path("${meta.id}.metabat_bins.tsv"), emit: bins
    tuple val(meta), path("bins/"),                       emit: fastas

    script:
    """
    mkdir -p bins

    set +e
    metabat2 \\
        -i "${assembly}" \\
        -a "${depths}" \\
        -o metabat_out/bin \\
        --saveCls \\
        -m ${params.metabat_min_cls} \\
        -t ${task.cpus}
    metabat_exit=\$?
    set -e

    > "${meta.id}.metabat_bins.tsv"

    if [ \$metabat_exit -ne 0 ]; then
        echo "[WARNING] MetaBAT2 exited with code \$metabat_exit for ${meta.id}" >&2
    else
        bin_num=0
        for bin_file in metabat_out/bin*.fa; do
            [ -e "\$bin_file" ] || continue
            case "\$(basename "\$bin_file")" in
                *unbinned*|*tooShort*|*lowDepth*) continue ;;
            esac
            bin_num=\$((bin_num + 1))
            bin_name=\$(printf 'metabat_%03d' \$bin_num)
            cp "\$bin_file" "bins/\${bin_name}.fa"
            grep '>' "\$bin_file" | tr -d '>' | cut -f1 -d' ' | while read -r contig; do
                printf '%s\\t%s\\n' "\$contig" "\$bin_name"
            done >> "${meta.id}.metabat_bins.tsv"
        done
    fi

    if [ ! -s "${meta.id}.metabat_bins.tsv" ]; then
        echo "[WARNING] MetaBAT2 produced no bins for ${meta.id}" >&2
    fi
    """
}
