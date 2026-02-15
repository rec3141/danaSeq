// Eukaryotic contig classification: Tiara (deep learning k-mer NN) + Whokaryote (gene structure RF)
// Classify contigs; contigs still flow through existing binning pipeline unchanged

process TIARA_CLASSIFY {
    tag "tiara"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-tiara"
    publishDir "${params.outdir}/eukaryotic/tiara", mode: 'copy'

    input:
    path(contigs)

    output:
    path("tiara_output.tsv"),         emit: classifications
    path("tiara_probabilities.tsv"),  emit: probabilities, optional: true

    script:
    def min_len = params.tiara_min_len ?: 3000
    """
    set +e
    tiara \\
        -i "${contigs}" \\
        -o tiara_output.tsv \\
        --probabilities \\
        -t ${task.cpus} \\
        -m ${min_len}
    tiara_exit=\$?
    set -e

    if [ \$tiara_exit -ne 0 ]; then
        echo "[WARNING] Tiara exited with code \$tiara_exit" >&2
        printf 'sequence_id\\tclass_fst_stage\\tclass_snd_stage\\n' > tiara_output.tsv
    fi

    # Tiara writes probabilities to a file named after the output with .prob suffix
    # or --probabilities flag creates a separate file
    if [ -f "tiara_output.tsv.prob" ]; then
        mv tiara_output.tsv.prob tiara_probabilities.tsv
    elif [ -f "probabilities.tsv" ]; then
        mv probabilities.tsv tiara_probabilities.tsv
    fi

    n_classified=\$(tail -n +2 tiara_output.tsv | wc -l)
    echo "[INFO] Tiara: \${n_classified} contigs classified (min length: ${min_len} bp)" >&2
    """
}

process WHOKARYOTE_CLASSIFY {
    tag "whokaryote"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-whokaryote"
    publishDir "${params.outdir}/eukaryotic/whokaryote", mode: 'copy'

    input:
    path(contigs)
    path(gff)

    output:
    path("whokaryote_classifications.tsv"), emit: classifications

    script:
    def min_len = params.whokaryote_min_len ?: 5000
    def gff_arg = gff.name != 'NO_GFF' ? "--gff ${gff}" : ''
    """
    set +e
    whokaryote.py \\
        --contigs "${contigs}" \\
        --outdir whokaryote_out/ \\
        ${gff_arg} \\
        --minsize ${min_len} \\
        --threads ${task.cpus}
    whokaryote_exit=\$?
    set -e

    if [ \$whokaryote_exit -ne 0 ]; then
        echo "[WARNING] Whokaryote exited with code \$whokaryote_exit" >&2
        printf 'contig\\tprediction\\n' > whokaryote_classifications.tsv
    elif [ -f whokaryote_out/featuretable_predictions_T.tsv ]; then
        cp whokaryote_out/featuretable_predictions_T.tsv whokaryote_classifications.tsv
    else
        echo "[WARNING] Whokaryote produced no predictions" >&2
        printf 'contig\\tprediction\\n' > whokaryote_classifications.tsv
    fi

    n_classified=\$(tail -n +2 whokaryote_classifications.tsv | wc -l)
    echo "[INFO] Whokaryote: \${n_classified} contigs classified (min length: ${min_len} bp)" >&2
    """
}
