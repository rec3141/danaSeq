// Eukaryotic contig classification: Tiara (deep learning k-mer NN) + Whokaryote (gene structure RF)
// Phase 1: classify contigs and tag bins; contigs still flow through existing binning pipeline

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
    path("featuretable_predictions_T.tsv"), emit: classifications

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
        printf 'contig\\tprediction\\n' > featuretable_predictions_T.tsv
    elif [ -f whokaryote_out/featuretable_predictions_T.tsv ]; then
        cp whokaryote_out/featuretable_predictions_T.tsv .
    else
        echo "[WARNING] Whokaryote produced no featuretable_predictions_T.tsv" >&2
        printf 'contig\\tprediction\\n' > featuretable_predictions_T.tsv
    fi

    n_classified=\$(tail -n +2 featuretable_predictions_T.tsv | wc -l)
    echo "[INFO] Whokaryote: \${n_classified} contigs classified (min length: ${min_len} bp)" >&2
    """
}

process CLASSIFY_CONSENSUS {
    tag "euk_consensus"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-tiara"
    publishDir "${params.outdir}/eukaryotic/consensus", mode: 'copy'

    input:
    path(tiara_tsv)
    path(whokaryote_tsv)

    output:
    path("contig_classifications.tsv"), emit: classifications
    path("prokaryotic_contigs.txt"),    emit: prokaryotic_ids
    path("eukaryotic_contigs.txt"),     emit: eukaryotic_ids
    path("organellar_contigs.txt"),     emit: organellar_ids

    script:
    """
    classify_consensus.py \\
        --tiara "${tiara_tsv}" \\
        --whokaryote "${whokaryote_tsv}" \\
        --outdir .
    """
}

process TAG_BINS {
    tag "tag_bins"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-tiara"
    publishDir "${params.outdir}/eukaryotic/bin_tags", mode: 'copy'

    input:
    path(classifications)
    path(contigs)
    path(bin_files)
    val(bin_labels)

    output:
    path("all_bin_domain_tags.tsv"),  emit: all_tags
    path("*_bin_domain_tags.tsv"),    emit: per_binner_tags

    script:
    // Build binner:path pairs from collected inputs
    def files_list = bin_files instanceof List ? bin_files : [bin_files]
    def labels_list = bin_labels instanceof List ? bin_labels : [bin_labels]
    def pairs = []
    for (int i = 0; i < labels_list.size(); i++) {
        pairs.add("${labels_list[i]}:${files_list[i]}")
    }
    def bins_arg = pairs.join(',')
    """
    tag_bins.py \\
        --classifications "${classifications}" \\
        --contigs "${contigs}" \\
        --bins "${bins_arg}" \\
        --outdir .
    """
}
