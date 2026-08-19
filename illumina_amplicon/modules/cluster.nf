// t-SNE ordination of samples and ASVs via Bray-Curtis distances.
// Uses scipy + scikit-learn.


process CLUSTER_TSNE {
    tag "tsne"
    label 'process_high'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/clustering", mode: 'copy'

    input:
    path(seqtab)
    path(primer_assignment)

    output:
    path("sample_bray_tsne.pkl"), emit: sample_tsne
    path("seq_bray_tsne.pkl"), emit: seq_tsne
    path("sample_bray_dist.pkl"), emit: sample_dist
    path("seq_bray_dist.pkl"), emit: seq_dist

    script:
    // Which assay each sample belongs to. Samples amplified for different genes
    // share no sequence, so every distance between them is the maximum and an
    // ordination across them arranges the run on a circle; given this, each
    // assay is ordinated on its own. Absent for a run that removed no primers,
    // which is then ordinated as one.
    def assay_arg = primer_assignment.name != 'NO_PRIMER_ASSIGNMENT' ? "\"${primer_assignment}\"" : ''
    """
    cluster_tsne.py "${seqtab}" ${task.cpus} ${assay_arg}
    """
}
