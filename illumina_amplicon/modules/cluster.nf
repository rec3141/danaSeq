// t-SNE ordination of samples and ASVs via Bray-Curtis distances.
// Uses scipy + scikit-learn.


process CLUSTER_TSNE {
    tag "tsne"
    label 'process_high'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/clustering", mode: 'copy'

    input:
    path(seqtab)

    output:
    path("sample_bray_tsne.pkl"), emit: sample_tsne
    path("seq_bray_tsne.pkl"), emit: seq_tsne
    path("sample_bray_dist.pkl"), emit: sample_dist
    path("seq_bray_dist.pkl"), emit: seq_dist

    script:
    """
    cluster_tsne.py "${seqtab}" ${task.cpus}
    """
}
