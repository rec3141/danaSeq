// t-SNE ordination of samples and ASVs via Bray-Curtis distances.
// Supports --lang R (parallelDist + Rtsne) or --lang python (scipy + sklearn).

def clusterExt() { return params.lang == 'python' ? 'pkl' : 'rds' }

process CLUSTER_TSNE {
    tag "tsne"
    label 'process_high'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    publishDir "${params.outdir}/clustering", mode: 'copy'

    input:
    path(seqtab)

    output:
    path("sample_bray_tsne.${clusterExt()}"), emit: sample_tsne
    path("seq_bray_tsne.${clusterExt()}"), emit: seq_tsne
    path("sample_bray_dist.${clusterExt()}"), emit: sample_dist
    path("seq_bray_dist.${clusterExt()}"), emit: seq_dist

    script:
    if (params.lang == 'python')
    """
    cluster_tsne.py "${seqtab}" ${task.cpus}
    """
    else
    """
    cluster_tsne.R "${seqtab}" ${task.cpus}
    """
}
