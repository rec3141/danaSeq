// Co-occurrence network: CLR transform (pseudocount 0.5) + Pearson correlation
// between ASVs, keeping |r| > 0.1. Not the SparCC algorithm — the process and
// file names are kept for compatibility with existing results and the viz.


process NETWORK_SPARCC {
    tag "sparcc"
    label 'process_high'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/network", mode: 'copy'

    input:
    path(renorm_merged)
    val(min_prevalence)

    output:
    path("sparcc_correlations.pkl"), emit: correlations
    path("sparcc_stats.tsv"), emit: stats

    script:
    """
    network_sparcc.py "${renorm_merged}" ${min_prevalence} ${task.cpus}
    """
}
