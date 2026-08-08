// SparCC / CLR co-occurrence network analysis.
// R uses SpiecEasi::sparcc(); Python uses CLR + Pearson correlation.

def netExt() { return params.lang == 'python' ? 'pkl' : 'rds' }

process NETWORK_SPARCC {
    tag "sparcc"
    label 'process_high'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    publishDir "${params.outdir}/network", mode: 'copy'

    input:
    path(renorm_merged)
    val(min_prevalence)

    output:
    path("sparcc_correlations.${netExt()}"), emit: correlations
    path("sparcc_stats.tsv"), emit: stats

    script:
    if (params.lang == 'python')
    """
    network_sparcc.py "${renorm_merged}" ${min_prevalence} ${task.cpus}
    """
    else
    """
    network_sparcc.R "${renorm_merged}" ${min_prevalence} ${task.cpus}
    """
}
