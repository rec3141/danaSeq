// Load sample metadata and merge with sequence data.

def metaExt() { return params.lang == 'python' ? 'pkl' : 'rds' }

process LOAD_METADATA {
    tag "metadata"
    label 'process_low'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"

    input:
    path(seqtab)
    path(metadata_file)
    val(sample_id_column)

    output:
    path("metadata.${metaExt()}"), emit: metadata
    path("match_stats.tsv"), emit: stats

    script:
    if (params.lang == 'python')
    """
    load_metadata.py "${seqtab}" "${metadata_file}" "${sample_id_column}"
    """
    else
    """
    load_metadata.R "${seqtab}" "${metadata_file}" "${sample_id_column}"
    """
}
