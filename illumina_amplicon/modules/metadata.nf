// Load sample metadata and merge with sequence data.


process LOAD_METADATA {
    tag "metadata"
    label 'process_low'
    conda "${projectDir}/envs/python.yml"

    input:
    path(seqtab)
    path(metadata_file)
    val(sample_id_column)

    output:
    path("metadata.pkl"), emit: metadata
    path("match_stats.tsv"), emit: stats

    script:
    """
    load_metadata.py "${seqtab}" "${metadata_file}" "${sample_id_column}"
    """
}
