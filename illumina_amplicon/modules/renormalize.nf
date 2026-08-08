// Taxonomy-based renormalization.
//
// Separates ASVs into biologically meaningful groups (prokaryotes,
// eukaryotes, chloroplasts, mitochondria, unknown) based on their
// taxonomic assignment, then normalizes counts to proportions within
// each group.

def renormExt() { return 'pkl' }

process RENORMALIZE {
    tag "renormalize"
    label 'process_medium'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/renormalized", mode: 'copy'

    input:
    path(seqtab)
    tuple val(db_name), path(taxonomy), path(bootstrap)

    output:
    path("renorm_by_group.${renormExt()}"), emit: by_group
    path("renorm_merged.${renormExt()}"), emit: merged
    path("renorm_table_list.${renormExt()}"), emit: table_list
    path("renorm_stats.tsv"), emit: stats

    script:
    """
    renormalize.py "${seqtab}" "${taxonomy}" "${db_name}"
    """
}
