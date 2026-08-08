// Taxonomy assignment — naive Bayesian kmer classifier.
//
// Each reference database runs as an independent task, enabling parallel
// classification across databases (custom Wang et al. 2007 implementation
// with multiprocessing).

def taxExt() { return 'pkl' }

process ASSIGN_TAXONOMY {
    tag "${db_name}"
    label 'process_high'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/taxonomy", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/taxonomy" : null

    input:
    tuple val(db_name), path(ref_db), val(tax_levels)
    path(seqtab)

    output:
    tuple val(db_name), path("${db_name}_taxonomy.${taxExt()}"), path("${db_name}_bootstrap.${taxExt()}"), emit: taxonomy
    path("${db_name}_taxonomy.tsv"), emit: taxonomy_tsv

    script:
    def levels_arg = tax_levels ? "\"${tax_levels}\"" : "null"
    def levels_cli = tax_levels ? "--tax-levels " + tax_levels.split(",").collect{ it.trim() }.join(" ") : ""
    """
    papa2 assign-taxonomy \
        "${seqtab}" "${ref_db}" "${db_name}" \
        ${levels_cli} \
        --threads ${task.cpus} \
        --verbose
    """
}
