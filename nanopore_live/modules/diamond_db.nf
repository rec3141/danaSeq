// DIAMOND blastp of annotation proteins (Bakta/Prokka) against a curated
// protein reference DB (e.g. HADEG — hydrocarbon/plastic degradation and
// biosurfactant genes). One task per [protein file x .dmnd database], stored
// beside the HMM results under <store>/<FC>/<BC>/protdb/.
//
// Cutoffs come from --protein_db_opts (default id=50,qcov=50,scov=50,
// evalue=1e-5: HADEG's recommended Proteinortho settings — identity 50%,
// Proteinortho's default 50% coverage and E 1e-5). Best hit per protein.
// Search + formatting is bin/protein_db_search.sh (same output as HMM_SEARCH).

process DIAMOND_DB_SEARCH {
    tag "${meta.id}_${dbname}"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-bakta"
    publishDir { "${params.outdir}/${meta.flowcell}/${meta.barcode}/protdb" }, mode: 'copy', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/${meta.flowcell}/${meta.barcode}/protdb" : null }

    input:
    tuple val(meta), path(faa), val(dbname), path(dmnd), val(annot), val(opts)

    output:
    tuple val(meta), path("${meta.id}.${dbname}.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}.${dbname}.tbl"), emit: tbl

    script:
    """
    protein_db_search.sh diamond "${dmnd}" "${dbname}" "${meta.id}.${dbname}" ${task.cpus} \\
        "${opts}" "${annot}" ${faa}
    """
}
