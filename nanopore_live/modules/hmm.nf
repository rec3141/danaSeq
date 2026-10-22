// Profile-HMM search of annotation proteins (Bakta/Prokka) against one HMM DB.
// Runs once per [protein file x HMM database]; the search + TSV formatting
// lives in bin/protein_db_search.sh so out-of-pipeline backfills produce
// identical output. Default cutoff is the models' trusted cutoff (--cut_tc),
// as recommended for CANT-HYD (Khot et al. 2022).
//
// Output <id>.<db>.tsv: '#'-provenance line (incl. n_proteins searched), header,
// then one row per hit annotated from <db>.genes.tsv (gene, category, pathway)
// when that catalogue sits next to the .hmm file. Loaded by import_all.py into
// the protein_hits / protein_db_annot / protein_db_searches tables.

process HMM_SEARCH {
    tag "${meta.id}_${dbname}"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-tools"
    publishDir { "${params.outdir}/${meta.flowcell}/${meta.barcode}/hmm" }, mode: 'copy', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/${meta.flowcell}/${meta.barcode}/hmm" : null }

    input:
    tuple val(meta), path(faa), val(dbname), path(hmm_db), val(annot), val(hmm_cut)

    output:
    tuple val(meta), path("${meta.id}.${dbname}.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}.${dbname}.tbl"), emit: tbl

    script:
    """
    protein_db_search.sh hmm "${hmm_db}" "${dbname}" "${meta.id}.${dbname}" ${task.cpus} \\
        "${hmm_cut}" "${annot}" ${faa}
    """
}
