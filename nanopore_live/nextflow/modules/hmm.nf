// HMM search on Prokka-predicted proteins
// Runs hmmsearch with trusted cutoffs for each [sample, hmm_db] combination

process HMM_SEARCH {
    tag "${meta.id}_${dbname}"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-tools"
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/hmm", mode: 'copy'

    input:
    tuple val(meta), path(faa), val(dbname), path(hmm_db)

    output:
    tuple val(meta), path("${meta.id}.${dbname}.tsv"), emit: tsv
    tuple val(meta), path("${meta.id}.${dbname}.tbl"), emit: tbl

    script:
    """
    hmmsearch \
        --cut_tc \
        --cpu 1 \
        --tblout "${meta.id}.${dbname}.tbl" \
        "${hmm_db}" \
        "${faa}" \
        > /dev/null

    # Parse tblout to TSV
    awk 'BEGIN {OFS="\\t"; print "gene_id", "hmm_name", "evalue", "score"}
         !/^#/ && NF >= 5 {print \$1, \$3, \$5, \$6}' \
         "${meta.id}.${dbname}.tbl" > "${meta.id}.${dbname}.tsv"
    """
}
