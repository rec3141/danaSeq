// Sendsketch taxonomic profiling via k-mer sketching against NCBI nt

process SENDSKETCH {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-bbmap"
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/sketch", mode: 'copy'

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.txt"), emit: sketch

    script:
    def mem_mb = task.memory ? (task.memory.toMega() * 85 / 100).intValue() : 3400
    """
    sendsketch.sh \
        -Xmx${mem_mb}m \
        in="${fasta}" \
        address=nt \
        out="${meta.id}.txt" \
        format=3
    """
}
