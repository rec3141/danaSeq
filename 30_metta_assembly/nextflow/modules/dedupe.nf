// Deduplication: cascade deduplication of multi-assembler contigs

process DEDUPE_ASSEMBLIES {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/assembly/${meta.id}", mode: 'copy', pattern: '*.{fasta,txt}'
    storeDir params.store_dir ? "${params.store_dir}/assembly/${meta.id}" : null

    input:
    tuple val(meta), path(contigs)

    output:
    tuple val(meta), path("${meta.id}.dedupe.fasta"), emit: assembly
    path("${meta.id}.assembly_stats.txt"),            emit: stats

    script:
    """
    # Combine all assembly FASTAs (skip empty files)
    > combined.fasta
    for f in ${contigs}; do
        if [ -s "\$f" ]; then
            cat "\$f" >> combined.fasta
        fi
    done

    if [ ! -s combined.fasta ]; then
        echo "[ERROR] No contigs from any assembler for ${meta.id}" >&2
        touch "${meta.id}.dedupe.fasta"
        echo "No assemblies produced" > "${meta.id}.assembly_stats.txt"
        exit 0
    fi

    # Cascade deduplication: 100% -> 99% -> final threshold
    dedupe.sh \\
        in=combined.fasta \\
        out=d100.fasta.gz \\
        sort=length uniquenames=t minidentity=100 \\
        t=${task.cpus} ow=t

    dedupe.sh \\
        in=d100.fasta.gz \\
        out=d99.fasta.gz \\
        sort=length uniquenames=t minidentity=99 \\
        t=${task.cpus} ow=t

    dedupe.sh \\
        in=d99.fasta.gz \\
        out="${meta.id}.dedupe.fasta" \\
        sort=length uniquenames=t minidentity=${params.dedupe_identity} \\
        t=${task.cpus} ow=t

    # Assembly statistics for all individual + deduplicated assemblies
    set +e
    statswrapper.sh \\
        ${contigs} "${meta.id}.dedupe.fasta" \\
        format=3 \\
        out="${meta.id}.assembly_stats.txt" \\
        t=${task.cpus} ow=t
    set -e

    if [ ! -s "${meta.id}.assembly_stats.txt" ]; then
        echo "Stats unavailable" > "${meta.id}.assembly_stats.txt"
    fi
    """
}
