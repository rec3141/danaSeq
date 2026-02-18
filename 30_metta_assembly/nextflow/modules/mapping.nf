// Read mapping: map preprocessed reads back to assembly, calculate depths

process MAP_READS_BBMAP {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-metta-bbmap"
    publishDir "${params.outdir}/mapping/${meta.id}", mode: 'copy', pattern: '*.{bam,bai,txt}'
    storeDir params.store_dir ? "${params.store_dir}/mapping/${meta.id}" : null

    input:
    tuple val(meta), path(reads), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: bam
    path("${meta.id}.covhist.txt"),  emit: covhist
    path("${meta.id}.covstats.txt"), emit: covstats

    script:
    """
    bbmap.sh \\
        in="${reads}" \\
        ref="${assembly}" \\
        out=stdout.sam \\
        nodisk \\
        maxindel=200 minid=90 \\
        qtrim=10 untrim \\
        ambig=all \\
        covhist="${meta.id}.covhist.txt" \\
        covstats="${meta.id}.covstats.txt" \\
        ow=t \\
        t=${task.cpus} \\
        | samtools sort -@ ${task.cpus} -o "${meta.id}.sorted.bam" -

    samtools index -@ ${task.cpus} "${meta.id}.sorted.bam"

    if [ ! -s "${meta.id}.sorted.bam" ]; then
        echo "[ERROR] Mapping produced empty BAM for ${meta.id}" >&2
        exit 1
    fi
    """
}

process CALCULATE_DEPTHS {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-metta-binning"
    publishDir "${params.outdir}/mapping/${meta.id}", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/mapping/${meta.id}" : null

    input:
    tuple val(meta), path(bams), path(assembly)

    output:
    tuple val(meta), path("${meta.id}.depths.txt"), emit: depths

    script:
    """
    jgi_summarize_bam_contig_depths \\
        --outputDepth "${meta.id}.depths.txt" \\
        *.sorted.bam

    if [ ! -s "${meta.id}.depths.txt" ]; then
        echo "[ERROR] Depth calculation produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}
