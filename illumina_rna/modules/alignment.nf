// Reference-based mapping with BBmap. Each (sample, reference) pair is one task.
// References come from a user-supplied directory containing <name>.fasta files,
// typically the output of an upstream assembly or MAG run.

process BBMAP_INDEX {
    tag "${ref.name}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    storeDir params.store_dir ? "${params.store_dir}/references/${ref.name}" : null

    input:
    tuple val(ref), path(fasta)

    output:
    tuple val(ref), path("${ref.name}.fasta"), path("ref_${ref.name}"), emit: index

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    # Stage the FASTA with a canonical name so downstream tasks can find it
    cp "${fasta}" "${ref.name}.fasta"

    bbmap.sh ${xmx} \\
        ref="${ref.name}.fasta" \\
        path="ref_${ref.name}" \\
        t=${task.cpus}

    if [ ! -d "ref_${ref.name}/ref" ]; then
        echo "[ERROR] BBmap index build failed for ${ref.name}" >&2
        exit 1
    fi
    """
}

process MAP_READS_BBMAP {
    tag "${meta.id}_vs_${ref.name}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir "${params.outdir}/mapping/${meta.id}/${ref.name}", mode: 'copy',
        pattern: '*.{bam,bai,covhist,covstats,flagstat,idxstats}'
    storeDir params.store_dir ? "${params.store_dir}/mapping/${meta.id}/${ref.name}" : null

    input:
    tuple val(meta), path(reads), val(ref), path(ref_fasta), path(ref_index)

    output:
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.sorted.bam"),
                                path("${meta.id}_vs_${ref.name}.sorted.bam.bai"), emit: bam
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.idxstats.tsv"),    emit: idxstats
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.flagstat.txt"),    emit: flagstat
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.covstats.txt"),    emit: covstats
    path("${meta.id}_vs_${ref.name}.covhist.txt"),                                emit: covhist

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    def tag = "${meta.id}_vs_${ref.name}"
    """
    # Skip mapping if reference is empty
    if [ ! -s "${ref_fasta}" ] || ! grep -q '>' "${ref_fasta}"; then
        echo "[WARNING] Empty reference ${ref.name} — creating empty BAM" >&2
        echo -e '@HD\tVN:1.6\tSO:coordinate' | samtools view -b -o "${tag}.sorted.bam" -
        samtools index "${tag}.sorted.bam"
        touch "${tag}.idxstats.tsv" "${tag}.flagstat.txt" "${tag}.covstats.txt" "${tag}.covhist.txt"
        exit 0
    fi

    # OpenJDK 23 C2 JIT bug workaround (same as illumina_assembly mapping)
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -XX:-UseLoopPredicate"

    bbmap.sh ${xmx} \\
        in="${reads}" \\
        path="${ref_index}" \\
        out=stdout.sam \\
        maxindel=200 minid=${params.min_identity} \\
        qtrim=10 untrim \\
        ambig=all \\
        covhist="${tag}.covhist.txt" \\
        covstats="${tag}.covstats.txt" \\
        ow=t \\
        t=${task.cpus} \\
        | samtools sort -@ ${task.cpus} -o "${tag}.sorted.bam" -

    samtools index -@ ${task.cpus} "${tag}.sorted.bam"
    samtools flagstat "${tag}.sorted.bam" > "${tag}.flagstat.txt"
    samtools idxstats "${tag}.sorted.bam" > "${tag}.idxstats.tsv"

    if [ ! -s "${tag}.sorted.bam" ]; then
        echo "[ERROR] Mapping produced empty BAM for ${tag}" >&2
        exit 1
    fi
    """
}
