// Reference-based mapping. Each (sample, reference) pair is one task.
//
// strobealign (default) indexes even a multi-Gbp co-assembly in minutes and
// maps short reads far faster than BBmap; the index is built once per
// reference (storeDir) and reused by every sample. BBmap stays available via
// --mapper bbmap.

process STROBEALIGN_INDEX {
    tag "${ref.name}"
    label 'process_map'
    label 'env_align'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    storeDir { params.store_dir ? "${params.store_dir}/references/${ref.name}" : null }

    input:
    tuple val(ref), path(fasta, stageAs: 'ref/*')

    output:
    tuple val(ref), path("${ref.name}.fasta"), path("${ref.name}.fasta.r${params.read_len}.sti"), emit: index

    script:
    """
    # Canonical name so the .sti sits next to the FASTA it was built from.
    # Keep only the first header word: BAM @SQ names must match the GFF seqids.
    awk '/^>/ { print \$1; next } { print }' ${fasta} > "${ref.name}.fasta"

    strobealign --create-index -t ${task.cpus} -r ${params.read_len} "${ref.name}.fasta"

    if [ ! -s "${ref.name}.fasta.r${params.read_len}.sti" ]; then
        echo "[ERROR] strobealign index build failed for ${ref.name}" >&2
        exit 1
    fi
    """
}

process MAP_READS_STROBEALIGN {
    tag "${meta.id}_vs_${ref.name}"
    label 'process_map'
    label 'env_align'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir { "${params.outdir}/mapping/${meta.id}/${ref.name}" }, mode: 'symlink',
        pattern: '*.{bam,flagstat.txt,log}', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/mapping/${meta.id}/${ref.name}" : null }

    input:
    tuple val(meta), path(reads), val(ref), path(ref_fasta), path(ref_index)

    output:
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.bam"),          emit: bam
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.flagstat.txt"), emit: flagstat
    path("${meta.id}_vs_${ref.name}.strobealign.log"),                         emit: log

    script:
    def tag = "${meta.id}_vs_${ref.name}"
    def sort_threads = Math.max(1, (task.cpus / 8) as int)
    """
    set -o pipefail

    # Name-grouped (unsorted) BAM: featureCounts -p reads mates adjacently,
    # and skipping the coordinate sort saves a full rewrite of every BAM.
    # -U drops unmapped reads; the QC counts already record the input.
    strobealign --use-index -t ${task.cpus} -U \\
        --rg-id="${meta.id}" --rg="SM:${meta.id}" \\
        --interleaved "${ref_fasta}" "${reads}" 2> "${tag}.strobealign.log" \\
    | samtools view -@ ${sort_threads} -b -o "${tag}.bam" -

    samtools flagstat -@ ${sort_threads} "${tag}.bam" > "${tag}.flagstat.txt"

    if [ ! -s "${tag}.bam" ]; then
        cat "${tag}.strobealign.log" >&2
        echo "[ERROR] Mapping produced empty BAM for ${tag}" >&2
        exit 1
    fi
    """
}

process BBMAP_INDEX {
    tag "${ref.name}"
    label 'process_map'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    storeDir { params.store_dir ? "${params.store_dir}/references/${ref.name}" : null }

    input:
    tuple val(ref), path(fasta, stageAs: 'ref/*')

    output:
    tuple val(ref), path("${ref.name}.fasta"), path("ref_${ref.name}"), emit: index

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    awk '/^>/ { print \$1; next } { print }' ${fasta} > "${ref.name}.fasta"

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
    label 'process_map'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    // Symlink works across filesystems (hardlinks fail between devices).
    publishDir { "${params.outdir}/mapping/${meta.id}/${ref.name}" }, mode: 'symlink',
        pattern: '*.{bam,bai,covhist.txt,covstats.txt,flagstat.txt,idxstats.tsv}', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/mapping/${meta.id}/${ref.name}" : null }

    input:
    tuple val(meta), path(reads), val(ref), path(ref_fasta), path(ref_index)

    output:
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.sorted.bam"),
                                path("${meta.id}_vs_${ref.name}.sorted.bam.bai"), emit: bam
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.flagstat.txt"),    emit: flagstat
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.covstats.txt"),    emit: covstats

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    def tag = "${meta.id}_vs_${ref.name}"
    """
    set -o pipefail
    # OpenJDK 23 C2 JIT bug workaround (same as illumina_assembly mapping)
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -XX:-UseLoopPredicate"

    bbmap.sh ${xmx} \\
        in="${reads}" interleaved=t \\
        path="${ref_index}" \\
        out=stdout.sam \\
        maxindel=200 minid=${params.min_identity} \\
        qtrim=10 untrim \\
        ambig=random \\
        covstats="${tag}.covstats.txt" \\
        ow=t \\
        t=${task.cpus} \\
    | samtools sort -@ ${task.cpus} -o "${tag}.sorted.bam" -

    samtools index -@ ${task.cpus} "${tag}.sorted.bam"
    samtools flagstat "${tag}.sorted.bam" > "${tag}.flagstat.txt"

    if [ ! -s "${tag}.sorted.bam" ]; then
        echo "[ERROR] Mapping produced empty BAM for ${tag}" >&2
        exit 1
    fi
    """
}
