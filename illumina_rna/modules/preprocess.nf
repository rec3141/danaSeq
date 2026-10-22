// Preprocessing: adapter trimming, artifact removal, rRNA split, human
// decontamination and QC reporting.
//
// QC_READS streams each sample through every bbduk stage in one pipe, so a
// sample is read once and written once however many filters run. Each stage's
// stderr summary supplies the read counts for the viz read-flow funnel.

process QC_READS {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir { "${params.outdir}/preprocess/${meta.id}" }, mode: 'symlink', pattern: '*.{fq.gz,log}', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null }

    input:
    tuple val(meta), path(r1), path(r2)
    path(rrna_ref)

    output:
    tuple val(meta), path("${meta.id}.qc.fq.gz"),         emit: reads
    tuple val(meta), path("${meta.id}.qc.counts.tsv"),    emit: counts
    path("${meta.id}.qc.log"),                            emit: log

    script:
    def mem_gb   = task.memory ? task.memory.toGiga() : 16
    def do_rrna  = params.rrna_method == 'bbduk'
    // Trim and filter use small references (adapters, artifacts, phiX); the
    // rRNA stage holds the SILVA k-mer table and gets the rest of the memory.
    def small_gb = do_rrna ? 4 : Math.max(4, (mem_gb * 0.4) as int)
    def rrna_gb  = Math.max(4, (mem_gb * 0.85 - 2 * small_gb) as int)
    def t_small  = Math.max(1, (task.cpus / 4) as int)
    def t_rrna   = Math.max(1, task.cpus - 2 * t_small)
    def cap      = params.max_read_pairs ? "reads=${params.max_read_pairs}" : ""
    def id       = meta.id
    """
    set -o pipefail
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -XX:-UseLoopPredicate"

    if [ "${params.run_filterbytile}" = "true" ]; then
        filterbytile.sh -Xmx${small_gb}g in1="${r1}" in2="${r2}" out=tiles.fq.gz ow=t t=${task.cpus} 2> filterbytile.log \\
            || { echo "[WARNING] filterbytile failed — using raw reads" >&2; rm -f tiles.fq.gz; }
    fi
    if [ -s tiles.fq.gz ]; then in_args="in=tiles.fq.gz interleaved=t"; else in_args="in1=${r1} in2=${r2}"; fi

    bbduk.sh -Xmx${small_gb}g \$in_args out=stdout.fq \\
        ktrim=r k=23 mink=11 hdist=1 tbo tpe \\
        minlen=${params.min_readlen} ref=adapters \\
        ftm=5 ordered ${cap} t=${t_small} 2> trim.log \\
    | bbduk.sh -Xmx${small_gb}g in=stdin.fq interleaved=t \\
        out=${do_rrna ? 'stdout.fq' : "${id}.qc.fq.gz"} \\
        k=31 ref=artifacts,phix entropy=0.95 \\
        ordered t=${t_small} 2> filter.log \\
    ${do_rrna ? """| bbduk.sh -Xmx${rrna_gb}g in=stdin.fq interleaved=t \\
        ref="${rrna_ref}" k=31 \\
        outu="${id}.qc.fq.gz" \\
        ordered t=${t_rrna} 2> rrna.log""" : ''}

    for f in trim.log filter.log ${do_rrna ? 'rrna.log' : ''}; do
        echo "==> \$f" >> "${id}.qc.log"; cat "\$f" >> "${id}.qc.log"
    done

    # bbduk prints "Input: N reads" and "Result: N reads" (reads, not pairs)
    input_n() { awk '/^Input:/ {print \$2; exit}' "\$1"; }
    result_n() { awk '/^Result:/ {print \$2; exit}' "\$1"; }
    {
        printf '%s\\traw\\t%s\\n'      "${id}" "\$(input_n trim.log)"
        printf '%s\\ttrimmed\\t%s\\n'  "${id}" "\$(result_n trim.log)"
        printf '%s\\tfiltered\\t%s\\n' "${id}" "\$(result_n filter.log)"
        if [ -f rrna.log ]; then
            printf '%s\\tnorrna\\t%s\\n' "${id}" "\$(result_n rrna.log)"
        fi
    } > "${id}.qc.counts.tsv"

    if ! awk -F'\\t' '\$3 == "" {bad=1} END {exit bad}' "${id}.qc.counts.tsv"; then
        echo "[ERROR] could not parse read counts for ${id}:" >&2; cat "${id}.qc.counts.tsv" >&2
        exit 1
    fi
    if [ ! -s "${id}.qc.fq.gz" ]; then
        echo "[ERROR] QC produced empty output for ${id}" >&2
        exit 1
    fi
    """
}

process HUMAN_INDEX {
    // removehuman.sh needs a BBTools index; build it once from a FASTA.
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    storeDir { params.store_dir ? "${params.store_dir}/references/human" : null }

    input:
    path(fasta)

    output:
    path("human_index"), emit: index

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    # k=14 + usemodulo match the index removehuman.sh maps with; any other
    # index makes every sample rebuild it (into this shared directory).
    bbmap.sh ${xmx} ref="${fasta}" path=human_index k=14 usemodulo t=${task.cpus} 2> human_index.log
    if [ ! -d human_index/ref ]; then
        cat human_index.log >&2
        echo "[ERROR] human index build failed" >&2
        exit 1
    fi
    """
}

process REMOVE_HUMAN {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir { "${params.outdir}/preprocess/${meta.id}" }, mode: 'symlink', pattern: '*.fq.gz', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null }

    input:
    tuple val(meta), path(reads)
    path(human_index)

    output:
    tuple val(meta), path("${meta.id}.nohuman.fq.gz"),       emit: reads
    tuple val(meta), path("${meta.id}.nohuman.counts.tsv"),  emit: counts

    script:
    def xmx = task.memory ? "-Xmx${(task.memory.toGiga() * 0.85).intValue()}g" : ""
    """
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -XX:-UseLoopPredicate"
    removehuman.sh ${xmx} \\
        in="${reads}" interleaved=t \\
        outu="${meta.id}.nohuman.fq.gz" \\
        path="${human_index}" \\
        t=${task.cpus} 2> removehuman.log

    # bbmap reports per-read mapping; unmapped reads (outu) are the non-human ones
    n=\$(reformat.sh in="${meta.id}.nohuman.fq.gz" interleaved=t 2>&1 | awk '/^Input:/ {print \$2; exit}')
    printf '%s\\tnohuman\\t%s\\n' "${meta.id}" "\${n:-0}" > "${meta.id}.nohuman.counts.tsv"

    if [ ! -s "${meta.id}.nohuman.fq.gz" ]; then
        cat removehuman.log >&2
        echo "[ERROR] Human removal produced empty output for ${meta.id}" >&2
        exit 1
    fi
    """
}

process FASTQC {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir { "${params.outdir}/preprocess/${meta.id}/fastqc" }, mode: 'copy', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}/fastqc" : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("*.html"), path("*.zip"), emit: reports

    script:
    """
    reformat.sh in="${reads}" interleaved=t out1=r1.fq.gz out2=r2.fq.gz ow=t

    fastqc -t ${task.cpus} --noextract -o . r1.fq.gz r2.fq.gz

    for f in *.html *.zip; do
        mv "\$f" "${meta.id}.\$f"
    done
    """
}
