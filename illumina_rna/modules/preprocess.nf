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
    // Required whenever bbduk splits rRNA, so stored runs from before it existed redo QC
    tuple val(meta), path("${meta.id}.rrna_sub.fq.gz"),   emit: rrna, optional: params.rrna_method != 'bbduk'

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
        outu="${id}.qc.fq.gz" outm=rrna_all.fq.gz \\
        ordered t=${t_rrna} 2> rrna.log""" : ''}

    # Keep the first rRNA pairs for rRNA-read taxonomy (RRNA_TAXONOMY)
    if [ -s rrna_all.fq.gz ]; then
        reformat.sh in=rrna_all.fq.gz interleaved=t out="${id}.rrna_sub.fq.gz" reads=${params.rrna_tax_pairs} ow=t 2> /dev/null
        rm -f rrna_all.fq.gz
    fi

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
    // strobealign index of the (masked) human reference, built once.
    // (BBMap 40.05's index builder livelocks in IndexMaker4 in this image,
    // which took removehuman.sh down with it.)
    label 'process_map'
    label 'env_align'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    storeDir { params.store_dir ? "${params.store_dir}/references/human" : null }

    input:
    path(fasta, stageAs: 'src/*')

    output:
    tuple path("human.fa"), path("human.fa.r${params.read_len}.sti"), emit: index

    script:
    """
    case "${fasta}" in *.gz) gzip -dc "${fasta}" ;; *) cat "${fasta}" ;; esac \
        | awk '/^>/ { print \$1; next } { print }' > human.fa
    strobealign --create-index -t ${task.cpus} -r ${params.read_len} human.fa
    [ -s "human.fa.r${params.read_len}.sti" ] || { echo "[ERROR] human index build failed" >&2; exit 1; }
    """
}

process REMOVE_HUMAN {
    tag "${meta.id}"
    label 'process_map'
    label 'env_align'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir { "${params.outdir}/preprocess/${meta.id}" }, mode: 'symlink', pattern: '*.fq.gz', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null }

    input:
    tuple val(meta), path(reads)
    tuple path(human_fa), path(human_index)

    output:
    tuple val(meta), path("${meta.id}.nohuman.fq.gz"),       emit: reads
    tuple val(meta), path("${meta.id}.nohuman.counts.tsv"),  emit: counts

    script:
    """
    set -o pipefail
    # A mate is human at >=95% identity over >=80% of its length (the
    # stringency of removehuman.sh); a pair goes if either mate is.
    strobealign --use-index -t ${task.cpus} --interleaved "${human_fa}" "${reads}" 2> strobealign.log \
    | python3 "${projectDir}/bin/filter_human_sam.py" --counts human.counts \
    | gzip -1 > "${meta.id}.nohuman.fq.gz"

    read kept removed < human.counts
    printf '%s\tnohuman\t%s\n' "${meta.id}" "\$kept" > "${meta.id}.nohuman.counts.tsv"
    echo "[INFO] ${meta.id}: removed \$removed human reads, kept \$kept" >&2

    if [ "\$kept" -eq 0 ]; then
        cat strobealign.log >&2
        echo "[ERROR] Human removal kept no reads for ${meta.id}" >&2
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
