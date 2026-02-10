// Quality control: BBDuk adapter/quality trimming, Filtlong length filtering, FASTA conversion

process QC_BBDUK {
    tag "${meta.id}"
    label 'process_medium'
    conda 'bioconda::bbmap'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.bbduk.fastq.gz"), emit: trimmed

    script:
    """
    bbduk.sh \
        in="${fastq}" \
        out="${meta.id}.bbduk.fastq.gz" \
        ref=adapters,artifacts,phix,lambda \
        qtrim=rl trimq=15 entropy=0.75 qin=33 \
        minlength=${params.min_readlen}
    """
}

process QC_FILTLONG {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::filtlong'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.filtlong.fastq"), emit: filtered

    script:
    """
    filtlong \
        --min_length ${params.min_readlen} \
        --keep_percent ${params.keep_percent} \
        "${fastq}" > "${meta.id}.filtlong.fastq"

    # Fail if output is empty (all reads filtered out)
    if [ ! -s "${meta.id}.filtlong.fastq" ]; then
        echo "[WARNING] All reads filtered out for ${meta.id}" >&2
        exit 1
    fi
    """
}

process CONVERT_TO_FASTA {
    tag "${meta.id}"
    label 'process_low'
    conda 'bioconda::bbmap'
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/fa", mode: 'copy'

    input:
    tuple val(meta), path(fastq)

    output:
    tuple val(meta), path("${meta.id}.fa"), emit: fasta

    script:
    """
    reformat.sh \
        in="${fastq}" \
        out=stdout.fa \
        fastawrap=0 \
    | cut -f1 -d' ' > "${meta.id}.fa"

    if [ ! -s "${meta.id}.fa" ]; then
        echo "[WARNING] Empty FASTA output for ${meta.id}" >&2
        exit 1
    fi
    """
}
