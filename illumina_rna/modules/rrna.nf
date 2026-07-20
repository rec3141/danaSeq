// rRNA depletion via SortMeRNA. Default-on for metatranscriptomic libraries where
// 50–90% of reads can be ribosomal. Falls back to passthrough if the SortMeRNA
// reference dir is missing.

process REMOVE_RRNA {
    tag "${meta.id}"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-illumina-rna-rnaseq"
    publishDir { "${params.outdir}/preprocess/${meta.id}" }, mode: 'copy', pattern: '*.{fq.gz,log}'
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.norrna.fq.gz"), emit: reads
    path("${meta.id}.rrna.fq.gz"),                    emit: rrna, optional: true
    path("${meta.id}.sortmerna.log"),                 emit: log

    script:
    """
    # Deinterleave for SortMeRNA paired-end mode
    reformat.sh in="${reads}" out1=r1.fq.gz out2=r2.fq.gz ow=t

    # Build the --ref args from rRNA FASTAs in the references dir
    ref_args=""
    for fa in "${params.sortmerna_refs}"/*.fa "${params.sortmerna_refs}"/*.fasta; do
        [ -e "\$fa" ] || continue
        ref_args="\$ref_args --ref \$fa"
    done

    if [ -z "\$ref_args" ]; then
        echo "[WARNING] No SortMeRNA references found in ${params.sortmerna_refs} — passing through" >&2
        cp "${reads}" "${meta.id}.norrna.fq.gz"
        echo "SortMeRNA skipped: no references found at ${params.sortmerna_refs}" > "${meta.id}.sortmerna.log"
        exit 0
    fi

    mkdir -p workdir
    sortmerna \$ref_args \\
        --reads r1.fq.gz --reads r2.fq.gz \\
        --workdir workdir \\
        --paired_in \\
        --fastx \\
        --aligned aligned \\
        --other other \\
        --threads ${task.cpus} \\
        --out2

    # Reinterleave the non-rRNA output
    reformat.sh in1=other_fwd.fq.gz in2=other_rev.fq.gz out="${meta.id}.norrna.fq.gz" ow=t
    if [ -s aligned_fwd.fq.gz ]; then
        reformat.sh in1=aligned_fwd.fq.gz in2=aligned_rev.fq.gz out="${meta.id}.rrna.fq.gz" ow=t
    fi
    cp workdir/out/aligned.log "${meta.id}.sortmerna.log"

    if [ ! -s "${meta.id}.norrna.fq.gz" ]; then
        echo "[ERROR] SortMeRNA produced empty non-rRNA output for ${meta.id}" >&2
        exit 1
    fi
    """
}
