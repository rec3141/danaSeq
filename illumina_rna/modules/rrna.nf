// rRNA depletion. Metatranscriptomic libraries are often 50–90% ribosomal.
//
// Default (--rrna_method bbduk): QC_READS splits rRNA off with bbduk k=31
// against the reference prepared here, in the same pipe as trimming.
// Alternative (--rrna_method sortmerna): alignment-based SortMeRNA, much slower.

process PREPARE_RRNA_REF {
    // One gzipped DNA-alphabet FASTA from a file or a directory of FASTAs.
    // SILVA ships RNA sequences (U); bbduk k-mers need T.
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    storeDir { params.store_dir ? "${params.store_dir}/references/rrna" : null }

    input:
    path(refs)

    output:
    path("rrna_ref.fa.gz"), emit: ref

    script:
    """
    if [ -d "${refs}" ]; then
        find -L "${refs}/" -maxdepth 1 -type f \\( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' -o -name '*.fa.gz' -o -name '*.fasta.gz' \\) | sort > ref_list.txt
    else
        echo "${refs}" > ref_list.txt
    fi
    if [ ! -s ref_list.txt ]; then
        echo "[ERROR] no rRNA FASTAs found in ${refs}" >&2
        exit 1
    fi

    while read -r fa; do
        case "\$fa" in *.gz) gzip -dc "\$fa" ;; *) cat "\$fa" ;; esac
    done < ref_list.txt \\
    | awk '/^>/ { print \$1; next } { gsub(/[Uu]/, "T"); print toupper(\$0) }' \\
    | gzip -1 > rrna_ref.fa.gz

    n=\$(gzip -dc rrna_ref.fa.gz | grep -c '^>' || true)
    echo "[INFO] rRNA reference: \$n sequences from \$(wc -l < ref_list.txt) file(s)" >&2
    [ "\$n" -gt 0 ] || { echo "[ERROR] empty rRNA reference" >&2; exit 1; }
    """
}

process REMOVE_RRNA_SORTMERNA {
    tag "${meta.id}"
    label 'process_high'
    label 'env_rnaseq'
    conda "${projectDir}/conda-envs/dana-illumina-rna-rnaseq"
    publishDir { "${params.outdir}/preprocess/${meta.id}" }, mode: 'symlink', pattern: '*.{fq.gz,log}', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/preprocess/${meta.id}" : null }

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.norrna.fq.gz"),       emit: reads
    tuple val(meta), path("${meta.id}.norrna.counts.tsv"),  emit: counts
    path("${meta.id}.sortmerna.log"),                       emit: log

    script:
    """
    # SortMeRNA paired mode wants separate mates
    reformat.sh in="${reads}" interleaved=t out1=r1.fq.gz out2=r2.fq.gz ow=t

    ref_args=""
    if [ -d "${params.rrna_refs}" ]; then
        for fa in "${params.rrna_refs}"/*.fa "${params.rrna_refs}"/*.fasta; do
            [ -e "\$fa" ] || continue
            ref_args="\$ref_args --ref \$fa"
        done
    else
        ref_args="--ref ${params.rrna_refs}"
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

    reformat.sh in1=other_fwd.fq.gz in2=other_rev.fq.gz out="${meta.id}.norrna.fq.gz" ow=t 2> reformat.log
    n=\$(awk '/^Input:/ {print \$2; exit}' reformat.log)
    printf '%s\\tnorrna\\t%s\\n' "${meta.id}" "\${n:-0}" > "${meta.id}.norrna.counts.tsv"
    cp aligned.log "${meta.id}.sortmerna.log" 2>/dev/null || cp workdir/out/aligned.log "${meta.id}.sortmerna.log"

    if [ ! -s "${meta.id}.norrna.fq.gz" ]; then
        echo "[ERROR] SortMeRNA produced empty non-rRNA output for ${meta.id}" >&2
        exit 1
    fi
    """
}
