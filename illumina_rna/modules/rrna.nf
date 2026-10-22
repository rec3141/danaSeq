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

process PREPARE_RRNA_TAXDB {
    // SILVA SSU+LSU as a DNA FASTA keyed by accession, its lineages, and a
    // strobealign index, for classifying the rRNA reads QC_READS sets aside.
    label 'process_map'
    label 'env_align'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    storeDir { params.store_dir ? "${params.store_dir}/references/rrna_tax" : null }

    input:
    path(refs)

    output:
    tuple path("rrna_tax.fa"), path("rrna_tax.fa.r${params.read_len}.sti"), path("rrna_tax.lineages.tsv"), emit: db

    script:
    """
    if [ -d "${refs}" ]; then
        find -L "${refs}/" -maxdepth 1 -type f \\( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' -o -name '*.fa.gz' -o -name '*.fasta.gz' \\) | sort > ref_list.txt
    else
        echo "${refs}" > ref_list.txt
    fi
    while read -r fa; do
        case "\$fa" in *.gz) gzip -dc "\$fa" ;; *) cat "\$fa" ;; esac
    done < ref_list.txt \\
    | awk -v tax=rrna_tax.lineages.tsv '
        /^>/ { acc = substr(\$1, 2); lin = substr(\$0, length(\$1) + 2)
               if (seen[acc]++) { skip = 1; next } skip = 0
               print acc "\\t" lin > tax; print ">" acc; next }
        !skip { gsub(/[Uu]/, "T"); print toupper(\$0) }' > rrna_tax.fa

    n=\$(grep -c '^>' rrna_tax.fa || true)
    if ! awk -F'\\t' 'NF == 2 && \$2 ~ /;/ { ok++ } END { exit !(ok > 0) }' rrna_tax.lineages.tsv; then
        echo "[ERROR] no lineages in the rRNA FASTA headers; rRNA taxonomy needs SILVA-style headers (>ACC Domain;Phylum;...)" >&2
        exit 1
    fi
    echo "[INFO] rRNA taxonomy DB: \$n sequences" >&2
    strobealign --create-index -t ${task.cpus} -r ${params.read_len} rrna_tax.fa
    """
}

process RRNA_TAXONOMY {
    tag "${meta.id}"
    label 'process_map'
    label 'env_align'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir { "${params.outdir}/rrna_taxonomy/${meta.id}" }, mode: 'copy', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/rrna_taxonomy/${meta.id}" : null }

    input:
    tuple val(meta), path(reads), path(db_fa), path(db_index), path(lineages)

    output:
    tuple val(meta), path("${meta.id}.rrna_taxonomy.tsv"), emit: taxonomy

    script:
    """
    set -o pipefail
    # Best SILVA hit per read (both mates counted as reads); the classifier
    # applies the per-rank identity cutoffs.
    strobealign --use-index -t ${task.cpus} "${db_fa}" "${reads}" 2> strobealign.log \\
    | python3 "${projectDir}/bin/classify_rrna.py" --taxonomy "${lineages}" \\
        --sample "${meta.id}" --out "${meta.id}.rrna_taxonomy.tsv"
    """
}
