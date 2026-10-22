// Quantification with featureCounts: genes (from the reference GFF) and whole
// contigs (always). Both SAFs are built once per reference by REF_FEATURES, so
// each sample's counts share one row order, which lets the merge stream them.

process REF_FEATURES {
    tag "${ref.name}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-illumina-rna-bbmap"
    publishDir "${params.outdir}/references/${ref.name}", mode: 'copy', pattern: '*.{tsv.gz,tsv}', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/references/${ref.name}" : null }

    input:
    tuple val(ref), path(fasta, stageAs: 'ref/*'), path(gff, stageAs: 'gff/*')

    output:
    tuple val(ref), path("${ref.name}.genes.saf"), path("${ref.name}.contigs.saf"),            emit: saf
    tuple val(ref), path("${ref.name}.genes.tsv.gz"), path("${ref.name}.contig_lengths.tsv"), emit: annotation

    script:
    """
    # Contig lengths straight from the FASTA (first header word, as in the BAM)
    awk '/^>/ { if (n) print n "\\t" l; n = substr(\$1, 2); l = 0; next } { l += length(\$0) }
         END { if (n) print n "\\t" l }' ${fasta} > ref.fai

    gff_arg=""
    [ -s "${gff}" ] && gff_arg="--gff ${gff}"
    python3 "${projectDir}/bin/gff_to_saf.py" \\
        --fai ref.fai \$gff_arg \\
        --types "${params.feature_types}" \\
        --attr "${params.attr_type}" \\
        --out "${ref.name}"
    """
}

process FEATURECOUNTS {
    tag "${meta.id}_vs_${ref.name}"
    label 'process_medium'
    label 'env_rnaseq'
    conda "${projectDir}/conda-envs/dana-illumina-rna-rnaseq"
    publishDir { "${params.outdir}/quantify/${meta.id}/${ref.name}" }, mode: 'copy',
        pattern: '*.{tsv.gz,summary}', enabled: !params.store_dir
    storeDir { params.store_dir ? "${params.store_dir}/quantify/${meta.id}/${ref.name}" : null }

    input:
    tuple val(meta), val(ref), path(bam), path(genes_saf), path(contigs_saf)

    output:
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.genes.tsv.gz"),
                               path("${meta.id}_vs_${ref.name}.contigs.tsv.gz"),
                               path("${meta.id}_vs_${ref.name}.genes.summary"), emit: counts

    script:
    def tag = "${meta.id}_vs_${ref.name}"
    def strand_flag = params.strandedness == 'reverse' ? 2 : params.strandedness == 'forward' ? 1 : 0
    // Pairs count once; a pair overlapping two genes goes to the larger
    // overlap (prokaryotic genes often overlap by a few bp). Primary only.
    def common = "-F SAF -p --countReadPairs --primary -Q ${params.min_mapq} -s ${strand_flag} -T ${task.cpus}"
    """
    set -o pipefail
    slim() {  # Geneid + count column, renamed to the sample id
        awk -v sid="${meta.id}" 'BEGIN { OFS = "\\t" } /^#/ { next } !h { print "Geneid", sid; h = 1; next } { print \$1, \$NF }' "\$1" \\
        | gzip -1 > "\$2"
    }

    if [ \$(wc -l < "${genes_saf}") -gt 1 ]; then
        featureCounts ${common} -O --largestOverlap -a "${genes_saf}" -o genes.raw "${bam}"
        slim genes.raw "${tag}.genes.tsv.gz"
        mv genes.raw.summary "${tag}.genes.summary"
    else
        printf 'Geneid\\t%s\\n' "${meta.id}" | gzip > "${tag}.genes.tsv.gz"
        printf 'Status\\t%s\\n' "${meta.id}" > "${tag}.genes.summary"
    fi

    featureCounts ${common} -a "${contigs_saf}" -o contigs.raw "${bam}"
    slim contigs.raw "${tag}.contigs.tsv.gz"
    """
}
