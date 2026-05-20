// Quantification: featureCounts when GFF is provided, samtools idxstats otherwise.
// Per-(sample, reference) gene counts feed into the merged matrix in summarize.nf.

process FEATURECOUNTS {
    tag "${meta.id}_vs_${ref.name}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-illumina-rna-rnaseq"
    publishDir "${params.outdir}/quantify/${meta.id}/${ref.name}", mode: 'copy',
        pattern: '*.{tsv,summary}'
    storeDir params.store_dir ? "${params.store_dir}/quantify/${meta.id}/${ref.name}" : null

    input:
    tuple val(meta), val(ref), path(bam), path(bai), path(gff)

    output:
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.counts.tsv"),         emit: counts
    tuple val(meta), val(ref), path("${meta.id}_vs_${ref.name}.counts.tsv.summary"), emit: summary

    script:
    def tag = "${meta.id}_vs_${ref.name}"
    def strand_flag = params.strandedness == 'reverse' ? 2 : params.strandedness == 'forward' ? 1 : 0
    def pair_flag = '-p --countReadPairs'
    """
    if [ ! -s "${gff}" ]; then
        echo "[WARNING] Empty GFF for ${ref.name} — emitting empty counts" >&2
        printf "Geneid\\tChr\\tStart\\tEnd\\tStrand\\tLength\\t${meta.id}\\n" > "${tag}.counts.tsv"
        touch "${tag}.counts.tsv.summary"
        exit 0
    fi

    # featureCounts wants -t (feature) and -g (group attribute). Bakta/Prokka emit
    # CDS features with locus_tag; fall back to 'gene' / 'ID' for more general GFFs.
    feature_type="${params.feature_type}"
    attr_type="${params.attr_type}"

    featureCounts \\
        -a "${gff}" \\
        -F GFF \\
        -t "\${feature_type}" \\
        -g "\${attr_type}" \\
        ${pair_flag} \\
        -s ${strand_flag} \\
        -T ${task.cpus} \\
        -o "${tag}.counts.raw.tsv" \\
        "${bam}"

    # Slim down: keep Geneid…Length columns + one count column renamed to sample id
    awk -v sid="${meta.id}" '
        BEGIN { OFS="\\t" }
        /^#/ { print; next }
        NR==2 { \$NF = sid; print; next }
        { print }
    ' "${tag}.counts.raw.tsv" > "${tag}.counts.tsv"
    mv "${tag}.counts.raw.tsv.summary" "${tag}.counts.tsv.summary"
    """
}
