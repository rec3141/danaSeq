// Cross-sample summaries: merged matrices, MAG/function aggregates and the
// viz JSONs.

process SUMMARIZE_EXPRESSION {
    tag "${ref_name}"
    label 'process_high'
    label 'env_rnaseq'
    conda "${projectDir}/conda-envs/dana-illumina-rna-rnaseq"
    publishDir "${params.outdir}/expression/${ref_name}", mode: 'copy', pattern: "${ref_name}.*"

    input:
    tuple val(ref_name),
          path(gene_counts,   stageAs: 'genes/*'),
          path(contig_counts, stageAs: 'contigs/*'),
          path(fc_summaries,  stageAs: 'fcsummary/*'),
          path(flagstats,     stageAs: 'flagstat/*'),
          path(gene_annot,    stageAs: 'ref/genes.tsv.gz'),
          path(contig_lengths, stageAs: 'ref/contig_lengths.tsv'),
          path(bins,          stageAs: 'ref/bins.tsv'),
          path(checkm2,       stageAs: 'ref/checkm2.tsv'),
          path(gtdbtk,        stageAs: 'ref/gtdbtk.tsv'),
          path(annotations,   stageAs: 'ref/annotations.tsv'),
          path(contig_explorer, stageAs: 'ref/contig_explorer'),
          path(kofam,         stageAs: 'ref/kofam.tsv'),
          path(read_counts,   stageAs: 'readcounts/*')

    output:
    tuple val(ref_name), path("${ref_name}.*"),       emit: tables
    tuple val(ref_name), path("viz/${ref_name}.*"),   emit: viz_tables

    script:
    """
    # KO definitions from KofamScan's detail table (one line per hit; the
    # marine run's is 20 GB), first definition per KO
    ko_arg=""
    if [ -s ref/kofam.tsv ]; then
        awk -F'\\t' '\$3 ~ /^K[0-9]+\$/ && !seen[\$3]++ { d = \$7; gsub(/^"|"\$/, "", d); print \$3 "\\t" d }' ref/kofam.tsv > ko_names.tsv
        ko_arg="--ko-names ko_names.tsv"
    fi

    python3 "${projectDir}/bin/summarize_expression.py" \\
        --ref "${ref_name}" \\
        --genes-dir genes --contigs-dir contigs \\
        --summary-dir fcsummary --flagstat-dir flagstat \\
        --readcounts-dir readcounts \\
        --gene-annot ref/genes.tsv.gz --contig-lengths ref/contig_lengths.tsv \\
        --bins ref/bins.tsv --checkm2 ref/checkm2.tsv --gtdbtk ref/gtdbtk.tsv \\
        --annotations ref/annotations.tsv \\
        --contig-explorer ref/contig_explorer \$ko_arg \\
        --top-genes ${params.viz_top_genes} --top-functions ${params.viz_top_functions} \\
        --out .
    """
}

process VIZ_PREPROCESS {
    label 'process_low'
    publishDir "${params.outdir}/viz", mode: 'copy'

    input:
    path(viz_tables, stageAs: 'tables/*')

    output:
    path("*.json*"), emit: json

    script:
    """
    python3 "${projectDir}/bin/build_viz_jsons.py" \\
        --tables tables \\
        --title "${params.viz_title ?: ''}" \\
        --mag-url "${params.viz_mag_url ?: ''}" \\
        --out .
    """
}
