// DuckDB integration: load pipeline outputs into dana.duckdb
// Runs R scripts that expect a specific directory layout with publishDir outputs.
// Uses val(barcode_dir) since R scripts operate on the published output directory
// (host filesystem), not the Nextflow work directory.
// maxForks=1 to prevent concurrent DuckDB access

process DB_INTEGRATION {
    tag "${barcode_dir}"
    label 'process_medium'
    conda 'conda-forge::r-base conda-forge::r-dbi conda-forge::r-duckdb conda-forge::r-readr'
    maxForks 1
    executor 'local'

    input:
    val barcode_dir

    output:
    val barcode_dir

    script:
    """
    # R scripts do setwd(args[1]) then open dana.duckdb in that directory
    # They expect kraken/, prokka/, sketch/, tetra/ subdirectories

    if [ -d "${barcode_dir}/kraken" ] && ls "${barcode_dir}"/kraken/*.tsv >/dev/null 2>&1; then
        Rscript ${params.danadir}/kraken-db.r "${barcode_dir}" || true
        Rscript ${params.danadir}/krakenreport-db.r "${barcode_dir}" || true
    fi

    if [ -d "${barcode_dir}/prokka" ]; then
        Rscript ${params.danadir}/prokka-db.r "${barcode_dir}" || true
    fi

    if [ -d "${barcode_dir}/sketch" ] && ls "${barcode_dir}"/sketch/*.txt >/dev/null 2>&1; then
        Rscript ${params.danadir}/sketch-db.r "${barcode_dir}" || true
    fi

    if [ -d "${barcode_dir}/tetra" ] && ls "${barcode_dir}"/tetra/*.lrn >/dev/null 2>&1; then
        Rscript ${params.danadir}/tetra-db.r "${barcode_dir}" || true
    fi
    """
}
