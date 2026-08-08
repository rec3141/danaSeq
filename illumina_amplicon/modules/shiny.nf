// Build visualization data package from all pipeline outputs.
//
// R mode: builds app_data.rds + app.R (Shiny)
// Python mode: builds JSON files for Svelte app (samples.json, asvs.json.gz, etc.)
//
// When --build_viz_site is true, BUNDLE_VIZ_SITE builds a static Svelte app
// that can be opened directly in a browser (no server needed).

process BUILD_VIZ {
    tag "build_viz"
    label 'process_low'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/viz", mode: 'copy'

    input:
    path(seqtab)
    path(renorm)
    path(taxonomy_files)
    path(metadata)
    path(sample_tsne)
    path(seq_tsne)
    path(network)
    path(primer_assignment)

    output:
    path("*.json"),    emit: json, optional: true
    path("*.json.gz"), emit: json_gz, optional: true
    path("app_data.rds"), emit: app_data, optional: true
    path("app.R"),        emit: app_script, optional: true

    script:
    """
    build_viz.py \
        "${seqtab}" \
        "${renorm}" \
        "." \
        "${metadata}" \
        "${sample_tsne}" \
        "${seq_tsne}" \
        "${network}" \
        "${primer_assignment}"
    """
}

process BUNDLE_VIZ_SITE {
    tag "bundle_site"
    label 'process_low'
    publishDir "${params.outdir}/site", mode: 'copy'

    input:
    path(viz_json)
    path(viz_json_gz)

    output:
    path("dist/**"), emit: site

    script:
    """
    # Copy the viz app source
    cp -r ${projectDir}/viz viz_build
    cd viz_build
    npm install --prefer-offline 2>/dev/null || npm install
    # Copy data files into public/data for the static build
    mkdir -p public/data
    cp ${viz_json} public/data/ 2>/dev/null || true
    cp ${viz_json_gz} public/data/ 2>/dev/null || true
    npx vite build
    mv dist ..
    """
}
