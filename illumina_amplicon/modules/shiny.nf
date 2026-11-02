// Build visualization data package from all pipeline outputs.
//
// BUILD_VIZ writes the JSON files the Svelte app loads (samples.json,
// asvs.json.gz, etc.). The file name is historical: the Shiny/R app it once
// built is gone, and nothing writes app_data.rds or app.R any more.
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

    script:
    // The primary reference — first named in --ref_databases, as RENORMALIZE and
    // the app take it — so the per-ASV lineage strings come from the same
    // database and not from whichever taxonomy file sorts first.
    def primary_db = (params.ref_databases
        ? params.ref_databases.tokenize(';')[0].tokenize(':')[0].trim() : '')
    """
    PRIMARY_DB="${primary_db}" build_viz.py \
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
    path(tree)          // phylo_tree.nwk from BUILD_PHYLOGENY, or the NO_TREE sentinel

    output:
    path("dist/**"), emit: site

    script:
    """
    # Copy the viz app source
    cp -r ${projectDir}/viz viz_build
    cd viz_build
    npm install --prefer-offline 2>/dev/null || npm install
    # Copy data files into public/data for the static build. The inputs are
    # staged in the task dir, one level up from viz_build.
    mkdir -p public/data
    for f in ${viz_json} ${viz_json_gz}; do cp "../\$f" public/data/; done
    # The app fetches data/tree.nwk for its Phylogeny view.
    [ "${tree}" = NO_TREE ] || cp ../${tree} public/data/tree.nwk
    npx vite build
    mv dist ..
    """
}
