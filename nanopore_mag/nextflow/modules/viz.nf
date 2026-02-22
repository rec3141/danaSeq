// Viz: preprocess pipeline results into JSON + optional static site build

process VIZ_PREPROCESS {
    tag "viz_preprocess"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-viz"
    publishDir "${params.outdir}/viz", mode: 'link'
    // No storeDir — always regenerate when pipeline runs

    input:
    val(ready)           // barrier signal, not used in script

    output:
    path("data"),  emit: json_data
    path("site"),  emit: static_site, optional: true

    script:
    def results = params.store_dir
        ? "${params.store_dir}/results_${file(params.outdir).name}"
        : "${params.outdir}"
    """
    # Generate JSON data files
    mkdir -p data
    python3 ${projectDir}/viz/preprocess/preprocess.py \
        --results "${results}" \
        --output data/

    # Build static site (node/npm provided by conda env)
    cd ${projectDir}/viz
    npm ci --prefer-offline 2>/dev/null || npm install --no-audit --no-fund
    # Copy preprocessed JSON into public/data for build
    cp \${OLDPWD}/data/* public/data/
    npm run build
    cp -r dist \${OLDPWD}/site
    """
}
