// Viz: preprocess pipeline results into JSON + optional static site build

process VIZ_PREPROCESS {
    tag "viz_preprocess"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-viz"
    // No publishDir or storeDir — writes directly to outdir/viz/

    input:
    val(ready)           // barrier signal, not used in script
    val(skip_tsne)       // true = skip t-SNE (Stage 2+, embeddings already computed)
    val(skip_umap)       // true = skip UMAP (Stage 2+, embeddings already computed)

    output:
    val(true), emit: done

    script:
    // --results = publishDir root; --store-dir = storeDir root (resolve_path checks storeDir first)
    def storeRoot = params.store_dir ?: ''
    def storeFlag = storeRoot ? "--store-dir ${storeRoot}" : ""
    def tsne_flag = skip_tsne ? '--skip-tsne' : ''
    def umap_flag = skip_umap ? '--skip-umap' : ''
    def vizDir = "${params.outdir}/viz"
    """
    # Write JSON data directly to outdir
    VIZ_DIR="${vizDir}"
    mkdir -p "\${VIZ_DIR}/data"
    python3 ${projectDir}/viz/preprocess/preprocess.py \
        --results "${params.outdir}" \
        --output "\${VIZ_DIR}/data/" \
        ${storeFlag} \
        ${tsne_flag} \
        ${umap_flag}

    # Generate genes.json: prefer BAKTA_EXTRA > BAKTA_BASIC > PROKKA annotation,
    # augmented with barrnap rRNA and Aragorn tRNA gene calls.
    STORE="${storeRoot}"
    OUT="${params.outdir}"
    find_first() { for f in "\$@"; do [ -f "\${f}" ] && echo "\${f}" && return; done; }
    ANNOT_TSV=\$(find_first \
        "\${STORE:+\${STORE}/annotation/bakta/extra/annotation.tsv}" \
        "\${OUT}/annotation/bakta/extra/annotation.tsv" \
        "\${STORE:+\${STORE}/annotation/bakta/basic/annotation.tsv}" \
        "\${OUT}/annotation/bakta/basic/annotation.tsv" \
        "\${STORE:+\${STORE}/annotation/prokka/annotation.tsv}" \
        "\${OUT}/annotation/prokka/annotation.tsv")
    RRNA_TSV=\$(find_first \
        "\${STORE:+\${STORE}/taxonomy/rrna/rrna_genes.tsv}" \
        "\${OUT}/taxonomy/rrna/rrna_genes.tsv")
    TRNA_TSV=\$(find_first \
        "\${STORE:+\${STORE}/taxonomy/rrna/trna_genes.tsv}" \
        "\${OUT}/taxonomy/rrna/trna_genes.tsv")
    if [ -n "\${ANNOT_TSV}" ]; then
        python3 ${projectDir}/viz/preprocess/genes_to_json.py \
            "\${ANNOT_TSV}" "\${VIZ_DIR}/data/genes.json" "\${RRNA_TSV}" "\${TRNA_TSV}"
    else
        echo '{}' > "\${VIZ_DIR}/data/genes.json"
        echo '{}' | gzip > "\${VIZ_DIR}/data/genes.json.gz"
    fi

    # Build static site (node/npm provided by conda env)
    cd ${projectDir}/viz
    npm ci --prefer-offline 2>/dev/null || npm install --no-audit --no-fund
    # Copy preprocessed JSON into public/data for build
    cp "\${VIZ_DIR}"/data/* public/data/
    npm run build
    cp -r dist "\${VIZ_DIR}/site"

    # Start/restart the preview server (port 5174, all interfaces)
    pkill -f "vite preview" 2>/dev/null || true
    sleep 1
    nohup npm run serve > /tmp/vite_preview.log 2>&1 &
    disown \$!
    _server_ip=\$(hostname -I | awk '{print \$1}')
    echo "Viz dashboard: http://\${_server_ip}:5174/"
    """
}
