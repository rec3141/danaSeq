// Viz: preprocess pipeline results into JSON + optional static site build

process VIZ_PREPROCESS {
    tag "viz_preprocess"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-pathviz"
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
    def vizDir = params.store_dir ? "${params.store_dir}/viz" : "${params.outdir}/viz"
    def vizPort = params.viz_port ?: 5174
    """
    # Ensure the viz conda env's Python is used (has numpy, pandas, scipy, etc.)
    # In container mode, conda is disabled and the wrapper may not be on PATH
    if [ -d /opt/conda/envs/dana-mag-pathviz/bin ]; then
        export PATH="/opt/conda/envs/dana-mag-pathviz/bin:\$PATH"
    fi

    # Copy Nextflow log for pipeline status parsing
    cp "${projectDir}/.nextflow.log" "${vizDir}/../pipeline_info/nextflow.log" 2>/dev/null || true

    # Write JSON data directly to outdir
    VIZ_DIR="${vizDir}"
    mkdir -p "\${VIZ_DIR}/data"
    python3 ${projectDir}/viz/preprocess/preprocess.py \
        --results "${params.store_dir ?: params.outdir}" \
        --output "\${VIZ_DIR}/data/" \
        ${storeFlag} \
        ${tsne_flag} \
        ${umap_flag}

    # Generate genes.json: prefer BAKTA_EXTRA > BAKTA_BASIC > PROKKA annotation,
    # augmented with barrnap rRNA and Aragorn tRNA gene calls.
    STORE="${storeRoot}"
    OUT="${params.outdir}"
    find_first() { for f in "\$@"; do [ -f "\${f}" ] && echo "\${f}" && return 0; done; return 0; }
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
    GENE_DEPTHS=\$(find_first \
        "\${STORE:+\${STORE}/mapping/gene_depths.tsv}" \
        "\${OUT}/mapping/gene_depths.tsv")
    ASSEMBLY=\$(find_first \
        "\${STORE:+\${STORE}/assembly/assembly.fasta}" \
        "\${OUT}/assembly/assembly.fasta")
    if [ -n "\${ANNOT_TSV}" ]; then
        python3 ${projectDir}/viz/preprocess/genes_to_json.py \
            "\${ANNOT_TSV}" "\${VIZ_DIR}/data/genes.json" "\${RRNA_TSV}" "\${TRNA_TSV}" "\${GENE_DEPTHS}" "\${ASSEMBLY}"
    else
        echo '{}' > "\${VIZ_DIR}/data/genes.json"
        echo '{}' | gzip > "\${VIZ_DIR}/data/genes.json.gz"
    fi

    # Copy ECOSSDB ecosystem services data if available
    ES_JSON=\$(find_first \
        "\${STORE:+\${STORE}/metabolism/ecossdb/ecosystem_services.json}" \
        "\${OUT}/metabolism/ecossdb/ecosystem_services.json")
    ES_SDG=\$(find_first \
        "\${STORE:+\${STORE}/metabolism/ecossdb_sdg/es_sdg.json}" \
        "\${OUT}/metabolism/ecossdb/sdg/es_sdg.json")
    [ -n "\${ES_JSON}" ] && cp "\${ES_JSON}" "\${VIZ_DIR}/data/ecosystem_services.json" && echo "Copied ecosystem_services.json"
    [ -n "\${ES_SDG}" ] && cp "\${ES_SDG}" "\${VIZ_DIR}/data/es_sdg.json" && echo "Copied es_sdg.json"

    # Build static site (node/npm provided by conda env)
    # Copy viz source to a writable location (container filesystem may be read-only)
    VIZ_BUILD="\${VIZ_DIR}/.build"
    # Kill any running vite server before cleaning .build — it holds node_modules open
    VIZ_PORT=${vizPort}
    pkill -f "vite preview.*--port \${VIZ_PORT}" 2>/dev/null || true
    sleep 1
    rm -rf "\${VIZ_BUILD}"
    cp -r ${projectDir}/viz "\${VIZ_BUILD}"
    cd "\${VIZ_BUILD}"
    npm ci --prefer-offline 2>/dev/null || npm install --no-audit --no-fund
    # Copy preprocessed JSON into public/data for build
    mkdir -p public/data
    cp "\${VIZ_DIR}"/data/* public/data/
    npm run build
    cp -r dist "\${VIZ_DIR}/site"
    # Copy pipeline_info (report, timeline, trace) into site for dashboard linking
    cp -r "${params.outdir}/pipeline_info" "\${VIZ_DIR}/site/pipeline_info" 2>/dev/null || true

    # Start/restart the preview server (all interfaces, configurable port).
    # setsid creates a new session so the child is not in the Nextflow process group.
    # Fds 0-2 are redirected to /dev/null and a log file.  Fds 3-9 are explicitly
    # closed because Nextflow's .command.run wraps every task in a tee pipeline:
    #   (nxf_launch | tee .command.out) 3>&1 1>&2 2>&3 | tee .command.err
    # That creates fd 3 holding the outer tee pipe's write end.  Without closing
    # it, tee never gets EOF → the subshell never exits → .exitcode is never
    # written → the Nextflow task hangs indefinitely.
    cd "\${VIZ_BUILD}"
    setsid nohup npm run serve -- --host 0.0.0.0 --port \${VIZ_PORT} \
        </dev/null >/tmp/vite_preview_\${VIZ_PORT}.log 2>&1 \
        3>&- 4>&- 5>&- 6>&- 7>&- 8>&- 9>&- &
    _server_ip=\$(hostname -I | awk '{print \$1}')
    echo "Viz dashboard: http://\${_server_ip}:\${VIZ_PORT}/"

    # Write a convenience script to restart the server without re-running the pipeline
    cat > "\${VIZ_DIR}/start-server.sh" << 'LAUNCHER'
#!/bin/bash
PORT="\${1:-${vizPort}}"
DIR="\$(cd "\$(dirname "\$0")/.build" && pwd)"
cd "\$DIR" || { echo "Build dir not found: \$DIR"; exit 1; }
pkill -f "vite preview.*--port \${PORT}" 2>/dev/null || true
sleep 1
echo "Starting viz server on port \${PORT}..."
npm run serve -- --host 0.0.0.0 --port "\${PORT}"
LAUNCHER
    chmod +x "\${VIZ_DIR}/start-server.sh"
    echo "Restart server later: \${VIZ_DIR}/start-server.sh [\${VIZ_PORT}]"
    """
}
