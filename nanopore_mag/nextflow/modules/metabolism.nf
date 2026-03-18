// Metabolic profiling: KofamScan (KEGG Orthology), eggNOG-mapper (COG/GO/EC/Pfam),
// dbCAN3 (CAZyme detection), annotation merging, bin mapping, and KEGG module scoring.
//
// Architecture: "annotate once, map to bins" — all tools run on the full .faa,
// then annotations are partitioned by contig2bin assignments (DAS_Tool + all binners).

process KOFAMSCAN {
    tag "kofamscan"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-classify"
    publishDir "${params.outdir}/metabolism/kofamscan", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/metabolism/kofamscan" : null

    input:
    path(proteins)
    path(kofam_profiles)
    path(ko_list)

    output:
    path("kofamscan_results.tsv"), emit: ko_assignments

    script:
    """
    # KofamScan: assign KEGG Orthology (KO) numbers using adaptive score thresholds
    # --format detail-tsv: output with significance mark, gene, KO, threshold, score, E-value
    # Profiles are HMMs from KEGG, ko_list provides per-KO adaptive thresholds

    if [ ! -s "${proteins}" ]; then
        echo "[WARNING] No protein sequences — skipping KofamScan" >&2
        printf '# gene name\\tKO\\tthrshld\\tscore\\tE-value\\tKO definition\\n' > kofamscan_results.tsv
        exit 0
    fi

    set +e
    exec_annotation \\
        --profile "${kofam_profiles}" \\
        --ko-list "${ko_list}" \\
        --cpu ${task.cpus} \\
        --format detail-tsv \\
        -o kofamscan_results.tsv \\
        "${proteins}"
    kofam_exit=\$?
    set -e

    if [ \$kofam_exit -ne 0 ] || [ ! -f kofamscan_results.tsv ]; then
        echo "[WARNING] KofamScan exited with code \$kofam_exit" >&2
        printf '# gene name\\tKO\\tthrshld\\tscore\\tE-value\\tKO definition\\n' > kofamscan_results.tsv
        exit 0
    fi
    """
}

process EMAPPER {
    tag "emapper"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-annotate"
    publishDir "${params.outdir}/metabolism/emapper", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/metabolism/emapper" : null

    input:
    path(proteins)
    path(eggnog_db)

    output:
    path("emapper_results.emapper.annotations"), emit: annotations

    script:
    def cpus = task.cpus
    def batch_size = params.emapper_batch_size
    """
    # eggNOG-mapper: functional annotation via DIAMOND search against eggNOG database
    # Assigns COG categories, GO terms, EC numbers, KEGG KOs, Pfam domains, descriptions
    #
    # Two-stage batched mode to prevent /tmp overflow from DIAMOND temp files:
    #   Stage 1: DIAMOND search per chunk (--no_annot), no --dbmem needed
    #   Stage 2: Annotate all hits at once (-m no_search --dbmem)
    #
    # DIAMOND writes unlinked temp files (~5 MB/s). Batching bounds temp size per chunk.
    # Chunk size controlled by params.emapper_batch_size (default 50000 proteins).

    HEADER='## eggNOG-mapper\\n#query\\tseed_ortholog\\tevalue\\tscore\\teggNOG_OGs\\tmax_annot_lvl\\tCOG_category\\tDescription\\tPreferred_name\\tGOs\\tEC\\tKEGG_ko\\tKEGG_Pathway\\tKEGG_Module\\tKEGG_Reaction\\tKEGG_rclass\\tBRITE\\tKEGG_TC\\tCAZy\\tBiGG_Reaction\\tPFAMs'

    if [ ! -s "${proteins}" ]; then
        echo "[WARNING] No protein sequences — skipping eggNOG-mapper" >&2
        printf "\$HEADER\\n" > emapper_results.emapper.annotations
        exit 0
    fi

    # Count proteins and split into chunks (same pattern as SENDSKETCH_CLASSIFY)
    N_SEQS=\$(grep -c '^>' "${proteins}")
    BATCH_SIZE=${batch_size}
    echo "[INFO] EMAPPER: \$N_SEQS proteins, batch size \$BATCH_SIZE" >&2

    mkdir -p chunks
    if [ "\$N_SEQS" -le "\$BATCH_SIZE" ]; then
        ln -s "\$(readlink -f ${proteins})" chunks/chunk_000.faa
    else
        awk -v bs="\$BATCH_SIZE" -v dir="chunks" '
            BEGIN { n=0; chunk=0; fn=sprintf("%s/chunk_%03d.faa", dir, chunk) }
            /^>/ { if (n >= bs) { close(fn); chunk++; fn=sprintf("%s/chunk_%03d.faa", dir, chunk); n=0 } n++ }
            { print > fn }
        ' "${proteins}"
    fi

    N_CHUNKS=\$(ls chunks/chunk_*.faa | wc -l)
    echo "[INFO] Split into \$N_CHUNKS chunk(s)" >&2

    # Stage 1: DIAMOND search per chunk (no --dbmem — search needs ~16 GB, not 48 GB)
    # Each chunk's temp files are cleaned before the next starts.
    > merged.seed_orthologs
    search_ok=true
    chunk_i=0
    for chunk in chunks/chunk_*.faa; do
        chunk_i=\$((chunk_i + 1))
        name=\$(basename "\$chunk" .faa)
        echo "[INFO] Stage 1: DIAMOND search chunk \$chunk_i/\$N_CHUNKS (\$name)" >&2

        set +e
        emapper.py \\
            -i "\$chunk" \\
            --data_dir "${eggnog_db}" \\
            -m diamond \\
            --no_annot \\
            --cpu ${cpus} \\
            --output "\$name" \\
            --override \\
            --dmnd_iterate no \\
            --temp_dir .
        chunk_exit=\$?
        set -e

        if [ \$chunk_exit -ne 0 ]; then
            echo "[WARNING] DIAMOND search failed for \$name (exit \$chunk_exit)" >&2
            search_ok=false
            break
        fi

        if [ -f "\${name}.emapper.seed_orthologs" ]; then
            cat "\${name}.emapper.seed_orthologs" >> merged.seed_orthologs
        fi

        # Clean up chunk temp files to free disk space
        rm -f "\${name}".emapper.* 2>/dev/null || true
    done

    if [ "\$search_ok" != "true" ] || [ ! -s merged.seed_orthologs ]; then
        echo "[WARNING] eggNOG-mapper search stage failed or produced no hits" >&2
        printf "\$HEADER\\n" > emapper_results.emapper.annotations
        exit 0
    fi

    echo "[INFO] Stage 1 complete: \$(wc -l < merged.seed_orthologs) seed orthologs" >&2

    # Stage 2: Annotate all hits at once (--dbmem loads 44 GB sqlite into RAM)
    echo "[INFO] Stage 2: Annotating merged seed orthologs" >&2
    set +e
    emapper.py \\
        -m no_search \\
        --annotate_hits_table merged.seed_orthologs \\
        --data_dir "${eggnog_db}" \\
        --dbmem \\
        --cpu ${cpus} \\
        --output emapper_results \\
        --override
    annot_exit=\$?
    set -e

    if [ \$annot_exit -ne 0 ] || [ ! -f emapper_results.emapper.annotations ]; then
        echo "[WARNING] eggNOG-mapper annotation stage exited with code \$annot_exit" >&2
        printf "\$HEADER\\n" > emapper_results.emapper.annotations
        exit 0
    fi
    """
}

process DBCAN {
    tag "dbcan"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-genomic"
    publishDir "${params.outdir}/metabolism/dbcan", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/metabolism/dbcan" : null

    input:
    path(proteins)
    path(dbcan_db)

    output:
    path("overview.tsv"),    emit: overview
    path("hmmer.out"),       emit: hmmer, optional: true
    path("diamond.out"),     emit: diamond, optional: true
    path("substrate.out"),   emit: substrate, optional: true

    script:
    """
    # dbCAN3 (run_dbcan): CAZyme annotation using three methods:
    #   1. HMMER against dbCAN HMM database
    #   2. DIAMOND against CAZy sequence database
    #   3. dbCAN-sub (HMM of CAZyme subfamilies)
    # Consensus: CAZyme families agreed upon by >= 2/3 methods

    if [ ! -s "${proteins}" ]; then
        echo "[WARNING] No protein sequences — skipping dbCAN" >&2
        printf 'Gene ID\\tEC#\\tdbCAN_hmm\\tdbCAN_sub\\tDIAMOND\\t#ofTools\\tRecommend Results\\tSubstrate\\n' > overview.tsv
        exit 0
    fi

    set +e
    run_dbcan CAZyme_annotation \\
        --input_raw_data "${proteins}" \\
        --mode protein \\
        --db_dir "${dbcan_db}" \\
        --output_dir dbcan_out \\
        --threads ${task.cpus}
    dbcan_exit=\$?
    set -e

    if [ \$dbcan_exit -ne 0 ] || [ ! -d dbcan_out ]; then
        echo "[WARNING] dbCAN exited with code \$dbcan_exit" >&2
        printf 'Gene ID\\tEC#\\tdbCAN_hmm\\tdbCAN_sub\\tDIAMOND\\t#ofTools\\tRecommend Results\\tSubstrate\\n' > overview.tsv
        exit 0
    fi

    # Copy overview (run_dbcan v4+ outputs .tsv, older versions output .txt)
    if [ -f dbcan_out/overview.tsv ]; then
        cp dbcan_out/overview.tsv .
    elif [ -f dbcan_out/overview.txt ]; then
        cp dbcan_out/overview.txt overview.tsv
    else
        printf 'Gene ID\\tEC#\\tdbCAN_hmm\\tdbCAN_sub\\tDIAMOND\\t#ofTools\\tRecommend Results\\tSubstrate\\n' > overview.tsv
    fi

    # Optional detailed outputs
    [ -f dbcan_out/hmmer.out ] && cp dbcan_out/hmmer.out . || true
    [ -f dbcan_out/diamond.out ] && cp dbcan_out/diamond.out . || true
    [ -f dbcan_out/substrate.out ] && cp dbcan_out/substrate.out . || true
    """
}

process MERGE_ANNOTATIONS {
    tag "merge"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-classify"
    publishDir "${params.outdir}/metabolism/merged", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/metabolism/merged" : null

    input:
    val(labels)
    path(annotation_files)

    output:
    path("merged_annotations.tsv"), emit: merged

    script:
    // Build flags dynamically based on which annotation tools ran
    def kofam_file   = labels.contains('kofamscan') ? annotation_files.find { it.name.contains('kofamscan') } : null
    def emapper_file = labels.contains('emapper')   ? annotation_files.find { it.name.contains('emapper') }   : null
    def dbcan_file   = labels.contains('dbcan')      ? annotation_files.find { it.name.contains('overview') }  : null
    def flags = []
    if (kofam_file)   flags << "--kofamscan ${kofam_file}"
    if (emapper_file) flags << "--emapper ${emapper_file}"
    if (dbcan_file)   flags << "--dbcan ${dbcan_file}"
    """
    # Merge annotation sources into unified per-protein TSV
    # Tools that ran: ${labels.join(', ')}

    merge_annotations.py ${flags.join(' ')} -o merged_annotations.tsv
    """
}

process MAP_TO_BINS {
    tag "map_to_bins"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-classify"
    publishDir "${params.outdir}/metabolism", mode: 'copy', enabled: !params.store_dir, pattern: 'per_mag'
    publishDir "${params.outdir}/metabolism/community", mode: 'copy', enabled: !params.store_dir, pattern: 'community_annotations.tsv'
    storeDir params.store_dir ? "${params.store_dir}/metabolism/per_mag" : null

    input:
    path(merged)
    path(contig2bin)
    path(binner_bins)
    path(gff)

    output:
    path("per_mag/"),                     emit: per_mag
    path("community_annotations.tsv"),    emit: community

    script:
    """
    # Map merged annotations to bins via combined contig2bin (DAS_Tool + all binners)
    # Uses GFF for accurate gene->contig mapping, then contig->bin assignments
    # Outputs: per-bin TSVs + community-wide table with bin_id column

    cat ${contig2bin} ${binner_bins} > combined_contig2bin.tsv

    map_annotations_to_bins.py \\
        --annotations "${merged}" \\
        --contig2bin combined_contig2bin.tsv \\
        --gff "${gff}" \\
        --outdir per_mag \\
        --community community_annotations.tsv
    """
}

process KEGG_MODULES {
    tag "kegg_modules"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-classify"
    publishDir "${params.outdir}/metabolism/modules", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/metabolism/modules" : null

    input:
    path(per_mag_dir)

    output:
    path("module_completeness.tsv"), emit: modules
    path("module_heatmap.svg"),      emit: heatmap, optional: true

    script:
    """
    # Evaluate ~150 KEGG module definitions against KO presence per MAG
    # Outputs completeness matrix (MAG x module) + clustered SVG heatmap
    # Covers: carbon fixation, nitrogen, sulfur, methane, ETC, vitamins, central carbon

    kegg_module_completeness.py \\
        --input "${per_mag_dir}" \\
        --output module_completeness.tsv \\
        --heatmap module_heatmap.svg
    """
}

process MINPATH {
    tag "minpath"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-pathviz"
    publishDir "${params.outdir}/metabolism/minpath", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/metabolism/minpath" : null

    input:
    path(per_mag_dir)

    output:
    path("minpath_pathways.tsv"),  emit: pathways
    path("details/"),              emit: details, optional: true

    script:
    def minpath_dir = "${projectDir}/conda-envs/dana-mag-pathviz/share/minpath"
    """
    # MinPath (Ye & Doak, 2009): parsimony pathway reconstruction
    # Finds the minimum set of KEGG pathways consistent with observed KOs,
    # preventing pathway inflation in draft MAGs.
    # Uses GLPK integer programming solver.

    run_minpath_per_mag.py \\
        --input "${per_mag_dir}" \\
        --minpath_dir "${minpath_dir}" \\
        --output minpath_pathways.tsv \\
        --details details
    """
}

process KEGG_DECODER {
    tag "kegg_decoder"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-pathviz"
    publishDir "${params.outdir}/metabolism/kegg_decoder", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/metabolism/kegg_decoder" : null

    input:
    path(per_mag_dir)

    output:
    path("kegg_decoder_output.tsv"),  emit: output
    path("function_heatmap.svg"),     emit: heatmap, optional: true

    script:
    """
    # KEGG-Decoder (Graham et al. 2018): scores ~80 environmentally-curated
    # biogeochemical modules and generates publication-quality heatmaps.
    # Covers: carbon, nitrogen, sulfur, photosynthesis, transporters, etc.

    # Step 1: Convert per-MAG TSVs to KEGG-Decoder input format
    prepare_keggdecoder_input.py \\
        --input "${per_mag_dir}" \\
        --output keggdecoder_input.tsv

    # Step 2: Run KEGG-Decoder (skip if no KO entries; requires >= 3 genomes)
    if [ -s keggdecoder_input.tsv ]; then
        # Count distinct genomes (prefix before first '_')
        n_genomes=\$(cut -f1 keggdecoder_input.tsv | cut -d'_' -f1 | sort -u | wc -l)
        if [ "\$n_genomes" -lt 3 ]; then
            echo "[WARNING] Only \$n_genomes genomes — KEGG-Decoder requires >= 3, skipping" >&2
            printf 'Function\\n' > kegg_decoder_output.tsv
        else
            set +e
            KEGG-decoder -i keggdecoder_input.tsv -o kegg_decoder_output.tsv -v static
            decoder_exit=\$?
            set -e

            if [ \$decoder_exit -ne 0 ] || [ ! -f kegg_decoder_output.tsv ]; then
                echo "[WARNING] KEGG-Decoder exited with code \$decoder_exit" >&2
                printf 'Function\\n' > kegg_decoder_output.tsv
            fi

            # KEGG-Decoder names the SVG after the -o file (e.g. kegg_decoder_output.svg)
            if [ ! -f "function_heatmap.svg" ] && ls *.svg 1>/dev/null 2>&1; then
                mv *.svg function_heatmap.svg 2>/dev/null || true
            fi
        fi
    else
        echo "[WARNING] No KO entries found — skipping KEGG-Decoder" >&2
        printf 'Function\\n' > kegg_decoder_output.tsv
    fi
    """
}

process ANTISMASH {
    tag "antismash"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-antismash"
    storeDir params.store_dir ? "${params.store_dir}/metabolism/antismash" : null

    input:
    path(assembly)
    path(gff)

    output:
    path("antismash_summary.tsv"),   emit: summary
    path("antismash_geneclusters/"), emit: geneclusters
    path("antismash_json/"),         emit: json, optional: true

    script:
    def db_arg = params.antismash_db ? "--databases ${params.antismash_db}" : ""
    def genefinding = gff.name != 'NO_GFF' ? "--genefinding-gff3 ${gff}" : "--genefinding-tool prodigal-m"
    """
    set +e
    antismash \\
        --cpus ${task.cpus} \\
        --output-dir antismash_out \\
        --output-basename antismash \\
        ${db_arg} \\
        ${genefinding} \\
        --cb-general \\
        --cb-knownclusters \\
        --cc-mibig \\
        --minlength 1000 \\
        --allow-long-headers \\
        --taxon bacteria \\
        "${assembly}"
    antismash_exit=\$?
    set -e

    mkdir -p antismash_geneclusters antismash_json

    if [ \$antismash_exit -ne 0 ]; then
        echo "[WARNING] antiSMASH exited with code \$antismash_exit" >&2
        printf 'region\\tcontig\\tstart\\tend\\ttype\\tmost_similar_known_cluster\\tsimilarity\\n' > antismash_summary.tsv
        exit 0
    fi

    # Collect region GenBank files
    find antismash_out -name "*.region*.gbk" -exec cp {} antismash_geneclusters/ \\; 2>/dev/null || true

    # Collect JSON output
    [ -f antismash_out/antismash.json ] && cp antismash_out/antismash.json antismash_json/ || true

    # Parse JSON to summary TSV
    python3 -c "
import json, os, sys

header = 'region\\tcontig\\tstart\\tend\\ttype\\tmost_similar_known_cluster\\tsimilarity'
rows = []
jf = 'antismash_out/antismash.json'
if os.path.isfile(jf):
    with open(jf) as f:
        data = json.load(f)
    for rec in data.get('records', []):
        contig = rec.get('id', 'unknown')
        for feat in rec.get('areas', rec.get('features', [])):
            if feat.get('type') == 'region':
                loc = feat.get('location', '')
                start = feat.get('start', '')
                end = feat.get('end', '')
                products = ','.join(feat.get('products', feat.get('product', ['unknown'])))
                knowncluster = feat.get('knownclusterblast', [{}])
                if knowncluster:
                    kc = knowncluster[0] if isinstance(knowncluster, list) else knowncluster
                    similar = kc.get('description', 'NA')
                    sim_pct = kc.get('similarity', 'NA')
                else:
                    similar = 'NA'
                    sim_pct = 'NA'
                rows.append(f'{contig}_region\\t{contig}\\t{start}\\t{end}\\t{products}\\t{similar}\\t{sim_pct}')

# Fallback: parse region GBKs if JSON parsing yielded nothing
if not rows:
    for gbk in sorted(os.listdir('antismash_geneclusters')):
        if gbk.endswith('.gbk'):
            name = gbk.replace('.gbk','')
            rows.append(f'{name}\\tNA\\tNA\\tNA\\tNA\\tNA\\tNA')

print(header)
for r in rows:
    print(r)
" > antismash_summary.tsv
    """
}

// ============================================================================
// ECOSSDB — Ecosystem Services Profiling
// ============================================================================

def ecossdb_bin = "${projectDir}/ecossdb/bin"

process ECOSSDB_MAP {
    tag "ecossdb_map"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-classify"
    beforeScript "export PATH=${ecossdb_bin}:\$PATH"
    publishDir "${params.outdir}/metabolism/ecossdb", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/metabolism/ecossdb" : null

    input:
    path community_annotations
    path contig2bin
    path es_mapping
    path es_ontology

    output:
    path 'es_gene_catalog.tsv',    emit: catalog
    path 'es_per_contig.tsv',      emit: contigs
    path 'es_per_mag.tsv',         emit: mag_profiles
    path 'es_mapping_stats.tsv',   emit: stats

    script:
    """
    set -euo pipefail

    # 1. Normalize annotations (strip ko: prefix, unify columns)
    normalize_annotations.py \
        --input ${community_annotations} \
        --format danaseq \
        --output normalized.tsv

    # 2. Map annotations to ecosystem services
    map_to_es.py \
        --annotations normalized.tsv \
        --mapping ${es_mapping} \
        --ontology ${es_ontology} \
        --output es_gene_hits.tsv \
        --stats es_mapping_stats.tsv

    # 3. Aggregate to contigs (merge + dedup)
    aggregate_contigs.py \
        --gene-hits es_gene_hits.tsv \
        --output-catalog es_gene_catalog.tsv \
        --output-contigs es_per_contig.tsv

    # 4. Aggregate to MAGs
    aggregate_bins.py \
        --catalog es_gene_catalog.tsv \
        --contig2bin ${contig2bin} \
        --output es_per_mag.tsv
    """
}

process ECOSSDB_SCORE {
    tag "ecossdb_score"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-classify"
    beforeScript "export PATH=${ecossdb_bin}:\$PATH"
    publishDir "${params.outdir}/metabolism/ecossdb", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/metabolism/ecossdb_score" : null

    input:
    path catalog
    path es_mapping
    val  role_weights

    output:
    path 'es_scores.tsv',     emit: scores
    path 'es_confidence.tsv', emit: confidence

    script:
    """
    score_es.py \
        --catalog ${catalog} \
        --mapping ${es_mapping} \
        --role-weights '${role_weights}' \
        --output es_scores.tsv \
        --confidence-out es_confidence.tsv
    """
}

process ECOSSDB_SDG {
    tag "ecossdb_sdg"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-classify"
    beforeScript "export PATH=${ecossdb_bin}:\$PATH"
    publishDir "${params.outdir}/metabolism/ecossdb/sdg", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/metabolism/ecossdb_sdg" : null

    input:
    path scores
    path crosswalk
    path sdg_targets

    output:
    path 'es_sdg_targets.tsv', emit: targets
    path 'es_sdg_goals.tsv',   emit: goals
    path 'es_sdg.json',        emit: json

    script:
    """
    map_es_to_sdg.py \
        --scores ${scores} \
        --crosswalk ${crosswalk} \
        --sdg-targets ${sdg_targets} \
        --output-targets es_sdg_targets.tsv \
        --output-goals es_sdg_goals.tsv \
        --output-json es_sdg.json
    """
}

process ECOSSDB_VIZ {
    tag "ecossdb_viz"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-classify"
    beforeScript "export PATH=${ecossdb_bin}:\$PATH"
    publishDir "${params.outdir}/metabolism/ecossdb", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/metabolism/ecossdb_viz" : null

    input:
    path scores
    path catalog
    path mag_profiles
    path hierarchy_json

    output:
    path 'ecosystem_services.json',    emit: json
    path 'ecosystem_services.json.gz', emit: json_gz

    script:
    """
    es_to_json.py \
        --scores ${scores} \
        --catalog ${catalog} \
        --mag-profiles ${mag_profiles} \
        --hierarchy ${hierarchy_json} \
        --output ecosystem_services.json
    """
}
