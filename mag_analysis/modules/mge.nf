// Mobile genetic element detection.
//
// Processes (all optional, gated by params.run_* + database paths):
//   GENOMAD_CLASSIFY   — geNomad end-to-end: marker genes + neural network → viruses,
//                         plasmids, proviruses. Runs directly on the assembly FASTA.
//   CHECKV_QUALITY     — CheckV viral QA: AAI + HMM completeness, host trimming.
//                         Requires geNomad virus FASTA as input.
//   INTEGRONFINDER     — Integron detection: integrase + attC/attI + gene cassettes.
//                         Runs on assembly; no annotation dependency. A subworkflow:
//                         the assembly is split into chunks run as separate tasks.
//   ISLANDPATH_DIMOB   — Genomic island detection via dinucleotide bias + mobility
//                         gene HMMs. Requires annotation GFF + FAA.
//   MACSYFINDER        — Secretion systems (TXSScan) + conjugation (CONJScan).
//                         Requires protein FAA.
//   DEFENSEFINDER      — Anti-phage defense systems (CRISPR, R-M, BREX, Abi, etc.).
//                         Requires protein FAA + GFF. Parallelized per-contig.

process GENOMAD_CLASSIFY {
    tag "genomad"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-quality"
    publishDir "${params.outdir}/mge/genomad", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/mge/genomad" : null

    input:
    path(assembly)

    output:
    path("virus_summary.tsv"),     emit: virus_summary
    path("plasmid_summary.tsv"),   emit: plasmid_summary
    path("virus.fna"),             emit: virus_fasta
    path("plasmid.fna"),           emit: plasmid_fasta
    path("virus_proteins.faa"),    emit: virus_proteins
    path("plasmid_proteins.faa"),  emit: plasmid_proteins
    path("virus_genes.tsv"),       emit: virus_genes
    path("plasmid_genes.tsv"),     emit: plasmid_genes
    path("provirus.tsv"),          emit: provirus_coords
    path("provirus.fna"),          emit: provirus_fasta
    path("taxonomy.tsv"),          emit: taxonomy
    path("genomad_summary.tsv"),   emit: summary

    script:
    def db_path = params.genomad_db
    """
    # geNomad end-to-end: marker gene annotation → neural network classification
    # Detects viruses, plasmids, and proviruses in a single pass
    # Note: no --cleanup so intermediate files (annotate, find_proviruses) are preserved
    set +e
    genomad end-to-end \\
        --splits ${task.cpus} \\
        "${assembly}" \\
        genomad_out \\
        "${db_path}"
    genomad_exit=\$?
    set -e

    # geNomad names output files based on the input filename
    input_base=\$(basename "${assembly}" | sed 's/\\.[^.]*\$//')

    if [ \$genomad_exit -ne 0 ]; then
        echo "[WARNING] geNomad exited with code \$genomad_exit" >&2
        touch virus_summary.tsv plasmid_summary.tsv virus.fna plasmid.fna \\
              virus_proteins.faa plasmid_proteins.faa virus_genes.tsv plasmid_genes.tsv \\
              provirus.tsv provirus.fna taxonomy.tsv genomad_summary.tsv
        exit 0
    fi

    # Helper: copy file if exists, else touch empty
    copy_or_touch() {
        if [ -f "\$1" ]; then cp "\$1" "\$2"; else touch "\$2"; fi
    }

    # Summary outputs (virus/plasmid summaries, sequences, proteins, gene annotations)
    copy_or_touch "genomad_out/\${input_base}_summary/\${input_base}_virus_summary.tsv"     virus_summary.tsv
    copy_or_touch "genomad_out/\${input_base}_summary/\${input_base}_plasmid_summary.tsv"   plasmid_summary.tsv
    copy_or_touch "genomad_out/\${input_base}_summary/\${input_base}_virus.fna"             virus.fna
    copy_or_touch "genomad_out/\${input_base}_summary/\${input_base}_plasmid.fna"           plasmid.fna
    copy_or_touch "genomad_out/\${input_base}_summary/\${input_base}_virus_proteins.faa"    virus_proteins.faa
    copy_or_touch "genomad_out/\${input_base}_summary/\${input_base}_plasmid_proteins.faa"  plasmid_proteins.faa
    copy_or_touch "genomad_out/\${input_base}_summary/\${input_base}_virus_genes.tsv"       virus_genes.tsv
    copy_or_touch "genomad_out/\${input_base}_summary/\${input_base}_plasmid_genes.tsv"     plasmid_genes.tsv

    # Provirus detection results
    copy_or_touch "genomad_out/\${input_base}_find_proviruses/\${input_base}_provirus.tsv"  provirus.tsv
    copy_or_touch "genomad_out/\${input_base}_find_proviruses/\${input_base}_provirus.fna"  provirus.fna

    # Per-contig taxonomy from annotation step
    copy_or_touch "genomad_out/\${input_base}_annotate/\${input_base}_taxonomy.tsv"         taxonomy.tsv

    # Aggregated classification scores (all contigs)
    copy_or_touch "genomad_out/\${input_base}_aggregated_classification/\${input_base}_aggregated_classification.tsv" genomad_summary.tsv
    """
}

process CHECKV_QUALITY {
    tag "checkv"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-quality"
    publishDir "${params.outdir}/mge/checkv", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/mge/checkv" : null

    input:
    path(virus_fasta)

    output:
    path("quality_summary.tsv"), emit: quality
    path("viruses.fna"),         emit: viruses
    path("proviruses.fna"),      emit: proviruses

    script:
    def db_path = params.checkv_db
    """
    # CheckV: assess viral genome completeness and contamination
    # Trims host contamination from proviruses, estimates completeness via
    # AAI comparison to reference genomes or HMM-based gene density models

    if [ ! -s "${virus_fasta}" ]; then
        echo "[WARNING] No viral contigs to assess — skipping CheckV" >&2
        printf 'contig_id\\tcontig_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tcompleteness\\tcontamination\\n' > quality_summary.tsv
        touch viruses.fna proviruses.fna
        exit 0
    fi

    set +e
    checkv end_to_end \\
        "${virus_fasta}" \\
        checkv_out \\
        -d "${db_path}" \\
        -t ${task.cpus}
    checkv_exit=\$?
    set -e

    if [ \$checkv_exit -ne 0 ]; then
        echo "[WARNING] CheckV exited with code \$checkv_exit" >&2
        printf 'contig_id\\tcontig_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tcompleteness\\tcontamination\\n' > quality_summary.tsv
        touch viruses.fna proviruses.fna
        exit 0
    fi

    # Copy outputs
    if [ -f checkv_out/quality_summary.tsv ]; then
        cp checkv_out/quality_summary.tsv .
    else
        printf 'contig_id\\tcontig_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tcompleteness\\tcontamination\\n' > quality_summary.tsv
    fi

    if [ -f checkv_out/viruses.fna ]; then
        cp checkv_out/viruses.fna .
    else
        touch viruses.fna
    fi

    if [ -f checkv_out/proviruses.fna ]; then
        cp checkv_out/proviruses.fna .
    else
        touch proviruses.fna
    fi
    """
}

// IntegronFinder writes a .integrons and a .summary for every replicon into one
// directory and merges them only when the whole input is done, so a single run
// peaks at two files per contig: ~3.1M on a 1.55M-contig co-assembly (#69).
// The assembly is therefore split into chunks of params.integron_chunk_size
// contigs, each chunk is its own task that deletes its per-contig files before
// it ends, and at most params.integron_max_forks chunks run at once. The peak
// is about 2 x chunk size x max forks files, whatever the assembly's size.

process INTEGRONFINDER_SPLIT {
    tag "integronfinder"
    label 'process_low'

    input:
    path(assembly)

    output:
    path("chunks/chunk_*.fa"), emit: chunks

    script:
    """
    mkdir -p chunks
    awk -v n=${params.integron_chunk_size} '
        /^>/ { if (c % n == 0) { if (out) close(out); out = sprintf("chunks/chunk_%05d.fa", ++k) } c++ }
        out  { print > out }
    ' "${assembly}"
    n_contigs=\$(grep -c '^>' "${assembly}" || true)
    if [ "\${n_contigs:-0}" -lt 1 ]; then
        echo "[ERROR] ${assembly} contains no contigs" >&2
        exit 1
    fi
    echo "[INFO] \${n_contigs} contigs in \$(ls chunks | wc -l) chunks of up to ${params.integron_chunk_size}"
    """
}

process INTEGRONFINDER_CHUNK {
    tag "${chunk.baseName}"
    label 'process_low'
    maxForks params.integron_max_forks
    conda "${projectDir}/conda-envs/dana-mag-genomic"

    input:
    path(chunk)

    output:
    path("${chunk.baseName}.integrons.tsv"), emit: integrons
    path("${chunk.baseName}.summary.tsv"),   emit: summary

    script:
    def base = chunk.baseName
    """
    # IntegronFinder: detect integrons (integrase + attC/attI sites + gene cassettes)
    # --local-max:     thorough local detection of attC sites (more sensitive)
    # --func-annot:    annotate gene cassettes with Resfams HMM profiles (AMR)
    # --promoter-attI: also search for Pc promoter and attI recombination sites
    # --linear:        contigs from Flye assembly are linear, not circular replicons
    # --cpu:           threading for INFERNAL (cmsearch) and HMMER (hmmsearch)
    integron_finder \\
        --local-max \\
        --func-annot \\
        --promoter-attI \\
        --linear \\
        --cpu ${task.cpus} \\
        --outdir integron_out \\
        "${chunk}"

    results_dir="integron_out/Results_Integron_Finder_${base}"
    cp "\${results_dir}/${base}.integrons" ${base}.integrons.tsv
    cp "\${results_dir}/${base}.summary"   ${base}.summary.tsv

    # IntegronFinder writes one summary row per replicon it searched. It skips,
    # with no row, a sequence of 50 bp or less ("is too short") and a replicon
    # with no predicted proteins ("Skip replicon"); any other shortfall means it
    # did not get through the chunk.
    n_contigs=\$(grep -c '^>' "${chunk}")
    n_rows=\$(grep -v -e '^#' -e '^ID_replicon' ${base}.summary.tsv | grep -c . || true)
    n_short=\$(grep -c 'is too short' "\${results_dir}/integron_finder.out" || true)
    n_noprot=\$(grep -c 'Skip replicon' "\${results_dir}/integron_finder.out" || true)
    if [ "\$((n_rows + n_short + n_noprot))" -ne "\${n_contigs}" ]; then
        echo "[ERROR] ${base}: \${n_rows} summary rows + \${n_short} too short + \${n_noprot} without proteins, for \${n_contigs} contigs" >&2
        exit 1
    fi
    [ "\$((n_short + n_noprot))" -eq 0 ] || echo "[INFO] ${base}: skipped \${n_short} contigs of 50 bp or less and \${n_noprot} without predicted proteins"

    # The per-contig files this bounds (#69).
    rm -rf integron_out
    """
}

process INTEGRONFINDER_MERGE {
    tag "integronfinder"
    label 'process_low'
    publishDir "${params.outdir}/mge/integrons", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/mge/integrons" : null

    input:
    path(integrons, stageAs: 'chunks/*')
    path(summaries, stageAs: 'chunks/*')
    val(n_chunks)

    output:
    path("integrons.tsv"), emit: integrons
    path("summary.tsv"),   emit: summary

    script:
    """
    # A chunk that failed twice is ignored by the global errorStrategy, so a
    # missing one must fail here rather than drop its contigs from the result.
    n_int=\$(ls chunks/*.integrons.tsv | wc -l)
    n_sum=\$(ls chunks/*.summary.tsv | wc -l)
    if [ "\${n_int}" -ne ${n_chunks} ] || [ "\${n_sum}" -ne ${n_chunks} ]; then
        echo "[ERROR] expected ${n_chunks} chunks, got \${n_int} integron and \${n_sum} summary tables" >&2
        exit 1
    fi

    # One header per table, comment lines dropped ("# No Integron found" is the
    # whole file for a chunk without integrons). Chunk order is contig order.
    merge() {  # \$1=first header field  \$2=header to use if no chunk has one
        awk -v key="\$1" -v fallback="\$2" '
            /^#/ { next }
            \$1 == key { if (!hdr) { hdr = \$0; print } next }
            NF { if (!hdr) { hdr = fallback; print hdr } print }
            END { if (!hdr) print fallback }
        ' FS='\\t' "\${@:3}"
    }
    merge ID_integron \\
        "\$(printf 'ID_integron\\tID_replicon\\telement\\tpos_beg\\tpos_end\\tstrand\\tevalue\\ttype_elt\\tannotation\\tmodel\\ttype\\tdefault\\tdistance_2attC\\tconsidered_topology')" \\
        \$(ls chunks/*.integrons.tsv | sort) > integrons.tsv
    merge ID_replicon \\
        "\$(printf 'ID_replicon\\tCALIN\\tcomplete\\tIn0\\ttopology\\tsize')" \\
        \$(ls chunks/*.summary.tsv | sort) > summary.tsv

    echo "[INFO] \$(grep -vc '^ID_replicon' summary.tsv) contigs, \$(awk -F'\\t' 'NR > 1 { print \$2 }' integrons.tsv | sort -u | grep -c . || true) with integron elements"
    """
}

workflow INTEGRONFINDER {
    take:
    assembly

    main:
    def stored = params.store_dir ? file("${params.store_dir}/mge/integrons") : null
    if (stored && stored.resolve('integrons.tsv').exists() && stored.resolve('summary.tsv').exists()) {
        // Store mode skips a finished process by its outputs, but a subworkflow
        // has none of its own: without this, a resubmission would split and
        // rerun every chunk before the merge found its outputs stored.
        integrons = Channel.value(stored.resolve('integrons.tsv'))
        summary   = Channel.value(stored.resolve('summary.tsv'))
    } else {
        INTEGRONFINDER_SPLIT(assembly)
        chunks   = INTEGRONFINDER_SPLIT.out.chunks.flatten()
        n_chunks = INTEGRONFINDER_SPLIT.out.chunks.map { it instanceof List ? it.size() : 1 }
        INTEGRONFINDER_CHUNK(chunks)
        INTEGRONFINDER_MERGE(
            INTEGRONFINDER_CHUNK.out.integrons.collect(),
            INTEGRONFINDER_CHUNK.out.summary.collect(),
            n_chunks
        )
        integrons = INTEGRONFINDER_MERGE.out.integrons
        summary   = INTEGRONFINDER_MERGE.out.summary
    }

    emit:
    integrons
    summary
}

process ISLANDPATH_DIMOB {
    tag "islandpath"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-genomic"
    publishDir "${params.outdir}/mge/islandpath", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/mge/islandpath" : null

    input:
    path(assembly)
    path(gff)
    path(faa)

    output:
    path("genomic_islands.tsv"), emit: islands

    script:
    def hmm_db = params.islandpath_hmm_db ? "${params.islandpath_hmm_db}/Pfam-A_mobgenes_201512_prok" : "${projectDir}/conda-envs/dana-mag-genomic/opt/islandpath/hmmpfam/Pfam-A_mobgenes_201512_prok"
    """
    # IslandPath-DIMOB: detect genomic islands via dinucleotide bias + mobility genes
    # Python reimplementation — works directly with GFF + FASTA + FAA from Prokka
    # Reference-free method — HMM profiles for mobility genes (Pfam-A, 2015)
    # hmmscan benefits from multiple CPUs

    set +e
    islandpath_dimob.py \\
        --gff "${gff}" \\
        --fasta "${assembly}" \\
        --faa "${faa}" \\
        --hmm_db "${hmm_db}" \\
        --cpus ${task.cpus} \\
        -o genomic_islands.tsv
    dimob_exit=\$?
    set -e

    if [ \$dimob_exit -ne 0 ] || [ ! -f genomic_islands.tsv ]; then
        echo "[WARNING] IslandPath-DIMOB exited with code \$dimob_exit" >&2
        printf 'island_id\\tcontig\\tstart\\tend\\n' > genomic_islands.tsv
        exit 0
    fi
    """
}

process MACSYFINDER {
    tag "macsyfinder"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-genomic"
    publishDir "${params.outdir}/mge/macsyfinder", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/mge/macsyfinder" : null

    input:
    path(proteins)

    output:
    path("all_systems.tsv"),   emit: systems
    path("all_systems.txt"),   emit: systems_txt

    script:
    def models_dir = params.macsyfinder_models
    """
    # MacSyFinder v2: detect secretion systems + conjugation in metagenome proteins
    # --db-type unordered: no gene order considered (appropriate for fragmented contigs)
    # --replicon-topology linear: assembly contigs are linear fragments
    # --models TXSScan all: all 20 secretion/appendage systems (T1SS-T9SS, flagellum, pili)
    # --models CONJScan all: all 17 conjugation systems (8 conjugative + 8 decayed + MOB)
    # -w: parallel HMMER searches

    if [ ! -s "${proteins}" ]; then
        echo "[WARNING] No protein sequences — skipping MacSyFinder" >&2
        printf 'replicon\\thit_id\\tgene_name\\thit_pos\\tmodel_fqn\\tsys_id\\tsys_loci\\tlocus_num\\tsys_wholeness\\tsys_score\\tsys_occ\\thit_gene_ref\\thit_status\\thit_seq_len\\thit_i_eval\\thit_score\\thit_profile_cov\\thit_seq_cov\\thit_begin_match\\thit_end_match\\n' > all_systems.tsv
        echo "# No systems found (empty input)" > all_systems.txt
        exit 0
    fi

    set +e
    macsyfinder \\
        --db-type unordered \\
        --sequence-db "${proteins}" \\
        --replicon-topology linear \\
        --models-dir "${models_dir}" \\
        --models TXSScan all \\
        --models CONJScan all \\
        -w ${task.cpus} \\
        -o msf_out \\
        --mute
    msf_exit=\$?
    set -e

    if [ \$msf_exit -ne 0 ] || [ ! -d msf_out ]; then
        echo "[WARNING] MacSyFinder exited with code \$msf_exit" >&2
        printf 'replicon\\thit_id\\tgene_name\\thit_pos\\tmodel_fqn\\tsys_id\\tsys_loci\\tlocus_num\\tsys_wholeness\\tsys_score\\tsys_occ\\thit_gene_ref\\thit_status\\thit_seq_len\\thit_i_eval\\thit_score\\thit_profile_cov\\thit_seq_cov\\thit_begin_match\\thit_end_match\\n' > all_systems.tsv
        echo "# MacSyFinder failed" > all_systems.txt
        exit 0
    fi

    # Copy outputs (unordered mode produces all_systems.tsv and all_systems.txt)
    if [ -f msf_out/all_systems.tsv ]; then
        cp msf_out/all_systems.tsv .
    else
        printf 'replicon\\thit_id\\tgene_name\\thit_pos\\tmodel_fqn\\tsys_id\\tsys_loci\\tlocus_num\\tsys_wholeness\\tsys_score\\tsys_occ\\thit_gene_ref\\thit_status\\thit_seq_len\\thit_i_eval\\thit_score\\thit_profile_cov\\thit_seq_cov\\thit_begin_match\\thit_end_match\\n' > all_systems.tsv
    fi

    if [ -f msf_out/all_systems.txt ]; then
        cp msf_out/all_systems.txt .
    else
        echo "# No systems found" > all_systems.txt
    fi
    """
}

process DEFENSEFINDER {
    tag "defensefinder"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-genomic"
    publishDir "${params.outdir}/mge/defensefinder", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/mge/defensefinder" : null

    input:
    path(proteins)
    path(gff)

    output:
    path("systems.tsv"), emit: systems
    path("genes.tsv"),   emit: genes
    path("hmmer.tsv"),   emit: hmmer

    script:
    def models_opt = params.defensefinder_models ? "--models-dir ${params.defensefinder_models}" : ""
    """
    # DefenseFinder: detect anti-phage defense systems (CRISPR, R-M, BREX, Abi, etc.)
    # Uses HMMER searches across ~280 defense system HMM profiles
    # Parallelized by splitting proteins into per-contig chunks (gembase format)
    # so MacSyFinder's system detection runs on N smaller replicons concurrently

    if [ ! -s "${proteins}" ]; then
        echo "[WARNING] No protein sequences — skipping DefenseFinder" >&2
        printf 'sys_id\\ttype\\tsubtype\\tprotein_in_syst\\tgenes_count\\tspec\\n' > systems.tsv
        printf 'replicon\\thit_id\\tgene_name\\n' > genes.tsv
        printf 'hit_id\\treplicon\\tposition_hit\\thit_sequence_length\\n' > hmmer.tsv
        exit 0
    fi

    # If no pre-downloaded models, fetch them first
    if [ -z "${models_opt}" ]; then
        defense-finder update
    fi

    set +e
    parallel_defensefinder.py \\
        --proteins "${proteins}" \\
        --gff "${gff}" \\
        --output-dir df_out \\
        --workers ${task.cpus} \\
        ${models_opt}
    df_exit=\$?
    set -e

    if [ \$df_exit -ne 0 ] || [ ! -d df_out ]; then
        echo "[WARNING] parallel_defensefinder.py exited with code \$df_exit" >&2
        printf 'sys_id\\ttype\\tsubtype\\tprotein_in_syst\\tgenes_count\\tspec\\n' > systems.tsv
        printf 'replicon\\thit_id\\tgene_name\\n' > genes.tsv
        printf 'hit_id\\treplicon\\tposition_hit\\thit_sequence_length\\n' > hmmer.tsv
        exit 0
    fi

    # Copy merged outputs from wrapper
    cp df_out/systems.tsv .
    cp df_out/genes.tsv .
    cp df_out/hmmer.tsv .
    """
}
