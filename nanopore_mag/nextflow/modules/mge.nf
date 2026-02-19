// Mobile genetic element detection: geNomad (virus + plasmid), CheckV (viral QA),
// IntegronFinder (integron + gene cassette detection),
// IslandPath-DIMOB (genomic island detection via dinucleotide bias),
// MacSyFinder (secretion systems + conjugation detection),
// DefenseFinder (anti-phage defense system detection)

process GENOMAD_CLASSIFY {
    tag "genomad"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-genomad"
    publishDir "${params.outdir}/mge/genomad", mode: 'copy'
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
    conda "${projectDir}/conda-envs/dana-mag-checkv"
    publishDir "${params.outdir}/mge/checkv", mode: 'copy'
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

process INTEGRONFINDER {
    tag "integronfinder"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-integron"
    publishDir "${params.outdir}/mge/integrons", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/mge/integrons" : null

    input:
    path(assembly)

    output:
    path("integrons.tsv"),     emit: integrons
    path("summary.tsv"),       emit: summary

    script:
    """
    # IntegronFinder: detect integrons (integrase + attC/attI sites + gene cassettes)
    # --local-max:     thorough local detection of attC sites (more sensitive)
    # --func-annot:    annotate gene cassettes with Resfams HMM profiles (AMR)
    # --promoter-attI: also search for Pc promoter and attI recombination sites
    # --linear:        contigs from Flye assembly are linear, not circular replicons
    # --cpu:           threading for INFERNAL (cmsearch) and HMMER (hmmsearch)

    set +e
    integron_finder \\
        --local-max \\
        --func-annot \\
        --promoter-attI \\
        --linear \\
        --cpu ${task.cpus} \\
        --outdir integron_out \\
        "${assembly}"
    if_exit=\$?
    set -e

    # IntegronFinder output directory: integron_out/Results_Integron_Finder_<basename>/
    input_base=\$(basename "${assembly}" | sed 's/\\.[^.]*\$//')
    results_dir="integron_out/Results_Integron_Finder_\${input_base}"

    if [ \$if_exit -ne 0 ] || [ ! -d "\${results_dir}" ]; then
        echo "[WARNING] IntegronFinder exited with code \$if_exit" >&2
        printf 'ID_integron\\tID_replicon\\telement\\tpos_beg\\tpos_end\\tstrand\\tevalue\\ttype_elt\\tmodel\\ttype\\tannotation\\n' > integrons.tsv
        printf 'ID_replicon\\tComplete\\tIn0\\tCALIN\\n' > summary.tsv
        exit 0
    fi

    # Main results: per-element annotations (integrase, attC, attI, gene cassettes)
    if [ -f "\${results_dir}/\${input_base}.integrons" ]; then
        cp "\${results_dir}/\${input_base}.integrons" integrons.tsv
    else
        printf 'ID_integron\\tID_replicon\\telement\\tpos_beg\\tpos_end\\tstrand\\tevalue\\ttype_elt\\tmodel\\ttype\\tannotation\\n' > integrons.tsv
    fi

    # Summary: counts of complete integrons, In0, and CALIN per replicon
    if [ -f "\${results_dir}/\${input_base}.summary" ]; then
        cp "\${results_dir}/\${input_base}.summary" summary.tsv
    else
        printf 'ID_replicon\\tComplete\\tIn0\\tCALIN\\n' > summary.tsv
    fi
    """
}

process ISLANDPATH_DIMOB {
    tag "islandpath"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-islandpath"
    publishDir "${params.outdir}/mge/islandpath", mode: 'copy'
    storeDir params.store_dir ? "${params.store_dir}/mge/islandpath" : null

    input:
    path(assembly)
    path(gff)
    path(faa)

    output:
    path("genomic_islands.tsv"), emit: islands

    script:
    def hmm_db = params.islandpath_hmm_db ? "${params.islandpath_hmm_db}/Pfam-A_mobgenes_201512_prok" : "${projectDir}/conda-envs/dana-mag-islandpath/opt/islandpath/hmmpfam/Pfam-A_mobgenes_201512_prok"
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
    conda "${projectDir}/conda-envs/dana-mag-macsyfinder"
    publishDir "${params.outdir}/mge/macsyfinder", mode: 'copy'
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
    conda "${projectDir}/conda-envs/dana-mag-defensefinder"
    publishDir "${params.outdir}/mge/defensefinder", mode: 'copy'
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
