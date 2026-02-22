// Eukaryotic analysis:
//   Tiara (deep learning k-mer NN) + Whokaryote (gene structure RF) — contig classification
//   MetaEuk — eukaryotic gene prediction (multi-exon, intron-aware, homology-based)

process TIARA_CLASSIFY {
    tag "tiara"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-tiara"
    publishDir "${params.outdir}/eukaryotic/tiara", mode: 'link'
    storeDir params.store_dir ? "${params.store_dir}/eukaryotic/tiara" : null

    input:
    path(contigs)

    output:
    path("tiara_output.tsv"),         emit: classifications
    path("tiara_probabilities.tsv"),  emit: probabilities, optional: true

    script:
    def min_len = params.tiara_min_len ?: 3000
    """
    set +e
    tiara \\
        -i "${contigs}" \\
        -o tiara_output.tsv \\
        --probabilities \\
        -t ${task.cpus} \\
        -m ${min_len}
    tiara_exit=\$?
    set -e

    if [ \$tiara_exit -ne 0 ]; then
        echo "[WARNING] Tiara exited with code \$tiara_exit" >&2
        printf 'sequence_id\\tclass_fst_stage\\tclass_snd_stage\\n' > tiara_output.tsv
    fi

    # Tiara writes probabilities to a file named after the output with .prob suffix
    # or --probabilities flag creates a separate file
    if [ -f "tiara_output.tsv.prob" ]; then
        mv tiara_output.tsv.prob tiara_probabilities.tsv
    elif [ -f "probabilities.tsv" ]; then
        mv probabilities.tsv tiara_probabilities.tsv
    fi

    n_classified=\$(tail -n +2 tiara_output.tsv | wc -l)
    echo "[INFO] Tiara: \${n_classified} contigs classified (min length: ${min_len} bp)" >&2
    """
}

process WHOKARYOTE_CLASSIFY {
    tag "whokaryote"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-whokaryote"
    publishDir "${params.outdir}/eukaryotic/whokaryote", mode: 'link'
    storeDir params.store_dir ? "${params.store_dir}/eukaryotic/whokaryote" : null

    input:
    path(contigs)
    path(gff)

    output:
    path("whokaryote_classifications.tsv"), emit: classifications

    script:
    def min_len = params.whokaryote_min_len ?: 5000
    def gff_arg = gff.name != 'NO_GFF' ? "--gff ${gff}" : ''
    """
    set +e
    whokaryote.py \\
        --contigs "${contigs}" \\
        --outdir whokaryote_out/ \\
        ${gff_arg} \\
        --minsize ${min_len} \\
        --threads ${task.cpus}
    whokaryote_exit=\$?
    set -e

    if [ \$whokaryote_exit -ne 0 ]; then
        echo "[WARNING] Whokaryote exited with code \$whokaryote_exit" >&2
        printf 'contig\\tprediction\\n' > whokaryote_classifications.tsv
    elif [ -f whokaryote_out/featuretable_predictions_T.tsv ]; then
        cp whokaryote_out/featuretable_predictions_T.tsv whokaryote_classifications.tsv
    else
        echo "[WARNING] Whokaryote produced no predictions" >&2
        printf 'contig\\tprediction\\n' > whokaryote_classifications.tsv
    fi

    n_classified=\$(tail -n +2 whokaryote_classifications.tsv | wc -l)
    echo "[INFO] Whokaryote: \${n_classified} contigs classified (min length: ${min_len} bp)" >&2
    """
}

process METAEUK_PREDICT {
    tag "metaeuk"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-metaeuk"
    publishDir "${params.outdir}/eukaryotic/metaeuk", mode: 'link'
    storeDir params.store_dir ? "${params.store_dir}/eukaryotic/metaeuk" : null

    input:
    path(contigs)
    path(tiara_tsv)
    path(whokaryote_tsv)

    output:
    path("metaeuk_proteins.fas"),       emit: proteins
    path("metaeuk_codon.fas"),          emit: codons
    path("metaeuk.gff"),                emit: gff
    path("metaeuk_headers.tsv"),        emit: headers

    script:
    def mem_limit = params.metaeuk_mem_limit ?: '50G'
    def min_len   = params.metaeuk_min_length ?: 20
    def max_intron = params.metaeuk_max_intron ?: 10000
    """
    # Union of Tiara non-prokaryotic + Whokaryote eukaryotic contigs
    # Tiara: keep eukarya, organelle, unknown (exclude bacteria/archaea/prokarya)
    awk -F'\\t' 'NR>1 && \$2!="bacteria" && \$2!="archaea" && \$2!="prokarya" {print \$1}' \\
        "${tiara_tsv}" > tiara_keep.txt
    # Whokaryote: keep anything classified as eukaryote
    awk -F'\\t' 'NR>1 && \$NF=="eukaryote" {print \$1}' \\
        "${whokaryote_tsv}" > whok_keep.txt
    # Union (deduplicated)
    sort -u tiara_keep.txt whok_keep.txt > keep_ids.txt

    n_total=\$(tail -n +2 "${tiara_tsv}" | wc -l)
    n_tiara=\$(wc -l < tiara_keep.txt)
    n_whok=\$(wc -l < whok_keep.txt)
    n_keep=\$(wc -l < keep_ids.txt)
    echo "[INFO] MetaEuk: Tiara non-prok: \${n_tiara}, Whokaryote euk: \${n_whok}, union: \${n_keep} of \${n_total} contigs" >&2

    # Extract matching contigs from assembly FASTA
    awk 'BEGIN{while((getline line < "keep_ids.txt")>0) ids[line]=1}
         /^>/{keep=(substr(\$1,2) in ids)} keep' "${contigs}" > filtered_contigs.fasta

    if [ ! -s filtered_contigs.fasta ]; then
        echo "[WARNING] No non-prokaryotic contigs to predict — skipping MetaEuk" >&2
        touch metaeuk_proteins.fas metaeuk_codon.fas metaeuk.gff metaeuk_headers.tsv
    else
        set +e
        metaeuk easy-predict \\
            filtered_contigs.fasta \\
            "${params.metaeuk_db}" \\
            metaeuk_out \\
            tmp_metaeuk \\
            --threads ${task.cpus} \\
            --split-memory-limit ${mem_limit} \\
            -e 100 \\
            --metaeuk-eval 0.0001 \\
            --metaeuk-tcov 0.6 \\
            --min-length ${min_len} \\
            --max-intron ${max_intron} \\
            --remove-tmp-files 1
        metaeuk_exit=\$?
        set -e

        if [ \$metaeuk_exit -ne 0 ]; then
            echo "[WARNING] MetaEuk exited with code \$metaeuk_exit" >&2
            touch metaeuk_proteins.fas metaeuk_codon.fas metaeuk.gff metaeuk_headers.tsv
        else
            # MetaEuk outputs: prefix.fas, prefix.codon.fas, prefix.gff, prefix.headersMap.tsv
            mv metaeuk_out.fas          metaeuk_proteins.fas
            mv metaeuk_out.codon.fas    metaeuk_codon.fas
            mv metaeuk_out.gff          metaeuk.gff
            mv metaeuk_out.headersMap.tsv metaeuk_headers.tsv
        fi
    fi

    n_proteins=\$(grep -c '^>' metaeuk_proteins.fas 2>/dev/null || echo 0)
    echo "[INFO] MetaEuk: \${n_proteins} eukaryotic proteins predicted" >&2
    """
}

process MARFERRET_CLASSIFY {
    tag "marferret"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-marferret"
    publishDir "${params.outdir}/eukaryotic/marferret", mode: 'link'
    storeDir params.store_dir ? "${params.store_dir}/eukaryotic/marferret" : null

    input:
    path(proteins)
    path(marferret_db)

    output:
    path("marferret_proteins.tsv"), emit: proteins
    path("marferret_contigs.tsv"),  emit: contigs

    script:
    """
    # Find MarFERReT database files
    DMND=\$(ls "${marferret_db}"/*.dmnd 2>/dev/null | head -1)
    TAX=\$(ls "${marferret_db}"/*.taxonomies.tab.gz 2>/dev/null | head -1)
    PFAM=\$(ls "${marferret_db}"/*.best_pfam_annotations.csv.gz 2>/dev/null | head -1)

    if [ -z "\$DMND" ]; then
        echo "[ERROR] No .dmnd file found in ${marferret_db}" >&2
        printf 'protein_id\\tcontig_id\\tbest_hit\\tpident\\taln_length\\tevalue\\tbitscore\\tstitle\\ttaxon_id\\ttaxonomy\\tpfam\\n' > marferret_proteins.tsv
        printf 'contig_id\\tn_proteins\\tn_classified\\ttop_taxonomy\\tpfam_domains\\n' > marferret_contigs.tsv
        exit 0
    fi

    # Check for empty protein input
    n_input=\$(grep -c '^>' "${proteins}" 2>/dev/null || echo 0)
    if [ "\$n_input" -eq 0 ]; then
        echo "[WARNING] MarFERReT: no input proteins — producing empty output" >&2
        printf 'protein_id\\tcontig_id\\tbest_hit\\tpident\\taln_length\\tevalue\\tbitscore\\tstitle\\ttaxon_id\\ttaxonomy\\tpfam\\n' > marferret_proteins.tsv
        printf 'contig_id\\tn_proteins\\tn_classified\\ttop_taxonomy\\tpfam_domains\\n' > marferret_contigs.tsv
        exit 0
    fi

    # DIAMOND blastp
    set +e
    diamond blastp \\
        --db "\$DMND" \\
        --query "${proteins}" \\
        --out diamond_hits.tsv \\
        --outfmt 6 qseqid sseqid pident length evalue bitscore stitle \\
        --max-target-seqs 1 \\
        --evalue 1e-5 \\
        --threads ${task.cpus} \\
        --block-size 4 \\
        --index-chunks 1
    diamond_exit=\$?
    set -e

    if [ \$diamond_exit -ne 0 ]; then
        echo "[WARNING] DIAMOND exited with code \$diamond_exit" >&2
        touch diamond_hits.tsv
    fi

    # Parse results: join with taxonomy + Pfam, aggregate per-contig
    TAX_ARG=""
    PFAM_ARG=""
    [ -n "\$TAX" ] && TAX_ARG="--taxonomy \$TAX"
    [ -n "\$PFAM" ] && PFAM_ARG="--pfam \$PFAM"

    parse_marferret_results.py \\
        --diamond diamond_hits.tsv \\
        \$TAX_ARG \\
        \$PFAM_ARG \\
        --out-proteins marferret_proteins.tsv \\
        --out-contigs marferret_contigs.tsv

    n_hits=\$(tail -n +2 marferret_proteins.tsv | wc -l)
    echo "[INFO] MarFERReT: \${n_hits} proteins classified from \${n_input} input" >&2
    """
}
