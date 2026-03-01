// RNA gene detection and classification
// rRNA: barrnap + vsearch against SILVA
// tRNA + tmRNA: Aragorn (fast, no Infernal dependency)
// Parsing logic in bin/parse_rrna_results.py and bin/parse_aragorn_results.py

process RNA_CLASSIFY {
    tag "rna"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-rrna"
    publishDir "${params.outdir}/taxonomy/rrna", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/taxonomy/rrna" : null

    input:
    path(contigs)   // Assembly FASTA — no annotation needed

    output:
    path("rrna_genes.tsv"),      emit: gene_classifications
    path("rrna_contigs.tsv"),    emit: contig_classifications
    path("rrna_sequences.fasta"), emit: sequences
    path("trna_genes.tsv"),      emit: trna_genes

    script:
    def silva_ssu = params.silva_ssu_db
    def silva_lsu = params.silva_lsu_db ?: ''
    def min_id    = params.rrna_min_identity
    """
    set -euo pipefail

    # ── A. barrnap rRNA detection ──────────────────────────────────────────
    # Run three kingdoms separately; --outseq extracts sequences directly
    for kingdom in bac arc euk; do
        set +e
        barrnap --kingdom \$kingdom --threads ${task.cpus} --reject 0.01 \
            --outseq rrna_\${kingdom}.fasta \
            "${contigs}" > barrnap_\${kingdom}.gff 2>barrnap_\${kingdom}.log
        barrnap_exit=\$?
        set -e

        if [ \$barrnap_exit -ne 0 ]; then
            echo "[WARNING] barrnap \$kingdom exited with code \$barrnap_exit" >&2
            cat barrnap_\${kingdom}.log >&2
            touch rrna_\${kingdom}.fasta
            printf '##gff-version 3\\n' > barrnap_\${kingdom}.gff
        fi
    done

    # Concatenate all extracted rRNA sequences
    cat rrna_bac.fasta rrna_arc.fasta rrna_euk.fasta > rrna_sequences.fasta

    # Check if any rRNA genes were found
    n_seqs=\$(grep -c '^>' rrna_sequences.fasta 2>/dev/null || echo 0)
    if [ "\$n_seqs" -eq 0 ]; then
        echo "[INFO] No rRNA genes detected in assembly" >&2
        printf 'gene_id\\tcontig_id\\tstart\\tend\\tstrand\\trrna_type\\tkingdom\\tgene_length\\tgene_completeness\\tbarrnap_score\\tbest_match\\tvsearch_identity\\tvsearch_coverage\\ttaxonomy\\n' > rrna_genes.tsv
        printf 'contig_id\\tn_rrna_genes\\trrna_types\\tkingdoms\\tbest_ssu_taxonomy\\tbest_ssu_identity\\tbest_lsu_taxonomy\\tbest_lsu_identity\\n' > rrna_contigs.tsv
    else
        echo "[INFO] Detected \$n_seqs rRNA gene(s) across all kingdoms" >&2

        # ── B. vsearch classification against SILVA SSU ────────────────────────
        for kingdom in bac arc euk; do
            if [ ! -s rrna_\${kingdom}.fasta ]; then
                touch vsearch_ssu_\${kingdom}.tsv
                continue
            fi

            set +e
            vsearch --usearch_global rrna_\${kingdom}.fasta \
                --db "${silva_ssu}" \
                --id ${min_id} \
                --maxaccepts 1 \
                --strand both \
                --threads ${task.cpus} \
                --blast6out vsearch_ssu_\${kingdom}.tsv \
                --top_hits_only \
                --output_no_hits \
                2>vsearch_ssu_\${kingdom}.log
            vsearch_exit=\$?
            set -e

            if [ \$vsearch_exit -ne 0 ]; then
                echo "[WARNING] vsearch SSU \$kingdom exited with code \$vsearch_exit" >&2
                cat vsearch_ssu_\${kingdom}.log >&2
                touch vsearch_ssu_\${kingdom}.tsv
            fi
        done

        # ── C. Optional: vsearch classification against SILVA LSU ──────────────
        if [ -n "${silva_lsu}" ] && [ -f "${silva_lsu}" ]; then
            for kingdom in bac arc euk; do
                if [ ! -s rrna_\${kingdom}.fasta ]; then
                    touch vsearch_lsu_\${kingdom}.tsv
                    continue
                fi

                set +e
                vsearch --usearch_global rrna_\${kingdom}.fasta \
                    --db "${silva_lsu}" \
                    --id ${min_id} \
                    --maxaccepts 1 \
                    --strand both \
                    --threads ${task.cpus} \
                    --blast6out vsearch_lsu_\${kingdom}.tsv \
                    --top_hits_only \
                    --output_no_hits \
                    2>vsearch_lsu_\${kingdom}.log
                vsearch_exit=\$?
                set -e

                if [ \$vsearch_exit -ne 0 ]; then
                    echo "[WARNING] vsearch LSU \$kingdom exited with code \$vsearch_exit" >&2
                    cat vsearch_lsu_\${kingdom}.log >&2
                    touch vsearch_lsu_\${kingdom}.tsv
                fi
            done
        else
            for kingdom in bac arc euk; do
                touch vsearch_lsu_\${kingdom}.tsv
            done
        fi

        # ── D. Build accession→taxonomy map from SILVA headers ──────────────
        # blast6out only stores the first word (accession), so we need a lookup
        grep '^>' "${silva_ssu}" | sed 's/^>//' | awk '{acc=\$1; \$1=""; sub(/^ /,""); print acc"\\t"\$0}' > silva_taxonomy.tsv
        if [ -n "${silva_lsu}" ] && [ -f "${silva_lsu}" ]; then
            grep '^>' "${silva_lsu}" | sed 's/^>//' | awk '{acc=\$1; \$1=""; sub(/^ /,""); print acc"\\t"\$0}' >> silva_taxonomy.tsv
        fi

        # ── E. Parse barrnap GFF + vsearch hits into per-gene and per-contig TSVs
        parse_rrna_results.py --taxonomy-map silva_taxonomy.tsv
    fi

    # ── F. Aragorn tRNA + tmRNA detection ─────────────────────────────────
    echo "[INFO] Running Aragorn for tRNA/tmRNA detection..." >&2
    aragorn -t -m -gcbact -w "${contigs}" > aragorn_raw.txt 2>/dev/null

    # Parse Aragorn batch output into TSV
    parse_aragorn_results.py aragorn_raw.txt trna_genes.tsv
    n_trna=\$(tail -n +2 trna_genes.tsv | grep -c tRNA || echo 0)
    n_tmrna=\$(tail -n +2 trna_genes.tsv | grep -c tmRNA || echo 0)
    echo "[INFO] Aragorn found \$n_trna tRNA and \$n_tmrna tmRNA genes" >&2
    """
}
