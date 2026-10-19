// ============================================================================
// MAG dereplication
// ============================================================================
//
// DEREP_GALAH clusters the consensus bins (DAS Tool, plus Binette and MAGScoT
// when they ran) at species-level ANI with galah and keeps the best genome of
// each cluster, ranked by CheckM2 completeness and contamination.

process DEREP_GALAH {
    tag "galah"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-derep"
    publishDir "${params.outdir}/binning/derep", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/derep" : null

    input:
    path(bins_dirs, stageAs: 'bins_?')   // collected list of consensus bins/ directories
    path(checkm2_report)

    output:
    path("derep_clusters.tsv"), emit: clusters
    path("representatives/"),   emit: representatives

    script:
    """
    mkdir -p genomes representatives

    # One flat set of genomes, named by file stem (the CheckM2 'Name' column).
    for d in bins_*; do
        [ -d "\$d" ] || continue
        for f in "\$d"/*.fa; do
            [ -e "\$f" ] || continue
            n=\$(basename "\$f")
            if [ -e "genomes/\$n" ]; then
                echo "[WARNING] derep: \$n appears in more than one bin set; keeping the first" >&2
                continue
            fi
            ln -s "\$(readlink -f "\$f")" "genomes/\$n"
        done
    done

    # galah aborts on a genome missing from the quality report, so only bins
    # CheckM2 scored are clustered (it leaves out oversized bins).
    tail -n +2 ${checkm2_report} | cut -f1-3 | sort -u -k1,1 > quality.tsv
    : > genome_list.txt
    : > no_quality.txt
    ls genomes | awk '/\\.fa\$/ { print "genomes/" \$0 }' > all_genomes.txt
    awk -F'\\t' 'NR == FNR { q[\$1] = 1; next }
        { n = \$0; sub(/.*\\//, "", n); sub(/\\.fa\$/, "", n)
          if (n in q) print \$0 > "genome_list.txt"; else print n > "no_quality.txt" }' quality.tsv all_genomes.txt
    if [ -s no_quality.txt ]; then
        echo "[WARNING] derep: \$(wc -l < no_quality.txt) bins have no CheckM2 estimate and are left out" >&2
    fi

    # A fragment can still score as half-complete, and then be picked over a
    # whole genome as its cluster's representative, so bins shorter than
    # derep_min_size_bp are left out of clustering.
    : > too_small.txt
    if [ ${params.derep_min_size_bp} -gt 0 ] && [ -s genome_list.txt ]; then
        : > kept.txt
        while read -r g; do
            bp=\$(grep -v '>' "\$g" | tr -d '\\n\\r' | wc -c)
            if [ "\$bp" -lt ${params.derep_min_size_bp} ]; then
                n=\$(basename "\$g" .fa); echo "\$n" >> too_small.txt
            else
                echo "\$g" >> kept.txt
            fi
        done < genome_list.txt
        mv kept.txt genome_list.txt
        [ -s too_small.txt ] && echo "[INFO] derep: \$(wc -l < too_small.txt) bins under ${params.derep_min_size_bp} bp left out" >&2
    fi

    HEADER='genome\\trepresentative\\tani\\taf_genome\\taf_representative\\tcompleteness\\tcontamination\\tstatus'
    if [ ! -s genome_list.txt ]; then
        echo "[WARNING] derep: no consensus bins with CheckM2 estimates -- nothing to dereplicate" >&2
        printf "\$HEADER\\n" > derep_clusters.tsv
        awk -v OFS='\\t' '{ print \$1, "NA", "NA", "NA", "NA", "NA", "NA", "no_checkm2" }' no_quality.txt >> derep_clusters.tsv
        awk -F'\\t' -v OFS='\\t' 'NR == FNR { comp[\$1] = \$2; cont[\$1] = \$3; next }
            { print \$1, "NA", "NA", "NA", "NA", comp[\$1], cont[\$1], "too_small" }' quality.tsv too_small.txt >> derep_clusters.tsv
        exit 0
    fi
    echo "[INFO] derep: clustering \$(wc -l < genome_list.txt) bins at ${params.derep_ani}% ANI" >&2

    # galah ranks genomes within a cluster by completeness - 5 x contamination,
    # less small penalties for contig count and ambiguous bases
    # (Parks2020_reduced), and drops genomes below the quality thresholds.
    galah cluster \\
        --genome-fasta-list genome_list.txt \\
        --checkm2-quality-report ${checkm2_report} \\
        --min-completeness ${params.derep_min_completeness} \\
        --max-contamination ${params.derep_max_contamination} \\
        --ani ${params.derep_ani} \\
        --min-aligned-fraction ${params.derep_min_af} \\
        --threads ${task.cpus} \\
        --output-cluster-definition galah_clusters.tsv \\
        --output-representative-fasta-directory-copy galah_reps

    # galah writes no ANI values; skani measures each member against its
    # representative, one cluster per job.
    mkdir -p pairs
    awk -F'\\t' -v OFS='\\t' '\$1 != \$2' galah_clusters.tsv | sort -k1,1 | awk -F'\\t' '
        { r = \$1; sub(/.*\\//, "", r); f = "pairs/" r ".members"
          if (f != prev) { if (prev != "") close(prev); prev = f }
          print \$2 > f }'
    # xargs exits non-zero if any skani job fails, which fails the task.
    # skani writes no line for a pair below its 80% ANI screen; such members
    # get NA in the table.
    ls pairs | awk '/\\.members\$/ { print "pairs/" \$0 }' | xargs -r -P ${task.cpus} -I {} sh -c \\
        'm="\$1"; skani dist -t 1 -r "genomes/\$(basename "\$m" .members)" --ql "\$m" -o "\${m%.members}.ani" 2> "\${m%.members}.log" || { cat "\${m%.members}.log" >&2; exit 1; }' _ {}
    cat pairs/*.ani 2>/dev/null | awk -F'\\t' '\$1 != "Ref_file"' > member_ani.tsv || true

    # Cluster table: every consensus bin, with its representative and ANI.
    {
        printf "\$HEADER\\n"
        awk -F'\\t' -v OFS='\\t' '
            function stem(p) { sub(/.*\\//, "", p); sub(/\\.fa\$/, "", p); return p }
            FILENAME == "quality.tsv"      { comp[\$1] = \$2; cont[\$1] = \$3; next }
            FILENAME == "member_ani.tsv"   { k = stem(\$1) SUBSEP stem(\$2); ani[k] = \$3; afr[k] = \$4; afq[k] = \$5; next }
            FILENAME == "galah_clusters.tsv" { rep[stem(\$2)] = stem(\$1); next }
            FILENAME == "genome_list.txt"  {
                g = stem(\$1)
                if (!(g in rep)) { print g, "NA", "NA", "NA", "NA", comp[g], cont[g], "failed_quality"; next }
                r = rep[g]
                if (r == g) { print g, r, "100", "100", "100", comp[g], cont[g], "representative"; next }
                k = r SUBSEP g
                if (k in ani) print g, r, ani[k], afq[k], afr[k], comp[g], cont[g], "member"
                else          print g, r, "NA", "NA", "NA", comp[g], cont[g], "member"
            }
            FILENAME == "no_quality.txt"   { print \$1, "NA", "NA", "NA", "NA", "NA", "NA", "no_checkm2" }
            FILENAME == "too_small.txt"    { print \$1, "NA", "NA", "NA", "NA", comp[\$1], cont[\$1], "too_small" }
        ' quality.tsv member_ani.tsv galah_clusters.tsv genome_list.txt no_quality.txt too_small.txt
    } > derep_clusters.tsv.tmp

    if [ -d galah_reps ]; then
        find galah_reps -mindepth 1 -maxdepth 1 -exec mv {} representatives/ \\;
    fi
    mv derep_clusters.tsv.tmp derep_clusters.tsv

    awk -F'\\t' 'NR > 1 { s[\$8]++ } END {
        printf "[INFO] derep: %d representatives, %d members, %d below quality, %d too small, %d without CheckM2\\n",
            s["representative"], s["member"], s["failed_quality"], s["too_small"], s["no_checkm2"] }' derep_clusters.tsv >&2
    """
}
