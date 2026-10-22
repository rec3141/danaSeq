#!/usr/bin/env bash
# protein_db_search.sh — search annotation proteins against one protein-level
# reference DB (HMM profile DB via hmmsearch, or a DIAMOND .dmnd) and write a
# self-describing hit TSV that import_all.py loads into dana.duckdb.
#
# Shared by modules/hmm.nf (HMM_SEARCH), modules/diamond_db.nf (DIAMOND_DB_SEARCH)
# and any out-of-pipeline backfill driver, so every route produces byte-for-byte
# the same output for the same inputs.
#
# Usage:
#   protein_db_search.sh hmm     <db.hmm>  <dbname> <out_prefix> <cpus> <hmm_cut>  <annot.tsv|NONE> <faa>...
#   protein_db_search.sh diamond <db.dmnd> <dbname> <out_prefix> <cpus> <dmnd_opts> <annot.tsv|NONE> <faa>...
#
#   hmm_cut    : tc | ga | nc (HMMER --cut_<x>), an E-value like 1e-10 (-E), or
#                'catalogue': per-target thresholds from the annot.tsv columns
#                `threshold` + `score_type` (full|domain), applied exactly like
#                kofam_scan (full-sequence score, or best 1-domain score)
#   dmnd_opts  : "id=50,qcov=50,scov=50,evalue=1e-5" (any subset; these are defaults)
#   annot.tsv  : optional gene catalogue, header
#                  target gene category pathway subpathway description [threshold score_type motif]
#                (target = HMM NAME / first '|' field of the DIAMOND subject id)
#
# Outputs:
#   <out_prefix>.tsv  line 1: '#' key=value provenance incl. n_proteins searched
#                     line 2: header
#                       gene_id db target gene category pathway subpathway evalue score pident qcov scov
#                     then one row per (protein, target) hit; pident/qcov/scov
#                     are empty for hmmsearch.
#   <out_prefix>.tbl  raw tool output (hmmsearch --tblout / DIAMOND outfmt 6)
set -euo pipefail

mode=$1; db=$2; dbname=$3; out=$4; cpus=$5; opts=$6; annot=$7
shift 7
[ $# -ge 1 ] || { echo "protein_db_search.sh: no protein files given" >&2; exit 2; }

# Merge the chunk's protein files (Bakta .faa / Prokka PROKKA_*.faa; .gz ok)
: > combined_proteins.faa
for f in "$@"; do
    case "$f" in
        *.gz) gzip -dc "$f" >> combined_proteins.faa ;;
        *)    cat "$f"      >> combined_proteins.faa ;;
    esac
done
n_prot=$(grep -c '^>' combined_proteins.faa || true)

case "$mode" in
hmm)
    case "$opts" in
        tc|ga|nc)  cut_args=(--cut_"$opts") ;;
        catalogue)
            [ "$annot" != NONE ] && [ -s "$annot" ] || { echo "protein_db_search.sh: hmm_cut=catalogue needs annot.tsv" >&2; exit 2; }
            # Report everything at/above the lowest per-target threshold; the
            # per-target filter below then applies each target's own cutoff.
            tmin=$(awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="threshold") c=i; next}
                                 c && $c != "" && (m == "" || $c+0 < m) {m = $c+0} END {print (m == "" ? 0 : m)}' "$annot")
            cut_args=(-T "$tmin" --domT "$tmin") ;;
        *)         cut_args=(-E "$opts") ;;
    esac
    if [ "$n_prot" -gt 0 ]; then
        hmmsearch "${cut_args[@]}" --noali --cpu "$cpus" \
            --tblout "$out.tbl" -o /dev/null "$db" combined_proteins.faa
    else
        : > "$out.tbl"
    fi
    # tblout: target(1) query(3) full-seq E-value(5) full score(6) best-1-domain score(9)
    grep -v '^#' "$out.tbl" | awk -v OFS='\t' -v mode="$opts" -v annot="$annot" '
        BEGIN {
            if (mode == "catalogue") {
                while ((getline line < annot) > 0) {
                    n = split(line, a, "\t")
                    if (a[1] == "target") { for (i = 1; i <= n; i++) col[a[i]] = i; continue }
                    thr[a[1]] = a[col["threshold"]]; st[a[1]] = a[col["score_type"]]
                }
            }
        }
        NF >= 9 {
            if (mode == "catalogue") {
                if (!($3 in thr) || thr[$3] == "") next
                sc = (st[$3] == "domain") ? $9 : $6
                if (sc + 0 < thr[$3] + 0) next
            }
            print $1, $3, $5, $6, "", "", ""
        }' > hits.raw || true
    params="method=hmmsearch	cutoff=$opts"
    ;;
diamond)
    id=50; qcov=50; scov=50; evalue=1e-5
    IFS=',' read -ra kv <<< "$opts"
    for p in "${kv[@]}"; do
        case "$p" in
            id=*) id=${p#id=} ;; qcov=*) qcov=${p#qcov=} ;;
            scov=*) scov=${p#scov=} ;; evalue=*) evalue=${p#evalue=} ;;
            "") ;; *) echo "protein_db_search.sh: unknown diamond option '$p'" >&2; exit 2 ;;
        esac
    done
    if [ "$n_prot" -gt 0 ]; then
        mkdir -p dmnd_tmp
        diamond blastp --quiet --threads "$cpus" --tmpdir dmnd_tmp \
            -q combined_proteins.faa -d "$db" -o "$out.tbl" \
            --id "$id" --query-cover "$qcov" --subject-cover "$scov" -e "$evalue" \
            --max-target-seqs 1 --max-hsps 1 \
            -f 6 qseqid sseqid pident evalue bitscore qcovhsp scovhsp
        rm -rf dmnd_tmp
    else
        : > "$out.tbl"
    fi
    # subject ids are 'ACCESSION|gene ...' -> target = ACCESSION
    awk -F'\t' -v OFS='\t' '{split($2, s, "|"); print $1, s[1], $4, $5, $3, $6, $7}' "$out.tbl" > hits.raw
    params="method=diamond_blastp	id=$id	qcov=$qcov	scov=$scov	evalue=$evalue"
    ;;
*)
    echo "protein_db_search.sh: unknown mode '$mode'" >&2; exit 2 ;;
esac

[ "$annot" = "NONE" ] || [ -s "$annot" ] || annot=NONE

# Optional sequence-motif check (catalogue column `motif`, an ERE such as the
# HgcA cap-helix N[VI][WY]C[AS][AGST]G): a hit on a target with a motif must
# contain it; otherwise it is relabelled <target>_nomotif when the catalogue
# lists that pseudo-target (so near-misses stay visible), else dropped.
if [ "$annot" != NONE ] && head -1 "$annot" | tr '\t' '\n' | grep -qx motif; then
    awk -F'\t' -v OFS='\t' -v annot="$annot" '
        BEGIN {
            while ((getline line < annot) > 0) {
                n = split(line, a, "\t")
                if (a[1] == "target") { for (i = 1; i <= n; i++) col[a[i]] = i; continue }
                known[a[1]] = 1
                if (col["motif"] && a[col["motif"]] != "") motif[a[1]] = a[col["motif"]]
            }
        }
        FNR == NR { if ($2 in motif) need[$1] = 1; next }           # pass 1: hits.raw
        FILENAME != "hits.raw" && /^>/ {                               # pass 2: proteins
            id = substr($0, 2); sub(/[ \t].*/, "", id); keep = (id in need); next
        }
        FILENAME != "hits.raw" { if (keep) seq[id] = seq[id] $0; next }
        END {
            while ((getline line < "hits.raw") > 0) {
                nf = split(line, h, "\t")
                if (h[2] in motif && seq[h[1]] !~ motif[h[2]]) {
                    if (!((h[2] "_nomotif") in known)) continue
                    h[2] = h[2] "_nomotif"
                    line = h[1]; for (i = 2; i <= nf; i++) line = line "\t" h[i]
                }
                print line
            }
        }' hits.raw combined_proteins.faa > hits.motif
    mv hits.motif hits.raw
fi
{
    printf '#db=%s\t%s\tn_proteins=%s\n' "$dbname" "$params" "$n_prot"
    printf 'gene_id\tdb\ttarget\tgene\tcategory\tpathway\tsubpathway\tevalue\tscore\tpident\tqcov\tscov\n'
    awk -F'\t' -v OFS='\t' -v db="$dbname" -v annot="$annot" '
        BEGIN {
            if (annot != "NONE") {
                while ((getline line < annot) > 0) {
                    split(line, a, "\t")
                    if (a[1] == "target") continue
                    ann[a[1]] = a[2] "\t" a[3] "\t" a[4] "\t" a[5]
                }
            }
        }
        {
            g = ($2 in ann) ? ann[$2] : $2 "\t\t\t"
            print $1, db, $2, g, $3, $4, $5, $6, $7
        }' hits.raw
} > "$out.tsv"
rm -f hits.raw combined_proteins.faa
