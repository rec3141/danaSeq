# Shared helpers for ablate-depth.sh and ablate-pairs.sh. Callers set TYPE, DEP,
# A (assembly_gsa_<type>), O (output dir) and SIF, and call container_setup first.
# subset_depth <out> <sample names...>: keep the given samples' mean/variance column
# pairs, in source order, and recompute totalAvgDepth the way the source table does.
# jgi_summarize (short reads) sums the per-sample means; nanopore's
# jgi_depth_table.awk (long reads) averages them.
case "$TYPE" in long) TOTAL=mean ;; short) TOTAL=sum ;; *) echo "type must be long or short" >&2; exit 1 ;; esac
subset_depth() {
    local out=$1; shift
    awk -F'\t' -v keep="$*" -v total="$TOTAL" 'BEGIN{OFS="\t"; n=split(keep,k," "); for(i=1;i<=n;i++) want[k[i]]=1}
        NR==1{ printf "%s\t%s\t%s", $1,$2,$3
               for(i=4;i<=NF;i+=2){ s=$i; sub(/\.sorted(\.bam)?$/,"",s); if(s in want){cols[++m]=i; printf "\t%s\t%s",$i,$(i+1)} }
               print ""; if(m!=n){print "[ERROR] matched "m" of "n" samples" > "/dev/stderr"; exit 3}; next }
        { t=0; for(j=1;j<=m;j++) t+=$(cols[j]); if(total=="mean") t/=m; printf "%s\t%s\t%s", $1,$2,t
          for(j=1;j<=m;j++) printf "\t%s\t%s", $(cols[j]), $(cols[j]+1); print "" }' "$DEP" > "$out"
}

run_metabat() {  # run_metabat <config> [depth file]
    local cfg=$1 dep=${2:-}
    local d=$O/binning/$cfg; mkdir -p "$d"
    rm -f "$d"/bin.*.fa
    cexec "$SIF" metabat2 -i "$A/assembly/assembly.fasta" ${dep:+-a "$dep"} -m 1500 --seed 1 \
        -t "$SLURM_CPUS_PER_TASK" -o "$d/bin"
    # The bin FASTAs are the authoritative membership: this MetaBAT2 build's
    # --saveCls matrix holds internal cluster IDs (including clusters below
    # minClsSize), not bin numbers. Keep only the contig<TAB>bin table.
    local f n
    for f in "$d"/bin.*.fa; do
        [ -e "$f" ] || continue
        n=${f##*/bin.}; n=${n%.fa}
        awk -v b="$(printf '%s_%04d' "$cfg" "$n")" '/^>/{sub(/^>/,""); sub(/[ \t].*/,""); print $0"\t"b}' "$f"
    done > "$d/${cfg}_bins.tsv"
    rm -f "$d"/bin.*.fa
    [ -s "$d/${cfg}_bins.tsv" ] || { echo "[ERROR] $cfg: no bins written" >&2; exit 1; }
    echo "$cfg: $(cut -f2 "$d/${cfg}_bins.tsv" | sort -u | wc -l) bins, $(wc -l < "$d/${cfg}_bins.tsv") contigs binned"
}
