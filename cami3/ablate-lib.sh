# Shared helpers for ablate-depth.sh and ablate-pairs.sh. Callers set TYPE, DEP,
# A (assembly_gsa_<type>), O (output dir) and SIF, and call container_setup first.
# subset_depth <out> <sample names...>: keep the given samples' mean/variance column
# pairs, in source order, and recompute totalAvgDepth the way the source table does.
# jgi_summarize (short reads) sums the per-sample means; nanopore's
# jgi_depth_table.awk (long reads) averages them.
case "$TYPE" in long) TOTAL=mean ;; short) TOTAL=sum ;; *) echo "type must be long or short" >&2; exit 1 ;; esac
subset_depth() {
    local out=$1; shift
    local tmp
    tmp=$(mktemp -d -- "${out}.XXXXXX") || return 1
    if awk -F'\t' -v keep="$*" -v total="$TOTAL" '
        function fail(message) { print "[ERROR] " message > "/dev/stderr"; exit 3 }
        function depth(value) { return value ~ /^[+]?[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?$/ }
        BEGIN {
            OFS="\t"
            n=split(keep, samples, " ")
            if (!n) fail("no samples requested")
            for (i=1; i<=n; i++) {
                s=samples[i]
                if (s in want) fail("duplicate requested sample " s)
                want[s]=1
            }
        }
        NR==1 {
            width=NF
            if ($1!="contigName" || $2!="contigLen" || $3!="totalAvgDepth" || NF<5 || (NF-3)%2)
                fail("invalid depth table header")
            printf "%s\t%s\t%s", $1,$2,$3
            for (i=4; i<=NF; i+=2) {
                if ($(i+1) != $i "-var") fail("missing variance column for " $i)
                s=$i
                sub(/\.sorted(\.bam)?$/, "", s)
                if (s in seen) fail("duplicate source sample " s)
                seen[s]=1
                if (s in want) { cols[++m]=i; printf "\t%s\t%s", $i,$(i+1) }
            }
            for (s in want) if (!(s in seen)) fail("missing sample " s)
            print ""
            next
        }
        {
            if (NF != width) fail("incorrect field count at row " NR)
            t=0
            for (j=1; j<=m; j++) {
                i=cols[j]
                if (!depth($i) || !depth($(i+1))) fail("invalid selected depth at row " NR)
                t+=$i
            }
            if (total=="mean") t/=m
            printf "%s\t%s\t%s", $1,$2,t
            for (j=1; j<=m; j++) printf "\t%s\t%s", $(cols[j]),$(cols[j]+1)
            print ""
        }
        END { if (!NR) fail("empty depth table") }
    ' "$DEP" > "$tmp/table" && mv -f -- "$tmp/table" "$out"; then
        rmdir -- "$tmp"
        return 0
    else
        local status=$?
        rm -f -- "$tmp/table"
        rmdir -- "$tmp"
        return "$status"
    fi
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
