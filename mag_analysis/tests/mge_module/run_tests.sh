#!/usr/bin/env bash
# Module tests for the MGE steps with fake geNomad/CheckV. Needs nextflow on PATH (or NEXTFLOW=...).
# Usage: tests/mge_module/run_tests.sh   (from anywhere; uses a temporary work tree)
set -uo pipefail
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
NF="${NEXTFLOW:-nextflow}"
export PATH="$here/fakebin:$here/../../bin:$PATH"
tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' EXIT
cd "$tmp"
export NXF_HOME="$tmp/.nxf"
printf '>host1\n%s\n>host2\n%s\n>plas1\n%s\n' "$(head -c 10000 /dev/zero | tr '\0' A)" \
  "$(head -c 5000 /dev/zero | tr '\0' C)" "$(head -c 3000 /dev/zero | tr '\0' G)" | gzip > asm.fa.gz
printf 'host1\tbinA\nhost2\tbinA\nplas1\tbinA\n' > c2b.tsv
fails=0
check() { if eval "$2"; then echo "PASS $1"; else echo "FAIL $1"; fails=$((fails+1)); fi; }
run() {  # name, extra args...; returns nextflow exit status
  local name="$1"; shift
  "$NF" -q run "$here/test.nf" -work-dir "$tmp/work_$name" --assembly asm.fa.gz --contig2bin c2b.tsv \
      --outdir "$tmp/out_$name" "$@" > "$tmp/$name.log" 2>&1
}
st() { awk -F'\t' 'NR==2{print $6}' "$1"; }

FAKE_GENOMAD_MODE=ok run ok --mge_alt_export --mge_export_min_plasmid_score 0.9 --mge_export_min_virus_score 0.9; rc=$?
check "ok: workflow succeeds on .fa.gz" "[ $rc -eq 0 ]"
check "ok: geNomad status ok" "[ \"\$(st out_ok/mge/genomad/status.tsv)\" = ok ]"
check "ok: CheckV status ok" "[ \"\$(st out_ok/mge/checkv/status.tsv)\" = ok ]"
check "ok: plasmid called" "grep -q \$'^plas1\t3000\tok\t.*\tplasmid\t' out_ok/mge/annotate/contig_class.tsv"
check "ok: provirus region + integration link" "grep -q predicted_integration out_ok/mge/annotate/host_links.tsv"
check "ok: CheckV joined on region" "grep -q \$'provirus_101_2100.*\t91.0\tAAI-based' out_ok/mge/annotate/mge_groups.tsv"
check "ok: plasmid peeled, host1 stays" "grep -q \$'^plas1\tmge_plas1\tpeel' out_ok/mge/annotate/bins_with_mge_separate.tsv && grep -q \$'^host1\tbinA\tkept' out_ok/mge/annotate/bins_with_mge_separate.tsv"

FAKE_GENOMAD_MODE=nohits run nohits; rc=$?
check "nohits: workflow succeeds" "[ $rc -eq 0 ]"
check "nohits: status no_hits" "[ \"\$(st out_nohits/mge/genomad/status.tsv)\" = no_hits ]"
check "nohits: calls not_scored" "[ \"\$(cut -f7 out_nohits/mge/annotate/contig_class.tsv | tail -n +2 | sort -u)\" = not_scored ]"

FAKE_GENOMAD_MODE=fail run fail --run_checkv false; rc=$?
check "fail: workflow succeeds without export" "[ $rc -eq 0 ]"
check "fail: status failed" "[ \"\$(st out_fail/mge/genomad/status.tsv)\" = failed ]"
check "fail: calls unavailable" "[ \"\$(cut -f7 out_fail/mge/annotate/contig_class.tsv | tail -n +2 | sort -u)\" = unavailable ]"

FAKE_GENOMAD_MODE=fail run failexport --run_checkv false --mge_alt_export --mge_export_min_plasmid_score 0.9 --mge_export_min_virus_score 0.9; rc=$?
check "fail+export: workflow stops (terminate)" "[ $rc -ne 0 ]"

FAKE_GENOMAD_MODE=missing_provirus run missingpv --run_checkv false; rc=$?
check "missing provirus.tsv: status failed" "[ \"\$(st out_missingpv/mge/genomad/status.tsv)\" = failed ]"

FAKE_GENOMAD_MODE=ok FAKE_CHECKV_MODE=fail run cvfail; rc=$?
check "checkv fail: status failed, classification kept" "[ \"\$(st out_cvfail/mge/checkv/status.tsv)\" = failed ] && grep -q \$'\tplasmid\t' out_cvfail/mge/annotate/contig_class.tsv"

echo "module tests: $fails failure(s)"; [ $fails -eq 0 ]
