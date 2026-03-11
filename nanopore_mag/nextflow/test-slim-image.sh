#!/usr/bin/env bash
# ============================================================================
# Test suite for danaseq-mag-base-slim container image
# Usage: bash test-slim-image.sh [path-to-sif]
# ============================================================================
set -euo pipefail

SIF="${1:-danaseq-mag-base-slim.sif}"
if [[ ! -f "$SIF" ]]; then
    echo "ERROR: SIF not found: $SIF"
    exit 1
fi

PASS=0
FAIL=0
SKIP=0
FAILURES=""

run_test() {
    local env="$1" tool="$2" cmd="$3"
    printf "  %-35s " "$tool"
    # Use the wrapper path so we test the dana-tools wrappers
    output=$(timeout 30 apptainer exec "$SIF" bash -c "$cmd" 2>&1) && rc=0 || rc=$?
    if [[ $rc -eq 124 ]]; then
        printf "TIMEOUT\n"
        FAILURES+="  $env/$tool: timed out after 30s\n"
        ((FAIL++))
        return
    elif [[ $rc -eq 0 ]]; then
        # Extract version from first line of output
        ver=$(echo "$output" | head -1 | grep -oP '[\d]+\.[\d]+[.\d]*' | head -1)
        printf "PASS  %s\n" "${ver:-ok}"
        ((PASS++))
    else
        printf "FAIL  (exit %d)\n" "$rc"
        FAILURES+="  $env/$tool: $(echo "$output" | tail -1)\n"
        ((FAIL++))
    fi
}

run_test_stderr() {
    # Some tools print version to stderr or exit nonzero for --version
    local env="$1" tool="$2" cmd="$3"
    printf "  %-35s " "$tool"
    output=$(timeout 30 apptainer exec "$SIF" bash -c "$cmd" 2>&1) || true
    if echo "$output" | grep -qiP '(version|v?\d+\.\d+|usage|help)'; then
        ver=$(echo "$output" | grep -oP '[\d]+\.[\d]+[.\d]*' | head -1)
        printf "PASS  %s\n" "${ver:-ok}"
        ((PASS++))
    else
        printf "FAIL\n"
        FAILURES+="  $env/$tool: $(echo "$output" | tail -1)\n"
        ((FAIL++))
    fi
}

echo "============================================"
echo "Testing slim image: $SIF"
echo "Size: $(du -h "$SIF" | cut -f1)"
echo "============================================"
echo

# ── 1. Assembly env ──────────────────────────────────────────────────────
echo "── Assembly (dana-mag-assembly) ──"
run_test assembly flye             "flye --version"
run_test assembly filtlong         "filtlong --version"
run_test assembly minimap2         "minimap2 --version"
run_test assembly samtools         "samtools --version | head -1"
run_test assembly coverm           "coverm --version"
run_test assembly bbmap.sh         "bbmap.sh --version 2>&1 | head -2"
run_test assembly bbduk.sh         "bbduk.sh --version 2>&1 | head -2"
run_test assembly dedupe.sh        "dedupe.sh --version 2>&1 | head -2"
run_test assembly reformat.sh      "reformat.sh --version 2>&1 | head -2"
run_test assembly pigz             "pigz --version"
run_test_stderr assembly metaMDBG  "metaMDBG --help | head -1"
run_test_stderr assembly myloasm   "myloasm --help | head -1"
run_test assembly nextflow         "nextflow -version 2>&1 | grep version"
run_test assembly java             "java -version 2>&1 | head -1"
echo

# ── 2. Binning env ──────────────────────────────────────────────────────
echo "── Binning (dana-mag-binning) ──"
run_test_stderr binning metabat2              "metabat2 --help 2>&1 | head -3"
run_test_stderr binning "run_MaxBin.pl"       "run_MaxBin.pl -v 2>&1 | head -1"
run_test_stderr binning DAS_Tool              "DAS_Tool --version 2>&1 | head -1"
run_test_stderr binning jgi_summarize         "jgi_summarize_bam_contig_depths --help 2>&1 | head -3"
run_test binning SemiBin2                     "SemiBin2 --version"
run_test_stderr binning LorBin                "python3 -c 'import lorbin; print(\"lorbin ok\")' 2>&1"
echo

# ── 3. Quality env ──────────────────────────────────────────────────────
echo "── Quality (dana-mag-quality) ──"
run_test quality genomad    "PATH=/opt/conda/envs/dana-mag-quality/bin:\$PATH genomad --version"
run_test quality checkm2    "PATH=/opt/conda/envs/dana-mag-quality/bin:\$PATH checkm2 --version"
run_test_stderr quality checkv     "PATH=/opt/conda/envs/dana-mag-quality/bin:\$PATH checkv -h 2>&1 | head -3"
run_test quality prodigal   "PATH=/opt/conda/envs/dana-mag-quality/bin:\$PATH prodigal -v 2>&1 | head -1"
run_test quality whokaryote "PATH=/opt/conda/envs/dana-mag-quality/bin:\$PATH python3 -c 'import whokaryote; print(\"whokaryote ok\")'"
run_test quality tiara      "PATH=/opt/conda/envs/dana-mag-quality/bin:\$PATH tiara --help 2>&1 | head -1"
echo

# ── 4. Annotate env ─────────────────────────────────────────────────────
echo "── Annotate (dana-mag-annotate) ──"
run_test annotate bakta     "PATH=/opt/conda/envs/dana-mag-annotate/bin:\$PATH bakta --version"
run_test annotate emapper   "PATH=/opt/conda/envs/dana-mag-annotate/bin:\$PATH emapper.py --version"
run_test annotate diamond   "PATH=/opt/conda/envs/dana-mag-annotate/bin:\$PATH diamond version"
echo

# ── 5. Classify env ─────────────────────────────────────────────────────
echo "── Classify (dana-mag-classify) ──"
run_test_stderr classify kaiju     "kaiju 2>&1 | head -3"
run_test_stderr classify kraken2   "kraken2 --version"
run_test classify barrnap          "barrnap --version"
run_test classify vsearch          "vsearch --version 2>&1 | head -1"
run_test_stderr classify aragorn   "aragorn -h 2>&1 | head -3"
run_test_stderr classify exec_annotation  "PATH=/opt/conda/envs/dana-mag-classify/bin:\$PATH exec_annotation --help 2>&1 | head -3"
run_test classify hmmpress         "hmmpress -h 2>&1 | head -3"
run_test_stderr classify metaeuk   "metaeuk -h 2>&1 | head -3"
echo

# ── 6. Genomic env ──────────────────────────────────────────────────────
echo "── Genomic (dana-mag-genomic) ──"
run_test genomic defense-finder    "PATH=/opt/conda/envs/dana-mag-genomic/bin:\$PATH defense-finder --version 2>&1"
run_test_stderr genomic macsyfinder "PATH=/opt/conda/envs/dana-mag-genomic/bin:\$PATH macsyfinder --version 2>&1"
run_test_stderr genomic integron_finder "PATH=/opt/conda/envs/dana-mag-genomic/bin:\$PATH integron_finder --version 2>&1"
run_test_stderr genomic run_dbcan  "PATH=/opt/conda/envs/dana-mag-genomic/bin:\$PATH run_dbcan --help 2>&1 | head -3"
run_test genomic hmmscan           "PATH=/opt/conda/envs/dana-mag-genomic/bin:\$PATH hmmscan -h 2>&1 | head -3"
run_test genomic hmmsearch         "PATH=/opt/conda/envs/dana-mag-genomic/bin:\$PATH hmmsearch -h 2>&1 | head -3"
echo

# ── 7. Derep env ────────────────────────────────────────────────────────
echo "── Derep (dana-mag-derep) ──"
run_test derep galah       "PATH=/opt/conda/envs/dana-mag-derep/bin:\$PATH galah --version"
run_test derep skani       "PATH=/opt/conda/envs/dana-mag-derep/bin:\$PATH skani --version"
run_test derep sourmash    "PATH=/opt/conda/envs/dana-mag-derep/bin:\$PATH sourmash --version"
run_test_stderr derep skder "PATH=/opt/conda/envs/dana-mag-derep/bin:\$PATH skder --help 2>&1 | head -3"
run_test derep dRep        "PATH=/opt/conda/envs/dana-mag-derep/bin:\$PATH dRep --version"
echo

# ── 8. Strain env ───────────────────────────────────────────────────────
echo "── Strain (dana-mag-strain) ──"
run_test_stderr strain strainy     "PATH=/opt/conda/envs/dana-mag-strain/bin:\$PATH strainy --help 2>&1 | head -3"
run_test_stderr strain floria      "PATH=/opt/conda/envs/dana-mag-strain/bin:\$PATH floria --help 2>&1 | head -3"
run_test strain inStrain           "PATH=/opt/conda/envs/dana-mag-strain/bin:\$PATH inStrain --version 2>&1"
echo

# ── 9. PathViz env ──────────────────────────────────────────────────────
echo "── PathViz (dana-mag-pathviz) ──"
run_test pathviz KEGG-decoder "PATH=/opt/conda/envs/dana-mag-pathviz/bin:\$PATH KEGG-decoder --help 2>&1 | head -3"
run_test pathviz node         "PATH=/opt/conda/envs/dana-mag-pathviz/bin:\$PATH node --version"
run_test pathviz python3      "PATH=/opt/conda/envs/dana-mag-pathviz/bin:\$PATH python3 --version"
echo

# ── Summary ─────────────────────────────────────────────────────────────
echo "============================================"
echo "Results: $PASS passed, $FAIL failed, $SKIP skipped"
echo "============================================"
if [[ $FAIL -gt 0 ]]; then
    echo
    echo "Failures:"
    printf "$FAILURES"
    echo
    exit 1
fi
