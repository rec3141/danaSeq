#!/usr/bin/env bash
# Emit the ACTUAL version of every tool this pipeline runs, as JSON, by asking each
# binary rather than reading the environment spec.
#
# The spec says what was requested ("cutadapt", "r-base >=4.1"); only the binary knows
# what was resolved. A dataset analysed months later needs the second one — a claim that
# rests on DADA2's error model is not reproducible against "bioconductor-dada2".
#
# Never fails: a tool that is absent or does not answer is reported as null, because a
# manifest that is missing is worse than one with a hole in it.
set -uo pipefail

version_line() {
    # First line containing something version-shaped, not simply the first line.
    # This container's nextflow prints "Illegal option --" ahead of its real output —
    # for BOTH -v and -version, so it is not the flag — and that error was being
    # recorded as the version. `-v` gives the terser single line of the two.
    tr -d '\r' | grep -m1 -E '[0-9]+\.[0-9]+' | sed 's/^[[:space:]]*//;s/[[:space:]]*$//'
}
# Most tools print to stdout; a few (FastTree, mafft) use stderr.
probe() {
    # `bin` is the first word of the command; if the shell cannot find it the tool is
    # ABSENT, which is a different fact from "present but does not report a version".
    # Merging stderr made "command not found" look like a version banner.
    local bin=${1%% *}
    command -v "$bin" >/dev/null 2>&1 || { printf '__ABSENT__'; return; }
    { eval "$1"; } 2>&1 | version_line
}

# tool -> command, and a regex to pull the bare version out of a chatty banner
emit() {
    local name="$1" raw="$2" ver
    ver=$(printf '%s' "$raw" | grep -oE '[0-9]+\.[0-9]+(\.[0-9]+)?([a-z0-9._-]*)?' | head -n1)
    if [ -z "$raw" ] || [ "$raw" = "__ABSENT__" ]; then
        printf '  "%s": null' "$name"
    else
        printf '  "%s": {"version": %s, "reported": "%s"}' "$name" \
            "$([ -n "$ver" ] && printf '"%s"' "$ver" || printf 'null')" \
            "$(printf '%s' "$raw" | sed 's/\\/\\\\/g; s/"/\\"/g' | cut -c1-160)"
    fi
}

# Python packages installed from a git URL report a version string that never changes —
# papa2 has said "0.1.0" across every commit — but PEP 610 writes the actual commit
# into direct_url.json. That commit is the only real answer for papa2, which is the
# tool most likely to change the numbers.
python_packages() {
    python - <<'PYEOF' 2>/dev/null || printf 'null'
import json
try:
    import importlib.metadata as md
except ImportError:
    print("null"); raise SystemExit
out = {}
for name in ("papa2", "biopython", "scikit-learn", "scipy",
             "numpy", "pandas", "cutadapt", "matplotlib"):
    try:
        d = md.distribution(name)
    except Exception:
        continue
    rec = {"version": d.version, "commit": None, "url": None}
    try:
        o = json.loads(d.read_text("direct_url.json") or "{}")
        rec["url"] = o.get("url")
        rec["commit"] = (o.get("vcs_info") or {}).get("commit_id")
    except Exception:
        pass
    out[name] = rec
print(json.dumps(out))
PYEOF
}

printf '{\n'
printf '  "python_packages": %s,\n' "$(python_packages)"
sep=""
for spec in \
    "cutadapt|cutadapt --version" \
    "mafft|mafft --version" \
    "fasttree|FastTree -expert" \
    "python|python --version" \
    "R|R --version" \
    "nextflow|nextflow -v" \
    "dada2|Rscript -e 'cat(as.character(packageVersion(\"dada2\")))'" \
    "decipher|Rscript -e 'cat(as.character(packageVersion(\"DECIPHER\")))'" \
    "ape|Rscript -e 'cat(as.character(packageVersion(\"ape\")))'" \
    "biopython|python -c 'import Bio;print(Bio.__version__)'" \
    "scikit_learn|python -c 'import sklearn;print(sklearn.__version__)'" \
    "scipy|python -c 'import scipy;print(scipy.__version__)'" \
    "papa2|python -c 'import papa2;print(getattr(papa2,\"__version__\",\"unknown\"))'" \
; do
    name="${spec%%|*}"; cmd="${spec#*|}"
    printf '%s' "$sep"; emit "$name" "$(probe "$cmd")"; sep=$',\n'
done
printf '\n}\n'
