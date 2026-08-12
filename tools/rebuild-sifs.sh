#!/usr/bin/env bash
set -uo pipefail

# ============================================================================
# rebuild-sifs.sh — pull every danaSeq pipeline image and place its .sif
# ============================================================================
#
# Each wrapper (run-<pipeline>.sh) looks for a dot-file beside itself:
#   <repo>/<component>/.danaseq-<name>.sif
# This script is what writes those. Rebuild all six, or name the ones you want:
#
#   ./tools/rebuild-sifs.sh                      # all six
#   ./tools/rebuild-sifs.sh mag_analysis         # just one
#   ./tools/rebuild-sifs.sh --local              # never submit to SLURM
#   ./tools/rebuild-sifs.sh --slurm              # submit everything as one job
#   ./tools/rebuild-sifs.sh --list               # show what would be built
#
# Why the SLURM split (the reason this script exists):
#
# On an HPC login node the obvious staging directory, /tmp, is tmpfs — it is
# RAM, and it is charged against your per-user memory cgroup. mksquashfs on a
# multi-GB image then trips that cap and the pull dies with an unhelpful
# "create command failed: signal: killed". The two largest images (live and
# mag_analysis) fail that way every time on a login node while the smaller
# four build there in ~2 minutes each.
#
# So by default the large ones are submitted to SLURM, where the job holds real
# RAM and staging goes to node-local disk. Override with --local or --slurm.
#
# Staging: $SLURM_TMPDIR when the site defines it (Compute Canada does; plain
# SLURM may not), otherwise a private dir under /tmp — which is correct on a
# compute node and is exactly the trap described above on a login node.
#
# Portability: needs bash 4+, apptainer or singularity, and — only for the
# SLURM path — sbatch. A host with just Docker needs none of this; the wrappers
# take --runtime docker and pull the image themselves.
#
# SLURM account:  --account, or $DANASEQ_SLURM_ACCOUNT; otherwise asked of
#                 SLURM (sacctmgr, then sshare). Never guessed.
# Job resources:  --time / --mem / --cpus, or $DANASEQ_SIFBUILD_{TIME,MEM,CPUS}.
# SLURM job log:  --log-dir, or $DANASEQ_SIFBUILD_LOGDIR, default $HOME. It must
# be on a shared filesystem — SLURM resolves --output on the compute node, so a
# log under /tmp lands on node-local disk and cannot be read afterwards.
# ============================================================================

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REGISTRY="ghcr.io/rec3141"

# Left empty on purpose: hardcoding one person's allocation into a shared repo
# means it silently charges the wrong account, or fails, for everyone else.
# Resolved from SLURM itself at submit time unless given.
ACCOUNT="${DANASEQ_SLURM_ACCOUNT:-}"
JOB_TIME="${DANASEQ_SIFBUILD_TIME:-03:00:00}"
JOB_MEM="${DANASEQ_SIFBUILD_MEM:-32G}"
JOB_CPUS="${DANASEQ_SIFBUILD_CPUS:-8}"

# name | image | component dir | sif basename | needs_slurm
IMAGES=(
    "nanopore_live|danaseq-live|nanopore_live|.danaseq-live.sif|1"
    "mag_analysis|danaseq-mag-analysis|mag_analysis|.danaseq-mag-analysis.sif|1"
    "nanopore_assembly|danaseq-nanopore-assembly|nanopore_assembly|.danaseq-nanopore-assembly.sif|0"
    "illumina_assembly|danaseq-illumina-assembly|illumina_assembly|.danaseq-illumina-assembly.sif|0"
    "illumina_rna|danaseq-illumina-rna|illumina_rna|.danaseq-illumina-rna.sif|0"
    "illumina_amplicon|danaseq-illumina-amplicon|illumina_amplicon|.danaseq-illumina-amplicon.sif|0"
)

MODE="auto"       # auto | local | slurm
LOG_DIR="${DANASEQ_SIFBUILD_LOGDIR:-$HOME}"
SELECTED=()

while (( $# )); do
    case "$1" in
        --local)   MODE="local"; shift ;;
        --slurm)   MODE="slurm"; shift ;;
        --account) ACCOUNT="$2"; shift 2 ;;
        --log-dir) LOG_DIR="$2"; shift 2 ;;
        --time)    JOB_TIME="$2"; shift 2 ;;
        --mem)     JOB_MEM="$2"; shift 2 ;;
        --cpus)    JOB_CPUS="$2"; shift 2 ;;
        --list)    MODE="list"; shift ;;
        -h|--help) sed -n '/^# ===/,/^# ===$/p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*)        echo "[ERROR] Unknown option: $1" >&2; exit 1 ;;
        *)         SELECTED+=("$1"); shift ;;
    esac
done

# Resolve which entries to build; an unknown name is an error, not a silent no-op.
entries=()
if (( ${#SELECTED[@]} == 0 )); then
    entries=("${IMAGES[@]}")
else
    for want in "${SELECTED[@]}"; do
        found=""
        for e in "${IMAGES[@]}"; do
            [[ "${e%%|*}" == "$want" ]] && { entries+=("$e"); found=1; break; }
        done
        [[ -n "$found" ]] || { echo "[ERROR] Unknown pipeline: $want" >&2
                               echo "        Known: $(printf '%s ' "${IMAGES[@]%%|*}")" >&2; exit 1; }
    done
fi

if [[ "$MODE" == "list" ]]; then
    printf "%-20s %-34s %s\n" PIPELINE IMAGE "DESTINATION"
    for e in "${entries[@]}"; do
        IFS='|' read -r name img comp sif _ <<< "$e"
        printf "%-20s %-34s %s\n" "$name" "${REGISTRY}/${img}:latest" "${comp}/${sif}"
    done
    exit 0
fi

RUNTIME=""
command -v apptainer   >/dev/null 2>&1 && RUNTIME=apptainer
[[ -z "$RUNTIME" ]] && command -v singularity >/dev/null 2>&1 && RUNTIME=singularity
if [[ -z "$RUNTIME" ]]; then
    module load apptainer >/dev/null 2>&1 || true
    command -v apptainer >/dev/null 2>&1 && RUNTIME=apptainer
fi
[[ -n "$RUNTIME" ]] || { echo "[ERROR] Neither apptainer nor singularity found." >&2; exit 1; }

# Inside a SLURM job, stage on the node's local disk. Otherwise /tmp, which is
# fine for the small images (see the header for why it is not, for the large).
STAGE="${SLURM_TMPDIR:-/tmp/$USER-sifbuild-$$}"
mkdir -p "$STAGE" || { echo "[ERROR] Cannot create staging dir $STAGE" >&2; exit 1; }
export APPTAINER_CACHEDIR="$STAGE/.cache" APPTAINER_TMPDIR="$STAGE/.tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
# Only clean up a directory we made ourselves; never $SLURM_TMPDIR.
[[ -n "${SLURM_TMPDIR:-}" ]] || trap 'rm -rf "$STAGE"' EXIT

rc=0

build_one() {
    local name="$1" img="$2" comp="$3" sif="$4"
    local staged="$STAGE/$sif" dest="$REPO/$comp/$sif"
    local t0=$SECONDS

    echo
    echo "================================================================"
    echo ">>> $(date -u +%H:%M:%SZ) [$(hostname -s)] $name"

    if [[ ! -d "$REPO/$comp" ]]; then
        echo "    [SKIP] no such component directory: $REPO/$comp" >&2
        rc=1; return 1
    fi

    rm -f "$staged"
    if ! "$RUNTIME" pull --force "$staged" "docker://${REGISTRY}/${img}:latest"; then
        echo "    [FAIL] pull failed: $img" >&2
        echo "           If this died with 'signal: killed', it ran out of memory —" >&2
        echo "           re-run with --slurm." >&2
        rc=1
        rm -rf "${APPTAINER_CACHEDIR:?}"/* "${APPTAINER_TMPDIR:?}"/* 2>/dev/null
        return 1
    fi
    if [[ ! -s "$staged" ]]; then
        echo "    [FAIL] empty image: $img" >&2; rc=1; return 1
    fi

    echo "    built in $((SECONDS-t0))s, $(du -h "$staged" | cut -f1)"

    # cp onto the existing file rather than mv/rename: that reuses the inode,
    # so this still works when the destination filesystem is over its inode
    # quota, which /project regularly is.
    if cp -f "$staged" "$dest"; then
        echo "    -> $dest"
    else
        echo "    [FAIL] could not write $dest" >&2; rc=1
    fi

    rm -f "$staged"
    rm -rf "${APPTAINER_CACHEDIR:?}"/* "${APPTAINER_TMPDIR:?}"/* 2>/dev/null
}

submit_slurm() {
    local names=("$@")
    command -v sbatch >/dev/null 2>&1 || {
        echo "[WARN] sbatch not found — building ${names[*]} locally instead." >&2
        echo "       Expect an out-of-memory failure on a login node." >&2
        for n in "${names[@]}"; do
            for e in "${IMAGES[@]}"; do
                IFS='|' read -r nm img comp sif _ <<< "$e"
                [[ "$nm" == "$n" ]] && build_one "$nm" "$img" "$comp" "$sif"
            done
        done
        return
    }

    # --output is resolved on the COMPUTE node, where /tmp is node-local — a log
    # written there is unreachable from the login node you submitted from. Use a
    # shared filesystem so the job's output can actually be read.
    local log="$LOG_DIR/danaseq-sifbuild-%j.log"

    # Ask SLURM for the account rather than assuming one.
    if [[ -z "$ACCOUNT" ]]; then
        ACCOUNT=$(sacctmgr -nP show user "$USER" format=DefaultAccount 2>/dev/null | head -1)
        [[ -n "$ACCOUNT" ]] || ACCOUNT=$(sshare -nP -U -o Account 2>/dev/null | head -1)
    fi
    if [[ -z "$ACCOUNT" ]]; then
        echo "[ERROR] No SLURM account found. Pass --account, or set \$DANASEQ_SLURM_ACCOUNT." >&2
        rc=1; return 1
    fi

    local jid
    jid=$(sbatch --parsable \
            --job-name=danaseq-sifbuild \
            --account="$ACCOUNT" \
            --time="$JOB_TIME" \
            --cpus-per-task="$JOB_CPUS" \
            --mem="$JOB_MEM" \
            --output="$log" \
            --wrap="bash '$REPO/tools/rebuild-sifs.sh' --local ${names[*]}")
    if [[ -n "$jid" ]]; then
        echo "[INFO] Submitted ${names[*]} as SLURM job $jid (account $ACCOUNT, $JOB_MEM)"
        echo "[INFO] Log: ${log/\%j/$jid}"
    else
        echo "[ERROR] sbatch failed for: ${names[*]}" >&2; rc=1
    fi
}

echo "=== danaSeq SIF rebuild — $(date -u +%FT%TZ) on $(hostname -s) ==="
echo "    repo:    $REPO"
echo "    runtime: $RUNTIME ($("$RUNTIME" --version 2>&1 | head -1))"
echo "    staging: $STAGE"

defer=()
for e in "${entries[@]}"; do
    IFS='|' read -r name img comp sif large <<< "$e"
    if [[ "$MODE" == "slurm" || ( "$MODE" == "auto" && "$large" == "1" && -z "${SLURM_JOB_ID:-}" ) ]]; then
        defer+=("$name")
        continue
    fi
    build_one "$name" "$img" "$comp" "$sif"
done

(( ${#defer[@]} )) && { echo; submit_slurm "${defer[@]}"; }

echo
echo "=== done rc=$rc — $(date -u +%FT%TZ) ==="
(( ${#defer[@]} )) && echo "    (deferred to SLURM: ${defer[*]} — check the job log)"
exit $rc
