#!/usr/bin/env bash
set -euo pipefail

# ============================================================================
# Dana MAG Analysis Pipeline Launcher
# ============================================================================
#
# Runs the technology-agnostic downstream analysis pipeline on pre-computed
# assembly outputs. Accepts assembly FASTA + depth table from any assembler
# (nanopore_assembly, illumina_assembly, or external).
#
# Usage:
#   ./run-mag-analysis.sh --assembly /path/to/assembly.fasta \
#       --depths /path/to/depths.txt --outdir /path/to/output [options]
#
# ============================================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTAINER_IMAGE="ghcr.io/rec3141/danaseq-mag-analysis:latest"
USE_CONTAINER=false
CONTAINER_RUNTIME=""
SIF_PATH=""
PULL_SIF=false
NF_ARGS=()
DB_ARGS=()
# Squashed database images (<db_dir>/<name>.sqsh) mounted read-only into the
# container in place of the plain directory <db_dir>/<name with - -> _>, e.g.
# bakta.sqsh -> bakta/, gtdbtk-db.sqsh -> gtdbtk_db/. One file instead of
# tens of thousands keeps the scratch inode quota in check. Apptainer/
# Singularity only (bind ... :image-src=/); Docker falls back to the plain dir.
SQSH_IMAGES=()   # entries: "<host path to .sqsh>|<subdir name>"
MOUNTS=()
ORIGINAL_ARGS=("$@")

die() { echo "[ERROR] $1" >&2; exit 1; }

# Resolve database paths from a standard base directory
resolve_db_dir() {
    local base="$1"
    local path

    local -A fixed_dbs=(
        [--bakta_db]="bakta/db"
        [--bakta_light_db]="bakta/db-light"
        [--genomad_db]="genomad_db"
        [--checkv_db]="checkv_db"
        [--kofam_db]="kofam_db"
        [--eggnog_db]="eggnog_db"
        [--dbcan_db]="dbcan_db"
        [--macsyfinder_models]="macsyfinder_models"
        [--defensefinder_models]="defensefinder_models"
        [--marferret_db]="marferret_db"
        [--gtdbtk_db]="gtdbtk_db"
        [--metaeuk_db]="metaeuk_db/metaeuk_db"
        [--antismash_db]="antismash_db"
        [--magscot_hmm_dir]="magscot_hmm"
    )
    # Squashed images take precedence over plain dirs when a container runtime
    # can mount them (decided later; Docker gets the plain dir if present).
    local img sub
    for img in "${base}"/*.sqsh; do
        [[ -f "$img" ]] || continue
        sub="$(basename "${img%.sqsh}")"; sub="${sub//-/_}"
        SQSH_IMAGES+=("${img}|${sub}")
    done
    sqsh_covers() {   # $1 = path under $base; true if a squashed image provides it
        local rel="${1#${base}/}" e
        for e in "${SQSH_IMAGES[@]}"; do
            [[ "$rel" == "${e#*|}" || "$rel" == "${e#*|}/"* ]] && return 0
        done
        return 1
    }

    for flag in "${!fixed_dbs[@]}"; do
        path="${base}/${fixed_dbs[$flag]}"
        if [[ -e "$path" ]] || sqsh_covers "$path"; then DB_ARGS+=("$flag" "$path"); fi
    done

    # Names written by download-databases.sh first, older layouts after.
    # CheckM2
    for path in "${base}/checkm2_db" "${base}/checkm2"; do
        [[ -e "$path" ]] && { DB_ARGS+=(--checkm2_db "$path"); break; }
    done

    # Kaiju
    if [[ -f "${base}/kaiju_db/nodes.dmp" ]]; then
        DB_ARGS+=(--kaiju_db "${base}/kaiju_db")
    elif [[ -d "${base}/kaiju/refseq_ref" ]]; then
        DB_ARGS+=(--kaiju_db "${base}/kaiju/refseq_ref")
    elif [[ -d "${base}/kaiju" ]]; then
        path=$(ls -d "${base}/kaiju"/*/2>/dev/null | head -1 || true)
        [[ -n "$path" ]] && DB_ARGS+=(--kaiju_db "$path")
    fi

    # Kraken2
    if [[ -f "${base}/kraken2_db/hash.k2d" ]]; then
        DB_ARGS+=(--kraken2_db "${base}/kraken2_db")
    elif [[ -d "${base}/krakendb/pluspfp_08gb" ]]; then
        DB_ARGS+=(--kraken2_db "${base}/krakendb/pluspfp_08gb")
    elif [[ -f "${base}/krakendb/hash.k2d" ]]; then
        DB_ARGS+=(--kraken2_db "${base}/krakendb")
    elif [[ -d "${base}/krakendb" ]]; then
        path=$(ls -d "${base}/krakendb"/*/2>/dev/null | head -1 || true)
        [[ -n "$path" ]] && DB_ARGS+=(--kraken2_db "$path")
    fi

    # SILVA
    path=$(ls "${base}/silva_db/SILVA_"*"_SSURef_NR99.fasta" 2>/dev/null | head -1 || true)
    [[ -n "$path" ]] && DB_ARGS+=(--silva_ssu_db "$path")
    path=$(ls "${base}/silva_db/SILVA_"*"_LSURef_NR99.fasta" 2>/dev/null | head -1 || true)
    [[ -n "$path" ]] && DB_ARGS+=(--silva_lsu_db "$path")
}

save_run_command() {
    local outdir="$1" session="$2"
    local save_args=()
    local skip_next=false has_session=false
    for arg in "${ORIGINAL_ARGS[@]}"; do
        if $skip_next; then skip_next=false; continue; fi
        if [[ "$arg" == "--resume" ]]; then
            has_session=true
            continue
        fi
        if $has_session && [[ "$arg" =~ ^[0-9a-f-]{36}$ ]]; then
            has_session=false
            continue
        fi
        has_session=false
        save_args+=("$arg")
    done
    if [[ -n "$session" ]]; then
        save_args+=("--resume" "$session")
    fi
    printf '%s\n' "$(realpath "$0") ${save_args[*]}" >> "${outdir}/pipeline_info/run_command.txt"
}

usage() {
    echo "Usage: $0 --assembly FILE --depths FILE --outdir DIR [options]"
    echo ""
    echo "Required:"
    echo "  --assembly FILE    Path to assembly FASTA"
    echo "  --depths FILE      Path to MetaBAT2-format depth table"
    echo "  --outdir DIR       Output directory"
    echo ""
    echo "Optional inputs:"
    echo "  --bam_dir DIR      Directory with *.sorted.bam + .bai (for SemiBin2, LorBin, COMEBin)"
    echo "  --tnf FILE         Tetranucleotide frequency table (for viz)"
    echo ""
    echo "Mode:"
    echo "  --docker           Run in Docker container"
    echo "  --apptainer        Run in Apptainer/Singularity container"
    echo "  --container        Auto-detect container runtime"
    echo ""
    echo "Shortcuts:"
    echo "  --all              Enable all optional modules (still requires DB paths)"
    echo "  --db_dir DIR       Auto-resolve database paths from download-databases.sh layout"
    echo ""
    echo "Caching & Resume:"
    echo "  --workdir DIR      Nextflow work directory (default: \$SLURM_TMPDIR/mag_analysis_work, else /tmp/...)"
    echo "  --store_dir DIR    Persistent cache directory (storeDir) [default: --outdir]"
    echo "  --publish          Copy outputs instead of storing them; a resubmit then redoes everything"
    echo "  --resume [ID]      Resume a previous run"
    echo ""
    echo "Hybrid Slurm execution:"
    echo "  --slurm                   Run Nextflow on this host and send the heavy steps"
    echo "                            (Bakta, the big binners, GTDB-Tk, Kaiju, eggNOG,"
    echo "                            MetaEuk, geNomad) out as their own Slurm jobs; the"
    echo "                            rest run here. Every task runs in the pipeline image."
    echo "  --slurm_account ACCT      sbatch --account for those jobs (required)"
    echo "  --slurm_queue P, --slurm_time T, --slurm_extra 'OPTS'  passed through"
    echo "  Needs nextflow on PATH and a --workdir on shared storage [default: <outdir>/work]."
    echo ""
    echo "microscape.app deploy (needs --run_viz or --all):"
    echo "  --deploy_slug SLUG        URL slug for the run on microscape.app"
    echo "  --deploy_name NAME        Display name [default: slug]"
    echo "  --deploy_visibility V     private (default) | shared | public"
    echo "  Writes <outdir>/deploy.sh, which each viz stage runs, so the published"
    echo "  run updates as results arrive. Needs \$MICROSCAPE_API_KEY or"
    echo "  ~/.config/microscape/api-key."
    echo ""
    echo "Pipeline flags: --annotator, --run_metabolism, --run_genomad, --run_gtdbtk,"
    echo "  --run_kaiju, --run_kraken2, --run_sendsketch, --run_rrna, --run_eukaryotic,"
    echo "  --run_antismash, --run_viz, --run_semibin, --run_maxbin, --run_lorbin,"
    echo "  --run_comebin, --run_vamb, --run_binette, --run_magscot, etc."
    echo ""
    echo "Run '$0 --help-pipeline' for full Nextflow help."
    exit 0
}

# ============================================================================
# Parse arguments
# ============================================================================

ASSEMBLY_HOST=""
DEPTHS_HOST=""
BAM_DIR_HOST=""
OUTDIR_HOST=""
WORKDIR_HOST="${SLURM_TMPDIR:-/tmp}/mag_analysis_work"
RESUME_SESSION=""
DO_RESUME=false
DB_DIR_HOST=""
STORE_DIR_HOST=""

# Store mode keeps every finished stage in $OUTDIR (storeDir), so a timeout or
# resubmit skips whatever already completed. --publish restores plain copies.
STORE_MODE=true
SLURM_MODE=false
SLURM_ACCOUNT=""
DEPLOY_SLUG=""
DEPLOY_NAME=""
DEPLOY_VISIBILITY="private"

while (( $# )); do
    case "$1" in
        -h|--help)
            usage ;;
        --docker)
            USE_CONTAINER=true; CONTAINER_RUNTIME=docker
            shift ;;
        --apptainer|--singularity)
            USE_CONTAINER=true; CONTAINER_RUNTIME="${1#--}"
            shift ;;
        --container)
            USE_CONTAINER=true; CONTAINER_RUNTIME=auto
            shift ;;
        --image)
            [[ -z "${2:-}" ]] && die "--image requires an image name"
            CONTAINER_IMAGE="$2"
            shift 2 ;;
        --sif)
            [[ -z "${2:-}" ]] && die "--sif requires a path"
            SIF_PATH="$2"
            shift 2 ;;
        --pull)
            PULL_SIF=true
            shift ;;
        --assembly)
            [[ -z "${2:-}" ]] && die "--assembly requires a file path"
            ASSEMBLY_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            NF_ARGS+=("--assembly" "$ASSEMBLY_HOST")
            shift 2 ;;
        --depths)
            [[ -z "${2:-}" ]] && die "--depths requires a file path"
            DEPTHS_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            NF_ARGS+=("--depths" "$DEPTHS_HOST")
            shift 2 ;;
        --bam_dir)
            [[ -z "${2:-}" ]] && die "--bam_dir requires a directory path"
            BAM_DIR_HOST="$(realpath "$2" 2>/dev/null || echo "$2")"
            NF_ARGS+=("--bam_dir" "$BAM_DIR_HOST")
            shift 2 ;;
        --outdir)
            [[ -z "${2:-}" ]] && die "--outdir requires a directory path"
            OUTDIR_HOST="$(realpath -m "$2")"
            shift 2 ;;
        --workdir|-w)
            [[ -z "${2:-}" ]] && die "--workdir requires a directory path"
            WORKDIR_HOST="$(realpath -m "$2")"
            shift 2 ;;
        --db_dir)
            [[ -z "${2:-}" ]] && die "--db_dir requires a directory path"
            [[ -d "$2" ]] || die "--db_dir directory does not exist: $2"
            DB_DIR_HOST="$(realpath "$2")"
            resolve_db_dir "$DB_DIR_HOST"
            shift 2 ;;
        --store_dir)
            [[ -z "${2:-}" ]] && die "--store_dir requires a directory path"
            STORE_DIR_HOST="$(realpath -m "$2")"
            NF_ARGS+=("--store_dir" "$STORE_DIR_HOST")
            shift 2 ;;
        --publish)
            STORE_MODE=false
            shift ;;
        --tnf|--assembly_info)
            # Resolve to canonical path so container bind-path rewriting catches them.
            # Without realpath, /home/rec3141/scratch/X stays unresolved while
            # STORE_DIR_HOST gets canonicalized to /scratch/rec3141/X — the substring
            # substitution then misses these args and Nextflow tries to read the host
            # path inside the container where it isn't bound.
            [[ -z "${2:-}" ]] && die "$1 requires a file path"
            NF_ARGS+=("$1" "$(realpath "$2" 2>/dev/null || echo "$2")")
            shift 2 ;;
        --all)
            NF_ARGS+=(
                --bakta_extra
                --run_semibin true
                --run_maxbin true
                --run_lorbin true
                --run_comebin true
                --run_vamb true
                --run_vamb_tax true
                --run_binette true
                --run_magscot true
                --run_kraken2 true
                --run_sendsketch true
                --run_rrna true
                --run_metabolism true
                --run_antismash true
                --run_eukaryotic true
                --run_metaeuk true
                --run_marferret true
                --run_gtdbtk true
                --run_kaiju true
                --run_genomad true
                --run_checkv true
                --run_defensefinder true
                --run_integronfinder true
                --run_islandpath true
                --run_macsyfinder true
                --run_viz true
            )
            shift ;;
        --slurm)
            # Implies the container: every task runs in the pipeline image.
            SLURM_MODE=true; USE_CONTAINER=true
            [[ "$CONTAINER_RUNTIME" == "docker" ]] && die "--slurm needs apptainer or singularity"
            [[ -z "${CONTAINER_RUNTIME:-}" || "$CONTAINER_RUNTIME" == "auto" ]] && CONTAINER_RUNTIME=auto
            shift ;;
        --slurm_account)
            [[ -z "${2:-}" ]] && die "--slurm_account requires an account"
            SLURM_ACCOUNT="$2"; NF_ARGS+=("--slurm_account" "$2"); shift 2 ;;
        --deploy_slug)
            [[ -z "${2:-}" ]] && die "--deploy_slug requires a slug"
            DEPLOY_SLUG="$2"; shift 2 ;;
        --deploy_name)
            [[ -z "${2:-}" ]] && die "--deploy_name requires a name"
            DEPLOY_NAME="$2"; shift 2 ;;
        --deploy_visibility)
            case "${2:-}" in private|shared|public) ;; *) die "--deploy_visibility must be private|shared|public" ;; esac
            DEPLOY_VISIBILITY="$2"; shift 2 ;;
        --resume)
            DO_RESUME=true
            if [[ -n "${2:-}" && "${2}" != --* ]]; then
                RESUME_SESSION="$2"
                shift
            fi
            shift ;;
        *)
            NF_ARGS+=("$1")
            shift ;;
    esac
done

# Prepend DB_ARGS so explicit user flags override auto-detected paths
NF_ARGS=("${DB_ARGS[@]}" "${NF_ARGS[@]}")

# ============================================================================
# Validate
# ============================================================================

[[ -z "$ASSEMBLY_HOST" ]] && die "--assembly is required."
[[ -z "$DEPTHS_HOST" ]] && die "--depths is required."

if [[ -n "$OUTDIR_HOST" && -n "$STORE_DIR_HOST" ]]; then
    die "--outdir and --store_dir are mutually exclusive."
fi

if [[ -z "$OUTDIR_HOST" && -n "$STORE_DIR_HOST" ]]; then
    OUTDIR_HOST="$STORE_DIR_HOST"
    NF_ARGS+=("--outdir" "$OUTDIR_HOST")
fi
[[ -z "$OUTDIR_HOST" ]] && die "--outdir (or --store_dir) is required."

# Default to store mode. storeDir writes each module to the same place publishDir
# did ($OUTDIR/annotation, $OUTDIR/taxonomy, ...), so the layout is unchanged --
# what it adds is that a resubmit reuses whatever is already there instead of
# redoing it. That is what makes a node-local work directory safe: $SLURM_TMPDIR
# is erased when the job ends, and without this the finished bins, annotations
# and GTDB-Tk placements would have nothing to resume from.
if [[ "$STORE_MODE" == true && -z "$STORE_DIR_HOST" ]]; then
    STORE_DIR_HOST="$OUTDIR_HOST"
    NF_ARGS+=("--store_dir" "$STORE_DIR_HOST")
fi

[[ -f "$ASSEMBLY_HOST" ]] || die "Assembly file not found: $ASSEMBLY_HOST"
[[ -f "$DEPTHS_HOST" ]] || die "Depths file not found: $DEPTHS_HOST"

if [[ ! -d "$OUTDIR_HOST" ]]; then
    echo "[INFO] Creating output directory: $OUTDIR_HOST"
    mkdir -p "$OUTDIR_HOST" || die "Cannot create output directory: $OUTDIR_HOST"
fi

# microscape.app deploy hook: VIZ_PREPROCESS runs <outdir>/deploy.sh after each
# viz snapshot. The hook runs where the pipeline runs, so in a container it
# names the image's copy of viz/deploy.sh and the API key bound in below.
if [[ -n "$DEPLOY_SLUG" ]]; then
    [[ -z "$DEPLOY_NAME" ]] && DEPLOY_NAME="$DEPLOY_SLUG"
    if [[ "$USE_CONTAINER" == true ]]; then
        hook_deployer=/pipeline/viz/deploy.sh
        hook_key=(--api-key-file /data/microscape/api-key)
    else
        hook_deployer="${SCRIPT_DIR}/viz/deploy.sh"
        hook_key=()
    fi
    {
        echo '#!/usr/bin/env bash'
        echo "# Written by run-mag-analysis.sh; run by each viz stage with the viz dir as \$1."
        printf 'exec %q --preprocess-dir "$1/data" --slug %q --name %q --visibility %q' \
            "$hook_deployer" "$DEPLOY_SLUG" "$DEPLOY_NAME" "$DEPLOY_VISIBILITY"
        [[ ${#hook_key[@]} -gt 0 ]] && printf ' %q' "${hook_key[@]}"
        echo
    } > "${OUTDIR_HOST}/deploy.sh"
    chmod +x "${OUTDIR_HOST}/deploy.sh"
    echo "[INFO] Wrote deploy hook: ${OUTDIR_HOST}/deploy.sh (slug=$DEPLOY_SLUG, visibility=$DEPLOY_VISIBILITY)"
    if [[ -z "${MICROSCAPE_API_KEY:-}" && ! -f "${HOME}/.config/microscape/api-key" ]]; then
        echo "[WARNING] No \$MICROSCAPE_API_KEY and no ~/.config/microscape/api-key; deploys will fail until one is set." >&2
    fi
fi

# ============================================================================
# Auto-detect session
# ============================================================================

if [[ "$DO_RESUME" == true && -z "$RESUME_SESSION" ]]; then
    _base="${STORE_DIR_HOST:-$OUTDIR_HOST}"
    for _hist in "${_base}/.nextflow-cache/dotdir/history" \
                 "${_base}/.nextflow/history" \
                 "${OUTDIR_HOST}/.nextflow-cache/dotdir/history" \
                 "${OUTDIR_HOST}/.nextflow/history"; do
        if [[ -f "$_hist" ]]; then
            RESUME_SESSION=$(awk '{print $6}' "$_hist" | tail -1)
            break
        fi
    done
    if [[ -n "$RESUME_SESSION" ]]; then
        echo "[INFO] Auto-detected session: $RESUME_SESSION"
    fi
fi

# ============================================================================
# Container runtime
# ============================================================================

if [[ "$USE_CONTAINER" == true ]]; then
    if [[ "$CONTAINER_RUNTIME" == "auto" ]]; then
        if command -v apptainer &>/dev/null; then CONTAINER_RUNTIME=apptainer
        elif command -v singularity &>/dev/null; then CONTAINER_RUNTIME=singularity
        elif command -v docker &>/dev/null; then CONTAINER_RUNTIME=docker
        else die "No container runtime found"; fi
    fi

    if [[ "$CONTAINER_RUNTIME" == "apptainer" || "$CONTAINER_RUNTIME" == "singularity" ]]; then
        [[ -z "$SIF_PATH" ]] && SIF_PATH="${SCRIPT_DIR}/.danaseq-mag-analysis.sif"
        if [[ ! -f "$SIF_PATH" ]]; then
            if [[ "$PULL_SIF" == true ]]; then
                "$CONTAINER_RUNTIME" pull "$SIF_PATH" "docker://${CONTAINER_IMAGE}"
            else
                die "SIF not found: $SIF_PATH"
            fi
        fi
    fi
fi

# ============================================================================
# Hybrid Slurm mode (--slurm)
# ============================================================================
# Nextflow runs on this host (sbatch is not available inside the image) and
# every task runs in the image; -profile slurm_hybrid sends the heavy steps out
# as their own Slurm jobs. The pipeline code is copied out of the image, so the
# code and the tools always come from the same build.

if [[ "$SLURM_MODE" == true ]]; then
    [[ "$CONTAINER_RUNTIME" == "apptainer" || "$CONTAINER_RUNTIME" == "singularity" ]] \
        || die "--slurm needs apptainer or singularity (found: ${CONTAINER_RUNTIME})"
    command -v nextflow >/dev/null 2>&1 \
        || die "--slurm runs Nextflow on the host; put nextflow on PATH (e.g. module load nextflow)"
    [[ -n "$SLURM_ACCOUNT" ]] || die "--slurm needs --slurm_account"

    # Slurm-run steps land on other nodes, so node-local work is invisible to them.
    case "$WORKDIR_HOST" in
        "${SLURM_TMPDIR:-/nonexistent}"*|/tmp/*) WORKDIR_HOST="${OUTDIR_HOST}/work" ;;
    esac
    mkdir -p "$WORKDIR_HOST" || die "Cannot create work directory: $WORKDIR_HOST"

    img_sha=$("$CONTAINER_RUNTIME" exec "$SIF_PATH" printenv DANASEQ_GIT_SHA 2>/dev/null | cut -c1-12)
    PIPE_DIR="${OUTDIR_HOST}/.pipeline-${img_sha:-unknown}"
    if [[ ! -f "$PIPE_DIR/main.nf" ]]; then
        mkdir -p "$PIPE_DIR"
        "$CONTAINER_RUNTIME" exec -B "$PIPE_DIR" "$SIF_PATH" cp -r /pipeline/. "$PIPE_DIR"/ \
            || die "Could not copy the pipeline out of $SIF_PATH"
    fi

    # Tasks see everything at its host path, plus the squashed DB images at
    # /data/dbimg/<sub> with their --*_db arguments rewritten, as in container mode.
    RUNOPTS=()
    for d in "$(dirname "$ASSEMBLY_HOST")" "$(dirname "$DEPTHS_HOST")" "$OUTDIR_HOST" \
             "$WORKDIR_HOST" "$PIPE_DIR" ${BAM_DIR_HOST:+"$BAM_DIR_HOST"} \
             ${DB_DIR_HOST:+"$DB_DIR_HOST"} ${STORE_DIR_HOST:+"$STORE_DIR_HOST"}; do
        RUNOPTS+=(-B "$d")
    done
    for e in "${SQSH_IMAGES[@]}"; do
        img="${e%%|*}"; sub="${e#*|}"
        RUNOPTS+=(-B "${img}:/data/dbimg/${sub}:image-src=/")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            case "${NF_ARGS[$i]}" in
                "${DB_DIR_HOST}/${sub}"|"${DB_DIR_HOST}/${sub}"/*)
                    NF_ARGS[$i]="/data/dbimg/${sub}${NF_ARGS[$i]#${DB_DIR_HOST}/${sub}}" ;;
            esac
        done
    done
    if [[ -n "$DEPLOY_SLUG" && -f "${HOME}/.config/microscape/api-key" ]]; then
        RUNOPTS+=(-B "${HOME}/.config/microscape:/data/microscape:ro")
    fi
    # Slurm jobs inherit the submitting environment, and the container passes it on.
    if [[ -n "${DB_DIR_HOST:-}" && -f "${DB_DIR_HOST}/checkm_data/taxon_marker_sets.tsv" ]]; then
        export CHECKM_DATA_PATH="${DB_DIR_HOST}/checkm_data"
    fi

    mkdir -p "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info"
    RUNTIME_CFG="${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info/slurm_hybrid.runtime.config"
    {
        echo "// Written by run-mag-analysis.sh --slurm"
        printf "process.container = '%s'\n" "$SIF_PATH"
        printf "singularity.runOptions = '%s'\n" "${RUNOPTS[*]}"
    } > "$RUNTIME_CFG"

    NF_CMD=(nextflow -log "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info/nextflow.log"
            run "$PIPE_DIR/main.nf" -profile slurm_hybrid
            -c "$PIPE_DIR/docker.config" -c "$RUNTIME_CFG"
            -w "$WORKDIR_HOST" --outdir "$OUTDIR_HOST" "${NF_ARGS[@]}")
    [[ "$DO_RESUME" == true ]] && NF_CMD+=(-resume ${RESUME_SESSION})

    echo "[INFO] Mode:     Hybrid Slurm (image ${img_sha:-unknown}, account ${SLURM_ACCOUNT})"
    echo "[INFO] Work:     $WORKDIR_HOST"
    echo "[INFO] Running: ${NF_CMD[*]}"
    "${NF_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?
    NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info/nextflow.log" 2>/dev/null | tail -1) || true
    save_run_command "${STORE_DIR_HOST:-$OUTDIR_HOST}" "$NF_SESSION" || echo "[WARN] Could not record run_command.txt"
    exit $NF_EXIT
fi

# ============================================================================
# Container mode
# ============================================================================

if [[ "$USE_CONTAINER" == true ]]; then
    BINDS=()
    BINDS+=("$(dirname "$ASSEMBLY_HOST"):/data/assembly:ro")
    BINDS+=("$(dirname "$DEPTHS_HOST"):/data/depths:ro")
    BINDS+=("${OUTDIR_HOST}:/data/output")
    if [[ -n "$DEPLOY_SLUG" && -f "${HOME}/.config/microscape/api-key" ]]; then
        BINDS+=("${HOME}/.config/microscape:/data/microscape:ro")
    fi

    # Rewrite assembly/depths paths
    for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
        NF_ARGS[$i]="${NF_ARGS[$i]//${ASSEMBLY_HOST}/\/data\/assembly\/$(basename "$ASSEMBLY_HOST")}"
        NF_ARGS[$i]="${NF_ARGS[$i]//${DEPTHS_HOST}/\/data\/depths\/$(basename "$DEPTHS_HOST")}"
    done
    NF_ARGS=("--outdir" "/data/output" "${NF_ARGS[@]}")

    if [[ -n "${BAM_DIR_HOST:-}" ]]; then
        BINDS+=("${BAM_DIR_HOST}:/data/bams:ro")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${BAM_DIR_HOST}/\/data\/bams}"
        done
    fi

    if [[ -n "${DB_DIR_HOST:-}" ]]; then
        BINDS+=("${DB_DIR_HOST}:/data/db:ro")
        # Squashed images: mount each at /data/dbimg/<sub> and point the
        # matching --*_db paths there (before the generic /data/db rewrite).
        if [[ "$CONTAINER_RUNTIME" == "apptainer" || "$CONTAINER_RUNTIME" == "singularity" ]]; then
            for e in "${SQSH_IMAGES[@]}"; do
                img="${e%%|*}"; sub="${e#*|}"
                BINDS+=("${img}:/data/dbimg/${sub}:image-src=/")
                for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
                    case "${NF_ARGS[$i]}" in
                        "${DB_DIR_HOST}/${sub}"|"${DB_DIR_HOST}/${sub}"/*)
                            NF_ARGS[$i]="/data/dbimg/${sub}${NF_ARGS[$i]#${DB_DIR_HOST}/${sub}}" ;;
                    esac
                done
                echo "[INFO] Using squashed database image $(basename "$img") for ${sub}/"
            done
        else
            for e in "${SQSH_IMAGES[@]}"; do
                sub="${e#*|}"
                [[ -e "${DB_DIR_HOST}/${sub}" ]] || echo "[WARNING] ${sub}: only a .sqsh image is available; ${CONTAINER_RUNTIME} cannot mount it, database skipped"
            done
            for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
                for e in "${SQSH_IMAGES[@]}"; do
                    sub="${e#*|}"
                    case "${NF_ARGS[$i]}" in
                        "${DB_DIR_HOST}/${sub}"|"${DB_DIR_HOST}/${sub}"/*) [[ -e "${NF_ARGS[$i]}" ]] || { unset "NF_ARGS[$((i-1))]" "NF_ARGS[$i]"; } ;;
                    esac
                done
            done
            NF_ARGS=("${NF_ARGS[@]}")
        fi
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${DB_DIR_HOST}/\/data\/db}"
        done
    fi

    if [[ -n "${STORE_DIR_HOST:-}" ]]; then
        BINDS+=("${STORE_DIR_HOST}:/data/store")
        for (( i=0; i<${#NF_ARGS[@]}; i++ )); do
            NF_ARGS[$i]="${NF_ARGS[$i]//${STORE_DIR_HOST}/\/data\/store}"
        done
    fi

    mkdir -p "${WORKDIR_HOST}" 2>/dev/null || true
    BINDS+=("${WORKDIR_HOST}:/data/work")

    NF_CACHE="${OUTDIR_HOST}/.nextflow-cache"
    mkdir -p "${NF_CACHE}/dotdir" 2>/dev/null || true
    BINDS+=("${NF_CACHE}/dotdir:/home/dana/.nextflow")

    # Node-local scratch for the tools' own temporary files: binning, GTDB-Tk
    # and antiSMASH each write sizeable scratch, which otherwise lands in the
    # container overlay. SLURM_TMPDIR where a scheduler provides it, else
    # TMPDIR. Bound at the same path inside and out, and exported, because
    # Docker does not inherit the host environment.
    SCRATCH_HOST="${SLURM_TMPDIR:-${TMPDIR:-}}"
    if [[ -n "$SCRATCH_HOST" && -d "$SCRATCH_HOST" && -w "$SCRATCH_HOST" ]]; then
        BINDS+=("${SCRATCH_HOST}:${SCRATCH_HOST}")
    else
        SCRATCH_HOST=""
    fi

    # CheckM v1 marker sets for COMEBin, from the bound database directory.
    CONTAINER_ENV=()
    if [[ -n "${DB_DIR_HOST:-}" && -f "${DB_DIR_HOST}/checkm_data/taxon_marker_sets.tsv" ]]; then
        CONTAINER_ENV+=("CHECKM_DATA_PATH=/data/db/checkm_data")
    fi

    CONTAINER_CMD=()
    case "$CONTAINER_RUNTIME" in
        docker)
            CONTAINER_CMD+=(docker run --user "$(id -u):$(id -g)")
            CONTAINER_CMD+=("-e" "NXF_HOME=/home/dana/.nextflow")
            for e in "${CONTAINER_ENV[@]}"; do CONTAINER_CMD+=("-e" "$e"); done
            [[ -n "$SCRATCH_HOST" ]] && CONTAINER_CMD+=("-e" "TMPDIR=${SCRATCH_HOST}")
            for bind in "${BINDS[@]}"; do CONTAINER_CMD+=("-v" "$bind"); done
            CONTAINER_CMD+=("$CONTAINER_IMAGE" -log /data/output/pipeline_info/nextflow.log run /pipeline/main.nf)
            ;;
        apptainer|singularity)
            container_ca="/etc/ssl/certs/ca-certificates.crt"
            CONTAINER_CMD+=("$CONTAINER_RUNTIME" run)
            CONTAINER_CMD+=("--env" "NXF_HOME=/home/dana/.nextflow")
            for e in "${CONTAINER_ENV[@]}"; do CONTAINER_CMD+=("--env" "$e"); done
            [[ -n "$SCRATCH_HOST" ]] && CONTAINER_CMD+=("--env" "TMPDIR=${SCRATCH_HOST}")
            CONTAINER_CMD+=("--env" "REQUESTS_CA_BUNDLE=${container_ca}")
            CONTAINER_CMD+=("--env" "SSL_CERT_FILE=${container_ca}")
            CONTAINER_CMD+=("--env" "CURL_CA_BUNDLE=${container_ca}")
            for bind in "${BINDS[@]}"; do CONTAINER_CMD+=("--bind" "$bind"); done
            CONTAINER_CMD+=("$SIF_PATH" -log /data/output/pipeline_info/nextflow.log run /pipeline/main.nf)
            ;;
    esac
    CONTAINER_CMD+=(-w /data/work "${NF_ARGS[@]}")
    if [[ "$DO_RESUME" == true ]]; then
        CONTAINER_CMD+=(-resume ${RESUME_SESSION})
    fi

    echo "[INFO] Mode:   Container (${CONTAINER_RUNTIME})"
    echo "[INFO] Assembly: $ASSEMBLY_HOST"
    echo "[INFO] Depths:   $DEPTHS_HOST"
    echo "[INFO] Output:   $OUTDIR_HOST"
    echo ""

    # Also holds nextflow.log: without -log, Nextflow writes .nextflow.log
    # relative to the launch directory, so concurrent runs sharing a submit
    # directory rotate each other out and the one that failed is gone.
    mkdir -p "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info" 2>/dev/null || true
    "${CONTAINER_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

    NF_SESSION=$(awk '{print $6}' "${NF_CACHE}/dotdir/history" 2>/dev/null | tail -1)
    [[ -z "$NF_SESSION" ]] && NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info/nextflow.log" 2>/dev/null | tail -1)
    save_run_command "${STORE_DIR_HOST:-$OUTDIR_HOST}" "$NF_SESSION"
    exit $NF_EXIT
fi

# ============================================================================
# Local mode (default)
# ============================================================================

CONDA_BASE_BIN="$(dirname "$(which mamba 2>/dev/null || which conda)")"
export PATH="${CONDA_BASE_BIN}:${PATH}"

WORKDIR_FLAG=()
if [[ -n "$WORKDIR_HOST" ]]; then
    mkdir -p "$WORKDIR_HOST" || die "Cannot create work directory: $WORKDIR_HOST"
    WORKDIR_FLAG=(-w "$WORKDIR_HOST")
fi

RESUME_FLAG=()
if [[ "$DO_RESUME" == true ]]; then
    RESUME_FLAG=(-resume ${RESUME_SESSION})
fi

# Use the flye env (has Nextflow + Java) or any available env with Nextflow
NF_ENV="${SCRIPT_DIR}/conda-envs/dana-mag-assembly"
if [[ ! -d "$NF_ENV" ]]; then
    # Fallback: look for any env with nextflow
    for env in "${SCRIPT_DIR}"/conda-envs/dana-mag-*; do
        if [[ -x "${env}/bin/nextflow" ]]; then
            NF_ENV="$env"
            break
        fi
    done
fi

LOCAL_CMD=(
    mamba run -p "$NF_ENV"
    nextflow -log "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info/nextflow.log" run "${SCRIPT_DIR}/main.nf"
    --outdir "$OUTDIR_HOST"
    "${WORKDIR_FLAG[@]}"
    "${NF_ARGS[@]}"
    "${RESUME_FLAG[@]}"
)

echo "[INFO] Mode:     Local (conda)"
echo "[INFO] Assembly: $ASSEMBLY_HOST"
echo "[INFO] Depths:   $DEPTHS_HOST"
echo "[INFO] Output:   $OUTDIR_HOST"
echo ""

mkdir -p "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info" 2>/dev/null || true

"${LOCAL_CMD[@]}" && NF_EXIT=0 || NF_EXIT=$?

NF_SESSION=$(awk '{print $6}' .nextflow/history 2>/dev/null | tail -1)
[[ -z "$NF_SESSION" ]] && NF_SESSION=$(grep -oP 'Session UUID: \K[0-9a-f-]{36}' "${STORE_DIR_HOST:-$OUTDIR_HOST}/pipeline_info/nextflow.log" 2>/dev/null | tail -1)
save_run_command "${STORE_DIR_HOST:-$OUTDIR_HOST}" "$NF_SESSION"

exit $NF_EXIT
