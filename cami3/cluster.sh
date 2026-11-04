# Cluster settings for the cami3 scripts. Every script sources this file, so the
# same scripts run on fir and on grex. It detects the cluster (CAMI3_CLUSTER
# overrides) and sets:
#   CLUSTER, BASE, GENICE, DB, HUMAN_FA, HUMAN_BBMAP, ASM_MEM, BIG_PARTITION
#   container_setup   exports the apptainer/singularity environment
#   cexec SIF ARGS    runs a command inside an image, with the data paths bound
# fir:  data on /scratch (1M-file quota) and work dirs on node-local TMP where the
#       launcher allows. fir /project is near its file quota, so nothing goes there.
# grex: data and work dirs must live on /project (6043541), whose inode quota is
#       tight (4.0M soft), so check `lfs quota -p 6043541 /project` before big runs.
#       Big jobs go to genoa (192 cores, 750 GB), so the assembly memory is capped below that.
# fir sets CC_CLUSTER (as do other Alliance clusters) and its nodes are *.fir.* or fc*;
# grex has neither, but has singularity-ce under /opt. /project/6043541 exists on both.
if [ -z "${CAMI3_CLUSTER:-}" ]; then
    case "${CC_CLUSTER:-}:$(hostname -f 2>/dev/null)" in
        fir:*|*.fir.*|*:fc[0-9]*) CAMI3_CLUSTER=fir ;;
        *) if [ -d /opt/singularity-ce-4.4.1 ]; then CAMI3_CLUSTER=grex
           else echo "[ERROR] cannot tell fir from grex; set CAMI3_CLUSTER" >&2; exit 1; fi ;;
    esac
fi
CLUSTER=$CAMI3_CLUSTER
GENICE=/home/rec3141/GENICE
DB=/home/rec3141/scratch/databases
HUMAN_FA=$DB/human_ref/hg19_main_mask_ribo_animal_allplant_allfungus.fa.gz
HUMAN_BBMAP=$DB/human_ref
case "$CLUSTER" in
    fir)
        BASE=/scratch/rec3141/cami3
        ASM_MEM='700 GB'
        BIG_PARTITION=
        container_setup() {
            module load apptainer
            export APPTAINER_CACHEDIR=/home/rec3141/scratch/apptainer_cache
            export APPTAINER_TMPDIR=${SLURM_TMPDIR:-/tmp}
            # Tools inside the image use TMPDIR, so the node-local dir is also bound into
            # every container (APPTAINER_BINDPATH), including ones the launchers start.
            if [ -n "${SLURM_TMPDIR:-}" ]; then
                export TMPDIR=$SLURM_TMPDIR APPTAINERENV_TMPDIR=$SLURM_TMPDIR
                export APPTAINER_BINDPATH=${APPTAINER_BINDPATH:+$APPTAINER_BINDPATH,}$SLURM_TMPDIR
            fi
        }
        cexec() { local sif=$1; shift; apptainer exec --bind /scratch --bind "${SLURM_TMPDIR:-/tmp}" "$sif" "$@"; } ;;
    grex)
        BASE=/project/6043541/rec3141/scratch/cami3
        # Slurm lives outside the default PATH on grex.
        case ":$PATH:" in *:/opt/slurm/bin:*) ;; *) export PATH=/opt/slurm/bin:$PATH ;; esac
        ASM_MEM='650 GB'
        BIG_PARTITION=genoa
        # grex jobs get no SLURM_TMPDIR, and node /tmp may be too small for
        # decompressed reads, so each job gets a temp dir on /project that is
        # removed when the job exits.
        if [ -z "${SLURM_TMPDIR:-}" ] && [ -n "${SLURM_JOB_ID:-}" ]; then
            export SLURM_TMPDIR=$BASE/tmp/job_$SLURM_JOB_ID
            mkdir -p "$SLURM_TMPDIR"
            trap 'rm -rf "$SLURM_TMPDIR"' EXIT
        fi
        container_setup() {
            export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
            export APPTAINER_CACHEDIR=/home/rec3141/scratch/apptainer_cache
            export APPTAINER_TMPDIR=/home/rec3141/scratch/apptainer_cache/tmp
            # Singularity CE reads its own variable names, so set both spellings.
            export SINGULARITY_CACHEDIR=$APPTAINER_CACHEDIR SINGULARITY_TMPDIR=$APPTAINER_TMPDIR
            mkdir -p "$APPTAINER_TMPDIR"
            # Tools inside the image use TMPDIR; the per-job dir is on /project,
            # which the ~/scratch bind below covers.
            if [ -n "${SLURM_TMPDIR:-}" ]; then
                export TMPDIR=$SLURM_TMPDIR SINGULARITYENV_TMPDIR=$SLURM_TMPDIR APPTAINERENV_TMPDIR=$SLURM_TMPDIR
            fi
            # ~/scratch is a symlink into /project; bind the resolved path, or links through it dangle.
            local real; real=$(readlink -f "$HOME/scratch")
            export SINGULARITY_BIND="$real" APPTAINER_BINDPATH="$real"
        }
        cexec() { local sif=$1; shift; singularity exec --bind "$(readlink -f "$HOME/scratch")" --bind "${SLURM_TMPDIR:-/tmp}" "$sif" "$@"; } ;;
    *) echo "[ERROR] unknown cluster $CLUSTER" >&2; exit 1 ;;
esac
