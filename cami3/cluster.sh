# Cluster settings for the cami3 scripts. Every script sources this file, so the
# same scripts run on fir and on grex. It detects the cluster (grex has project
# 6043541; CAMI3_CLUSTER overrides) and sets:
#   CLUSTER, BASE, GENICE, DB, HUMAN_FA, HUMAN_BBMAP, ASM_MEM, BIG_PARTITION
#   container_setup   exports the apptainer/singularity environment
#   cexec SIF ARGS    runs a command inside an image, with the data paths bound
# fir:  data on /scratch (1M-file quota) and work dirs on node-local TMP where the
#       launcher allows. fir /project is near its file quota, so nothing goes there.
# grex: data and work dirs must live on /project (6043541), whose inode quota is
#       tight (4.0M soft), so check `lfs quota -p 6043541 /project` before big runs.
#       Big jobs go to genoa (192 cores, 750 GB), so the assembly memory is capped below that.
if [ -z "${CAMI3_CLUSTER:-}" ]; then
    if [ -d /project/6043541 ]; then CAMI3_CLUSTER=grex; else CAMI3_CLUSTER=fir; fi
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
        }
        cexec() { local sif=$1; shift; apptainer exec --bind /scratch --bind "${SLURM_TMPDIR:-/tmp}" "$sif" "$@"; } ;;
    grex)
        BASE=/project/6043541/rec3141/scratch/cami3
        # Slurm lives outside the default PATH on grex.
        case ":$PATH:" in *:/opt/slurm/bin:*) ;; *) export PATH=/opt/slurm/bin:$PATH ;; esac
        ASM_MEM='650 GB'
        BIG_PARTITION=genoa
        container_setup() {
            export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
            export APPTAINER_CACHEDIR=/home/rec3141/scratch/apptainer_cache
            export APPTAINER_TMPDIR=/home/rec3141/scratch/apptainer_cache/tmp
            mkdir -p "$APPTAINER_TMPDIR"
            # ~/scratch is a symlink into /project; bind the resolved path, or links through it dangle.
            local real; real=$(readlink -f "$HOME/scratch")
            export SINGULARITY_BIND="$real" APPTAINER_BINDPATH="$real"
        }
        cexec() { local sif=$1; shift; singularity exec --bind "$(readlink -f "$HOME/scratch")" --bind "${SLURM_TMPDIR:-/tmp}" "$sif" "$@"; } ;;
    *) echo "[ERROR] unknown cluster $CLUSTER" >&2; exit 1 ;;
esac
