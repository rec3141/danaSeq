// Read-level t-SNE, decoupled from DB_SYNC.
//
// DB_SYNC ticks every db_sync_minutes (default 10) but t-SNE on ~1M tetra_data
// rows takes longer than a tick. Running it inline inside DB_SYNC's while-loop
// stampedes — a slow t-SNE blocks the next sync cycle and the read-explorer
// JSON falls minutes-to-hours behind reality.
//
// Design: DB_SYNC drops a single trigger file (.tsne_pending) when (a) no
// trigger is already pending and (b) rows grew since the last trigger.
// This process consumes the trigger via Channel.watchPath. maxForks=1
// provides back-pressure — a new trigger arriving during a run queues
// until the current run exits. The script removes the pending file as its
// first act so DB_SYNC may queue a fresh trigger during the run; queue
// depth is bounded at 1 running + 1 pending. The "interval" is emergent:
// (t-SNE runtime) + (next DB_SYNC tick). No time knob.
//
// Failure mode: if compute_read_tsne.py crashes, the pending file has
// already been deleted, so we don't busy-loop on a poison input — the
// next DB_SYNC tick that sees row growth will queue again.

process READ_TSNE {
    tag "read-tsne"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-tools"
    maxForks 1
    executor 'local'

    input:
    path trigger
    val  outdir

    output:
    path 'tsne.done'

    script:
    """
    # Clear the pending slot first so DB_SYNC can post a fresh trigger
    # while we work. Queue depth is bounded at 1 running + 1 pending.
    rm -f ${outdir}/.tsne_pending

    rows=\$(cat ${trigger} 2>/dev/null || echo 0)
    echo "[INFO] READ_TSNE: starting on \${rows} rows"

    mkdir -p ${outdir}/viz
    python3 ${projectDir}/viz/preprocess/compute_read_tsne.py \\
        --input ${outdir} \\
        --output ${outdir}/viz 2>&1 | sed 's/^/  [TSNE] /'

    date -Iseconds > tsne.done
    echo "[INFO] READ_TSNE: complete"
    """
}
