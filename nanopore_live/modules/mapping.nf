// Reference mapping with minimap2 (ONT preset).
//
// One [sample, reference] task per invocation — main.nf produces the
// Cartesian product. Output is a filtered SAM-ish stream (headers
// stripped, mapq>=1 and aligned_len>=10 — see filter_minimap2.awk).
//
// Conventions:
//   - Reference index basename is the canonical reference name (matches the
//     `mapping.reference` column in dana.duckdb and the ./data/ais_<name>.json
//     filename consumed by the SPA's AIS view).
//   - Output filename is `<refname>.txt`, dropped under map/ inside the
//     barcode dir so multiple references can co-exist.

process MAP_REFERENCE {
    tag "${meta.id}_${refname}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-tools"
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/map", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/${meta.flowcell}/${meta.barcode}/map" : null

    input:
    tuple val(meta), path(fastq), val(refname), path(ref_idx)

    output:
    tuple val(meta), val(refname), path("${refname}.txt"), emit: hits

    script:
    """
    # pipefail is essential here — minimap2 failing (e.g. wrong env, no
    # binary on PATH) must propagate past awk so Nextflow marks the task
    # failed. Without it, the awk leg returns 0 even when minimap2 wrote
    # nothing, and the empty .txt gets cached as "successful".
    set -o pipefail

    # Stream the QC'd fastq through minimap2 ONT preset, then the awk filter
    # (drops headers, keeps mapq>=1 + aligned_len>=10). --secondary=no keeps
    # one primary alignment per read, matching the standalone protocol.
    minimap2 \\
        -ax map-ont \\
        -t ${task.cpus} \\
        --secondary=no \\
        "${ref_idx}" \\
        "${fastq}" \\
      | awk -f ${projectDir}/bin/filter_minimap2.awk \\
      > "${refname}.txt"
    """
}
