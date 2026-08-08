// Merge sequence tables and quality control.
//
// Three-step post-DADA2 processing:
//   1. MERGE_SEQTABS   — combine per-plate sequence tables (long-format output)
//   2. REMOVE_CHIMERAS — sparse consensus chimera removal (long-format I/O)
//   3. FILTER_SEQTAB   — length, prevalence, abundance filtering (long → long + wide)
//
// Supports --lang R or --lang python.

def seqExt() { return params.lang == 'python' ? 'pkl' : 'rds' }

process MERGE_SEQTABS {
    tag "merge-all"
    label 'process_medium'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    // Publish the stats: without them the read/sample attrition between
    // denoising and the final table is invisible in the results, which made
    // a 218k -> 51k read drop impossible to attribute after the fact.
    publishDir "${params.outdir}/seqtab_final", mode: 'copy', pattern: '*.tsv'

    input:
    path(seqtab_files)

    output:
    path("seqtab_merged.${seqExt()}"), emit: seqtab
    path("merge_stats.tsv"), emit: stats
    path("merge_sample_reads.tsv"), emit: sample_reads, optional: true

    script:
    if (params.lang == 'python')
    """
    merge_seqtabs.py
    """
    else
    """
    merge_seqtabs.R
    """
}

// Sparse consensus chimera removal.
process REMOVE_CHIMERAS {
    tag "chimera-removal"
    label 'process_high'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    // Publish the stats: without them the read/sample attrition between
    // denoising and the final table is invisible in the results, which made
    // a 218k -> 51k read drop impossible to attribute after the fact.
    publishDir "${params.outdir}/seqtab_final", mode: 'copy', pattern: '*.tsv'

    input:
    path(seqtab)

    output:
    path("seqtab_nochim.${seqExt()}"), emit: seqtab
    path("chimera_stats.tsv"), emit: stats
    path("chimera_sample_reads.tsv"), emit: sample_reads, optional: true

    script:
    if (params.lang == 'python')
    """
    remove_chimeras.py "${seqtab}" ${task.cpus}
    """
    else
    """
    remove_chimeras.R "${seqtab}" ${task.cpus}
    """
}

// Length, prevalence, and abundance filtering.
process FILTER_SEQTAB {
    tag "filter-qc"
    label 'process_medium'
    conda params.lang == 'python' ? "${projectDir}/envs/python.yml" : "${projectDir}/envs/r.yml"
    publishDir "${params.outdir}/seqtab_final", mode: 'copy'

    input:
    path(seqtab)

    output:
    path("seqtab_final.${seqExt()}"), emit: seqtab
    path("seqtab_final_wide.${seqExt()}"), emit: seqtab_wide
    path("seqtab_orphans.${seqExt()}"), emit: orphans
    path("seqtab_small.${seqExt()}"), emit: small_samples
    path("filter_stats.tsv"), emit: stats
    path("final_sample_reads.tsv"), emit: sample_reads, optional: true
    path("sequence_summaries.pdf"), emit: plots

    script:
    if (params.lang == 'python')
    """
    filter_seqtab.py \
        "${seqtab}" \
        ${params.min_seq_length} ${params.min_samples} \
        ${params.min_seqs} ${params.min_reads}
    """
    else
    """
    filter_seqtab.R \
        "${seqtab}" \
        ${params.min_seq_length} ${params.min_samples} \
        ${params.min_seqs} ${params.min_reads}
    """
}
