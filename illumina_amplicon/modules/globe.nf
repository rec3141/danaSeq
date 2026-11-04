// EarthDNA: place each ASV on the globe of full-length SSU references (see bin/place_globe.py).
// Runs only with --earthdna_db; ITS, LSU and anything not SSU are counted as unplaced.

process PLACE_GLOBE {
    tag "earthdna"
    label 'process_medium'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/globe", mode: 'copy'

    input:
    path(seqtab)
    path(earthdna_db)

    output:
    path("globe_placements.tsv"), emit: placements
    path("globe_summary.json"),   emit: summary

    script:
    """
    place_globe.py "${seqtab}" "${earthdna_db}" ${task.cpus} "${params.fastsne_bin}"
    """
}
