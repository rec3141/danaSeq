// NCLB (No Contig Left Behind) bin refinement
//
// Four-step refinement process:
//   1. GATHER  — build identity cards for all contigs (Python only)
//   2. CONVERSE — LLM tool-use conversations for placement (requires LLM server)
//   3. ELDERS   — investigate SCG redundancy: ecotype vs contamination
//   4. INTEGRATE — apply proposals, extract FASTAs, generate chronicle
//
// Requires: --nclb_dir pointing to the NCLB repository
// Optional: --nclb_base_url for LLM server (default: local LM Studio)


process NCLB_GATHER {
    tag "nclb_gather"
    label 'process_medium'
    conda "${params.nclb_dir}/envs/nclb.yml"
    publishDir "${params.outdir}/binning/nclb", mode: 'link', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/nclb" : null

    input:
    val(ready)    // upstream completion signal (collected outputs)

    output:
    path("gathering.json"), emit: gathering

    script:
    def results_dir = file(params.outdir).toAbsolutePath()
    def nclb_dir    = file(params.nclb_dir).toAbsolutePath()
    """
    python "${nclb_dir}/bin/nclb_gather.py" \\
        --results "${results_dir}" \\
        --output gathering.json
    """
}


process NCLB_CONVERSE {
    tag "nclb_converse"
    label 'process_medium'
    conda "${params.nclb_dir}/envs/nclb.yml"
    publishDir "${params.outdir}/binning/nclb", mode: 'link', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/nclb" : null

    input:
    path(gathering)   // from NCLB_GATHER

    output:
    path("proposals.json"), emit: proposals

    script:
    def results_dir = file(params.outdir).toAbsolutePath()
    def nclb_dir    = file(params.nclb_dir).toAbsolutePath()
    def base_url    = params.nclb_base_url ?: "http://localhost:1234/v1"
    def model_arg   = params.nclb_model ? "--model ${params.nclb_model}" : ""
    """
    python "${nclb_dir}/bin/nclb_converse.py" \\
        --results "${results_dir}" \\
        --gathering "${gathering}" \\
        --base-url "${base_url}" \\
        ${model_arg} \\
        --output proposals.json
    """
}


process NCLB_ELDERS {
    tag "nclb_elders"
    label 'process_medium'
    conda "${params.nclb_dir}/envs/nclb.yml"
    publishDir "${params.outdir}/binning/nclb", mode: 'link', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/nclb" : null

    input:
    path(gathering)   // from NCLB_GATHER (DAG ordering)

    output:
    path("elder_reports.json"), emit: reports

    script:
    def results_dir = file(params.outdir).toAbsolutePath()
    def nclb_dir    = file(params.nclb_dir).toAbsolutePath()
    def ani_flag    = params.nclb_with_ani ? "--with-ani" : ""
    """
    python "${nclb_dir}/bin/nclb_elders.py" \\
        --results "${results_dir}" \\
        ${ani_flag} \\
        --output elder_reports.json
    """
}


process NCLB_INTEGRATE {
    tag "nclb_integrate"
    label 'process_medium'
    conda "${params.nclb_dir}/envs/nclb.yml"
    publishDir "${params.outdir}/binning/nclb", mode: 'link', enabled: !params.store_dir,
        saveAs: { fn -> fn }
    storeDir params.store_dir ? "${params.store_dir}/binning/nclb" : null

    input:
    path(proposals)       // from NCLB_CONVERSE

    output:
    path("communities/*.fa"),     emit: fastas, optional: true
    path("chronicle.json"),       emit: chronicle_json
    path("chronicle.md"),         emit: chronicle_md
    path("contig2community.tsv"), emit: membership
    path("quality_report.tsv"),   emit: quality
    path("valence_report.tsv"),   emit: valence

    script:
    def results_dir = file(params.outdir).toAbsolutePath()
    def nclb_dir    = file(params.nclb_dir).toAbsolutePath()
    """
    python "${nclb_dir}/bin/nclb_integrate.py" \\
        --proposals "${proposals}" \\
        --results "${results_dir}" \\
        --output .
    """
}
