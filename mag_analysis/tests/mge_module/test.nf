// Module test for GENOMAD_CLASSIFY -> CHECKV_QUALITY -> MGE_BIN_ANNOTATE with fake tools.
// Run through run_tests.sh, which puts fakebin/ and ../../bin on PATH.
nextflow.enable.dsl = 2
params.outdir = "results"
params.store_dir = null
params.genomad_db = "fake_db"
params.checkv_db = "fake_db"
params.run_checkv = true
params.mge_binner_label = 'test'
params.mge_alt_export = false
params.mge_rescue_unbinned = false
params.mge_export_min_plasmid_score = null
params.mge_export_min_virus_score = null
params.assembly = null
params.contig2bin = null

include { GENOMAD_CLASSIFY; CHECKV_QUALITY; MGE_BIN_ANNOTATE } from '../../modules/mge'

workflow {
    asm = Channel.value(file(params.assembly))
    GENOMAD_CLASSIFY(asm)
    def no_file = file("${projectDir}/../../assets/NO_FILE")
    cv_status = no_file; cv_quality = no_file
    if (params.run_checkv) {
        CHECKV_QUALITY(GENOMAD_CLASSIFY.out.virus_fasta)
        cv_status = CHECKV_QUALITY.out.status; cv_quality = CHECKV_QUALITY.out.quality
    }
    MGE_BIN_ANNOTATE(asm, Channel.value(file(params.contig2bin)), GENOMAD_CLASSIFY.out.status,
        GENOMAD_CLASSIFY.out.summary, GENOMAD_CLASSIFY.out.virus_summary, GENOMAD_CLASSIFY.out.plasmid_summary,
        GENOMAD_CLASSIFY.out.provirus_coords, cv_status, cv_quality)
}
