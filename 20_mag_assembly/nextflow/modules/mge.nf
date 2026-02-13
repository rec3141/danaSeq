// Mobile genetic element detection: geNomad (virus + plasmid), CheckV (viral QA)

process GENOMAD_CLASSIFY {
    tag "genomad"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-genomad"
    publishDir "${params.outdir}/mge/genomad", mode: 'copy'

    input:
    path(assembly)

    output:
    path("virus_summary.tsv"),   emit: virus_summary
    path("plasmid_summary.tsv"), emit: plasmid_summary
    path("virus.fna"),           emit: virus_fasta
    path("plasmid.fna"),         emit: plasmid_fasta
    path("genomad_summary.tsv"), emit: summary

    script:
    def db_path = params.genomad_db
    """
    # geNomad end-to-end: marker gene annotation → neural network classification
    # Detects viruses, plasmids, and proviruses in a single pass
    set +e
    genomad end-to-end \\
        --cleanup \\
        --splits ${task.cpus} \\
        "${assembly}" \\
        genomad_out \\
        "${db_path}"
    genomad_exit=\$?
    set -e

    if [ \$genomad_exit -ne 0 ]; then
        echo "[WARNING] geNomad exited with code \$genomad_exit" >&2
        touch virus_summary.tsv plasmid_summary.tsv virus.fna plasmid.fna genomad_summary.tsv
        exit 0
    fi

    # geNomad names output files based on the input filename
    # Find the actual output directory
    input_base=\$(basename "${assembly}" | sed 's/\\.[^.]*\$//')

    # Copy virus results
    if [ -f "genomad_out/\${input_base}_summary/\${input_base}_virus_summary.tsv" ]; then
        cp "genomad_out/\${input_base}_summary/\${input_base}_virus_summary.tsv" virus_summary.tsv
    else
        echo "[WARNING] No virus summary found" >&2
        touch virus_summary.tsv
    fi

    # Copy plasmid results
    if [ -f "genomad_out/\${input_base}_summary/\${input_base}_plasmid_summary.tsv" ]; then
        cp "genomad_out/\${input_base}_summary/\${input_base}_plasmid_summary.tsv" plasmid_summary.tsv
    else
        echo "[WARNING] No plasmid summary found" >&2
        touch plasmid_summary.tsv
    fi

    # Copy virus FASTA
    if [ -f "genomad_out/\${input_base}_summary/\${input_base}_virus.fna" ]; then
        cp "genomad_out/\${input_base}_summary/\${input_base}_virus.fna" virus.fna
    else
        touch virus.fna
    fi

    # Copy plasmid FASTA
    if [ -f "genomad_out/\${input_base}_summary/\${input_base}_plasmid.fna" ]; then
        cp "genomad_out/\${input_base}_summary/\${input_base}_plasmid.fna" plasmid.fna
    else
        touch plasmid.fna
    fi

    # Copy main summary (aggregated scores per contig)
    if [ -f "genomad_out/\${input_base}_aggregated_classification/\${input_base}_aggregated_classification.tsv" ]; then
        cp "genomad_out/\${input_base}_aggregated_classification/\${input_base}_aggregated_classification.tsv" genomad_summary.tsv
    else
        touch genomad_summary.tsv
    fi
    """
}

process CHECKV_QUALITY {
    tag "checkv"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-checkv"
    publishDir "${params.outdir}/mge/checkv", mode: 'copy'

    input:
    path(virus_fasta)

    output:
    path("quality_summary.tsv"), emit: quality
    path("viruses.fna"),         emit: viruses
    path("proviruses.fna"),      emit: proviruses

    script:
    def db_path = params.checkv_db
    """
    # CheckV: assess viral genome completeness and contamination
    # Trims host contamination from proviruses, estimates completeness via
    # AAI comparison to reference genomes or HMM-based gene density models

    if [ ! -s "${virus_fasta}" ]; then
        echo "[WARNING] No viral contigs to assess — skipping CheckV" >&2
        printf 'contig_id\\tcontig_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tcompleteness\\tcontamination\\n' > quality_summary.tsv
        touch viruses.fna proviruses.fna
        exit 0
    fi

    set +e
    checkv end_to_end \\
        "${virus_fasta}" \\
        checkv_out \\
        -d "${db_path}" \\
        -t ${task.cpus}
    checkv_exit=\$?
    set -e

    if [ \$checkv_exit -ne 0 ]; then
        echo "[WARNING] CheckV exited with code \$checkv_exit" >&2
        printf 'contig_id\\tcontig_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tcompleteness\\tcontamination\\n' > quality_summary.tsv
        touch viruses.fna proviruses.fna
        exit 0
    fi

    # Copy outputs
    if [ -f checkv_out/quality_summary.tsv ]; then
        cp checkv_out/quality_summary.tsv .
    else
        printf 'contig_id\\tcontig_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tcompleteness\\tcontamination\\n' > quality_summary.tsv
    fi

    if [ -f checkv_out/viruses.fna ]; then
        cp checkv_out/viruses.fna .
    else
        touch viruses.fna
    fi

    if [ -f checkv_out/proviruses.fna ]; then
        cp checkv_out/proviruses.fna .
    else
        touch proviruses.fna
    fi
    """
}
