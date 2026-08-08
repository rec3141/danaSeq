// Demultiplexing with Mr_Demuxy.
//
// Takes raw paired-end FASTQ files grouped by lane/flowcell and demultiplexes
// using inner barcodes. Produces per-sample paired-end FASTQ files.
// This step is optional — most users will start with pre-demultiplexed reads.

process DEMULTIPLEX {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/envs/python.yml"
    publishDir "${params.outdir}/demultiplexed", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/demultiplexed" : null

    input:
    tuple val(meta), path(fastqs)
    path(forward_bcs)
    path(reverse_bcs)

    output:
    tuple val(meta), path("demuxed/*.fastq.gz"), emit: reads

    script:
    """
    mkdir -p demuxed

    # Find R1 and R2 files
    R1=\$(ls *_R1_*.fastq.gz 2>/dev/null || ls *_R1.fastq.gz 2>/dev/null || true)
    R2=\$(ls *_R2_*.fastq.gz 2>/dev/null || ls *_R2.fastq.gz 2>/dev/null || true)

    if [ -z "\$R1" ] || [ -z "\$R2" ]; then
        echo "[ERROR] ${meta.id}: Could not find paired R1/R2 files" >&2
        exit 1
    fi

    # Decompress for Mr_Demuxy (requires uncompressed FASTQ)
    gunzip -c \$R1 > input_R1.fastq
    gunzip -c \$R2 > input_R2.fastq

    # Run demultiplexer
    pe_demuxer.py \\
        -r1 input_R1.fastq \\
        -r2 input_R2.fastq \\
        -r1_bc ${forward_bcs} \\
        -r2_bc ${reverse_bcs}

    # Rename and compress outputs
    if [ -d pe_demuxer_output ]; then
        for DIR in R1 R2; do
            if [ -d pe_demuxer_output/\$DIR ]; then
                for FILE in pe_demuxer_output/\$DIR/*.fastq; do
                    BASENAME=\$(basename \$FILE)
                    NEWNAME="${meta.id}_\${BASENAME}"
                    gzip -c \$FILE > demuxed/\${NEWNAME}.gz
                done
            fi
        done
    fi

    # Cleanup
    rm -f input_R1.fastq input_R2.fastq
    rm -rf pe_demuxer_output

    echo "[INFO] ${meta.id}: demultiplexed \$(ls demuxed/*.fastq.gz 2>/dev/null | wc -l) files" >&2
    """
}
