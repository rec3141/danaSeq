// Co-assembly: concatenate all reads, optional filtlong, Flye metagenomic assembly

process ASSEMBLY_FLYE {
    tag "co-assembly"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-flye"
    publishDir "${params.outdir}/assembly", mode: 'copy'

    input:
    path(reads)

    output:
    path("assembly.fasta"),      emit: assembly
    path("assembly_info.txt"),   emit: info
    path("assembly_graph.gfa"),  emit: graph

    script:
    def filtlong_cmd = ""
    if (params.filtlong_size) {
        filtlong_cmd = """
        filtlong -t ${params.filtlong_size} all_reads.fastq.gz | gzip > all_reads_filt.fastq.gz
        mv all_reads_filt.fastq.gz all_reads.fastq.gz
        """
    }

    """
    # Concatenate all input reads
    cat ${reads} > all_reads.fastq.gz

    ${filtlong_cmd}

    # Run Flye metagenomic assembly
    flye \\
        --meta \\
        --min-overlap ${params.min_overlap} \\
        --nano-hq all_reads.fastq.gz \\
        --out-dir flye_out \\
        --threads ${task.cpus}

    # Validate assembly
    if [ ! -s flye_out/assembly.fasta ]; then
        echo "[ERROR] Flye produced no assembly output" >&2
        exit 1
    fi

    # Check minimum contig size
    max_len=\$(awk '/^>/{if(len>0) print len; len=0; next} {len+=length(\$0)} END{print len}' flye_out/assembly.fasta | sort -rn | head -1)
    if [ "\$max_len" -lt 1000 ]; then
        echo "[ERROR] Assembly contains no contigs >= 1000bp (max: \${max_len}bp)" >&2
        exit 1
    fi

    cp flye_out/assembly.fasta assembly.fasta
    cp flye_out/assembly_info.txt assembly_info.txt
    cp flye_out/assembly_graph.gfa assembly_graph.gfa

    # Cleanup intermediate files
    rm -f all_reads.fastq.gz
    """
}
