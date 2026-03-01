// Co-assembly: concatenate all reads, optional filtlong, Flye metagenomic assembly

process ASSEMBLY_FLYE {
    tag "co-assembly"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-flye"
    publishDir "${params.outdir}/assembly", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/assembly" : null

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
        ${params.polish ? '' : '--iterations 0'} \\
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

    # Sort contigs by length (longest first) and rename with zero-padded IDs
    NCONTIGS=\$(grep -c '^>' flye_out/assembly.fasta)
    PAD=\${#NCONTIGS}
    TAB=\$(printf '\\t')

    # Build sorted FASTA + name map (old_name → new_name)
    awk '/^>/{if(h) print h "\\t" length(s) "\\t" s; h=substr(\$0,2); s=""; next} {s=s\$0} END{if(h) print h "\\t" length(s) "\\t" s}' flye_out/assembly.fasta \\
    | sort -t"\$TAB" -k2,2rn \\
    | awk -F'\\t' -v pad="\$PAD" '{n++; new=sprintf("contig_%0*d",pad,n); printf ">%s\\n",new>"assembly.fasta"; s=\$3; for(i=1;i<=length(s);i+=80) print substr(s,i,80)>"assembly.fasta"; print \$1"\\t"new>"name_map.tsv"}'

    # Update assembly_info.txt: rename first column, sort by length descending
    head -1 flye_out/assembly_info.txt > assembly_info.txt
    tail -n+2 flye_out/assembly_info.txt \\
    | awk -F'\\t' 'BEGIN{OFS="\\t"} NR==FNR{map[\$1]=\$2;next} {if(\$1 in map) \$1=map[\$1]; print}' name_map.tsv - \\
    | sort -t"\$TAB" -k2,2rn >> assembly_info.txt

    # Update GFA: rename contig names in P (path) lines only
    # S and L lines use edge_N names (graph structure) — unchanged
    awk -F'\\t' 'BEGIN{OFS="\\t"} NR==FNR{map[\$1]=\$2;next} \$1=="P"&&(\$2 in map){\$2=map[\$2]} {print}' name_map.tsv flye_out/assembly_graph.gfa > assembly_graph.gfa

    rm -f name_map.tsv

    # Cleanup intermediate files
    rm -f all_reads.fastq.gz
    """
}

process ASSEMBLY_METAMDBG {
    tag "co-assembly-metamdbg"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-metamdbg"
    publishDir "${params.outdir}/assembly", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/assembly" : null

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

    # Run metaMDBG assembly
    metaMDBG asm \\
        --out-dir metamdbg_out \\
        --in-ont all_reads.fastq.gz \\
        --threads ${task.cpus}

    # Generate GFA (separate metaMDBG command)
    # Find the highest available k value for the most resolved graph
    MAX_K=\$(metaMDBG gfa --assembly-dir metamdbg_out --k 0 2>&1 | awk '/^\\t- /{k=\$2} END{print k}')
    if [ -n "\$MAX_K" ]; then
        metaMDBG gfa --assembly-dir metamdbg_out --k "\$MAX_K" --threads ${task.cpus} || true
    fi

    # Validate assembly — metaMDBG outputs contigs.fasta.gz
    if [ ! -s metamdbg_out/contigs.fasta.gz ]; then
        echo "[ERROR] metaMDBG produced no assembly output" >&2
        exit 1
    fi

    # Decompress contigs
    gunzip -c metamdbg_out/contigs.fasta.gz > metamdbg_raw.fasta

    # Check minimum contig size
    max_len=\$(awk '/^>/{if(len>0) print len; len=0; next} {len+=length(\$0)} END{print len}' metamdbg_raw.fasta | sort -rn | head -1)
    if [ "\$max_len" -lt 1000 ]; then
        echo "[ERROR] Assembly contains no contigs >= 1000bp (max: \${max_len}bp)" >&2
        exit 1
    fi

    # Sort contigs by length (longest first) and rename with zero-padded IDs
    NCONTIGS=\$(grep -c '^>' metamdbg_raw.fasta)
    PAD=\${#NCONTIGS}
    TAB=\$(printf '\\t')

    # Parse metaMDBG headers: >ctgN length=LEN circular={yes,no}
    # Extract old_name (first word only), length, coverage, circularity
    # Build assembly_info.txt with Flye-compatible columns
    awk '/^>/{if(h) print h "\\t" length(s) "\\t" s; h=substr(\$0,2); s=""; next} {s=s\$0} END{if(h) print h "\\t" length(s) "\\t" s}' metamdbg_raw.fasta \\
    | sort -t"\$TAB" -k2,2rn \\
    | awk -F'\\t' -v pad="\$PAD" '{n++; new=sprintf("contig_%0*d",pad,n); printf ">%s\\n",new>"assembly.fasta"; s=\$3; for(i=1;i<=length(s);i+=80) print substr(s,i,80)>"assembly.fasta"; split(\$1,w," "); print w[1]"\\t"new>"name_map.tsv"}'

    # Build assembly_info.txt by parsing metaMDBG FASTA headers
    # Header format: >ctgN length=LEN circular={yes,no}
    echo -e "#seq_name\\tlength\\tcov.\\tcirc." > assembly_info.txt
    awk '/^>/{
        hdr = substr(\$0, 2)
        split(hdr, parts, " ")
        name = parts[1]
        len = 0; cov = 0; circ = "N"
        for (i = 2; i <= length(parts); i++) {
            if (parts[i] ~ /^length=/) { split(parts[i], kv, "="); len = kv[2] }
            if (parts[i] ~ /^coverage=/) { split(parts[i], kv, "="); cov = kv[2] }
            if (parts[i] ~ /^circular=/) { split(parts[i], kv, "="); circ = (kv[2] == "yes") ? "Y" : "N" }
        }
        print name "\\t" len "\\t" cov "\\t" circ
    }' metamdbg_raw.fasta \\
    | awk -F'\\t' 'BEGIN{OFS="\\t"} NR==FNR{map[\$1]=\$2;next} {if(\$1 in map) \$1=map[\$1]; print}' name_map.tsv - \\
    | sort -t"\$TAB" -k2,2rn >> assembly_info.txt

    # GFA: use metaMDBG graph if available (in tmp/pass_kN/assembly_graph.gfa), otherwise stub
    METAMDBG_GFA=\$(ls metamdbg_out/tmp/pass_k*/assembly_graph.gfa 2>/dev/null | sort -t'k' -k2 -rn | head -1)
    if [ -n "\$METAMDBG_GFA" ] && [ -f "\$METAMDBG_GFA" ]; then
        awk -F'\\t' 'BEGIN{OFS="\\t"} NR==FNR{map[\$1]=\$2;next} \$1=="S"&&(\$2 in map){\$2=map[\$2]} \$1=="P"&&(\$2 in map){\$2=map[\$2]} {print}' name_map.tsv "\$METAMDBG_GFA" > assembly_graph.gfa
    else
        # Stub GFA with S lines only (no graph topology)
        echo "H\\tVN:Z:1.0" > assembly_graph.gfa
        awk '/^>/{name=substr(\$0,2); next} {seq=seq\$0} /^>/{if(name) print "S\\t"name"\\t"seq; seq=""} END{if(name) print "S\\t"name"\\t"seq}' assembly.fasta \\
        | awk -F'\\t' 'NF==3' >> assembly_graph.gfa
    fi

    rm -f name_map.tsv metamdbg_raw.fasta
    rm -f all_reads.fastq.gz
    """
}

process ASSEMBLY_MYLOASM {
    tag "co-assembly-myloasm"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-myloasm"
    publishDir "${params.outdir}/assembly", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/assembly" : null

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

    # Run myloasm assembly
    myloasm all_reads.fastq.gz \\
        -o myloasm_out \\
        -t ${task.cpus} \\
        --clean-dir

    # Validate assembly — myloasm outputs assembly_primary.fa
    if [ ! -s myloasm_out/assembly_primary.fa ]; then
        echo "[ERROR] myloasm produced no assembly output" >&2
        exit 1
    fi

    # Check minimum contig size
    max_len=\$(awk '/^>/{if(len>0) print len; len=0; next} {len+=length(\$0)} END{print len}' myloasm_out/assembly_primary.fa | sort -rn | head -1)
    if [ "\$max_len" -lt 1000 ]; then
        echo "[ERROR] Assembly contains no contigs >= 1000bp (max: \${max_len}bp)" >&2
        exit 1
    fi

    # Sort contigs by length (longest first) and rename with zero-padded IDs
    NCONTIGS=\$(grep -c '^>' myloasm_out/assembly_primary.fa)
    PAD=\${#NCONTIGS}
    TAB=\$(printf '\\t')

    # Build sorted FASTA + name map (old_name first word → new_name)
    awk '/^>/{if(h) print h "\\t" length(s) "\\t" s; h=substr(\$0,2); s=""; next} {s=s\$0} END{if(h) print h "\\t" length(s) "\\t" s}' myloasm_out/assembly_primary.fa \\
    | sort -t"\$TAB" -k2,2rn \\
    | awk -F'\\t' -v pad="\$PAD" '{n++; new=sprintf("contig_%0*d",pad,n); printf ">%s\\n",new>"assembly.fasta"; s=\$3; for(i=1;i<=length(s);i+=80) print substr(s,i,80)>"assembly.fasta"; split(\$1,w," "); print w[1]"\\t"new>"name_map.tsv"}'

    # Build assembly_info.txt by parsing myloasm FASTA headers
    # Header format: >NAME_len-LEN_circular-{yes,no,possibly}_depth-D1-D2-D3_duplicated-{yes,no} mult=M
    echo -e "#seq_name\\tlength\\tcov.\\tcirc." > assembly_info.txt
    awk '/^>/{
        hdr = substr(\$0, 2)
        split(hdr, parts, " ")
        name = parts[1]
        len = 0; cov = 0; circ = "N"
        # Parse underscore-delimited fields in the contig name
        n_fields = split(name, fields, "_")
        for (i = 1; i <= n_fields; i++) {
            if (fields[i] ~ /^len-/) { split(fields[i], kv, "-"); len = kv[2] }
            if (fields[i] ~ /^circular-/) {
                split(fields[i], kv, "-")
                circ = (kv[2] == "yes" || kv[2] == "possibly") ? "Y" : "N"
            }
            if (fields[i] ~ /^depth-/) { split(fields[i], kv, "-"); cov = kv[2] }
        }
        # Also check for mult= in the space-separated part
        for (j = 2; j <= length(parts); j++) {
            if (parts[j] ~ /^mult=/) { split(parts[j], kv, "="); if (cov == 0) cov = kv[2] }
        }
        print name "\\t" len "\\t" cov "\\t" circ
    }' myloasm_out/assembly_primary.fa \\
    | awk -F'\\t' 'BEGIN{OFS="\\t"} NR==FNR{map[\$1]=\$2;next} {if(\$1 in map) \$1=map[\$1]; print}' name_map.tsv - \\
    | sort -t"\$TAB" -k2,2rn >> assembly_info.txt

    # GFA: use myloasm output if available, otherwise generate stub
    if [ -f myloasm_out/final_contig_graph.gfa ]; then
        awk -F'\\t' 'BEGIN{OFS="\\t"} NR==FNR{map[\$1]=\$2;next} \$1=="S"&&(\$2 in map){\$2=map[\$2]} \$1=="P"&&(\$2 in map){\$2=map[\$2]} {print}' name_map.tsv myloasm_out/final_contig_graph.gfa > assembly_graph.gfa
    else
        # Stub GFA with S lines only
        echo "H\\tVN:Z:1.0" > assembly_graph.gfa
        awk '/^>/{name=substr(\$0,2); next} {seq=seq\$0} /^>/{if(name) print "S\\t"name"\\t"seq; seq=""} END{if(name) print "S\\t"name"\\t"seq}' assembly.fasta \\
        | awk -F'\\t' 'NF==3' >> assembly_graph.gfa
    fi

    rm -f name_map.tsv
    rm -f all_reads.fastq.gz
    """
}

process CALCULATE_TNF {
    tag "tnf"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-flye"
    publishDir "${params.outdir}/assembly", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/assembly" : null

    input:
    path(assembly)

    output:
    path("tnf.tsv"), emit: tnf
    path("gc.tsv"),  emit: gc

    script:
    """
    tetramer_freqs.py -f ${assembly} -min 0 -max 999999999 -o tnf.tsv

    if [ ! -s tnf.tsv ]; then
        echo "[ERROR] TNF calculation produced empty output" >&2
        exit 1
    fi

    # Compute per-contig GC% from tetranucleotide frequencies
    # Each row: contig_id followed by 136 RC-collapsed 4-mer relative abundances.
    # GC% = weighted sum of each 4-mer frequency by its G+C base count / 4.
    python3 - tnf.tsv gc.tsv <<'PYEOF'
import sys
from itertools import product

BASES = "ACGT"
COMP = str.maketrans("ACGT", "TGCA")

# Build canonical (RC-collapsed) 4-mers in same order as tetramer_freqs.py
seen = set()
canonical = []
for kmer in product(BASES, repeat=4):
    kmer = "".join(kmer)
    rc = kmer.translate(COMP)[::-1]
    if rc not in seen:
        canonical.append(kmer)
        seen.add(kmer)

# Weight = fraction of G+C bases in each tetramer
gc_weights = [sum(1 for b in k if b in "GC") / 4.0 for k in canonical]

with open(sys.argv[1]) as fin, open(sys.argv[2], "w") as fout:
    fout.write("contig_id\\tgc_pct\\n")
    for line in fin:
        cols = line.rstrip("\\n").split("\\t")
        contig_id = cols[0]
        freqs = [float(x) for x in cols[1:]]
        gc = sum(f * w for f, w in zip(freqs, gc_weights))
        fout.write(f"{contig_id}\\t{gc * 100:.2f}\\n")
PYEOF
    """
}
