// Co-assembly: concatenate all reads, optional filtlong, Flye metagenomic assembly

process ASSEMBLY_FLYE {
    tag "co-assembly"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-flye"
    publishDir "${params.outdir}/assembly", mode: 'copy'
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

process CALCULATE_TNF {
    tag "tnf"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-flye"
    publishDir "${params.outdir}/assembly", mode: 'copy'
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
