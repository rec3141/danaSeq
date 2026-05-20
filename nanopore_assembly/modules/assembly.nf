// Flye metagenomic assembly without polishing (--iterations 0).
// Polishing is handled by the separate FLYE_POLISH process which can
// run on any assembler's output. The flye.yml conda env includes
// samtools>=1.17 to replace Flye's bundled samtools 1.9 which deadlocks
// on large BAMs (github.com/samtools/htslib/issues/831).

process FLYE_ASSEMBLE {
    tag "flye-assemble"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-assembly"
    publishDir "${params.outdir}/assembly", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/assembly" : null

    input:
    path("all_reads.fastq.gz")

    output:
    path("assembly.fasta"),      emit: assembly
    path("assembly_info.txt"),   emit: info
    path("assembly_graph.gfa"),  emit: graph

    script:
    """
    # Auto-detect read type from median quality scores
    if [ "${params.read_type}" = "auto" ]; then
        MEDIAN_Q=\$(zcat all_reads.fastq.gz | head -40000 | awk 'NR%4==0' | head -10000 | \
            python3 -c "
import sys
quals = []
for line in sys.stdin:
    line = line.strip()
    if line:
        avg = sum(ord(c)-33 for c in line) / len(line)
        quals.append(avg)
quals.sort()
print(int(quals[len(quals)//2]) if quals else 10)
")
        if [ -z "\$MEDIAN_Q" ] || [ "\$MEDIAN_Q" -lt 1 ]; then MEDIAN_Q=10; fi
        echo "[INFO] Median read quality: Q\${MEDIAN_Q}"
        if [ "\$MEDIAN_Q" -ge 20 ]; then
            FLYE_READ_TYPE="--nano-hq"
        else
            FLYE_READ_TYPE="--nano-raw"
        fi
    else
        FLYE_READ_TYPE="--${params.read_type}"
    fi
    echo "[INFO] Flye read type: \$FLYE_READ_TYPE"

    # Run Flye assembly without polishing (handled by FLYE_POLISH downstream)
    flye \\
        --meta \\
        --min-overlap ${params.min_overlap} \\
        --iterations 0 \\
        \$FLYE_READ_TYPE all_reads.fastq.gz \\
        --out-dir flye_out \\
        --threads ${task.cpus}

    # Validate assembly
    if [ ! -s flye_out/assembly.fasta ]; then
        echo "[ERROR] Flye produced no assembly output" >&2
        exit 1
    fi

    # Sort/rename contigs, remap assembly_info.txt + GFA (replaces fragile awk pipeline
    # that silently produced empty FASTA on multi-megabase contigs).
    normalize_assembly.py \\
        --input flye_out/assembly.fasta \\
        --output assembly.fasta \\
        --output-info assembly_info.txt \\
        --output-map name_map.tsv \\
        --header-format flye \\
        --source-info flye_out/assembly_info.txt \\
        --gfa-in flye_out/assembly_graph.gfa \\
        --gfa-out assembly_graph.gfa \\
        --min-len 1000

    rm -f name_map.tsv
    """
}

// Assembler-agnostic polishing using Flye's --polish-target.
// Runs on any assembly.fasta + reads — not Flye-specific despite using the Flye binary.
process FLYE_POLISH {
    tag "flye-polish"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-assembly"
    publishDir "${params.outdir}/assembly", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/assembly" : null

    input:
    path(assembly)
    path(info)
    path(graph)
    path(reads)

    output:
    path("assembly.fasta"),      emit: assembly
    path("assembly_info.txt"),   emit: info
    path("assembly_graph.gfa"),  emit: graph

    script:
    """
    # Determine read type
    if [ "${params.read_type}" = "auto" ]; then
        MEDIAN_Q=\$(zcat ${reads} | head -40000 | awk 'NR%4==0' | head -10000 | \
            python3 -c "
import sys
quals = []
for line in sys.stdin:
    line = line.strip()
    if line:
        avg = sum(ord(c)-33 for c in line) / len(line)
        quals.append(avg)
quals.sort()
print(int(quals[len(quals)//2]) if quals else 10)
")
        if [ -z "\$MEDIAN_Q" ] || [ "\$MEDIAN_Q" -lt 1 ]; then MEDIAN_Q=10; fi
        if [ "\$MEDIAN_Q" -ge 20 ]; then
            FLYE_READ_TYPE="--nano-hq"
        else
            FLYE_READ_TYPE="--nano-raw"
        fi
    else
        FLYE_READ_TYPE="--${params.read_type}"
    fi

    # Polish assembly with Flye's standalone polisher
    flye \\
        --polish-target ${assembly} \\
        \$FLYE_READ_TYPE ${reads} \\
        --out-dir polish_out \\
        --threads ${task.cpus}

    # Use polished output if it exists, otherwise keep original
    if [ -s polish_out/polished_1.fasta ]; then
        mv polish_out/polished_1.fasta assembly.fasta
        echo "[INFO] Polishing complete"
    else
        echo "[WARNING] Polishing produced no output, keeping unpolished assembly"
        cp ${assembly} assembly.fasta
    fi

    # Pass through info and graph unchanged (polishing doesn't alter these)
    # Use ln -f to handle case where staged input already has the target name
    [ "${info}" != "assembly_info.txt" ] && cp ${info} assembly_info.txt || true
    [ "${graph}" != "assembly_graph.gfa" ] && cp ${graph} assembly_graph.gfa || true
    """
}

process ASSEMBLY_METAMDBG {
    tag "co-assembly-metamdbg"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-assembly"
    publishDir "${params.outdir}/assembly", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/assembly" : null

    input:
    path("all_reads.fastq.gz")

    output:
    path("assembly.fasta"),      emit: assembly
    path("assembly_info.txt"),   emit: info
    path("assembly_graph.gfa"),  emit: graph

    script:
    """
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

    # Pick highest-k GFA if metaMDBG produced one
    METAMDBG_GFA=\$(ls metamdbg_out/tmp/pass_k*/assembly_graph.gfa 2>/dev/null | sort -t'k' -k2 -rn | head -1)

    # Sort/rename contigs, build Flye-style assembly_info.txt from headers,
    # remap GFA. Replaces fragile awk|sort|awk pipeline.
    normalize_assembly.py \\
        --input metamdbg_raw.fasta \\
        --output assembly.fasta \\
        --output-info assembly_info.txt \\
        --output-map name_map.tsv \\
        --header-format metamdbg \\
        \${METAMDBG_GFA:+--gfa-in "\$METAMDBG_GFA"} \\
        \${METAMDBG_GFA:+--gfa-out assembly_graph.gfa} \\
        --min-len 1000

    # If metaMDBG didn't emit a GFA, write a stub (S-lines only) from the renamed FASTA
    if [ ! -f assembly_graph.gfa ]; then
        python3 -c "
import sys
with open('assembly.fasta','rb') as fin, open('assembly_graph.gfa','wb') as fout:
    fout.write(b'H\\tVN:Z:1.0\\n')
    name = None
    seq = bytearray()
    for line in fin:
        if line.startswith(b'>'):
            if name:
                fout.write(b'S\\t' + name + b'\\t' + bytes(seq) + b'\\n')
            name = line.rstrip(b'\\r\\n')[1:]
            seq = bytearray()
        else:
            seq.extend(line.rstrip(b'\\r\\n'))
    if name:
        fout.write(b'S\\t' + name + b'\\t' + bytes(seq) + b'\\n')
"
    fi

    rm -f name_map.tsv metamdbg_raw.fasta
    """
}

process ASSEMBLY_MYLOASM {
    tag "co-assembly-myloasm"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-assembly"
    publishDir "${params.outdir}/assembly", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/assembly" : null

    input:
    path("all_reads.fastq.gz")

    output:
    path("assembly.fasta"),      emit: assembly
    path("assembly_info.txt"),   emit: info
    path("assembly_graph.gfa"),  emit: graph

    script:
    """
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

    # Sort/rename contigs, build assembly_info.txt from myloasm headers,
    # remap GFA. Replaces fragile awk|sort|awk that silently produced empty
    # outputs on multi-megabase contigs (which myloasm regularly emits).
    GFA_IN=myloasm_out/final_contig_graph.gfa
    normalize_assembly.py \\
        --input myloasm_out/assembly_primary.fa \\
        --output assembly.fasta \\
        --output-info assembly_info.txt \\
        --output-map name_map.tsv \\
        --header-format myloasm \\
        \$( [ -f "\$GFA_IN" ] && echo "--gfa-in \$GFA_IN --gfa-out assembly_graph.gfa" ) \\
        --min-len 1000

    # If myloasm didn't emit a GFA, write a stub (S-lines only) from the renamed FASTA
    if [ ! -f assembly_graph.gfa ]; then
        python3 -c "
import sys
with open('assembly.fasta','rb') as fin, open('assembly_graph.gfa','wb') as fout:
    fout.write(b'H\\tVN:Z:1.0\\n')
    name = None
    seq = bytearray()
    for line in fin:
        if line.startswith(b'>'):
            if name:
                fout.write(b'S\\t' + name + b'\\t' + bytes(seq) + b'\\n')
            name = line.rstrip(b'\\r\\n')[1:]
            seq = bytearray()
        else:
            seq.extend(line.rstrip(b'\\r\\n'))
    if name:
        fout.write(b'S\\t' + name + b'\\t' + bytes(seq) + b'\\n')
"
    fi

    rm -f name_map.tsv
    """
}

process CALCULATE_TNF {
    tag "tnf"
    label 'process_low'
    conda "${projectDir}/conda-envs/dana-mag-assembly"
    publishDir "${params.outdir}/assembly", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/assembly" : null

    input:
    path(assembly)

    output:
    path("tnf.tsv"), emit: tnf
    path("gc.tsv"),  emit: gc

    script:
    """
    tetramer_freqs -f ${assembly} -min 0 -max 999999999 -o tnf.tsv

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

# Build canonical (RC-collapsed) 4-mers in same order as tetramer_freqs
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
