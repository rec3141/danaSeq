// Binning: SemiBin2, MetaBAT2, MaxBin2, LorBin, COMEBin, VAMB,
//          DAS_Tool consensus, Binette consensus, MAGScoT consensus.
//
// Up to six binners run in parallel, each emitting a DAS_Tool-format TSV (contig\tbin).
// All outputs are mixed into ch_binner_results in main.nf and collected by
// DASTOOL_CONSENSUS for score-based consensus integration.
// BINETTE_CONSENSUS and MAGSCOT_CONSENSUS run in parallel with DAS_Tool as
// alternative consensus methods for quality comparison.
//
// Each binner catches failures and emits an empty TSV so the pipeline continues.
// DAS_Tool filters out empty inputs and handles the "no bins above threshold" case.
// CHECKM2 runs on all binner + consensus FASTAs for quality assessment.

process BIN_SEMIBIN2 {
    tag "semibin2"
    label 'process_gpu'
    conda "${projectDir}/conda-envs/dana-mag-binning"
    publishDir "${params.outdir}/binning/semibin", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/semibin" : null

    input:
    path(assembly)
    path(bams)

    output:
    path("semibin_bins.tsv"), emit: bins
    path("bins/"),            emit: fastas

    script:
    """
    mkdir -p bins

    # SemiBin2 can crash on very small datasets (0 bins → empty ORFs → hmmsearch fail)
    # Catch failures and produce an empty output so the pipeline can continue
    set +e
    SemiBin2 single_easy_bin \\
        -i "${assembly}" \\
        -b *.sorted.bam \\
        -o semibin_out \\
        --sequencing-type long_read
    semibin_exit=\$?
    set -e

    if [ \$semibin_exit -ne 0 ]; then
        echo "[WARNING] SemiBin2 exited with code \$semibin_exit (dataset may be too small)" >&2
        touch semibin_bins.tsv
    elif [ -d semibin_out/output_bins ] && ls semibin_out/output_bins/*.fa.gz 1>/dev/null 2>&1; then
        # Build TSV from FASTA headers (single source of truth)
        # Extract bin number from SemiBin2 filename (e.g. SemiBin2_42.fa.gz -> 43)
        > semibin_bins.tsv
        for bin_file in semibin_out/output_bins/*.fa.gz; do
            [ -e "\$bin_file" ] || continue
            bin_num=\$(basename "\$bin_file" .fa.gz | grep -oP '\\d+\$')
            bin_name=\$(printf 'semibin_%03d' \$((bin_num + 1)))
            zcat "\$bin_file" > "bins/\${bin_name}.fa"
            zcat "\$bin_file" | grep '>' | tr -d '>' | cut -f1 -d' ' | while read -r contig; do
                printf '%s\\t%s\\n' "\$contig" "\$bin_name"
            done >> semibin_bins.tsv
        done
    else
        echo "[WARNING] SemiBin2 produced no bins" >&2
        touch semibin_bins.tsv
    fi

    if [ ! -s semibin_bins.tsv ]; then
        echo "[WARNING] SemiBin2 produced no bins" >&2
    fi
    """
}

process BIN_METABAT2 {
    tag "metabat2"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-binning"
    publishDir "${params.outdir}/binning/metabat", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/metabat" : null

    input:
    path(assembly)
    path(jgi_depth)

    output:
    path("metabat_bins.tsv"), emit: bins
    path("bins/"),            emit: fastas

    script:
    """
    mkdir -p bins

    metabat2 \\
        -i "${assembly}" \\
        -o metabat_out/bin \\
        --saveCls \\
        --minClsSize ${params.metabat_min_cls} \\
        -a "${jgi_depth}"

    # Convert MetaBAT2 FASTA bins to contig_bins.tsv format
    # Skip unbinned/tooShort/lowDepth catch-all files from --saveCls
    > metabat_bins.tsv
    bin_num=0
    for bin_file in metabat_out/bin*.fa; do
        [ -e "\$bin_file" ] || continue
        case "\$(basename "\$bin_file")" in
            *unbinned*|*tooShort*|*lowDepth*) continue ;;
        esac
        bin_num=\$((bin_num + 1))
        bin_name=\$(printf 'metabat_%03d' \$bin_num)
        cp "\$bin_file" "bins/\${bin_name}.fa"
        grep '>' "\$bin_file" | tr -d '>' | cut -f1 -d' ' | while read -r contig; do
            printf '%s\\t%s\\n' "\$contig" "\$bin_name"
        done >> metabat_bins.tsv
    done

    if [ ! -s metabat_bins.tsv ]; then
        echo "[WARNING] MetaBAT2 produced no bins" >&2
    fi
    """
}

process BIN_MAXBIN2 {
    tag "maxbin2"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-binning"
    publishDir "${params.outdir}/binning/maxbin", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/maxbin" : null

    input:
    path(assembly)
    path(jgi_depth)

    output:
    path("maxbin_bins.tsv"), emit: bins
    path("bins/"),           emit: fastas

    script:
    """
    mkdir -p bins

    # Extract coverage column from JGI depth table for MaxBin2 format
    # JGI format: contigName contigLen totalAvgDepth sample1.var sample2 ...
    # MaxBin2 wants: contigName avgDepth
    cut -f1,3 "${jgi_depth}" | tail -n +2 > coverage.txt

    mkdir -p maxbin_out
    run_MaxBin.pl \\
        -contig "${assembly}" \\
        -abund coverage.txt \\
        -out maxbin_out/bin \\
        -thread ${task.cpus}

    # Convert MaxBin2 FASTA bins to contig_bins.tsv format
    # Extract bin number from filename (e.g. bin.042.fasta -> 42) for consistent naming
    > maxbin_bins.tsv
    for bin_file in maxbin_out/bin*.fasta; do
        [ -e "\$bin_file" ] || continue
        bin_num=\$(basename "\$bin_file" .fasta | grep -oP '\\d+\$')
        bin_name=\$(printf 'maxbin_%03d' "\$((10#\$bin_num))")
        cp "\$bin_file" "bins/\${bin_name}.fa"
        grep '>' "\$bin_file" | tr -d '>' | cut -f1 -d' ' | while read -r contig; do
            printf '%s\\t%s\\n' "\$contig" "\$bin_name"
        done >> maxbin_bins.tsv
    done

    if [ ! -s maxbin_bins.tsv ]; then
        echo "[WARNING] MaxBin2 produced no bins" >&2
    fi
    """
}

process BIN_LORBIN {
    tag "lorbin"
    label 'process_gpu'
    conda "${projectDir}/conda-envs/dana-mag-binning"
    publishDir "${params.outdir}/binning/lorbin", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/lorbin" : null

    input:
    path(assembly)
    path(bams)

    output:
    path("lorbin_bins.tsv"), emit: bins
    path("bins/"),            emit: fastas

    script:
    """
    mkdir -p bins

    # LorBin can fail on very small or unusual datasets
    # Note: --multi is only for concatenated per-sample assemblies (LorBin concat)
    # where contig names have a sample prefix delimited by '-'.
    # Co-assemblies (Flye) use plain contig names, so we omit --multi.
    set +e
    LorBin bin \\
        -o lorbin_out \\
        -fa "${assembly}" \\
        -b *.sorted.bam \\
        --num_process ${task.cpus} \\
        --bin_length ${params.lorbin_min_length}
    lorbin_exit=\$?
    set -e

    if [ \$lorbin_exit -ne 0 ]; then
        echo "[WARNING] LorBin exited with code \$lorbin_exit (dataset may be too small)" >&2
        touch lorbin_bins.tsv
    elif [ -d lorbin_out/output_bins ]; then
        # Convert LorBin FASTA bins to DAS_Tool format TSV (contig\\tbin)
        # Extract bin number from filename (e.g. bin.42.fa -> 42) for consistent naming
        > lorbin_bins.tsv
        for bin_file in lorbin_out/output_bins/bin.*.fa; do
            [ -e "\$bin_file" ] || continue
            bin_num=\$(basename "\$bin_file" .fa | grep -oP '\\d+\$')
            bin_name=\$(printf 'lorbin_%03d' "\$((10#\$bin_num))")
            cp "\$bin_file" "bins/\${bin_name}.fa"
            grep '>' "\$bin_file" | tr -d '>' | cut -f1 -d' ' | while read -r contig; do
                printf '%s\\t%s\\n' "\$contig" "\$bin_name"
            done >> lorbin_bins.tsv
        done
    else
        echo "[WARNING] LorBin produced no output_bins directory" >&2
        touch lorbin_bins.tsv
    fi

    if [ ! -s lorbin_bins.tsv ]; then
        echo "[WARNING] LorBin produced no bins" >&2
    fi
    """
}

process BIN_COMEBIN {
    tag "comebin"
    label 'process_gpu'
    conda "${projectDir}/conda-envs/dana-mag-binning"
    publishDir "${params.outdir}/binning/comebin", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/comebin" : null

    input:
    path(assembly)
    path(bams)

    output:
    path("comebin_bins.tsv"), emit: bins
    path("bins/"),            emit: fastas

    script:
    """
    mkdir -p bins

    # COMEBin (contrastive multi-view binning) — deep learning binner
    # run_comebin.sh wraps the two-step generate_coverage + run_comebin workflow
    # -p . because BAMs are staged in the working directory
    set +e
    run_comebin.sh \\
        -a "${assembly}" \\
        -o comebin_out \\
        -p . \\
        -t ${task.cpus}
    comebin_exit=\$?
    set -e

    if [ \$comebin_exit -ne 0 ]; then
        echo "[WARNING] COMEBin exited with code \$comebin_exit (dataset may be too small)" >&2
        touch comebin_bins.tsv
    elif [ -d comebin_out/comebin_res/comebin_res_bins ]; then
        # Convert COMEBin FASTA bins to DAS_Tool format TSV (contig\\tbin)
        # COMEBin uses sparse bin numbers (0, 11, 768, 4926, ...) — renumber sequentially
        > comebin_bins.tsv
        bin_num=0
        for bin_file in comebin_out/comebin_res/comebin_res_bins/*.fa; do
            [ -e "\$bin_file" ] || continue
            bin_num=\$((bin_num + 1))
            bin_name=\$(printf 'comebin_%03d' \$bin_num)
            cp "\$bin_file" "bins/\${bin_name}.fa"
            grep '>' "\$bin_file" | tr -d '>' | cut -f1 -d' ' | while read -r contig; do
                printf '%s\\t%s\\n' "\$contig" "\$bin_name"
            done >> comebin_bins.tsv
        done
    else
        echo "[WARNING] COMEBin produced no comebin_res_bins directory" >&2
        touch comebin_bins.tsv
    fi

    if [ ! -s comebin_bins.tsv ]; then
        echo "[WARNING] COMEBin produced no bins" >&2
    fi
    """
}

process BIN_VAMB {
    tag "vamb"
    label 'process_gpu'
    conda "${projectDir}/conda-envs/dana-mag-vamb"
    publishDir "${params.outdir}/binning/vamb", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/vamb" : null

    input:
    path(assembly)
    path(depths)

    output:
    path("vamb_bins.tsv"), emit: bins
    path("bins/"),         emit: fastas

    script:
    """
    mkdir -p bins

    # VAMB — variational autoencoder for metagenomic binning (TNF + coverage)
    # Convert CoverM MetaBAT2-format depth table to VAMB abundance TSV:
    #   - Rename contigName → contigname (VAMB requirement)
    #   - Drop contigLen, totalAvgDepth, and *-var columns
    #   - Drop samples (columns) where all depths are zero
    python3 -c "
import sys, csv

with open('${depths}') as f:
    reader = csv.reader(f, delimiter='\\t')
    header = next(reader)

# Identify depth columns (skip contigName, contigLen, totalAvgDepth, and *-var)
depth_cols = []
for i, h in enumerate(header):
    if i >= 3 and not h.endswith('-var'):
        depth_cols.append(i)

# Read all rows (VAMB does its own length filtering with -m 2000)
rows = []
with open('${depths}') as f:
    reader = csv.reader(f, delimiter='\\t')
    next(reader)  # skip header
    for row in reader:
        rows.append(row)

# Find columns with non-zero sum (VAMB crashes on all-zero samples)
keep = []
for j, ci in enumerate(depth_cols):
    col_sum = sum(float(rows[r][ci]) for r in range(len(rows)))
    if col_sum > 0:
        keep.append(ci)

skipped = len(depth_cols) - len(keep)
print(f'[INFO] VAMB abundance: {len(rows)} contigs, {len(keep)} samples ({skipped} zero-depth dropped)', file=sys.stderr)

# Write VAMB abundance TSV
with open('vamb_abundance.tsv', 'w') as out:
    # Header: contigname + sample names
    sample_names = [header[ci] for ci in keep]
    out.write('contigname\\t' + '\\t'.join(sample_names) + '\\n')
    for row in rows:
        out.write(row[0] + '\\t' + '\\t'.join(row[ci] for ci in keep) + '\\n')
"

    n_samples=\$(head -1 vamb_abundance.tsv | awk -F'\\t' '{print NF-1}')
    if [ "\$n_samples" -lt 1 ]; then
        echo "[WARNING] No samples with non-zero depth — skipping VAMB" >&2
        touch vamb_bins.tsv
    else

    # Auto-detect GPU: use --cuda if available, fall back to CPU
    CUDA_FLAG=""
    if python3 -c "import torch; assert torch.cuda.is_available()" 2>/dev/null; then
        CUDA_FLAG="--cuda"
        echo "[INFO] VAMB: GPU detected, using CUDA" >&2
    else
        echo "[INFO] VAMB: No GPU detected, using CPU" >&2
    fi

    set +e
    vamb bin default \
        --outdir vamb_out \
        --fasta "${assembly}" \
        --abundance_tsv vamb_abundance.tsv \
        -m 2000 \
        --minfasta 200000 \
        \$CUDA_FLAG \
        -p ${task.cpus}
    vamb_exit=\$?
    set -e

    if [ \$vamb_exit -ne 0 ]; then
        echo "[WARNING] VAMB exited with code \$vamb_exit" >&2
        touch vamb_bins.tsv
    elif [ -d vamb_out/bins ] && ls vamb_out/bins/*.fna 1>/dev/null 2>&1; then
        > vamb_bins.tsv
        bin_num=0
        for bin_file in vamb_out/bins/*.fna; do
            [ -e "\$bin_file" ] || continue
            bin_num=\$((bin_num + 1))
            bin_name=\$(printf 'vamb_%03d' \$bin_num)
            cp "\$bin_file" "bins/\${bin_name}.fa"
            grep '>' "\$bin_file" | tr -d '>' | cut -f1 -d' ' | while read -r contig; do
                printf '%s\\t%s\\n' "\$contig" "\$bin_name"
            done >> vamb_bins.tsv
        done
        echo "[INFO] VAMB produced \$bin_num bins" >&2
    elif [ -f vamb_out/vae_clusters_unsplit.tsv ]; then
        # No FASTA bins (all below --minfasta), build TSV from cluster file
        # VAMB format: clustername\tcontigname — swap to contig\tbin
        tail -n +2 vamb_out/vae_clusters_unsplit.tsv | \
            awk -F'\\t' '{print \$2 "\\t" \$1}' > vamb_bins.tsv
        echo "[INFO] VAMB produced clusters but no FASTA bins above size threshold" >&2
    else
        echo "[WARNING] VAMB produced no output" >&2
        touch vamb_bins.tsv
    fi

    if [ ! -s vamb_bins.tsv ]; then
        echo "[WARNING] VAMB produced no bins" >&2
    fi

    fi  # end: samples check
    """
}

process BIN_VAMB_TAX {
    tag "vamb_tax"
    label 'process_gpu'
    conda "${projectDir}/conda-envs/dana-mag-vamb"
    publishDir "${params.outdir}/binning/vamb_tax", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/vamb_tax" : null

    input:
    path(assembly)
    path(depths)
    path(sendsketch_contigs)

    output:
    path("vamb_tax_bins.tsv"), emit: bins
    path("bins/"),             emit: fastas

    script:
    """
    mkdir -p bins

    # Convert sendsketch per-contig GTDB taxonomy to VAMB format
    # sendsketch: contig_id\tstatus\tANI\tref_name\td__X;p__Y;...
    # VAMB wants: contigname\tDomain;Phylum;Class;... (no GTDB prefixes, no empty ranks)
    # VAMB requires ALL contigs ≥ -m length to be in the taxonomy file
    python3 -c "
import sys

# 1. Read sendsketch taxonomy
tax = {}
with open('${sendsketch_contigs}') as f:
    f.readline()  # skip header
    for line in f:
        parts = line.strip().split('\\t')
        if len(parts) < 5:
            continue
        contig = parts[0]
        lineage = parts[4]
        ranks = []
        for rank in lineage.split(';'):
            if '__' in rank:
                rank = rank.split('__', 1)[1]
            if rank:
                ranks.append(rank)
        if ranks:
            tax[contig] = ';'.join(ranks)

# 2. Read ALL contig names from FASTA to ensure completeness
all_contigs = []
with open('${assembly}') as f:
    for line in f:
        if line.startswith('>'):
            all_contigs.append(line[1:].strip().split()[0])

# 3. Write taxonomy for all contigs (unclassified get 'Bacteria' placeholder)
n_mapped = 0
with open('vamb_taxonomy.tsv', 'w') as out:
    out.write('contigs\\tpredictions\\n')
    for contig in all_contigs:
        if contig in tax:
            out.write(contig + '\\t' + tax[contig] + '\\n')
            n_mapped += 1
        else:
            out.write(contig + '\\tBacteria\\n')

print('[INFO] VAMB taxonomy: %d/%d contigs classified (%d unclassified → Bacteria)' % (n_mapped, len(all_contigs), len(all_contigs) - n_mapped), file=sys.stderr)
"

    # Reuse same depth conversion as BIN_VAMB
    python3 -c "
import sys, csv

with open('${depths}') as f:
    reader = csv.reader(f, delimiter='\\t')
    header = next(reader)

depth_cols = []
for i, h in enumerate(header):
    if i >= 3 and not h.endswith('-var'):
        depth_cols.append(i)

rows = []
with open('${depths}') as f:
    reader = csv.reader(f, delimiter='\\t')
    next(reader)
    for row in reader:
        rows.append(row)

keep = []
for j, ci in enumerate(depth_cols):
    col_sum = sum(float(rows[r][ci]) for r in range(len(rows)))
    if col_sum > 0:
        keep.append(ci)

skipped = len(depth_cols) - len(keep)
print('[INFO] VAMB-tax abundance: %d contigs, %d samples (%d zero-depth dropped)' % (len(rows), len(keep), skipped), file=sys.stderr)

with open('vamb_abundance.tsv', 'w') as out:
    sample_names = [header[ci] for ci in keep]
    out.write('contigname\\t' + '\\t'.join(sample_names) + '\\n')
    for row in rows:
        out.write(row[0] + '\\t' + '\\t'.join(row[ci] for ci in keep) + '\\n')
"

    n_tax=\$(tail -n +2 vamb_taxonomy.tsv | wc -l)
    n_samples=\$(head -1 vamb_abundance.tsv | awk -F'\\t' '{print NF-1}')
    if [ "\$n_tax" -lt 1 ] || [ "\$n_samples" -lt 1 ]; then
        echo "[WARNING] Insufficient data for taxvamb (tax=\$n_tax, samples=\$n_samples)" >&2
        touch vamb_tax_bins.tsv
    else

    CUDA_FLAG=""
    if python3 -c "import torch; assert torch.cuda.is_available()" 2>/dev/null; then
        CUDA_FLAG="--cuda"
        echo "[INFO] VAMB-tax: GPU detected, using CUDA" >&2
    else
        echo "[INFO] VAMB-tax: No GPU detected, using CPU" >&2
    fi

    set +e
    vamb bin taxvamb \
        --outdir vamb_out \
        --fasta "${assembly}" \
        --abundance_tsv vamb_abundance.tsv \
        --taxonomy vamb_taxonomy.tsv \
        -m 2000 \
        --minfasta 200000 \
        \$CUDA_FLAG \
        -p ${task.cpus}
    vamb_exit=\$?
    set -e

    if [ \$vamb_exit -ne 0 ]; then
        echo "[WARNING] VAMB taxvamb exited with code \$vamb_exit" >&2
        touch vamb_tax_bins.tsv
    elif [ -d vamb_out/bins ] && ls vamb_out/bins/*.fna 1>/dev/null 2>&1; then
        > vamb_tax_bins.tsv
        bin_num=0
        for bin_file in vamb_out/bins/*.fna; do
            [ -e "\$bin_file" ] || continue
            bin_num=\$((bin_num + 1))
            bin_name=\$(printf 'vamb_tax_%03d' \$bin_num)
            cp "\$bin_file" "bins/\${bin_name}.fa"
            grep '>' "\$bin_file" | tr -d '>' | cut -f1 -d' ' | while read -r contig; do
                printf '%s\\t%s\\n' "\$contig" "\$bin_name"
            done >> vamb_tax_bins.tsv
        done
        echo "[INFO] VAMB taxvamb produced \$bin_num bins" >&2
    elif [ -f vamb_out/vae_clusters_unsplit.tsv ]; then
        tail -n +2 vamb_out/vae_clusters_unsplit.tsv | \
            awk -F'\\t' '{print \$2 "\\t" \$1}' > vamb_tax_bins.tsv
        echo "[INFO] VAMB taxvamb produced clusters but no FASTA bins above size threshold" >&2
    else
        echo "[WARNING] VAMB taxvamb produced no output" >&2
        touch vamb_tax_bins.tsv
    fi

    if [ ! -s vamb_tax_bins.tsv ]; then
        echo "[WARNING] VAMB taxvamb produced no bins" >&2
    fi

    fi  # end: data check
    """
}

process BINETTE_CONSENSUS {
    tag "binette"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-binette"
    publishDir "${params.outdir}/binning/binette", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/binette" : null

    input:
    path(assembly)
    path(bin_dirs, stageAs: 'binner_?')

    output:
    path("binette_bins.tsv"), emit: bins
    path("bins/"),            emit: fastas

    script:
    // Build space-separated list of staged bin directories
    def dirs_list = bin_dirs instanceof List ? bin_dirs.join(' ') : bin_dirs
    """
    mkdir -p bins

    # Binette consensus refinement — set operations + internal CheckM2 evaluation
    # Resolve the DIAMOND .dmnd file from the CheckM2 database directory
    CHECKM2_DMND=\$(find ${params.checkm2_db} -name '*.dmnd' -type f | head -1)
    if [ -z "\$CHECKM2_DMND" ]; then
        echo "[ERROR] No .dmnd file found in ${params.checkm2_db}" >&2
        exit 1
    fi
    echo "[INFO] Using CheckM2 DIAMOND db: \$CHECKM2_DMND" >&2

    set +e
    binette \\
        --bin_dirs ${dirs_list} \\
        --contigs ${assembly} \\
        --checkm2_db "\$CHECKM2_DMND" \\
        --threads ${task.cpus} \\
        -o binette_out
    binette_exit=\$?
    set -e

    if [ \$binette_exit -ne 0 ]; then
        echo "[WARNING] Binette exited with code \$binette_exit" >&2
        touch binette_bins.tsv
    elif [ -d binette_out/final_bins ]; then
        > binette_bins.tsv
        bin_num=0
        for bin_file in binette_out/final_bins/*.fa; do
            [ -e "\$bin_file" ] || continue
            bin_num=\$((bin_num + 1))
            bin_name=\$(printf 'binette_%03d' \$bin_num)
            cp "\$bin_file" "bins/\${bin_name}.fa"
            grep '>' "\$bin_file" | tr -d '>' | cut -f1 -d' ' | while read -r contig; do
                printf '%s\\t%s\\n' "\$contig" "\$bin_name"
            done >> binette_bins.tsv
        done
    else
        echo "[WARNING] Binette produced no final_bins directory" >&2
        touch binette_bins.tsv
    fi

    if [ ! -s binette_bins.tsv ]; then
        echo "[WARNING] Binette produced no bins" >&2
    fi
    """
}

process MAGSCOT_CONSENSUS {
    tag "magscot"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-magscot"
    publishDir "${params.outdir}/binning/magscot", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/magscot" : null

    input:
    path(assembly)
    path(bin_files)
    val(bin_labels)

    output:
    path("magscot_bins.tsv"), emit: bins
    path("bins/"),            emit: fastas

    script:
    // Build comma-separated file and label lists from collected inputs
    def files_csv = bin_files instanceof List ? bin_files.join(',') : bin_files
    def labels_csv = bin_labels instanceof List ? bin_labels.join(',') : bin_labels
    """
    mkdir -p bins

    # 1. Predict proteins with Prodigal (meta mode, parallelized)
    # MAGScoT requires prodigal-style gene IDs (contig_N) to map genes → contigs
    mkdir -p prodigal_chunks prodigal_out
    awk 'BEGIN{n=0; f="prodigal_chunks/chunk_000.fa"} /^>/{n++; if(n%1000==1){f=sprintf("prodigal_chunks/chunk_%03d.fa",int(n/1000))}} {print > f}' ${assembly}
    ls prodigal_chunks/*.fa | parallel -j ${task.cpus} \
        'prodigal -i {} -a prodigal_out/{/.}.faa -o /dev/null -p meta -f gff 2>/dev/null'
    cat prodigal_out/*.faa > proteins.faa
    rm -rf prodigal_chunks prodigal_out
    echo "[INFO] Prodigal predicted \$(grep -c "^>" proteins.faa) proteins" >&2

    # 2. HMM search against GTDBtk rel207 markers (shipped in MAGScoT repo)
    MAGSCOT_DIR=\$(dirname \$(readlink -f \$(which MAGScoT.R)))
    echo "[INFO] MAGScoT dir: \$MAGSCOT_DIR" >&2

    # TIGRFAM and Pfam HMMs use different tblout column layouts
    hmmsearch --tblout tigr.tblout --noali --notextw --cut_nc \
        --cpu ${task.cpus} \$MAGSCOT_DIR/hmm/gtdbtk_rel207_tigrfam.hmm proteins.faa > /dev/null
    hmmsearch --tblout pfam.tblout --noali --notextw --cut_nc \
        --cpu ${task.cpus} \$MAGSCOT_DIR/hmm/gtdbtk_rel207_Pfam-A.hmm proteins.faa > /dev/null

    # 3. Parse HMM results: gene_id, marker_name, e-value
    #    TIGRFAM: target=\$1, query=\$3, evalue=\$5
    #    Pfam:    target=\$1, query=\$4, evalue=\$5
    grep -v '^#' tigr.tblout | awk '{print \$1 "\\t" \$3 "\\t" \$5}' > markers.tsv
    grep -v '^#' pfam.tblout | awk '{print \$1 "\\t" \$4 "\\t" \$5}' >> markers.tsv
    echo "[INFO] Found \$(wc -l < markers.tsv) marker gene hits" >&2

    # 4. Build combined contig-to-bin TSV: bin \\t contig \\t binner
    > contig2bin.tsv
    IFS=',' read -ra F_ARR <<< "${files_csv}"
    IFS=',' read -ra L_ARR <<< "${labels_csv}"
    for i in "\${!F_ARR[@]}"; do
        if [ -s "\${F_ARR[\$i]}" ]; then
            # Input TSVs are contig\\tbin; MAGScoT expects bin\\tcontig\\tbinner
            awk -v lbl="\${L_ARR[\$i]}" '{print \$2 "\\t" \$1 "\\t" lbl}' "\${F_ARR[\$i]}" >> contig2bin.tsv
        fi
    done

    if [ ! -s contig2bin.tsv ]; then
        echo "[WARNING] All binners produced empty output -- skipping MAGScoT" >&2
        touch magscot_bins.tsv
    else
        # 5. Run MAGScoT
        set +e
        Rscript \$(which MAGScoT.R) -i contig2bin.tsv --hmm markers.tsv -o magscot_out
        magscot_exit=\$?
        set -e

        if [ \$magscot_exit -ne 0 ]; then
            echo "[WARNING] MAGScoT exited with code \$magscot_exit" >&2
            touch magscot_bins.tsv
        elif [ -f magscot_out.refined.contig_to_bin.out ]; then
            # MAGScoT output: header + bin\\tcontig
            # Skip header, extract bin\\tcontig pairs
            tail -n +2 magscot_out.refined.contig_to_bin.out > magscot_raw.tsv

            # Get unique bin names (column 1) and renumber sequentially
            > magscot_bins.tsv
            cut -f1 magscot_raw.tsv | sort -u > magscot_bin_list.txt
            bin_num=0
            for old_bin in \$(cat magscot_bin_list.txt); do
                bin_num=\$((bin_num + 1))
                bin_name=\$(printf 'magscot_%03d' \$bin_num)
                # Extract contigs for this bin (column 2 where column 1 matches)
                awk -v b="\$old_bin" '\$1 == b {print \$2}' magscot_raw.tsv > tmp_contigs.txt
                while read -r contig; do
                    printf '%s\\t%s\\n' "\$contig" "\$bin_name"
                done < tmp_contigs.txt >> magscot_bins.tsv
                # Extract bin FASTA from assembly
                awk 'BEGIN{while((getline<"tmp_contigs.txt")>0) keep[\$1]=1}
                     /^>/{p=keep[substr(\$1,2)]}p' ${assembly} > "bins/\${bin_name}.fa"
                rm -f tmp_contigs.txt
            done
            echo "[INFO] MAGScoT produced \$bin_num bins" >&2
            rm -f magscot_bin_list.txt
        else
            echo "[WARNING] MAGScoT produced no output" >&2
            touch magscot_bins.tsv
        fi
    fi

    if [ ! -s magscot_bins.tsv ]; then
        echo "[WARNING] MAGScoT produced no bins" >&2
    fi
    """
}

process DASTOOL_CONSENSUS {
    tag "dastool"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-binning"
    publishDir "${params.outdir}/binning/dastool", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/dastool" : null

    input:
    path(assembly)
    path(bin_files)
    val(bin_labels)

    output:
    path("bins/"),             emit: bins
    path("contig2bin.tsv"),    emit: contig2bin
    path("allbins.fa"),       emit: allbins
    path("bin_quality.tsv"),  emit: bin_quality
    path("summary.tsv"),     emit: summary
    path("bacteria.scg"),    emit: bacteria_scg
    path("archaea.scg"),     emit: archaea_scg

    script:
    // Build comma-separated file and label lists from collected inputs
    // bin_files is a list of paths; bin_labels is a list of strings
    def files_csv = bin_files instanceof List ? bin_files.join(',') : bin_files
    def labels_csv = bin_labels instanceof List ? bin_labels.join(',') : bin_labels
    """
    # Filter out empty bin files (binners that produced 0 bins)
    FILTERED_FILES=""
    FILTERED_LABELS=""
    IFS=',' read -ra F_ARR <<< "${files_csv}"
    IFS=',' read -ra L_ARR <<< "${labels_csv}"
    for i in "\${!F_ARR[@]}"; do
        if [ -s "\${F_ARR[\$i]}" ]; then
            [ -n "\$FILTERED_FILES" ] && FILTERED_FILES+=","
            FILTERED_FILES+="\${F_ARR[\$i]}"
            [ -n "\$FILTERED_LABELS" ] && FILTERED_LABELS+=","
            FILTERED_LABELS+="\${L_ARR[\$i]}"
        else
            echo "[WARNING] Skipping empty binner output: \${L_ARR[\$i]}" >&2
        fi
    done

    if [ -z "\$FILTERED_FILES" ]; then
        echo "[WARNING] All binners produced empty output -- skipping DAS_Tool" >&2
        mkdir -p bins
        touch contig2bin.tsv allbins.fa bin_quality.tsv summary.tsv bacteria.scg archaea.scg
        exit 0
    fi

    # DAS_Tool exits 1 when no bins pass score_threshold (0.5 default)
    # This is expected with small datasets -- catch and produce empty output
    set +e
    DAS_Tool \\
        -i "\$FILTERED_FILES" \\
        -l "\$FILTERED_LABELS" \\
        -c "${assembly}" \\
        -o dastool_out/dastool \\
        --threads ${task.cpus} \\
        --write_bin_evals \\
        --write_bins
    dastool_exit=\$?
    set -e

    # Rename DAS_Tool verbose output to clean names
    mkdir -p bins

    if [ \$dastool_exit -ne 0 ]; then
        echo "[WARNING] DAS_Tool exited with code \$dastool_exit (no bins above score threshold)" >&2
    fi

    # ---- Collect DAS_Tool evaluation and summary files ----
    if [ -f dastool_out/dastool_allBins.eval ]; then
        cp dastool_out/dastool_allBins.eval bin_quality.tsv
    else
        touch bin_quality.tsv
    fi

    if [ -f dastool_out/dastool_DASTool_summary.tsv ]; then
        cp dastool_out/dastool_DASTool_summary.tsv summary.tsv
    else
        touch summary.tsv
    fi

    if [ -f dastool_out/dastool_DASTool_contig2bin.tsv ]; then
        cp dastool_out/dastool_DASTool_contig2bin.tsv contig2bin.tsv
    else
        touch contig2bin.tsv
    fi

    if [ -d dastool_out/dastool_DASTool_bins ]; then
        for f in dastool_out/dastool_DASTool_bins/*.fa; do
            [ -e "\$f" ] || continue
            cp "\$f" "bins/dastool-\$(basename "\$f")"
        done
    fi

    # Prefix bin names in contig2bin.tsv and summary.tsv so they match the
    # dastool- prefixed FASTA filenames
    if [ -s contig2bin.tsv ]; then
        sed -i 's/\\t/\\tdastool-/' contig2bin.tsv
    fi
    if [ -s summary.tsv ]; then
        # Skip header line, prefix bin name in first column
        sed -i '2,\$s/^/dastool-/' summary.tsv
    fi

    # Combine all bin FASTAs
    if ls bins/*.fa 1>/dev/null 2>&1; then
        cat bins/*.fa > allbins.fa
    else
        echo "[WARNING] DAS_Tool produced no consensus bins" >&2
        touch allbins.fa
    fi

    # Publish single-copy gene (SCG) assignments from DAS_Tool
    if [ -f dastool_out/dastool_proteins.faa.bacteria.scg ]; then
        cp dastool_out/dastool_proteins.faa.bacteria.scg bacteria.scg
    else
        touch bacteria.scg
    fi
    if [ -f dastool_out/dastool_proteins.faa.archaea.scg ]; then
        cp dastool_out/dastool_proteins.faa.archaea.scg archaea.scg
    else
        touch archaea.scg
    fi
    """
}

process CHECKM2 {
    tag "checkm2"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-mag-quality"
    publishDir "${params.outdir}/binning/checkm2", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/binning/checkm2" : null

    input:
    path(bins_dirs, stageAs: 'bins_?')   // collected list of bins/ directories

    output:
    path("quality_report.tsv"), emit: report

    script:
    """
    # Merge all bin directories into one flat directory
    # stageAs: 'bins_?' resolves the name collision (bins_1, bins_2, ...)
    mkdir -p all_bins
    for d in bins_*; do
        [ -d "\$d" ] || continue
        cp "\$d"/*.fa all_bins/ 2>/dev/null || true
    done

    if [ -z "\$(ls all_bins/*.fa 2>/dev/null)" ]; then
        echo "[WARNING] No bin FASTAs found -- skipping CheckM2" >&2
        printf 'Name\\tCompleteness\\tContamination\\n' > quality_report.tsv
        exit 0
    fi

    N_BINS=\$(ls all_bins/*.fa | wc -l)
    echo "[INFO] CheckM2: \$N_BINS bins" >&2

    # Resolve database path
    DB_PATH="${params.checkm2_db}"
    if [ -d "\$DB_PATH" ]; then
        DMND_FILE=\$(find "\$DB_PATH" -name '*.dmnd' -type f | head -1)
        if [ -z "\$DMND_FILE" ]; then
            echo "[ERROR] No .dmnd file found in \$DB_PATH" >&2
            exit 1
        fi
        DB_PATH="\$DMND_FILE"
    fi

    # Use /tmp for TMPDIR — CheckM2's Python multiprocessing creates AF_UNIX
    # sockets under TMPDIR, and Nextflow work dir paths can exceed the 108-char
    # socket path limit, causing "AF_UNIX path too long" errors.
    export TMPDIR="/tmp/checkm2_\$\$"
    mkdir -p "\$TMPDIR"
    trap 'rm -rf "\$TMPDIR"' EXIT

    # Shim prodigal → pyrodigal to avoid free(): invalid pointer crash.
    # Prodigal 2.6.3 (C) has a known memory-safety bug that triggers under
    # CheckM2's Python multiprocessing. Pyrodigal (Cython rewrite) is safe.
    # CheckM2 calls 'prodigal -p meta -q -m -f gff ...' but pyrodigal doesn't
    # accept -q (quiet), so the wrapper strips it.
    mkdir -p shims
    cat > shims/prodigal <<'SHIM'
#!/bin/bash
args=()
for a in "\$@"; do [ "\$a" != "-q" ] && args+=("\$a"); done
exec pyrodigal "\${args[@]}"
SHIM
    chmod +x shims/prodigal
    export PATH="\$PWD/shims:\$PATH"

    checkm2 predict \\
        --threads ${task.cpus} \\
        --input all_bins \\
        --output-directory checkm2_out \\
        -x fa \\
        --database_path "\$DB_PATH"

    cp checkm2_out/quality_report.tsv .
    """
}
