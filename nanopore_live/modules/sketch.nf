// Sendsketch per-read GTDB classification against a local sendsketch server.
//
// Ports the SENDSKETCH_CLASSIFY approach from mag_analysis/modules/taxonomy.nf:
// batches the input FASTA into chunks of up to 90k sequences (sendsketch
// per-request limit), sends each batch with persequence records=1, and parses
// the format=2 output into a per-read TSV.
//
// Output: sendsketch_reads.tsv (read_id, status, ANI, ref_name, lineage).
// Lineage is normalized from GTDB's d:/p:/c:/o:/f:/g:/s: tokens to the
// standard d__/p__/.../s__ form; all seven ranks are always emitted.

process SENDSKETCH {
    tag "${meta.id}"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-bbmap"
    publishDir "${params.outdir}/${meta.flowcell}/${meta.barcode}/sketch", mode: 'copy', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/${meta.flowcell}/${meta.barcode}/sketch" : null

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}.sendsketch_reads.tsv"), emit: sketch

    script:
    def address = params.sendsketch_address
    def mem_mb = task.memory ? (task.memory.toMega() * 85 / 100).intValue() : 3400
    """
    # SendSketch has a hard limit of 100k sequences per request.
    # Split input FASTA into chunks of 90k reads, send each batch, concatenate.
    BATCH_SIZE=90000
    N_SEQS=\$(grep -c '^>' "${fasta}")
    echo "[INFO] ${meta.id}: \$N_SEQS reads (batch size: \$BATCH_SIZE)" >&2

    mkdir -p chunks
    if [ "\$N_SEQS" -le "\$BATCH_SIZE" ]; then
        ln -s "\$(readlink -f ${fasta})" chunks/chunk_000.fa
    else
        awk -v bs="\$BATCH_SIZE" -v dir="chunks" '
            BEGIN { n=0; chunk=0; fn=sprintf("%s/chunk_%03d.fa", dir, chunk) }
            /^>/ { if (n >= bs) { close(fn); chunk++; fn=sprintf("%s/chunk_%03d.fa", dir, chunk); n=0 } n++ }
            { print > fn }
        ' "${fasta}"
    fi

    > sendsketch_raw.txt
    sketch_ok=true
    for chunk in chunks/chunk_*.fa; do
        echo "[INFO] Processing \$(basename \$chunk) ..." >&2
        set +e
        sendsketch.sh \\
            -Xmx${mem_mb}m \\
            in="\$chunk" \\
            address="${address}" \\
            k=31 \\
            format=2 \\
            persequence \\
            records=1 \\
            color=f \\
            printtaxa=t \\
            out=sendsketch_chunk.txt \\
            2>sendsketch_stderr.txt
        chunk_exit=\$?
        set -e

        if [ \$chunk_exit -ne 0 ]; then
            echo "[WARNING] sendsketch.sh exited with code \$chunk_exit for \$(basename \$chunk)" >&2
            cat sendsketch_stderr.txt >&2
            sketch_ok=false
            break
        fi
        cat sendsketch_chunk.txt >> sendsketch_raw.txt
    done

    if [ "\$sketch_ok" != "true" ]; then
        printf 'read_id\\tstatus\\tANI\\tref_name\\tlineage\\n' > "${meta.id}.sendsketch_reads.tsv"
        exit 0
    fi

    python3 - sendsketch_raw.txt "${meta.id}.sendsketch_reads.tsv" <<'PYEOF'
import sys

raw_path, out_path = sys.argv[1], sys.argv[2]

# GTDB lineage tokens are d:Foo;k:Bar;p:... — we drop k (kingdom) and emit
# the standard d__/p__/c__/o__/f__/g__/s__ order, always all seven ranks.
gtdb_to_std = {'d': 'd', 'p': 'p', 'c': 'c', 'o': 'o', 'f': 'f', 'g': 'g', 's': 's'}

def convert_lineage(gtdb_lin):
    rank_map = {}
    for token in gtdb_lin.split(';'):
        token = token.strip()
        if ':' not in token:
            continue
        prefix, name = token.split(':', 1)
        if prefix in gtdb_to_std:
            rank_map[gtdb_to_std[prefix]] = name
    if not rank_map:
        return 'Unclassified'
    all_ranks = ['d', 'p', 'c', 'o', 'f', 'g', 's']
    return ';'.join(f'{r}__{rank_map.get(r, "")}' for r in all_ranks)

n_total = 0
n_classified = 0
current_read = None
header_seen = False

with open(raw_path) as fin, open(out_path, 'w') as fout:
    fout.write('read_id\tstatus\tANI\tref_name\tlineage\n')
    for line in fin:
        line = line.rstrip('\n')

        if line.startswith('Query:'):
            # previous read had no hit data line
            if current_read is not None and not header_seen:
                fout.write(f'{current_read}\tU\t0\tUnclassified\tUnclassified\n')
                n_total += 1

            parts = line.split('\t')
            current_read = parts[0].replace('Query: ', '').strip()
            header_seen = False
            continue

        if line.startswith('WKID'):
            header_seen = True
            continue

        if line.strip() == 'No hits.':
            if current_read:
                fout.write(f'{current_read}\tU\t0\tUnclassified\tUnclassified\n')
                n_total += 1
                current_read = None
            continue

        # Data line (format=2): WKID KID ANI SSU Complt Contam Matches Unique TaxID gSize gSeqs taxName [seqName] taxonomy
        # 13 cols without seqName, 14 with. Taxonomy always last.
        if header_seen and current_read and '\t' in line:
            cols = line.split('\t')
            if len(cols) >= 12:
                ani_str = cols[2].rstrip('%')
                try:
                    ani = float(ani_str)
                except ValueError:
                    ani = 0.0
                ref_name = cols[11]
                taxonomy = cols[-1] if len(cols) >= 13 else ''
                lineage = convert_lineage(taxonomy) if taxonomy else 'Unclassified'
                fout.write(f'{current_read}\tC\t{ani:.2f}\t{ref_name}\t{lineage}\n')
                n_total += 1
                n_classified += 1
                current_read = None
            continue

    if current_read is not None:
        fout.write(f'{current_read}\tU\t0\tUnclassified\tUnclassified\n')
        n_total += 1

print(f'[INFO] SendSketch: {n_classified}/{n_total} reads classified', file=sys.stderr)
PYEOF
    """
}
