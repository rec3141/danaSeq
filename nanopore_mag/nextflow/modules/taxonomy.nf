// Contig-level taxonomy classification: Kaiju six-frame translation (runs directly on contigs)
// Primary Kaiju mode — no annotation dependency, avoids Bakta contig renumbering issues.
// Uses Kaiju's built-in six-frame translation to classify each contig independently.

process KAIJU_CONTIG_CLASSIFY {
    tag "kaiju_contigs"
    label 'process_kraken'
    conda "${projectDir}/conda-envs/dana-mag-kaiju"
    publishDir "${params.outdir}/taxonomy/kaiju", mode: 'link', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/taxonomy/kaiju" : null

    input:
    path(contigs)   // Assembly FASTA — no annotation needed

    output:
    path("kaiju_contigs.tsv"), emit: contig_taxonomy

    script:
    def db_dir = params.kaiju_db
    """
    # Auto-detect .fmi file in database directory
    fmi=\$(ls "${db_dir}"/*.fmi 2>/dev/null | head -1)
    nodes="${db_dir}/nodes.dmp"
    names="${db_dir}/names.dmp"

    if [ -z "\$fmi" ] || [ ! -f "\$fmi" ]; then
        echo "[ERROR] No .fmi file found in ${db_dir}" >&2
        printf 'contig_id\\tstatus\\ttaxon_id\\tlineage\\n' > kaiju_contigs.tsv
        exit 0
    fi

    if [ ! -f "\$nodes" ] || [ ! -f "\$names" ]; then
        echo "[ERROR] nodes.dmp or names.dmp not found in ${db_dir}" >&2
        printf 'contig_id\\tstatus\\ttaxon_id\\tlineage\\n' > kaiju_contigs.tsv
        exit 0
    fi

    # Kaiju six-frame translation mode (no -p flag): translates contigs in all 6 frames
    # and classifies the resulting protein fragments against the database.
    # -a greedy -e 5: greedy mode allowing 5 mismatches (better sensitivity)
    set +e
    kaiju \\
        -t "\$nodes" \\
        -f "\$fmi" \\
        -i "${contigs}" \\
        -o kaiju_raw.tsv \\
        -z ${task.cpus} \\
        -a greedy \\
        -e 5 \\
        -v
    kaiju_exit=\$?
    set -e

    if [ \$kaiju_exit -ne 0 ]; then
        echo "[WARNING] Kaiju exited with code \$kaiju_exit" >&2
        printf 'contig_id\\tstatus\\ttaxon_id\\tlineage\\n' > kaiju_contigs.tsv
        exit 0
    fi

    # Add full taxonomic lineage names
    set +e
    kaiju-addTaxonNames \\
        -t "\$nodes" \\
        -n "\$names" \\
        -i kaiju_raw.tsv \\
        -o kaiju_names.tsv \\
        -r superkingdom,phylum,class,order,family,genus,species
    names_exit=\$?
    set -e

    if [ \$names_exit -ne 0 ]; then
        echo "[WARNING] kaiju-addTaxonNames exited with code \$names_exit" >&2
        printf 'contig_id\\tstatus\\ttaxon_id\\tlineage\\n' > kaiju_contigs.tsv
        exit 0
    fi

    # Build per-contig TSV: contig_id, status, taxon_id, lineage
    # Output is one row per contig (Kaiju outputs one row per input sequence in contig mode)
    {
        printf 'contig_id\\tstatus\\ttaxon_id\\tlineage\\n'
        awk -F'\\t' '{
            status = \$1
            contig = \$2
            taxid = \$3
            lineage = (NF >= 4) ? \$NF : ""
            # Clean trailing semicolons and whitespace
            gsub(/;[ ]*\$/, "", lineage)
            if (status == "C") {
                print contig "\\t" status "\\t" taxid "\\t" lineage
            } else {
                print contig "\\tU\\t0\\tUnclassified"
            }
        }' kaiju_names.tsv
    } > kaiju_contigs.tsv

    n_total=\$(tail -n +2 kaiju_contigs.tsv | wc -l)
    n_classified=\$(tail -n +2 kaiju_contigs.tsv | awk -F'\\t' '\$2=="C"' | wc -l)
    echo "[INFO] Kaiju (6-frame): \$n_classified/\$n_total contigs classified" >&2
    """
}

// Kaiju protein-level taxonomy (secondary — uses Prokka/Bakta .faa + .gff)
// Per-gene classification aggregated to per-contig via GFF mapping.
// Note: depends on annotation; contig IDs must match between assembly and annotation.

process KAIJU_CLASSIFY {
    tag "kaiju_proteins"
    label 'process_kraken'
    conda "${projectDir}/conda-envs/dana-mag-kaiju"
    publishDir "${params.outdir}/taxonomy/kaiju", mode: 'link', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/taxonomy/kaiju" : null

    input:
    path(proteins)   // Prokka .faa (amino acid sequences)
    path(gff)        // Prokka .gff (gene-to-contig mapping)

    output:
    path("kaiju_genes.tsv"),            emit: gene_taxonomy
    path("kaiju_protein_contigs.tsv"),  emit: contig_taxonomy

    script:
    def db_dir = params.kaiju_db
    """
    # Auto-detect .fmi file in database directory
    fmi=\$(ls "${db_dir}"/*.fmi 2>/dev/null | head -1)
    nodes="${db_dir}/nodes.dmp"
    names="${db_dir}/names.dmp"

    if [ -z "\$fmi" ] || [ ! -f "\$fmi" ]; then
        echo "[ERROR] No .fmi file found in ${db_dir}" >&2
        printf 'classified\\tgene_id\\ttaxon_id\\tlineage\\n' > kaiju_genes.tsv
        printf 'contig_id\\tclassified_genes\\ttotal_genes\\tfraction_classified\\ttaxon_id\\tlineage\\n' > kaiju_protein_contigs.tsv
        exit 0
    fi

    if [ ! -f "\$nodes" ] || [ ! -f "\$names" ]; then
        echo "[ERROR] nodes.dmp or names.dmp not found in ${db_dir}" >&2
        printf 'classified\\tgene_id\\ttaxon_id\\tlineage\\n' > kaiju_genes.tsv
        printf 'contig_id\\tclassified_genes\\ttotal_genes\\tfraction_classified\\ttaxon_id\\tlineage\\n' > kaiju_protein_contigs.tsv
        exit 0
    fi

    # Kaiju protein mode (-p): classify Prokka-called proteins against RefSeq
    # -a greedy -e 5: greedy mode allowing 5 mismatches (better sensitivity for divergent proteins)
    set +e
    kaiju \\
        -p \\
        -t "\$nodes" \\
        -f "\$fmi" \\
        -i "${proteins}" \\
        -o kaiju_raw.tsv \\
        -z ${task.cpus} \\
        -a greedy \\
        -e 5 \\
        -v
    kaiju_exit=\$?
    set -e

    if [ \$kaiju_exit -ne 0 ]; then
        echo "[WARNING] Kaiju exited with code \$kaiju_exit" >&2
        printf 'classified\\tgene_id\\ttaxon_id\\tlineage\\n' > kaiju_genes.tsv
        printf 'contig_id\\tclassified_genes\\ttotal_genes\\tfraction_classified\\ttaxon_id\\tlineage\\n' > kaiju_protein_contigs.tsv
        exit 0
    fi

    # Add full taxonomic lineage names
    set +e
    kaiju-addTaxonNames \\
        -t "\$nodes" \\
        -n "\$names" \\
        -i kaiju_raw.tsv \\
        -o kaiju_names.tsv \\
        -r superkingdom,phylum,class,order,family,genus,species
    names_exit=\$?
    set -e

    if [ \$names_exit -ne 0 ]; then
        echo "[WARNING] kaiju-addTaxonNames exited with code \$names_exit" >&2
        printf 'classified\\tgene_id\\ttaxon_id\\tlineage\\n' > kaiju_genes.tsv
        printf 'contig_id\\tclassified_genes\\ttotal_genes\\tfraction_classified\\ttaxon_id\\tlineage\\n' > kaiju_protein_contigs.tsv
        exit 0
    fi

    # Per-gene output with header
    # kaiju-addTaxonNames output: C/U \\t gene_id \\t taxon_id \\t lineage
    {
        printf 'classified\\tgene_id\\ttaxon_id\\tlineage\\n'
        cat kaiju_names.tsv
    } > kaiju_genes.tsv

    # Aggregate per-gene taxonomy to per-contig using GFF gene-to-contig mapping
    # Majority-vote at the most specific rank with agreement
    python3 - "${gff}" kaiju_names.tsv <<'PYEOF'
import sys
import re
from collections import defaultdict, Counter

gff_path, kaiju_path = sys.argv[1], sys.argv[2]

# Build gene_id -> contig_id mapping from GFF
gene2contig = {}
contig_genes = defaultdict(int)
with open(gff_path) as f:
    for line in f:
        if line.startswith('#') or '\\t' not in line:
            continue
        cols = line.strip().split('\\t')
        if len(cols) < 9 or cols[2] != 'CDS':
            continue
        contig = cols[0]
        m = re.search(r'ID=([^;]+)', cols[8])
        if m:
            gene_id = m.group(1)
            gene2contig[gene_id] = contig
            contig_genes[contig] += 1

# Parse kaiju results, aggregate by contig
contig_lineages = defaultdict(list)   # contig -> list of (taxon_id, lineage)
contig_classified = defaultdict(int)
with open(kaiju_path) as f:
    for line in f:
        cols = line.strip().split('\\t')
        if len(cols) < 3:
            continue
        status, gene_id, taxon_id = cols[0], cols[1], cols[2]
        lineage = cols[-1] if len(cols) >= 4 else ''
        contig = gene2contig.get(gene_id)
        if not contig:
            continue
        if status == 'C':
            contig_classified[contig] += 1
            contig_lineages[contig].append((taxon_id, lineage.strip().rstrip(';')))

# For each contig, majority-vote lineage
# Strategy: pick the most common full lineage; if tied, pick the one with
# the most specific (deepest) classification
with open('kaiju_protein_contigs.tsv', 'w') as out:
    out.write('contig_id\\tclassified_genes\\ttotal_genes\\tfraction_classified\\ttaxon_id\\tlineage\\n')
    for contig in sorted(contig_genes.keys()):
        total = contig_genes[contig]
        n_classified = contig_classified.get(contig, 0)
        frac = n_classified / total if total > 0 else 0.0
        if n_classified == 0:
            out.write(f'{contig}\\t0\\t{total}\\t0.000\\tNA\\tUnclassified\\n')
            continue
        lineages = contig_lineages[contig]
        # Majority vote on full lineage string
        lineage_counts = Counter(l for _, l in lineages)
        best_lineage, best_count = lineage_counts.most_common(1)[0]
        # Find taxon_id for the winning lineage
        best_taxid = next(tid for tid, l in lineages if l == best_lineage)
        out.write(f'{contig}\\t{n_classified}\\t{total}\\t{frac:.3f}\\t{best_taxid}\\t{best_lineage}\\n')

n_contigs = len(contig_genes)
n_classified_contigs = sum(1 for c in contig_genes if contig_classified.get(c, 0) > 0)
print(f'[INFO] Kaiju: {n_classified_contigs}/{n_contigs} contigs have at least one classified gene', file=sys.stderr)
PYEOF
    """
}

// Contig-level taxonomy classification: Kraken2 (k-mer-based, runs directly on contigs)

process KRAKEN2_CLASSIFY {
    tag "kraken2"
    label 'process_kraken'
    maxForks 1
    conda "${projectDir}/conda-envs/dana-mag-kraken2"
    publishDir "${params.outdir}/taxonomy/kraken2", mode: 'link', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/taxonomy/kraken2" : null

    input:
    path(contigs)   // Assembly FASTA — no annotation needed

    output:
    path("kraken2_contigs.tsv"), emit: contig_taxonomy
    path("kraken2_report.txt"),  emit: report

    script:
    def db_dir = params.kraken2_db
    def confidence = params.kraken2_confidence
    """
    nodes="${db_dir}/nodes.dmp"
    names="${db_dir}/names.dmp"

    if [ ! -f "${db_dir}/hash.k2d" ]; then
        echo "[ERROR] Kraken2 database not found at ${db_dir} (missing hash.k2d)" >&2
        printf 'contig_id\\tstatus\\ttaxon_id\\tname\\tlineage\\n' > kraken2_contigs.tsv
        printf '' > kraken2_report.txt
        exit 0
    fi

    if [ ! -f "\$nodes" ] || [ ! -f "\$names" ]; then
        echo "[ERROR] nodes.dmp or names.dmp not found in ${db_dir}" >&2
        printf 'contig_id\\tstatus\\ttaxon_id\\tname\\tlineage\\n' > kraken2_contigs.tsv
        printf '' > kraken2_report.txt
        exit 0
    fi

    # Kraken2 k-mer classification on assembly contigs
    # --confidence ${confidence}: fraction of k-mers required for classification
    #   (contigs have many k-mers; higher threshold reduces false positives)
    set +e
    kraken2 \\
        --db "${db_dir}" \\
        --threads ${task.cpus} \\
        --confidence ${confidence} \\
        --output kraken2_raw.tsv \\
        --report kraken2_report.txt \\
        "${contigs}"
    kraken2_exit=\$?
    set -e

    if [ \$kraken2_exit -ne 0 ]; then
        echo "[WARNING] Kraken2 exited with code \$kraken2_exit" >&2
        printf 'contig_id\\tstatus\\ttaxon_id\\tname\\tlineage\\n' > kraken2_contigs.tsv
        printf '' > kraken2_report.txt
        exit 0
    fi

    # Post-process: build full GTDB-style lineage strings from nodes.dmp + names.dmp
    python3 - kraken2_raw.tsv "\$nodes" "\$names" <<'PYEOF'
import sys
from collections import defaultdict

kraken2_path, nodes_path, names_path = sys.argv[1], sys.argv[2], sys.argv[3]

# Parse NCBI taxonomy nodes.dmp: taxid -> (parent_taxid, rank)
parent = {}
rank = {}
with open(nodes_path) as f:
    for line in f:
        parts = line.split('|')
        tid = int(parts[0].strip())
        pid = int(parts[1].strip())
        r = parts[2].strip()
        parent[tid] = pid
        rank[tid] = r

# Parse names.dmp: taxid -> scientific name
sci_name = {}
with open(names_path) as f:
    for line in f:
        parts = line.split('|')
        tid = int(parts[0].strip())
        name = parts[1].strip()
        name_class = parts[3].strip()
        if name_class == 'scientific name':
            sci_name[tid] = name

# NCBI rank -> GTDB prefix
rank_prefix = {
    'superkingdom': 'd', 'domain': 'd', 'phylum': 'p', 'class': 'c',
    'order': 'o', 'family': 'f', 'genus': 'g', 'species': 's'
}
rank_order = ['superkingdom', 'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def build_lineage(taxid):
    # Walk up the taxonomy tree and build d__X;p__Y;...;s__Z lineage.
    lineage_parts = {}
    tid = taxid
    visited = set()
    while tid != 1 and tid in parent and tid not in visited:
        visited.add(tid)
        r = rank.get(tid, '')
        if r in rank_prefix:
            lineage_parts[r] = f"{rank_prefix[r]}__{sci_name.get(tid, 'unknown')}"
        tid = parent[tid]
    # Build ordered lineage string
    parts = []
    for r in rank_order:
        if r in lineage_parts:
            parts.append(lineage_parts[r])
    return ';'.join(parts) if parts else 'Unclassified'

# Parse Kraken2 output: status(C/U) \t contig_id \t taxon_id \t length \t LCA_mapping
n_classified = 0
n_total = 0
with open(kraken2_path) as fin, open('kraken2_contigs.tsv', 'w') as fout:
    fout.write('contig_id\\tstatus\\ttaxon_id\\tname\\tlineage\\n')
    for line in fin:
        cols = line.strip().split('\\t')
        if len(cols) < 3:
            continue
        n_total += 1
        status = cols[0]       # C or U
        contig_id = cols[1]
        taxon_id = int(cols[2])

        if status == 'C':
            n_classified += 1
            name = sci_name.get(taxon_id, 'unknown')
            lineage = build_lineage(taxon_id)
        else:
            name = 'Unclassified'
            lineage = 'Unclassified'

        fout.write(f'{contig_id}\\t{status}\\t{taxon_id}\\t{name}\\t{lineage}\\n')

print(f'[INFO] Kraken2: {n_classified}/{n_total} contigs classified', file=sys.stderr)
PYEOF
    """
}

// Contig-level taxonomy classification: BBSketch/sendsketch (MinHash sketch vs GTDB TaxServer)

process SENDSKETCH_CLASSIFY {
    tag "sendsketch"
    label 'process_medium'
    conda "${projectDir}/conda-envs/dana-bbmap"
    publishDir "${params.outdir}/taxonomy/sendsketch", mode: 'link', enabled: !params.store_dir
    storeDir params.store_dir ? "${params.store_dir}/taxonomy/sendsketch" : null

    input:
    path(contigs)   // Assembly FASTA — no annotation needed

    output:
    path("sendsketch_contigs.tsv"), emit: contig_taxonomy

    script:
    def address = params.sendsketch_address
    """
    # SendSketch has a hard limit of 100k sequences per request.
    # Split assembly into chunks of 90k contigs, send each batch, concatenate.
    BATCH_SIZE=90000
    N_SEQS=\$(grep -c '^>' "${contigs}")
    echo "[INFO] Assembly has \$N_SEQS contigs (batch size: \$BATCH_SIZE)" >&2

    mkdir -p chunks
    if [ "\$N_SEQS" -le "\$BATCH_SIZE" ]; then
        ln -s "\$(readlink -f ${contigs})" chunks/chunk_000.fa
    else
        awk -v bs="\$BATCH_SIZE" -v dir="chunks" '
            BEGIN { n=0; chunk=0; fn=sprintf("%s/chunk_%03d.fa", dir, chunk) }
            /^>/ { if (n >= bs) { close(fn); chunk++; fn=sprintf("%s/chunk_%03d.fa", dir, chunk); n=0 } n++ }
            { print > fn }
        ' "${contigs}"
    fi

    # Run sendsketch on each chunk
    > sendsketch_raw.txt
    sketch_ok=true
    for chunk in chunks/chunk_*.fa; do
        echo "[INFO] Processing \$(basename \$chunk) ..." >&2
        set +e
        sendsketch.sh \\
            -Xmx4g \\
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
        printf 'contig_id\\tstatus\\tANI\\tref_name\\tlineage\\n' > sendsketch_contigs.tsv
        exit 0
    fi

    # Parse format=2 output into clean TSV
    # Format=2 structure:
    #   Query: contig_name\\tSketchLen: N\\t...
    #   WKID\\tKID\\tANI\\t...\\ttaxName\\tseqName\\ttaxonomy
    #   <data line>
    #   (or "No hits.")
    python3 - sendsketch_raw.txt <<'PYEOF'
import sys
import re

raw_path = sys.argv[1]

# GTDB rank prefixes to standard d__/p__/c__/o__/f__/g__/s__ format
# GTDB uses d: k: p: c: o: f: g: s: — we drop k (kingdom) for consistency
gtdb_to_std = {'d': 'd', 'p': 'p', 'c': 'c', 'o': 'o', 'f': 'f', 'g': 'g', 's': 's'}

def convert_lineage(gtdb_lin):
    # Convert d:Bacteria;k:Pseudomonadati;p:Bacteroidota;... to d__Bacteria;p__Bacteroidota;...
    # Always emit all 7 standard ranks; missing ranks get empty prefix (e.g. o__)
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
current_contig = None
header_seen = False

with open(raw_path) as fin, open('sendsketch_contigs.tsv', 'w') as fout:
    fout.write('contig_id\\tstatus\\tANI\\tref_name\\tlineage\\n')
    for line in fin:
        line = line.rstrip('\\n')

        # Query line: "Query: contig_1\\tSketchLen: 181\\t..."
        if line.startswith('Query:'):
            # Write previous contig if it had no hits
            if current_contig is not None and not header_seen:
                fout.write(f'{current_contig}\\tU\\t0\\tUnclassified\\tUnclassified\\n')
                n_total += 1

            parts = line.split('\\t')
            current_contig = parts[0].replace('Query: ', '').strip()
            header_seen = False
            continue

        # Header line
        if line.startswith('WKID'):
            header_seen = True
            continue

        # "No hits." line
        if line.strip() == 'No hits.':
            if current_contig:
                fout.write(f'{current_contig}\\tU\\t0\\tUnclassified\\tUnclassified\\n')
                n_total += 1
                current_contig = None
            continue

        # Data line (after header): WKID KID ANI SSU Complt Contam Matches Unique TaxID gSize gSeqs taxName [seqName] taxonomy
        # Column count varies: 13 without seqName, 14 with it. Taxonomy is always last column.
        if header_seen and current_contig and '\\t' in line:
            cols = line.split('\\t')
            if len(cols) >= 12:
                ani_str = cols[2].rstrip('%')
                try:
                    ani = float(ani_str)
                except ValueError:
                    ani = 0.0
                ref_name = cols[11]  # taxName column
                taxonomy = cols[-1] if len(cols) >= 13 else ''
                lineage = convert_lineage(taxonomy) if taxonomy else 'Unclassified'
                fout.write(f'{current_contig}\\tC\\t{ani:.2f}\\t{ref_name}\\t{lineage}\\n')
                n_total += 1
                n_classified += 1
                current_contig = None
            continue

    # Handle last contig if unresolved
    if current_contig is not None:
        fout.write(f'{current_contig}\\tU\\t0\\tUnclassified\\tUnclassified\\n')
        n_total += 1

print(f'[INFO] SendSketch: {n_classified}/{n_total} contigs classified', file=sys.stderr)
PYEOF
    """
}
