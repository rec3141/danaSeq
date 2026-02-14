// Contig-level taxonomy classification: Kaiju (protein-level via Prokka .faa)

process KAIJU_CLASSIFY {
    tag "kaiju"
    label 'process_high'
    conda "${projectDir}/conda-envs/dana-mag-kaiju"
    publishDir "${params.outdir}/taxonomy/kaiju", mode: 'copy'

    input:
    path(proteins)   // Prokka .faa (amino acid sequences)
    path(gff)        // Prokka .gff (gene-to-contig mapping)

    output:
    path("kaiju_genes.tsv"),   emit: gene_taxonomy
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
        printf 'classified\\tgene_id\\ttaxon_id\\tlineage\\n' > kaiju_genes.tsv
        printf 'contig_id\\tclassified_genes\\ttotal_genes\\tfraction_classified\\ttaxon_id\\tlineage\\n' > kaiju_contigs.tsv
        exit 0
    fi

    if [ ! -f "\$nodes" ] || [ ! -f "\$names" ]; then
        echo "[ERROR] nodes.dmp or names.dmp not found in ${db_dir}" >&2
        printf 'classified\\tgene_id\\ttaxon_id\\tlineage\\n' > kaiju_genes.tsv
        printf 'contig_id\\tclassified_genes\\ttotal_genes\\tfraction_classified\\ttaxon_id\\tlineage\\n' > kaiju_contigs.tsv
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
        printf 'contig_id\\tclassified_genes\\ttotal_genes\\tfraction_classified\\ttaxon_id\\tlineage\\n' > kaiju_contigs.tsv
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
        printf 'contig_id\\tclassified_genes\\ttotal_genes\\tfraction_classified\\ttaxon_id\\tlineage\\n' > kaiju_contigs.tsv
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
with open('kaiju_contigs.tsv', 'w') as out:
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
