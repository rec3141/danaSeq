#!/usr/bin/env python3
"""MetaCyc pathway predictions per MAG and for the community, without running Pathway Tools.

Proteins get MetaCyc reactions through Rhea (DIAMOND hits to Rhea-annotated UniProtKB
proteins, Bakta/eggNOG ECs, KOs) and pathways are called with PathoLogic's published rules
(Karp et al. 2011, doi:10.4056/sigs.1794338), each call with an evidence-based confidence.
Standard library only; the package is bin/metacyc_scorer/.

  metacyc_scorer.py run --bakta-tsv annotation.tsv --merged merged_annotations.tsv \\
      --hits metacyc_hits.tsv.gz --contig2bin contig2bin.tsv \\
      --metacyc METACYC_DIR --refs REFS_DIR [--ko-list ko_list] \\
      [--gtdbtk gtdbtk_taxonomy.tsv | --sendsketch sendsketch_contigs.tsv --gff annotation.gff3] \\
      --outdir OUT
  metacyc_scorer.py build-refs --refs REFS_DIR [--threads N]      # open data; needs diamond, curl
  metacyc_scorer.py setup --db-dir DB_DIR (--flatfiles X | --installer X | --sif X)   # licensed MetaCyc

Outputs of `run` (OUT/):
  metacyc_reactions.tsv.gz  protein x MetaCyc reaction, with evidence class(es) and DIAMOND detail
  metacyc_pathways.tsv      bin (incl. _community) x pathway with any evidence: PathoLogic call,
                            deciding rule, confidence + tier, evidence counts, supporting hits, MinPath
  metacyc_bin_taxonomy.tsv  GTDB lineage per bin and the NCBI taxids used for taxonomic ranges
Exit status 3 when the MetaCyc reference is missing (after printing how to install it).
"""
import argparse
import gzip
import json
import os
import resource
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))
from metacyc_scorer import confidence, evidence, pathways, ptdat, taxonomy  # noqa: E402

COMMUNITY = '_community'
MISSING_METACYC = """
==============================================================================
MetaCyc reference not found: {path}
MetaCyc is licensed by SRI and is not downloaded automatically. Install it once with
  download-databases.sh --metacyc <flat-file tarball | Pathway Tools installer | .sif>
(or: metacyc_scorer.py setup --db-dir <db_dir> --flatfiles|--installer|--sif <path>)
and pass --metacyc_db <db_dir>/metacyc (found automatically under --db_dir).
==============================================================================
"""
PATHWAY_COLS = ['pathway', 'name', 'super', 'klass', 'call', 'reason', 'confidence', 'confidence_tier',
                'n_rxn', 'n_present', 'frac', 'n_enzymatic', 'n_spontaneous', 'n_missing',
                'n_seq', 'n_seq_id60', 'n_ec_only', 'n_ko_only', 'n_ec_direct_only', 'n_linked_only',
                'seq_id_median', 'seq_id_min', 'n_key', 'n_key_present', 'key_seq', 'n_unique',
                'n_unique_present', 'n_unique_seq', 'tax_range', 'in_tax_range',
                'supporting_hits', 'present_rxns', 'missing_rxns']


def log(msg):
    print(f'[metacyc] {msg}', file=sys.stderr, flush=True)


def read_contig2bin(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            f = line.rstrip('\n').split('\t')
            if len(f) >= 2 and f[0] not in ('contig', 'contig_id'):
                out[f[0]] = f[1]
    return out


def bin_lineages(a, c2b):
    """GTDB lineage per bin -> NCBI TAX frames ({bin: set}, {bin: info row})."""
    table_path = os.path.join(a.refs, 'gtdb', 'gtdb2ncbi_r226.tsv.gz')
    if not os.path.exists(table_path):
        log(f'no {table_path}; taxonomic-range rules off')
        return {}, {}
    if a.gtdbtk and os.path.getsize(a.gtdbtk) > 0:
        gl, src = taxonomy.read_gtdbtk(a.gtdbtk), 'gtdbtk'
    elif a.sendsketch and os.path.getsize(a.sendsketch) > 0:
        lengths = taxonomy.contig_lengths_gff(a.gff) if a.gff else {}
        gl, src = taxonomy.read_sendsketch_bins(a.sendsketch, c2b, lengths), 'sendsketch_contigs'
    else:
        log('no GTDB-Tk or SendSketch taxonomy; taxonomic-range rules off')
        return {}, {}
    table = taxonomy.load(table_path)
    lin, info = {}, {}
    for b in sorted(set(c2b.values())):
        if b not in gl:
            continue
        taxon, ncbi = taxonomy.ncbi_lineage(table, gl[b])
        lin[b] = ncbi
        info[b] = [src, gl[b], taxon or '', ';'.join(sorted(ncbi))]
    log(f'taxonomy from {src}: {len(lin)} of {len(set(c2b.values()))} bins')
    return lin, info


def cmd_run(a):
    t0 = time.time()
    mcdir = a.metacyc
    if not all(os.path.exists(os.path.join(mcdir, f)) for f in ('pathways.dat', 'reactions.dat', 'classes.dat')):
        print(MISSING_METACYC.format(path=mcdir), file=sys.stderr)
        sys.exit(3)
    for p in ('rhea/rhea2metacyc.tsv', 'seqref/seqref_rhea.tsv.gz'):
        if not os.path.exists(os.path.join(a.refs, p)):
            sys.exit(f'[metacyc] ERROR: {os.path.join(a.refs, p)} missing; build the open references with '
                     f'download-databases.sh --metacyc-refs')
    os.makedirs(a.outdir, exist_ok=True)
    stamp = os.path.join(mcdir, 'VERSION.json')
    release = json.load(open(stamp)).get('release', '?') if os.path.exists(stamp) else '?'
    mc = ptdat.MetaCyc(ptdat.load(os.path.join(mcdir, 'pathways.dat')),
                       ptdat.load(os.path.join(mcdir, 'reactions.dat')),
                       ptdat.load(os.path.join(mcdir, 'classes.dat')))
    refs = evidence.References(a.refs, mc, a.ko_list)
    enzrxns = os.path.join(mcdir, 'enzrxns.dat')
    links = confidence.reaction_links(mc, ptdat.load(enzrxns) if os.path.exists(enzrxns) else None)
    if not os.path.exists(enzrxns):
        log(f'no {enzrxns}; LINKED evidence from shared EC numbers only')
    c2b = read_contig2bin(a.contig2bin)
    with open(a.bakta_tsv) as fh:
        n_est = sum(1 for line in fh if not line.startswith('#'))
    k = max(1, -(-n_est // a.chunk))
    tmp = os.path.join(a.outdir, '.chunks')
    parts, subjects, merged_cols = evidence.partition(a.bakta_tsv, a.merged, a.hits, k, tmp)
    seqref = evidence.load_seqref(os.path.join(a.refs, 'seqref', 'seqref_rhea.tsv.gz'), subjects)
    del subjects
    log(f'MetaCyc {release}; {n_est:,} features in {k} chunk(s); {len(seqref):,} Rhea-annotated hit proteins')
    ev = confidence.BinEvidence()
    n_prot = n_with = n_rows = 0
    with gzip.open(os.path.join(a.outdir, 'metacyc_reactions.tsv.gz'), 'wt') as out:
        out.write('\t'.join(evidence.ASSIGN_HEADER) + '\n')
        for prot, assigns in evidence.iter_chunks(parts, merged_cols, seqref, refs,
                                                  dict(min_id=a.min_id, min_qcov=a.min_cov, min_scov=a.min_cov)):
            n_prot += len(prot)
            n_with += len({t['protein'] for t in assigns})
            for t in sorted(assigns, key=lambda t: (t['protein'], t['metacyc_rxn'])):
                out.write('\t'.join(str(x) for x in evidence.assign_row(t, prot, c2b)) + '\n')
                n_rows += 1
                b = c2b.get(prot[t['protein']]['contig'])
                if b:
                    ev.add(b, t)
                ev.add(COMMUNITY, t)
    os.rmdir(tmp)
    log(f'{n_prot:,} proteins; {n_with:,} with reactions; {n_rows:,} protein-reaction rows')
    lin, info = bin_lineages(a, c2b)
    with open(os.path.join(a.outdir, 'metacyc_bin_taxonomy.tsv'), 'w') as fh:
        fh.write('bin\tsource\tgtdb_lineage\tgtdb_taxon_used\tncbi_tax_frames\n')
        for b, v in sorted(info.items()):
            fh.write('\t'.join([b] + v) + '\n')
    n_called = 0
    with open(os.path.join(a.outdir, 'metacyc_pathways.tsv'), 'w') as fh:
        fh.write('\t'.join(['bin'] + PATHWAY_COLS + ['minpath', 'metacyc_release']) + '\n')
        for b in sorted(ev.bins, key=lambda x: (x != COMMUNITY, x)):
            bev = ev.bins[b]
            confidence.add_linked(bev, links)
            present = {r: 1 for r in bev}
            sc = pathways.score_pathways(mc, present, lin.get(b) if b != COMMUNITY else None)
            sel, _ = pathways.minpath(mc, present)
            for p, r in sorted(sc.items()):
                if not r['n_enzymatic']:
                    continue
                confidence.annotate(mc, r, bev)
                n_called += r['call'] and b != COMMUNITY
                row = [b] + [r.get(c, '') for c in PATHWAY_COLS] + [int(p in sel), release]
                fh.write('\t'.join(str(x) for x in row) + '\n')
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    log(f'{len(ev.bins) - (COMMUNITY in ev.bins)} bins, {n_called:,} MAG pathway calls; '
        f'{time.time() - t0:.0f} s, peak RSS {rss:.0f} MB')


def main():
    if sys.argv[1:2] == ['setup']:            # options belong to install.py's own parser
        from metacyc_scorer import install
        install.main(sys.argv[2:])
        return
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    # Python 3.6 (e.g. Grex's host python3) has no add_subparsers(required=...).
    sub = ap.add_subparsers(dest='cmd')
    sub.required = True
    r = sub.add_parser('run', help='reactions and pathways per bin')
    r.add_argument('--bakta-tsv', required=True, help='Bakta annotation.tsv (extra run, basic locus tags)')
    r.add_argument('--merged', help='merged_annotations.tsv (KofamScan/eggNOG KO and EC)')
    r.add_argument('--hits', required=True, help='DIAMOND hits vs REFS/seqref/seqref.dmnd (outfmt 6: '
                   + evidence.DIAMOND_FIELDS + ')')
    r.add_argument('--contig2bin', required=True)
    r.add_argument('--metacyc', required=True, help='licensed MetaCyc reference dir (setup)')
    r.add_argument('--refs', required=True, help='open reference dir (build-refs)')
    r.add_argument('--ko-list', help='KofamScan ko_list (KO definitions give KO -> EC)')
    r.add_argument('--gtdbtk', help='GTDB-Tk summary (gtdbtk_taxonomy.tsv)')
    r.add_argument('--sendsketch', help='sendsketch_contigs.tsv, used when there is no GTDB-Tk summary')
    r.add_argument('--gff', help='GFF3 with ##sequence-region lines (contig lengths for --sendsketch)')
    r.add_argument('--chunk', type=int, default=500000, help='proteins per chunk [500000]')
    r.add_argument('--min-id', type=float, default=35.0)
    r.add_argument('--min-cov', type=float, default=0.6)
    r.add_argument('--outdir', required=True)
    b = sub.add_parser('build-refs', help='fetch and build the open references')
    b.add_argument('--refs', required=True)
    b.add_argument('--threads', type=int, default=8)
    b.add_argument('--keep-intermediate', action='store_true')
    sub.add_parser('setup', help='install licensed MetaCyc (see: metacyc_scorer.py setup -h)')
    a = ap.parse_args()
    if a.cmd == 'run':
        cmd_run(a)
    else:
        from metacyc_scorer import refs
        refs.build(a.refs, a.threads, a.keep_intermediate)


if __name__ == '__main__':
    main()
