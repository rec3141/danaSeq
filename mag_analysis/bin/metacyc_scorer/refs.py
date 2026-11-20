"""Build the open references the scorer needs (no licensed data):

  <refs>/rhea/     Rhea TSVs (rhea2metacyc, rhea2ec, rhea2kegg_reaction, rhea-directions)
  <refs>/kegg/     kegg_ko2rn.tsv (KEGG REST link/reaction/ko)
  <refs>/gtdb/     gtdb2ncbi_r226.tsv.gz (GTDB taxon -> NCBI taxids, see taxonomy.py)
  <refs>/seqref/   seqref.dmnd + seqref_rhea.tsv.gz: the DIAMOND reference of Rhea-annotated
                   UniProtKB proteins and their Rhea reactions
  <refs>/VERSION.json

The sequence reference holds every Swiss-Prot entry with a Rhea annotation plus one
Rhea-annotated TrEMBL member per UniRef90 cluster that has no such Swiss-Prot member (the
cluster representative when it is annotated, else the first annotated member). UniRef90 members
are >= 90% identical, so one per cluster keeps hit identity within ~10 points of the full set
(19M sequences) at half the size. For each kept protein and reaction, the share of the cluster's
Rhea-annotated members that carry the same reaction is recorded as cluster agreement.
"""
import datetime
import gzip
import json
import os
import shutil
import subprocess
import sys
from collections import Counter, defaultdict

from . import taxonomy

RHEA = 'https://ftp.expasy.org/databases/rhea'
RHEA_FILES = ('rhea2metacyc.tsv', 'rhea2ec.tsv', 'rhea2kegg_reaction.tsv', 'rhea-directions.tsv')
RHEA_UNIPROT = ('rhea2uniprot_sprot.tsv', 'rhea2uniprot_trembl.tsv.gz')
UNIPROT = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase'
GTDB = 'https://data.gtdb.ecogenomic.org/releases/release226/226.0'
TAXDUMP = 'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'


def log(msg):
    print(f'[metacyc refs] {msg}', file=sys.stderr, flush=True)


def fetch(url, dest):
    if not os.path.exists(dest):
        subprocess.run(['curl', '-fsS', '--retry', '5', '--retry-delay', '30', '-o', dest + '.part', url], check=True)
        os.replace(dest + '.part', dest)


def _gunzip_cmd(threads=2):
    return f'pigz -p {threads} -dc' if shutil.which('pigz') else 'gzip -dc'


def _stream(cmd):
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, universal_newlines=True, bufsize=1 << 20)
    yield from p.stdout
    if p.wait() != 0:
        raise RuntimeError(f'failed: {cmd}')


def load_rhea_uniprot(rhea_dir):
    """-> ({accession: 'm1,m2'}, set(Swiss-Prot accessions)); Rhea IDs folded to master IDs."""
    acc = defaultdict(set)
    sp = set()
    for name in RHEA_UNIPROT:
        path = os.path.join(rhea_dir, name)
        with (gzip.open(path, 'rt') if name.endswith('.gz') else open(path)) as fh:
            next(fh)
            for line in fh:
                f = line.rstrip('\n').split('\t')
                acc[f[3]].add(f[2])
                if 'sprot' in name:
                    sp.add(f[3])
    intern = {}
    return {a: intern.setdefault(t, t) for a, s in acc.items() for t in [','.join(sorted(s))]}, sp


def select_reference(rhea_dir, idmapping, out_dir, tmp):
    """One pass over idmapping_selected (accession, UniRef90) plus an external sort by cluster."""
    acc2rhea, sprot = load_rhea_uniprot(rhea_dir)
    log(f'{len(acc2rhea):,} Rhea-annotated UniProtKB accessions ({len(sprot):,} Swiss-Prot)')
    raw = os.path.join(tmp, 'members.tsv')
    with open(raw, 'w') as out:
        for line in _stream(f'{_gunzip_cmd()} {idmapping} | cut -f1,9'):
            a, ur = line.rstrip('\n').split('\t')
            r = acc2rhea.get(a)
            if r is None:
                continue
            pr = 0 if a in sprot else (1 if ur.split('_', 1)[-1] == a else 2)
            out.write(f'{ur or "none:" + a}\t{pr}\t{a}\t{r}\n')
    del acc2rhea, sprot
    srt = raw + '.sorted'
    subprocess.run(f'LC_ALL=C sort -t "\t" -k1,1 -k2,2n -S 4G --parallel=4 -T {tmp} {raw} -o {srt}',
                   shell=True, check=True)
    os.remove(raw)
    chosen = set()
    n_cl = 0
    with open(srt) as fh, gzip.open(os.path.join(out_dir, 'seqref_rhea.tsv.gz'), 'wt') as out:
        out.write('accession\tuniref90\treviewed\trhea_master\tcluster_n_rhea_members\tcluster_n_with_rhea\n')

        def flush(ur, members):
            if not members:
                return 0
            cnt = Counter(m for _pr, _a, r in members for m in r.split(','))
            keep = [x for x in members if x[0] == 0] or [members[0]]
            for pr, a, r in keep:
                chosen.add(a)
                u = '' if ur.startswith('none:') else ur
                for m in r.split(','):
                    out.write(f'{a}\t{u}\t{int(pr == 0)}\t{m}\t{len(members)}\t{cnt[m]}\n')
            return 1

        cur, members = None, []
        for line in fh:
            ur, pr, a, r = line.rstrip('\n').split('\t')
            if ur != cur:
                n_cl += flush(cur, members)
                cur, members = ur, []
            members.append((int(pr), a, r))
        n_cl += flush(cur, members)
    os.remove(srt)
    log(f'{n_cl:,} UniRef90 clusters; reference = {len(chosen):,} sequences')
    return chosen


def fetch_sequences(chosen, out_path):
    """Stream Swiss-Prot and TrEMBL FASTA from the UniProt FTP, keeping only chosen accessions.
    TrEMBL is ~40 GB compressed; it is filtered in flight, never stored."""
    n = 0
    with gzip.open(out_path + '.part', 'wt', compresslevel=3) as gz:
        for name in ('uniprot_sprot.fasta.gz', 'uniprot_trembl.fasta.gz'):
            keep = False
            for line in _stream(f'curl -fsS --retry 5 --retry-delay 30 {UNIPROT}/complete/{name} | {_gunzip_cmd()}'):
                if line[0] == '>':
                    acc = line.split('|', 2)[1]
                    keep = acc in chosen
                    if keep:
                        gz.write(f'>{acc}\n')
                        n += 1
                elif keep:
                    gz.write(line)
            log(f'{name}: {n:,} sequences kept so far')
    os.replace(out_path + '.part', out_path)
    if n < len(chosen):
        log(f'{len(chosen) - n:,} chosen accessions were not in the FASTA release')
    return n


def build(refs_dir, threads=8, keep_intermediate=False):
    refs_dir = os.path.abspath(refs_dir)
    tmp = os.path.join(refs_dir, '.build_tmp')
    for d in ('rhea', 'kegg', 'gtdb', 'seqref'):
        os.makedirs(os.path.join(refs_dir, d), exist_ok=True)
    os.makedirs(tmp, exist_ok=True)
    rd = os.path.join(refs_dir, 'rhea')
    for f in RHEA_FILES:
        fetch(f'{RHEA}/tsv/{f}', os.path.join(rd, f))
    for f in RHEA_UNIPROT:
        fetch(f'{RHEA}/tsv/{f}', os.path.join(tmp, f))
    fetch(f'{RHEA}/rhea-release.properties', os.path.join(rd, 'rhea-release.properties'))
    fetch('https://rest.kegg.jp/link/reaction/ko', os.path.join(refs_dir, 'kegg', 'kegg_ko2rn.tsv'))
    log('GTDB r226 metadata + NCBI taxdump -> gtdb2ncbi')
    meta = []
    for f in ('bac120_metadata_r226.tsv.gz', 'ar53_metadata_r226.tsv.gz'):
        fetch(f'{GTDB}/{f}', os.path.join(tmp, f))
        meta.append(os.path.join(tmp, f))
    fetch(TAXDUMP, os.path.join(tmp, 'taxdump.tar.gz'))
    subprocess.run(['tar', 'xzf', os.path.join(tmp, 'taxdump.tar.gz'), '-C', tmp, 'nodes.dmp', 'merged.dmp'], check=True)
    taxonomy.build(meta, tmp, os.path.join(refs_dir, 'gtdb', 'gtdb2ncbi_r226.tsv.gz'))
    sq = os.path.join(refs_dir, 'seqref')
    idm = os.path.join(tmp, 'idmapping_selected.tab.gz')
    fetch(f'{UNIPROT}/idmapping/idmapping_selected.tab.gz', idm)
    chosen = select_reference(tmp, idm, sq, tmp)
    faa = os.path.join(sq, 'seqref.faa.gz')
    n = fetch_sequences(chosen, faa)
    subprocess.run(['diamond', 'makedb', '--in', faa, '-d', os.path.join(sq, 'seqref'), '-p', str(threads)], check=True)
    rel = {}
    with open(os.path.join(rd, 'rhea-release.properties')) as fh:
        for line in fh:
            if '=' in line:
                k, v = line.strip().split('=', 1)
                rel[k] = v
    with open(os.path.join(refs_dir, 'VERSION.json'), 'w') as fh:
        json.dump(dict(rhea_release=rel.get('rhea.release.number'), rhea_date=rel.get('rhea.release.date'),
                       seqref_sequences=n, gtdb='r226', built=datetime.date.today().isoformat()), fh, indent=1)
    if not keep_intermediate:
        shutil.rmtree(tmp)
        os.remove(faa)
    log(f'done: {refs_dir}')
