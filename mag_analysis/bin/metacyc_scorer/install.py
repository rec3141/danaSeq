"""Install a user's licensed MetaCyc data as the scorer's reference: <db_dir>/metacyc/.

Exactly one source:
  --flatfiles PATH   BioCyc/MetaCyc flat-file tarball (.tar.gz/.tgz/.tar) or a directory
                     containing pathways.dat, reactions.dat, classes.dat (any depth)
  --sif PATH         an existing Pathway Tools Apptainer image (e.g. MetaPathways build_pt)
  --installer PATH   a Pathway Tools Linux installer (pathway-tools-*-linux-64-tier1-install).
                     Installed unattended as the calling user into a scratch prefix, flat files
                     exported in batch Lisp mode, prefix deleted afterwards (--keep-install keeps it).

Pathway Tools runs need a runtime with Motif/X libraries and xvfb-run:
  --sif           the SIF itself is the runtime
  --installer     the 'ptools-base' root filesystem shipped in the MetaPathways image
                  (quay.io/hallamlab/metapathways): through docker by default; on Apptainer hosts
                  pass --runtime with that image as a .sif (its ptools-base is copied out and run
                  as a sandbox) or with an already extracted sandbox directory
Apptainer is called directly when it is on PATH; otherwise through the MetaPathways docker
image (which carries apptainer). No root needed in either case.

Output: <db_dir>/metacyc/{pathways,reactions,classes}.dat (+ compounds.dat), VERSION.json.
Files are written owner-only. They are licensed: never put them in images, git or public places.
"""
import argparse
import datetime
import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile

NEEDED = ('pathways.dat', 'reactions.dat', 'classes.dat')
OPTIONAL = ('compounds.dat', 'enzrxns.dat', 'proteins.dat')
MP_IMAGE = 'quay.io/hallamlab/metapathways:4.0.3'
PT_BASE = '/opt/metapathways/ptools-base'
# {org} is a PGDB org ID: meta for MetaCyc, ecoli for EcoCyc (tools/metacyc_benchmark.py)
EXPORT_LISP = ('(progn (so (quote {org})) (format t "~%MC-KB ~a~%" (current-kb)) '
               '(create-flat-files-for-current-kb :biopax? nil :sbml? nil :dump-fasta? nil) '
               '(format t "~%MC-EXPORT-DONE~%") (exit))')


def die(msg):
    sys.exit(f'[metacyc setup] ERROR: {msg}')


def version_of(pathways_dat):
    info = {}
    with open(pathways_dat, encoding='latin-1') as fh:
        for line in fh:
            if not line.startswith('#'):
                break
            m = re.match(r'#\s*(Version|Database|Date and time generated):\s*(.+)', line)
            if m:
                info[m.group(1)] = m.group(2).strip()
    return info


def install_files(src_dir, db_dir, source, log=sys.stderr):
    found = {}
    for name in NEEDED + OPTIONAL:
        hits = sorted(glob.glob(os.path.join(src_dir, '**', name), recursive=True), key=len)
        if hits:
            found[name] = hits[0]
    miss = [n for n in NEEDED if n not in found]
    if miss:
        die(f'no {", ".join(miss)} under {src_dir}')
    info = version_of(found['pathways.dat'])
    if 'MetaCyc' not in info.get('Database', 'MetaCyc'):
        die(f'{found["pathways.dat"]} is from {info.get("Database")}, not MetaCyc')
    out = os.path.join(db_dir, 'metacyc')
    os.makedirs(out, exist_ok=True)
    os.chmod(out, 0o700)
    for name, path in found.items():
        dst = os.path.join(out, name)
        shutil.copyfile(path, dst)
        os.chmod(dst, 0o600)
    stamp = dict(database='MetaCyc', release=info.get('Version', 'unknown'),
                 generated=info.get('Date and time generated', ''), source=source,
                 installed=datetime.datetime.now().isoformat(timespec='seconds'),
                 files=sorted(found), licensed=True)
    with open(os.path.join(out, 'VERSION.json'), 'w') as fh:
        json.dump(stamp, fh, indent=1)
    os.chmod(os.path.join(out, 'VERSION.json'), 0o600)
    print(f'[metacyc setup] MetaCyc {stamp["release"]} -> {out} ({len(found)} files)', file=log)
    return out


def from_flatfiles(path, db_dir, tmp):
    if os.path.isdir(path):
        return install_files(path, db_dir, f'flatfiles:{os.path.abspath(path)}')
    if not tarfile.is_tarfile(path):
        die(f'{path} is neither a directory nor a tar archive')
    want = set(NEEDED + OPTIONAL)
    with tarfile.open(path) as tf:
        members = [m for m in tf.getmembers() if m.isfile() and os.path.basename(m.name) in want]
        for m in members:                       # extract only the files we need, flattened
            m.name = os.path.basename(m.name)
            tf.extract(m, tmp)
    return install_files(tmp, db_dir, f'flatfiles:{os.path.abspath(path)}')


def _exec(container, binds, script, workdir):
    """Run `bash -c script` inside an Apptainer container/sandbox as the calling user."""
    bflags = sum((['-B', b] for b in binds), [])
    if shutil.which('apptainer') or shutil.which('singularity'):
        exe = shutil.which('apptainer') or shutil.which('singularity')
        cmd = [exe, 'exec', '--cleanenv'] + bflags + [container, 'bash', '-c', script]
    elif shutil.which('docker'):
        vols = []
        hosts = [b.split(':')[0] for b in binds]
        if os.path.exists(container) and not container.startswith('/opt/'):
            hosts.append(os.path.dirname(os.path.abspath(container)))
        for host in dict.fromkeys(hosts):
            # only the scratch work dir is writable; sources (installer, SIF) are mounted read-only
            ro = '' if host.startswith(workdir) else ':ro'
            vols += ['-v', f'{host}:{host}{ro}']
        cmd = ['docker', 'run', '--rm', '--device', '/dev/fuse', '--security-opt', 'seccomp=unconfined',
               '--security-opt', 'apparmor=unconfined', '--user', f'{os.getuid()}:{os.getgid()}',
               '-v', '/etc/passwd:/etc/passwd:ro', '-v', '/etc/group:/etc/group:ro',
               '-e', f'HOME={workdir}/home', '-w', workdir] + vols + \
              [MP_IMAGE, 'apptainer', 'exec'] + bflags + [container, 'bash', '-c', script]
    else:
        die('neither apptainer nor docker is available to run Pathway Tools')
    print('[metacyc setup] ' + ' '.join(cmd[:6]) + ' ...', file=sys.stderr)
    return subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)


def _export(container, binds, pt_exe, ptlocal, workdir, logpath, extra_env='', org='meta'):
    lisp = EXPORT_LISP.format(org=org)
    script = (f'unset DISPLAY; export HOME={workdir}/home; {extra_env} '
              f'xvfb-run -a {pt_exe} -no-patch-download -lisp -eval \'{lisp}\'')
    r = _exec(container, binds, script, workdir)
    with open(logpath, 'w') as fh:
        fh.write(r.stdout)
    # Pathway Tools echoes the command, so the marker counts only as a line of its own
    if 'MC-EXPORT-DONE' not in r.stdout.splitlines():
        die(f'Pathway Tools export did not finish; see {logpath}')


def from_sif(sif, db_dir, work):
    """The SIF's MetaCyc data dir is read-only; bind writable copies over it and over ptools-local."""
    os.makedirs(f'{work}/home', exist_ok=True)
    probe = _exec(sif, [work], 'ls -d /opt/pathway-tools/aic-export/pgdbs/biocyc/metacyc/*/data; '
                  'ls -d /opt/ptools-local-template /data/ptools-local 2>/dev/null', work)
    dirs = [x for x in probe.stdout.split() if x.startswith('/')]
    data = next((d for d in dirs if d.endswith('/data') and '/metacyc/' in d), None)
    tmpl = next((d for d in dirs if 'ptools-local' in d), None)
    if not data:
        die(f'no MetaCyc PGDB found inside {sif}')
    mc = f'{work}/mcdata'
    pl = f'{work}/ptools-local'
    os.makedirs(mc, exist_ok=True)
    os.makedirs(pl, exist_ok=True)
    if tmpl:
        _exec(sif, [work], f'cp -a {tmpl}/. {pl}/', work)
    binds = [work, f'{mc}:{data}', f'{pl}:/data/ptools-local']
    _export(sif, binds, '/opt/pathway-tools/pathway-tools', pl, work, f'{work}/export.log')
    return install_files(mc, db_dir, f'sif:{os.path.abspath(sif)}')


def _runtime_from_sif(sif, work):
    """The MetaPathways image carries the Pathway Tools runtime as a root filesystem at
    /opt/metapathways/ptools-base. Apptainer cannot start a container nested in another one
    without privileges, so the rootfs is copied out as the user and run as a sandbox."""
    sandbox = os.path.join(work, 'ptools-base')
    # A few root-only files in the image (dpkg locks, btmp, apt logs) can't be read as a user
    # and aren't needed, so copy errors are tolerated; the runtime's own dirs are checked instead.
    r = _exec(sif, [work], f'test -d {PT_BASE} && {{ cp -a {PT_BASE} {sandbox} 2>/dev/null; true; }}', work)
    if r.returncode or not all(os.path.isdir(os.path.join(sandbox, d)) for d in ('usr', 'lib', 'etc')):
        die(f'{sif} has no {PT_BASE}; pass a MetaPathways image or a Pathway Tools runtime sandbox\n{r.stdout}')
    return sandbox


def from_installer(installer, db_dir, work, runtime=None):
    """Unattended user-level install into <work>/pt, export, return reference dir."""
    os.makedirs(f'{work}/home', exist_ok=True)
    os.makedirs(f'{work}/bin', exist_ok=True)
    stub = f'{work}/bin/xterm'          # the installer refuses to start without a terminal emulator
    with open(stub, 'w') as fh:
        fh.write('#!/bin/sh\nwhile [ $# -gt 0 ]; do [ "$1" = "-e" ] && { shift; exec "$@"; }; shift; done\n')
    os.chmod(stub, 0o755)
    inst = os.path.realpath(installer)
    apptainer_host = bool(shutil.which('apptainer') or shutil.which('singularity'))
    if runtime is None and apptainer_host:
        die('--installer on an Apptainer host needs --runtime: the MetaPathways image '
            f'(apptainer pull mp.sif docker://{MP_IMAGE}) or a sandbox of its {PT_BASE}')
    container = runtime or PT_BASE
    if runtime and runtime.endswith('.sif'):
        container = _runtime_from_sif(os.path.realpath(runtime), work)
    binds = [work, f'{os.path.dirname(inst)}']
    script = (f'unset DISPLAY; export PATH={work}/bin:$PATH; '
              f'{inst} --mode unattended --InstallDir {work}/pt --PTOOLS_LOCAL_PATH {work}/ptl '
              f'--InstallDesktopShortcuts 0 --UninstallPrevious 0')
    r = _exec(container, binds, script, work)
    with open(f'{work}/install.log', 'w') as fh:
        fh.write(r.stdout)
    if r.returncode or not os.path.exists(f'{work}/pt/pathway-tools'):
        die(f'Pathway Tools installer failed; see {work}/install.log')
    _export(container, binds, f'{work}/pt/pathway-tools', f'{work}/ptl', work, f'{work}/export.log')
    data = glob.glob(f'{work}/pt/aic-export/pgdbs/biocyc/metacyc/*/data')
    if not data:
        die('MetaCyc PGDB not found in the installation')
    return install_files(data[0], db_dir, f'installer:{os.path.basename(inst)}')


def main(argv=None):
    ap = argparse.ArgumentParser(prog='metacyc_scorer.py setup', description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--db-dir', required=True, help='danaSeq database dir; MetaCyc goes to <db-dir>/metacyc')
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--flatfiles')
    g.add_argument('--sif')
    g.add_argument('--installer')
    ap.add_argument('--runtime', help='Apptainer SIF/sandbox with Pathway Tools runtime libs (installer mode)')
    ap.add_argument('--work-dir', help='scratch dir for the install/export (default: <db-dir>/.metacyc_setup)')
    ap.add_argument('--keep-install', action='store_true')
    a = ap.parse_args(argv)
    # Real paths throughout: a bind through a symlink (e.g. Grex's ~/scratch -> /project/...)
    # leaves the target missing inside the container.
    db = os.path.realpath(a.db_dir)
    os.makedirs(db, exist_ok=True)
    work = os.path.realpath(a.work_dir or os.path.join(db, '.metacyc_setup'))
    os.makedirs(work, exist_ok=True)
    try:
        if a.flatfiles:
            with tempfile.TemporaryDirectory(dir=work) as tmp:
                from_flatfiles(a.flatfiles, db, tmp)
        elif a.sif:
            from_sif(os.path.realpath(a.sif), db, work)
        else:
            from_installer(os.path.realpath(a.installer), db, work,
                           os.path.realpath(a.runtime) if a.runtime else None)
    finally:
        if not a.keep_install:
            shutil.rmtree(work, ignore_errors=True)


