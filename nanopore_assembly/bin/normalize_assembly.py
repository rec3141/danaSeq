#!/usr/bin/env python3
"""
Normalize a raw assembler FASTA: sort contigs by length, rename to contig_NNNNNN,
emit Flye-style assembly_info.txt, and (optionally) remap a GFA with the new names.

Replaces a fragile multi-stage awk|sort|awk pipeline that silently produced empty
outputs when contigs were multi-megabase (sequence-as-single-tab-field overflows
some awk/sort implementations).

Memory: streams sequences via a two-pass approach (index headers + offsets in
pass 1; seek and emit in pass 2). Peak memory is O(n_contigs * ~100B).
"""
import argparse
import os
import re
import sys


def index_fasta(path):
    """Pass 1: build (header_line, seq_offset, seq_length) for every contig."""
    entries = []
    cur_header = None
    cur_offset = None
    cur_length = 0
    with open(path, "rb") as fh:
        while True:
            pos = fh.tell()
            line = fh.readline()
            if not line:
                break
            if line.startswith(b">"):
                if cur_header is not None:
                    entries.append((cur_header, cur_offset, cur_length))
                cur_header = line.rstrip(b"\r\n").decode("utf-8", errors="replace")[1:]
                cur_offset = fh.tell()
                cur_length = 0
            else:
                cur_length += len(line.rstrip(b"\r\n"))
        if cur_header is not None:
            entries.append((cur_header, cur_offset, cur_length))
    return entries


def emit_wrapped(out_fh, seq_bytes, wrap):
    """Write sequence bytes with `wrap`-column line wrapping."""
    for i in range(0, len(seq_bytes), wrap):
        out_fh.write(seq_bytes[i:i + wrap])
        out_fh.write(b"\n")


def read_sequence(fh, offset, length):
    """Pass 2: seek and read sequence bytes, stripping newlines."""
    fh.seek(offset)
    buf = bytearray()
    remaining = length
    while remaining > 0:
        line = fh.readline()
        if not line:
            break
        if line.startswith(b">"):
            break
        stripped = line.rstrip(b"\r\n")
        buf.extend(stripped)
        remaining -= len(stripped)
    return bytes(buf)


def parse_header_flye(header):
    """Flye headers are just contig names — info comes from assembly_info.txt."""
    return {"old_name": header.split()[0], "len": None, "cov": None, "circ": None}


def parse_header_metamdbg(header):
    """metaMDBG: >ctgN length=LEN circular={yes,no} coverage=N"""
    parts = header.split()
    info = {"old_name": parts[0], "len": None, "cov": None, "circ": "N"}
    for p in parts[1:]:
        if "=" not in p:
            continue
        k, v = p.split("=", 1)
        if k == "length":
            info["len"] = int(v)
        elif k == "coverage":
            try:
                info["cov"] = float(v)
            except ValueError:
                pass
        elif k == "circular":
            info["circ"] = "Y" if v.lower() == "yes" else "N"
    return info


def parse_header_myloasm(header):
    """myloasm: >NAME_len-LEN_circular-{yes,no,possibly}_depth-D1-D2-D3_duplicated-{yes,no} mult=M"""
    parts = header.split()
    name = parts[0]
    info = {"old_name": name, "len": None, "cov": None, "circ": "N"}
    for field in name.split("_"):
        if field.startswith("len-"):
            try:
                info["len"] = int(field[4:])
            except ValueError:
                pass
        elif field.startswith("circular-"):
            v = field[len("circular-"):]
            info["circ"] = "Y" if v in ("yes", "possibly") else "N"
        elif field.startswith("depth-"):
            depth_parts = field[len("depth-"):].split("-")
            try:
                info["cov"] = float(depth_parts[0])
            except (ValueError, IndexError):
                pass
    if info["cov"] is None:
        for p in parts[1:]:
            if p.startswith("mult="):
                try:
                    info["cov"] = float(p[5:])
                except ValueError:
                    pass
                break
    return info


HEADER_PARSERS = {
    "flye": parse_header_flye,
    "metamdbg": parse_header_metamdbg,
    "myloasm": parse_header_myloasm,
}


def remap_gfa(in_path, out_path, name_map):
    """Rewrite S/P/L lines in a GFA with new contig names."""
    n_renamed = 0
    with open(in_path, "rb") as fin, open(out_path, "wb") as fout:
        for line in fin:
            if not line or line[:1] not in (b"S", b"L", b"P", b"H", b"C"):
                fout.write(line)
                continue
            fields = line.rstrip(b"\r\n").split(b"\t")
            if not fields:
                fout.write(line)
                continue
            rec = fields[0]
            if rec == b"S" and len(fields) >= 2:
                old = fields[1].decode()
                if old in name_map:
                    fields[1] = name_map[old].encode()
                    n_renamed += 1
            elif rec == b"L" and len(fields) >= 5:
                for idx in (1, 3):
                    old = fields[idx].decode()
                    if old in name_map:
                        fields[idx] = name_map[old].encode()
            elif rec == b"P" and len(fields) >= 3:
                # P lines have comma-separated segments with optional +/- orientation
                segs = fields[2].split(b",")
                new_segs = []
                for seg in segs:
                    m = re.match(rb"^([^+\-]+)([+\-]?)$", seg)
                    if m:
                        old = m.group(1).decode()
                        orient = m.group(2)
                        new_name = name_map.get(old, old).encode()
                        new_segs.append(new_name + orient)
                    else:
                        new_segs.append(seg)
                fields[2] = b",".join(new_segs)
            fout.write(b"\t".join(fields))
            fout.write(b"\n")
    return n_renamed


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, help="Raw assembler FASTA")
    ap.add_argument("--output", required=True, help="Renamed FASTA output")
    ap.add_argument("--output-info", required=True, help="assembly_info.txt output")
    ap.add_argument("--output-map", required=True, help="name_map.tsv output (old<TAB>new)")
    ap.add_argument("--header-format", required=True, choices=list(HEADER_PARSERS.keys()))
    ap.add_argument("--source-info", help="Existing assembly_info.txt to remap (Flye)")
    ap.add_argument("--gfa-in", help="Input GFA to remap")
    ap.add_argument("--gfa-out", help="Output GFA path (required if --gfa-in given)")
    ap.add_argument("--min-len", type=int, default=1000, help="Drop contigs shorter than this")
    ap.add_argument("--wrap", type=int, default=80, help="FASTA line width")
    args = ap.parse_args()

    if args.gfa_in and not args.gfa_out:
        ap.error("--gfa-out is required when --gfa-in is given")

    parser = HEADER_PARSERS[args.header_format]

    sys.stderr.write(f"[normalize_assembly] indexing {args.input}\n")
    entries = index_fasta(args.input)
    sys.stderr.write(f"[normalize_assembly] found {len(entries)} contigs\n")
    if not entries:
        sys.stderr.write("[ERROR] no contigs in input FASTA\n")
        sys.exit(1)

    # Filter by min length, sort by length descending (stable on header)
    filtered = [(h, off, ln) for (h, off, ln) in entries if ln >= args.min_len]
    if not filtered:
        max_len = max(ln for (_, _, ln) in entries)
        sys.stderr.write(
            f"[ERROR] no contigs >= {args.min_len}bp (max: {max_len}bp)\n"
        )
        sys.exit(1)
    filtered.sort(key=lambda x: (-x[2], x[0]))

    n_total = len(filtered)
    pad = len(str(n_total))

    # Parse optional source assembly_info.txt for Flye (keeps Flye's full metadata)
    source_info_by_name = None
    source_info_header = None
    if args.source_info:
        source_info_by_name = {}
        with open(args.source_info) as fh:
            source_info_header = fh.readline().rstrip("\n")
            for line in fh:
                if not line.strip():
                    continue
                fields = line.rstrip("\n").split("\t")
                if fields:
                    source_info_by_name[fields[0]] = fields[1:]

    # Pass 2: write outputs
    name_map = {}
    info_rows = []
    with open(args.input, "rb") as fin, \
         open(args.output, "wb") as fout, \
         open(args.output_map, "w") as fmap:
        for i, (header, offset, length) in enumerate(filtered, start=1):
            new_name = f"contig_{i:0{pad}d}"
            old_name = header.split()[0]
            name_map[old_name] = new_name
            fmap.write(f"{old_name}\t{new_name}\n")

            fout.write(f">{new_name}\n".encode())
            seq = read_sequence(fin, offset, length)
            if len(seq) != length:
                sys.stderr.write(
                    f"[ERROR] sequence length mismatch for {old_name}: "
                    f"expected {length}, got {len(seq)}\n"
                )
                sys.exit(1)
            emit_wrapped(fout, seq, args.wrap)

            if source_info_by_name is not None:
                # Flye path: re-emit Flye's row with renamed contig
                src = source_info_by_name.get(old_name)
                if src is not None:
                    info_rows.append([new_name] + src)
                else:
                    info_rows.append([new_name, str(length), "0", "N"])
            else:
                parsed = parser(header)
                length_v = parsed.get("len") or length
                cov_v = parsed.get("cov")
                circ_v = parsed.get("circ") or "N"
                info_rows.append([
                    new_name,
                    str(length_v),
                    f"{cov_v:.2f}" if isinstance(cov_v, float) else (str(cov_v) if cov_v is not None else "0"),
                    circ_v,
                ])

    # Write assembly_info.txt
    with open(args.output_info, "w") as fout:
        if source_info_header is not None:
            fout.write(source_info_header + "\n")
        else:
            fout.write("#seq_name\tlength\tcov.\tcirc.\n")
        for row in info_rows:
            fout.write("\t".join(row) + "\n")

    # Sanity-check: every contig made it through, FASTA non-empty
    out_size = os.path.getsize(args.output)
    if out_size == 0:
        sys.stderr.write("[ERROR] output FASTA is empty\n")
        sys.exit(1)

    sys.stderr.write(
        f"[normalize_assembly] wrote {n_total} contigs to {args.output} ({out_size} bytes)\n"
    )

    if args.gfa_in:
        if not os.path.exists(args.gfa_in):
            sys.stderr.write(
                f"[normalize_assembly] gfa-in {args.gfa_in} missing; writing stub\n"
            )
            with open(args.gfa_out, "w") as fout:
                fout.write("H\tVN:Z:1.0\n")
        else:
            sys.stderr.write(
                f"[normalize_assembly] remapping GFA {args.gfa_in} -> {args.gfa_out}\n"
            )
            n_renamed = remap_gfa(args.gfa_in, args.gfa_out, name_map)
            sys.stderr.write(f"[normalize_assembly] remapped {n_renamed} GFA S-lines\n")


if __name__ == "__main__":
    main()
