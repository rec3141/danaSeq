#!/usr/bin/env python3
"""Recover every intact FASTQ record from a damaged multi-member .fastq.gz.

MinKNOW writes one gzip member per chunk and concat_run_fastqs.sh joins them
byte for byte, so a damaged stretch only breaks the members it touches. A
plain decompressor (gzip, BBTools) stops at the first bad byte and loses
everything after it; this walks the file member by member and, when a member
fails, resynchronises on the next gzip header and carries on.

Records are written only when structurally whole (@header, sequence, +,
quality of the same length, sequence/quality alphabets as expected), so the
garbage a deflate stream can emit just before it detects damage never gets
through. A read ID is written once: damage seen so far has been the file's
own data written at the wrong offset, so the damaged stretch can hold intact
copies of members from elsewhere in the file. Writes uncompressed FASTQ to
stdout; a summary goes to stderr.

Usage: fastq_gz_salvage.py IN.fastq.gz | gzip > OUT.fastq.gz
"""
import mmap
import re
import sys
import zlib

MAGIC = b"\x1f\x8b\x08"
STEP = 1 << 20
SEQ_OK = re.compile(rb"^[ACGTNacgtnRYKMSWBDHVUryk]*$")


class Records:
    """Buffer decompressed text, emit whole valid records, drop the rest."""

    def __init__(self, out):
        self.out, self.buf = out, b""
        self.kept = self.dropped = self.dups = 0
        self.seen = set()

    def feed(self, data):
        self.buf += data
        lines = self.buf.split(b"\n")
        self.buf = lines.pop()
        i = 0
        while i + 4 <= len(lines):
            h, s, p, q = lines[i:i + 4]
            if (h.startswith(b"@") and p.startswith(b"+") and len(s) == len(q)
                    and SEQ_OK.match(s) and all(33 <= c <= 126 for c in q)):
                rid = h.split(None, 1)[0]
                if rid in self.seen:
                    self.dups += 1
                else:
                    self.seen.add(rid)
                    self.out.write(b"\n".join((h, s, p, q)) + b"\n")
                    self.kept += 1
                i += 4
            else:  # out of frame: slide to the next line that looks like a header
                self.dropped += 1
                i += 1
                while i < len(lines) and not lines[i].startswith(b"@"):
                    i += 1
        self.buf = b"\n".join(lines[i:] + [self.buf])

    def reset(self):  # a member died: whatever is buffered is untrustworthy
        if self.buf:
            self.dropped += 1
        self.buf = b""


def main():
    path = sys.argv[1]
    out = sys.stdout.buffer
    rec = Records(out)
    with open(path, "rb") as fh:
        mm = mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ)
        size, pos = len(mm), 0
        members = bad = 0
        skipped_bytes = 0
        while pos < size:
            d = zlib.decompressobj(31)
            p = pos
            try:
                while not d.eof and p < size:
                    chunk = mm[p:p + STEP]
                    rec.feed(d.decompress(chunk))
                    p += len(chunk)
                    if d.eof:
                        p -= len(d.unused_data)
                rec.feed(d.flush())
                if not d.eof:
                    raise zlib.error("truncated member")
                members += 1
                pos = p
            except zlib.error:
                bad += 1
                rec.reset()
                nxt = mm.find(MAGIC, pos + 1)
                nxt = size if nxt < 0 else nxt
                skipped_bytes += nxt - pos
                pos = nxt
        mm.close()
    rec.reset()
    out.flush()
    print(f"members_ok={members} members_bad={bad} bytes_skipped={skipped_bytes} "
          f"records_kept={rec.kept} duplicates_dropped={rec.dups} fragments_dropped={rec.dropped}", file=sys.stderr)
    sys.exit(0 if rec.kept else 1)


if __name__ == "__main__":
    main()
