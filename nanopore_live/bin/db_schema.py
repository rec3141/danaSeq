"""Shared DuckDB schema for the Dana real-time pipeline.

Usage:
    from db_schema import ensure_schema
    con = duckdb.connect('dana.duckdb')
    ensure_schema(con)
"""


def ensure_schema(con):
    """Create all pipeline tables if they don't exist."""

    con.execute("""
        CREATE TABLE IF NOT EXISTS import_log (
            filename    TEXT PRIMARY KEY,
            imported_at TIMESTAMP DEFAULT current_timestamp
        )
    """)

    con.execute("""
        CREATE TABLE IF NOT EXISTS kraken (
            seqid     TEXT PRIMARY KEY,
            taxid     INTEGER,
            taxa_name TEXT
        )
    """)

    con.execute("""
        CREATE TABLE IF NOT EXISTS krakenreport (
            percent   REAL,
            reads     INTEGER,
            direct    INTEGER,
            rank      TEXT,
            taxid     INTEGER,
            name      TEXT,
            level     INTEGER,
            flowcell  TEXT,
            barcode   TEXT,
            taxonomy  TEXT,
            rank_full TEXT
        )
    """)

    # Per-read GTDB classification from sendsketch against the local GTDB server.
    # Replaces the old per-fasta aggregate schema (fileid, ref_name, ani).
    con.execute("""
        CREATE TABLE IF NOT EXISTS sendsketch (
            read_id  TEXT,
            status   TEXT,
            ani      REAL,
            ref_name TEXT,
            lineage  TEXT
        )
    """)

    con.execute("""
        CREATE TABLE IF NOT EXISTS stats (
            seqid  TEXT,
            length INTEGER
        )
    """)

    con.execute("""
        CREATE TABLE IF NOT EXISTS sequence_index (
            seqid  TEXT PRIMARY KEY,
            fileid TEXT
        )
    """)

    con.execute("""
        CREATE TABLE IF NOT EXISTS locus_index (
            seqid     TEXT,
            locus_tag TEXT PRIMARY KEY
        )
    """)

    con.execute("""
        CREATE TABLE IF NOT EXISTS prokka_annotations (
            locus_tag TEXT PRIMARY KEY,
            ftype     TEXT,
            length_bp INTEGER,
            gene      TEXT,
            ec_number TEXT,
            cog       TEXT,
            product   TEXT
        )
    """)

    con.execute("""
        CREATE TABLE IF NOT EXISTS read_contig_map (
            read_id   TEXT,
            contig_id TEXT,
            fileid    TEXT
        )
    """)
