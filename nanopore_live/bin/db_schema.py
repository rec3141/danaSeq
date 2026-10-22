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
    # Drop legacy 3-col tables so the new schema applies; old aggregate data is unused.
    legacy = con.execute("""
        SELECT COUNT(*) FROM information_schema.columns
        WHERE table_name = 'sendsketch'
    """).fetchone()[0]
    if legacy and legacy != 5:
        con.execute("DROP TABLE sendsketch")
        con.execute("DELETE FROM import_log WHERE filename LIKE 'sketch/%'")
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

    # Per-read minimap2 alignments against arbitrary references (e.g. AIS
    # genomes). Loaded by import_mapping.py from map/<refname>.txt files.
    con.execute("""
        CREATE TABLE IF NOT EXISTS mapping (
            reference    TEXT,
            qname        TEXT,
            flag         INTEGER,
            rname        TEXT,
            pos          INTEGER,
            mapq         INTEGER,
            cigar        TEXT,
            nm           INTEGER,
            as_score     INTEGER,
            de           REAL,
            identity_pct REAL,
            aligned_len  INTEGER
        )
    """)

    # Protein-DB search results (HMM_SEARCH / DIAMOND_DB_SEARCH): pollutant
    # marker DBs such as CANT-HYD and HADEG, and any later HMM/protein DB.
    # Loaded by import_protein_hits.py from hmm/*.tsv and protdb/*.tsv.
    #   protein_hits        one row per (protein, target) hit; `db` = source DB,
    #                       `chunk` = annotation batch (meta.id) the protein came from
    #   protein_db_annot    target -> gene / category / pathway catalogue per DB
    #   protein_db_searches proteins searched per (db, chunk), so hit counts can
    #                       be normalised (hits per million proteins searched)
    con.execute("""
        CREATE TABLE IF NOT EXISTS protein_hits (
            db       TEXT,
            chunk    TEXT,
            gene_id  TEXT,
            target   TEXT,
            evalue   DOUBLE,
            score    DOUBLE,
            pident   REAL,
            qcov     REAL,
            scov     REAL
        )
    """)
    con.execute("""
        CREATE TABLE IF NOT EXISTS protein_db_annot (
            db         TEXT,
            target     TEXT,
            gene       TEXT,
            category   TEXT,
            pathway    TEXT,
            subpathway TEXT,
            PRIMARY KEY (db, target)
        )
    """)
    con.execute("""
        CREATE TABLE IF NOT EXISTS protein_db_searches (
            db         TEXT,
            chunk      TEXT,
            n_proteins BIGINT,
            n_hits     INTEGER,
            params     TEXT,
            PRIMARY KEY (db, chunk)
        )
    """)
