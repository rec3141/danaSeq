suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(readr)
})

# set up db schema

# clean up rogue connections
# try(DBI::dbDisconnect(DBI::dbConnect(duckdb::duckdb(), dbdir = "dana.duckdb"), shutdown = TRUE), silent = TRUE)

# dbExecute(con, "DROP TABLE stats;")

result = dbExecute(con,
                   "
CREATE TABLE IF NOT EXISTS sendsketch (
  fileid     TEXT,
  ref_name   TEXT,
  ani        REAL,
);")

result = dbExecute(con,
                   "
CREATE TABLE IF NOT EXISTS stats (
  seqid     TEXT PRIMARY KEY,
  length   INTEGER
);")



result = dbExecute(
  con,
  "
CREATE TABLE IF NOT EXISTS import_log (
  filename TEXT PRIMARY KEY,
  imported_at TIMESTAMP DEFAULT current_timestamp
);"
)



# make new table for the fileid and seqids
result = dbExecute(
  con,
  "
CREATE TABLE IF NOT EXISTS sequence_index (
  seqid TEXT PRIMARY KEY,  -- your standard sequence ID
  fileid TEXT              -- the batch / source file / barcode
);"
)


result = dbExecute(
  con,
  "
CREATE TABLE IF NOT EXISTS locus_index (
  seqid TEXT,  -- your standard sequence ID
  locus_tag TEXT PRIMARY KEY,
);"
)

result = dbExecute(
  con,
  "
CREATE TABLE IF NOT EXISTS prokka_annotations (
  locus_tag TEXT PRIMARY KEY,
  ftype     TEXT,
  length_bp INTEGER,
  gene      TEXT,
  ec_number TEXT,
  cog       TEXT,
  product   TEXT
);"
)

result = dbExecute(
  con,
  "
CREATE TABLE IF NOT EXISTS kraken (
  seqid     TEXT PRIMARY KEY,
  taxid   INTEGER,
  taxa_name        TEXT,
);"
)

result = dbExecute(
  con,
  "
CREATE TABLE IF NOT EXISTS krakenreport (
  percent REAL,
  reads INTEGER,
  direct INTEGER,
  rank TEXT,
  taxid INTEGER,
  name TEXT,
  level INTEGER,
  flowcell TEXT,
  barcode TEXT,
  taxonomy TEXT,
  rank_full TEXT,
);"
)


