suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
#setwd("/data/project_20250421X1/out_dana")

#open connection to db
con <- dbConnect(duckdb::duckdb(), dbdir = "dana.duckdb")

#setup schema
source("/work/apps/dana/log-db.r")

# make a temp staging table with kraken's schema
DBI::dbExecute(con, "CREATE TEMP TABLE IF NOT EXISTS kraken_stage AS SELECT * FROM kraken WHERE 0")

# Read the output, no header
kraken_files = list.files("kraken",pattern = "*.tsv",full.names = TRUE)

allfiles = dbGetQuery(con, "SELECT filename FROM import_log WHERE filename LIKE 'kraken%'")

if(nrow(allfiles) > 0) { kraken_files <- kraken_files[!(kraken_files %in% allfiles$filename)] }
for (file in kraken_files) {
  dbBegin(con)
  ok <- FALSE
  try({
    kraken_in <- readr::read_tsv(
      file,
      col_names = c("seqid", "taxa_w_id"),
      comment = "#", progress = FALSE, show_col_types = FALSE,
      col_select = c(1, 2)
    )
    
    if (nrow(kraken_in) > 0) {
      kraken_in$taxid    <- as.numeric(sub(".*\\(taxid ([0-9]+)\\)$", "\\1", kraken_in$taxa_w_id))
      kraken_in$taxa_name <- sub(" \\(taxid [0-9]+\\)$", "", kraken_in$taxa_w_id)
      kraken_in$taxa_w_id <- NULL
      
      # (optional) dedupe within the file to save work
      kraken_in <- kraken_in[!duplicated(kraken_in$seqid), , drop = FALSE]
      
      # stage → upsert (ignore existing seqid)
      DBI::dbAppendTable(con, "kraken_stage", kraken_in)
      DBI::dbExecute(con, "
        INSERT INTO kraken BY NAME
        SELECT * FROM kraken_stage
        ON CONFLICT (seqid) DO NOTHING
      ")
      DBI::dbExecute(con, "DELETE FROM kraken_stage")
    }
    
    # log file import; safe if already logged
    DBI::dbExecute(con,
                   "INSERT INTO import_log (filename) VALUES (?) ON CONFLICT DO NOTHING",
                   params = list(file)
    )
    
    ok <- TRUE
  }, silent = FALSE)
  
  if (ok) dbCommit(con) else { dbRollback(con); warning(sprintf("Import failed: %s", file)) }
}


#dbGetQuery(con, "SELECT * FROM sendsketch LIMIT 10")
#dbGetQuery(con, "SELECT COUNT(*) AS n FROM sendsketch")

dbDisconnect(con, shutdown = TRUE)

#source("/data/work/dana/prokka-db.r")

