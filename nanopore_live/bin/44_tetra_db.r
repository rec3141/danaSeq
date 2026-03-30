suppressPackageStartupMessages({
library(DBI); library(duckdb); library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

#setwd("/data/project_20250421X1/out_dana/")

#open connection to database
con <- dbConnect(duckdb::duckdb(), dbdir = "dana.duckdb")

# setup db schema
source("/work/apps/dana/log-db.r")

# Mapping: a data.frame with columns `old` and `new`
tnfs <- as.vector(t(read_tsv("tnfs.txt",col_names = FALSE, progress=FALSE, show_col_types = FALSE)))
type = c("TEXT",rep("REAL",length(tnfs)-1))
schema <- paste(tnfs, type, collapse = ",\n  ")

create_sql <- paste0(
  "CREATE TABLE IF NOT EXISTS tetra_data (\n ",
  schema,
  "\n);"
)

result = dbExecute(con, create_sql)

process_lrn_file <- function(filepath) {
  read_tsv(filepath, comment = "%", col_names = tnfs, progress = FALSE, show_col_types = FALSE)
}

lrnfiles <- list.files("tetra", pattern = "\\.lrn$", full.names = TRUE)
fileids = sub("\\..*", "", basename(lrnfiles))

allfiles = dbGetQuery(con, "SELECT filename FROM import_log WHERE filename LIKE 'tetra%'")


if(nrow(allfiles) > 0) { lrnfiles <- lrnfiles[!(lrnfiles %in% allfiles$filename)] }


# Process and write each one to DuckDB
for (i in seq_along(lrnfiles)) {
#  print(lrnfiles[i])
  
  df_lrn <- process_lrn_file(lrnfiles[i])
  if(nrow(df_lrn) > 0 ) {
  
  df_seqs = as.data.frame(cbind(seqid=df_lrn$seqid, fileid=fileids[i]))

  try(duckdb::duckdb_unregister(con, "new_seqs"), silent = TRUE)
  # df_seqs has columns: seqid, fileid
  duckdb::duckdb_register(con, "new_seqs", df_seqs)
  on.exit(duckdb::duckdb_unregister(con, "new_seqs"), add = TRUE)
  
  DBI::dbExecute(con, "
  INSERT INTO sequence_index BY NAME
  SELECT * FROM new_seqs
  ON CONFLICT DO NOTHING
")
  
#    result = dbAppendTable(con, "sequence_index", df_seqs)
    result = dbAppendTable(con, "tetra_data", df_lrn)
}
    result = dbExecute(con, "INSERT INTO import_log (filename) VALUES (?)", params = list(lrnfiles[i]))
}

#dbGetQuery(con, "SELECT * FROM tetra_data LIMIT 10")
#dbGetQuery(con, "SELECT COUNT(*) AS n FROM tetra_data")

dbDisconnect(con, shutdown = TRUE)

#source("/data/work/dana/sketch-db.r")






