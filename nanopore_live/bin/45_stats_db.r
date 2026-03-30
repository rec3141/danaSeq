suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(readr)
})


#setwd("/data/project_20250421X1/out_dana")

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

#open connection to db
con <- dbConnect(duckdb::duckdb(), dbdir = "dana.duckdb")

#setup schema
source("/work/apps/dana/log-db.r")

# Read the output, no header
ss_files = list.files("stats", pattern = "*.tsv", full.names = TRUE)

allfiles = dbGetQuery(con,
                      "SELECT filename FROM import_log WHERE filename LIKE 'stats%'")

if (nrow(allfiles) > 0) {
  ss_files <- ss_files[!(ss_files %in% allfiles$filename)]
}

for (file in ss_files) {
  #  print(file)
  #Name	Length	A	C	G	T	N	IUPAC	Other	GC
  
  ss <- read_tsv(
    file,
    comment = "",
    progress = FALSE,
    show_col_types = FALSE
  ) #
  if (nrow(ss) > 0) {
    ss = ss[, c(1, 2, 10)]
    colnames(ss) = c("seqid", "length", "gc")
    
    # Append to DuckDB
    result = dbAppendTable(con, "stats", ss)
  }
  result = dbExecute(con,
                     "INSERT INTO import_log (filename) VALUES (?)",
                     params = list(file))
}

#dbGetQuery(con, "SELECT * FROM stats LIMIT 10")
#dbGetQuery(con, "SELECT COUNT(*) AS n FROM stats")

dbDisconnect(con, shutdown = TRUE)
