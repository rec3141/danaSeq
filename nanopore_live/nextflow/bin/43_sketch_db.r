suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(readr)
})

#setwd("/data/project_boston/out_dana")

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

#open connection to db
con <- dbConnect(duckdb::duckdb(), dbdir = "dana.duckdb")

#setup schema
source("/work/apps/dana/log-db.r")

# Read the output, no header
ss_files = list.files("sketch",pattern = "*.txt",full.names = TRUE)

allfiles = dbGetQuery(con, "SELECT filename FROM import_log WHERE filename LIKE 'sketch%'")

if(nrow(allfiles) > 0) { ss_files <- ss_files[!(ss_files %in% allfiles$filename)] }

for(file in ss_files) {
#  print(file)
  ss <- read_tsv(file, col_names = c("fileid", "ref_name", "ani"), comment = "#", progress = FALSE, show_col_types = FALSE, col_select = c(1,2,3))
  if(nrow(ss) > 0) {

ss$fileid = sub(".fa","",ss$fileid, fixed = TRUE)

# Append to DuckDB
result = dbAppendTable(con, "sendsketch", ss)
}
result = dbExecute(con, "INSERT INTO import_log (filename) VALUES (?)", params = list(file))
}

#dbGetQuery(con, "SELECT * FROM sendsketch LIMIT 10")
#dbGetQuery(con, "SELECT COUNT(*) AS n FROM sendsketch")

dbDisconnect(con, shutdown = TRUE)

#source("/data/work/dana/prokka-db.r")

