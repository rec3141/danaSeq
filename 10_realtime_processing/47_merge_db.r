library(DBI)
library(duckdb)

setwd("/data/project_20250421X1/out_dana")

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

# Connect to the destination (merged) DB
con <- dbConnect(duckdb::duckdb(), dbdir = "merged.duckdb")

# List of source DBs
source_dbs <- list.files("/data/project_dbs", pattern = "\\.duckdb$", full.names = TRUE, recursive = TRUE)
tables_to_merge <- paste0("table", 1:5)

# For each table, create it once from the first DB
first_db <- source_dbs[[1]]
dbExecute(con, sprintf("ATTACH '%s' AS first;", first_db))

for (tbl in tables_to_merge) {
  dbExecute(con, sprintf("CREATE TABLE %s AS SELECT * FROM first.%s;", tbl, tbl))
}
dbExecute(con, "DETACH first;")

# For the rest of the DBs, attach and append
for (i in 2:length(source_dbs)) {
  alias <- paste0("src", i)
  db_path <- source_dbs[[i]]
  dbExecute(con, sprintf("ATTACH '%s' AS %s;", db_path, alias))
  
  for (tbl in tables_to_merge) {
    dbExecute(con, sprintf("INSERT INTO %s SELECT * FROM %s.%s;", tbl, alias, tbl))
  }
  
  dbExecute(con, sprintf("DETACH %s;", alias))
}
