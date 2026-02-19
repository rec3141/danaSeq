suppressPackageStartupMessages({
  library(DBI); library(duckdb); library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])
#setwd("/data/project_20250502/20250411/out_dana_bc/FAY98439/barcode11/")

#open connection to db
con <- dbConnect(duckdb::duckdb(), dbdir = "dana.duckdb")

#setup schema
source("/work/apps/dana/log-db.r")

#LOAD AND COMBINE KRAKEN REPORTS WITH INDENT-BASED TAXONOMY
read_kraken_report <- function(file) {
  lines <- read_lines(file)
  df <- tibble(raw = lines) %>%
    separate(raw, into = c("percent", "reads", "direct", "rank", "taxid", "name"), sep = "\t", fill = "right") %>%
    mutate(level = str_extract(name, "^ *") %>% str_length(),
           name = str_trim(name),
           flowcell = str_split(basename(file), "_", simplify = TRUE)[, 1],
           barcode = str_split(basename(file), "_", simplify = TRUE)[, 3])
  
  # Reconstruct full taxonomy
  taxonomy_stack <- character()
  taxonomy_full <- character(nrow(df))
  rank_stack <- character()
  rank_full <- character(nrow(df))
  
  for (i in seq_len(nrow(df))) {
    lvl <- df$level[i] %/% 2 + 1
    taxonomy_stack[lvl] <- df$name[i]
    taxonomy_stack <- taxonomy_stack[seq_len(lvl)]
    rank_stack[lvl] <- df$rank[i]
    rank_stack <- rank_stack[seq_len(lvl)]
    taxonomy_full[i] <- paste(taxonomy_stack, collapse = "; ")
    rank_full[i] <- paste(rank_stack, collapse = " ")
  }
  
  df$taxonomy <- taxonomy_full
  df$rank_full <- rank_full
  df$percent <- as.numeric(df$percent)
  return(df)
}

# Read the output, no header
kraken_files = list.files("kraken",pattern = "*.report",full.names = TRUE)

allfiles = dbGetQuery(con, "SELECT filename FROM import_log WHERE filename LIKE 'kraken%report'")

if(nrow(allfiles) > 0) { kraken_files <- kraken_files[!(kraken_files %in% allfiles$filename)] }

for(file in kraken_files) {
   # print(file)
  kraken_in <- read_kraken_report(file)

# would like to have each rank in it's own column but too complicated for now with all of the different ranks
#     kraken_in <- read_kraken_report(file) %>% mutate(
#     rank_values = map2(taxonomy, rank_full, ~ set_names(str_split(.x, " ")[[1]], .y))
#   )

#  kraken_in <- read_tsv(file, col_names = c("seqid", "taxa_w_id"), comment = "#", progress = FALSE, show_col_types = FALSE, col_select = c(1,2))
  if(nrow(kraken_in) > 0) { 

    # Append to DuckDB
  result = dbAppendTable(con, "krakenreport", kraken_in)
  }
  result = dbExecute(con, "INSERT INTO import_log (filename) VALUES (?)", params = list(file))
}

#dbGetQuery(con, "SELECT * FROM sendsketch LIMIT 10")
#dbGetQuery(con, "SELECT COUNT(*) AS n FROM sendsketch")

dbDisconnect(con, shutdown = TRUE)

#source("/data/work/dana/prokka-db.r")

