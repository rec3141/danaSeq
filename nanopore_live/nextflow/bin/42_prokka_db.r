library(DBI)
library(duckdb)
library(readr)
library(stringr)

args <- commandArgs(trailingOnly = TRUE)
setwd(args[1])

# set up db connection
con <- dbConnect(duckdb::duckdb(), dbdir = "dana.duckdb")

# set up db schema
source("/work/apps/dana/log-db.r")

read_prokka_tsv <- function(path) {
  df <- read_tsv(
    path,
    col_types = cols(),
    show_col_types = FALSE,
    progress = FALSE
  )
  colnames(df) = tolower(colnames(df))
  
  df = df[df$ftype %in% c("CDS", "rRNA", "tRNA"), ]
  df
}

read_bakta_tsv <- function(path) {
  df <- read_tsv(
    path,
    col_types = cols(),
    show_col_types = FALSE,
    progress = FALSE,
    comment = "#"
  )
  colnames(df) <- tolower(gsub(" ", "_", colnames(df)))

  df <- df[df$type %in% c("cds", "tRNA", "rRNA", "CDS"), ]
  if (nrow(df) == 0) return(df)

  # Map Bakta columns to Prokka DuckDB schema:
  # Prokka schema: locus_tag, ftype, length_bp, gene, ec_number, cog, product
  result <- data.frame(
    locus_tag = df$locus_tag,
    ftype     = toupper(df$type),
    length_bp = abs(as.integer(df$stop) - as.integer(df$start)) + 1L,
    gene      = ifelse(is.na(df$gene), "", df$gene),
    ec_number = "",
    cog       = "",
    product   = ifelse(is.na(df$product), "", df$product),
    stringsAsFactors = FALSE
  )

  # Extract EC numbers from DbXrefs if available
  if ("dbxrefs" %in% colnames(df)) {
    ec <- str_match(df$dbxrefs, "EC:(\\S+)")
    result$ec_number <- ifelse(is.na(ec[, 2]), "", ec[, 2])
    cog_match <- str_match(df$dbxrefs, "COG:(\\S+)")
    result$cog <- ifelse(is.na(cog_match[, 2]), "", cog_match[, 2])
  }

  result
}

read_bakta_gff <- function(path) {
  lines <- readLines(path)
  len_lines <- lines[grepl("##sequence-region", lines)]
  if (length(len_lines) == 1)
    len_lines = paste(len_lines, "\n")
  len <- suppressWarnings(read_table(
    I(len_lines),
    col_names = c("X1", "seqid", "X3", "length"),
    show_col_types = FALSE
  ))
  len = len[, c("seqid", "length")]

  # Bakta GFF3 has ##FASTA section at end
  fasta_idx <- match("##FASTA", lines, nomatch = length(lines) + 1)
  gff_lines <- lines[seq_len(fasta_idx - 1)]
  gff <- read_tsv(
    I(gff_lines),
    comment = "#",
    col_names = FALSE,
    show_col_types = FALSE
  )
  if (nrow(gff) > 0) {
    gff <- gff[gff$X3 %in% c("CDS", "rRNA", "tRNA"), ]
    # Bakta uses ID= in attributes; also extract locus_tag= if present
    gff$locus_tag <- str_match(gff$X9, "locus_tag=([^;]+)")[, 2]
    if (all(is.na(gff$locus_tag))) {
      gff$locus_tag <- str_match(gff$X9, "ID=([^;]+)")[, 2]
    }
    gff <- gff[, c("X1", "locus_tag")]
    names(gff) <- c("seqid", "locus_tag")
  } else {
    gff <- data.frame(seqid = character(), locus_tag = character())
  }
  list(gff = gff, len = len)
}

read_prokka_gff <- function(path) {
  lines <- readLines(path)
  len_lines <- lines[grepl("##sequence-region", lines)]
  if (length(len_lines) == 1)
    len_lines = paste(len_lines, "\n")
  len <- suppressWarnings(read_table(
    I(len_lines),
    col_names = c("X1", "seqid", "X3", "length"),
    show_col_types = FALSE
  ))
  len = len[,c("seqid", "length")]
  
  gff_lines <- lines[!grepl("^##FASTA", lines)]
  gff_lines <- gff_lines[seq_len(match("##FASTA", lines, nomatch = length(lines) + 1) - 1)]
  # Now parse only the valid annotation part
  gff <- read_tsv(
    I(gff_lines),
    comment = "#",
    col_names = FALSE,
    show_col_types = FALSE
  )
  if (nrow(gff) > 0) {
    gff = gff[gff$X3 %in% c("CDS", "rRNA", "tRNA"), ]
    
    ##cutting out CRISPR lines
    # 20f2a8f0-c20b-4194-b7f3-bbb26e802bb8    minced:0.4.2    repeat_region   768     1357    .       .       .       note=CRISPR with 9 repeat units;rpt_family=CRISPR;rpt_type=direct;rpt_unit_seq=GAACCAACGCCCTGATTAAGAAGGGATTAAGAC
    
    gff$locus_tag <- str_match(gff$X9, "locus_tag=([^;]+)")[, 2]
    gff <- gff[, c("X1", "locus_tag")]
    names(gff) <- c("seqid", "locus_tag")
    
  } else {
    gff = data.frame(seqid = character(), locus_tag = character())
  }
  list(gff = gff, len = len)
}


# Loop and append
files_tsv <- list.files(
  "prokka",
  pattern = "\\.tsv$",
  full.names = TRUE,
  recursive = TRUE
)

allfiles = dbGetQuery(con,
                      "SELECT filename FROM import_log WHERE filename LIKE 'prokka%tsv'")

if (nrow(allfiles) > 0) {
  files_tsv <- files_tsv[!(files_tsv %in% allfiles$filename)]
}

files_gff = sub(".tsv", ".gff", files_tsv, fixed = T)

for (i in 1:length(files_tsv)) {
  #  print(files_tsv[i])
  df_tsv <- read_prokka_tsv(files_tsv[i])
  if (nrow(df_tsv) > 0) {
    df_gff <- read_prokka_gff(files_gff[i])
    
    invisible(dbAppendTable(con, "prokka_annotations", df_tsv))
    invisible(dbAppendTable(con, "locus_index", df_gff$gff))
    invisible(dbAppendTable(con, "stats", df_gff$len))
    
  }
  invisible(dbExecute(
    con,
    "INSERT INTO import_log (filename) VALUES (?)",
    params = list(files_tsv[i])
  ))
  invisible(dbExecute(
    con,
    "INSERT INTO import_log (filename) VALUES (?)",
    params = list(files_gff[i])
  ))
  
}


# Also scan bakta/ directories and import into the same tables
files_bakta_tsv <- list.files(
  "bakta",
  pattern = "\\.tsv$",
  full.names = TRUE,
  recursive = TRUE
)

allfiles_bakta <- dbGetQuery(con,
                      "SELECT filename FROM import_log WHERE filename LIKE 'bakta%tsv'")

if (nrow(allfiles_bakta) > 0) {
  files_bakta_tsv <- files_bakta_tsv[!(files_bakta_tsv %in% allfiles_bakta$filename)]
}

files_bakta_gff <- sub("\\.tsv$", ".gff3", files_bakta_tsv)

if (length(files_bakta_tsv) > 0) {
  for (i in 1:length(files_bakta_tsv)) {
    df_tsv <- read_bakta_tsv(files_bakta_tsv[i])
    if (nrow(df_tsv) > 0) {
      df_gff <- read_bakta_gff(files_bakta_gff[i])

      invisible(dbAppendTable(con, "prokka_annotations", df_tsv))
      invisible(dbAppendTable(con, "locus_index", df_gff$gff))
      invisible(dbAppendTable(con, "stats", df_gff$len))
    }
    invisible(dbExecute(
      con,
      "INSERT INTO import_log (filename) VALUES (?)",
      params = list(files_bakta_tsv[i])
    ))
    invisible(dbExecute(
      con,
      "INSERT INTO import_log (filename) VALUES (?)",
      params = list(files_bakta_gff[i])
    ))
  }
}

#dbGetQuery(con, "SELECT * FROM prokka_annotations LIMIT 10")
#dbGetQuery(con, "SELECT * FROM locus_index LIMIT 10")
#dbGetQuery(con, "SELECT COUNT(*) AS n FROM prokka_annotations")

dbDisconnect(con, shutdown = TRUE)
