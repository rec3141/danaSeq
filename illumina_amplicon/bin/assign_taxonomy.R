#!/usr/bin/env Rscript
#
# assign_taxonomy.R — Assign taxonomy to ASVs using a single reference database
#
# Runs dada2::assignTaxonomy() on a set of unique ASV sequences against one
# reference database. Each database is processed as an independent Nextflow
# task, allowing parallel taxonomy assignment across databases.
#
# assignTaxonomy() uses a naive Bayesian classifier (Wang et al. 2007) that
# breaks each query into overlapping k-mers and compares against the
# reference training set. Bootstrap values indicate classification confidence
# at each taxonomic rank.
#
# Usage:
#   assign_taxonomy.R <seqtab.rds> <ref_db.fasta> <db_name> <cpus> [tax_levels]
#
# Arguments:
#   seqtab.rds     Long-format data.table (sample, sequence, count) or
#                  a character vector of unique sequences
#   ref_db.fasta   Reference database FASTA (dada2 training set format)
#   db_name        Short name for the database (used in output filenames)
#   cpus           Number of threads for assignTaxonomy
#   tax_levels     Optional comma-separated taxonomy level names
#                  (e.g., "Domain,Phylum,Class,Order,Family,Genus")
#
# Outputs:
#   <db_name>_taxonomy.rds     Taxonomy matrix (ASVs x ranks)
#   <db_name>_bootstrap.rds    Bootstrap confidence matrix (ASVs x ranks)
#   <db_name>_taxonomy.tsv     Tab-delimited taxonomy for inspection

library(dada2)
library(data.table)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: assign_taxonomy.R <seqtab.rds> <ref_db> <db_name> <cpus> [tax_levels]")
}

input_path  <- args[1]
ref_db      <- args[2]
db_name     <- args[3]
cpus        <- as.integer(args[4])
tax_levels  <- if (length(args) >= 5 && args[5] != "null") {
    strsplit(args[5], ",")[[1]]
} else {
    NULL
}

# ---------------------------------------------------------------------------
# Validate reference database
# ---------------------------------------------------------------------------
if (!file.exists(ref_db)) {
    stop("[ERROR] Reference database not found: ", ref_db)
}
cat("[INFO] Reference database:", ref_db, "\n")
cat("[INFO] Database name:", db_name, "\n")

# ---------------------------------------------------------------------------
# Extract unique sequences from input
#
# Input can be either:
#   - Long-format data.table (from merge/chimera/filter pipeline)
#   - Character vector of sequences (from sequences_unique.rds)
# ---------------------------------------------------------------------------
input <- readRDS(input_path)

if (is.data.table(input)) {
    seqs <- unique(input$sequence)
    cat("[INFO] Extracted", length(seqs), "unique sequences from long-format table\n")
} else if (is.character(input)) {
    seqs <- input
    cat("[INFO] Loaded", length(seqs), "sequences\n")
} else if (is.matrix(input)) {
    # Wide-format matrix — columns are sequences
    seqs <- colnames(input)
    cat("[INFO] Extracted", length(seqs), "sequences from wide matrix\n")
} else {
    stop("[ERROR] Unrecognized input format")
}

if (length(seqs) == 0) {
    stop("[ERROR] No sequences to classify")
}

# ---------------------------------------------------------------------------
# Assign taxonomy
#
# UNITE databases use a special internal format in dada2 and should not have
# custom taxLevels passed. All other databases accept a taxLevels vector.
# minBoot=0 retains all classifications; bootstrap values are saved
# separately so downstream filtering can apply any threshold.
# ---------------------------------------------------------------------------
cat("[INFO] Classifying", length(seqs), "sequences against", db_name, "\n")

if (grepl("unite", db_name, ignore.case = TRUE) || is.null(tax_levels)) {
    tax_result <- assignTaxonomy(seqs, ref_db,
                                  multithread = cpus,
                                  minBoot = 0,
                                  outputBootstraps = TRUE,
                                  verbose = TRUE)
} else {
    cat("[INFO] Using taxonomy levels:", paste(tax_levels, collapse = ", "), "\n")
    tax_result <- assignTaxonomy(seqs, ref_db,
                                  multithread = cpus,
                                  minBoot = 0,
                                  outputBootstraps = TRUE,
                                  verbose = TRUE,
                                  taxLevels = tax_levels)
}

tax_matrix  <- tax_result$tax
boot_matrix <- tax_result$boot

cat("[INFO] Classification complete:", nrow(tax_matrix), "sequences,",
    ncol(tax_matrix), "ranks\n")

# ---------------------------------------------------------------------------
# Summary: how many sequences classified at each rank?
# ---------------------------------------------------------------------------
for (rank in colnames(tax_matrix)) {
    n_classified <- sum(!is.na(tax_matrix[, rank]))
    pct <- round(n_classified / nrow(tax_matrix) * 100, 1)
    cat("[INFO]", rank, ":", n_classified, "/", nrow(tax_matrix),
        "(", pct, "%) classified\n")
}

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
saveRDS(tax_matrix,  paste0(db_name, "_taxonomy.rds"))
saveRDS(boot_matrix, paste0(db_name, "_bootstrap.rds"))

# TSV for inspection (sequences as rownames get their own column)
tax_df <- data.frame(sequence = rownames(tax_matrix), tax_matrix,
                     check.names = FALSE, stringsAsFactors = FALSE)
write.table(tax_df, paste0(db_name, "_taxonomy.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
