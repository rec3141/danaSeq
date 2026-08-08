#!/usr/bin/env Rscript
#
# load_metadata.R — Load sample metadata and merge with sequence data
#
# Loads a MIMARKS-compliant (or custom) TSV/CSV metadata file and matches
# its rows to sample IDs found in the pipeline's sequence table. Reports
# matching statistics so the user can catch mismatches early.
#
# Usage:
#   load_metadata.R <seqtab.rds> <metadata_file> [sample_id_column]
#
# Arguments:
#   seqtab.rds         Long-format data.table (sample, sequence, count)
#   metadata_file      TSV or CSV with sample metadata
#   sample_id_column   Column in metadata matching sample IDs (default: "sample_name")
#
# Outputs:
#   metadata.rds       data.frame with rownames = sample IDs
#   match_stats.tsv    Matching summary between metadata and seqtab

library(data.table)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: load_metadata.R <seqtab.rds> <metadata_file> [sample_id_column]")
}

seqtab_path   <- args[1]
meta_path     <- args[2]
id_column     <- if (length(args) >= 3) args[3] else "sample_name"

# ---------------------------------------------------------------------------
# Load sequence table to extract sample IDs
# ---------------------------------------------------------------------------
dt <- readRDS(seqtab_path)
if (!is.data.table(dt)) {
    stop("[ERROR] Expected long-format data.table with columns: sample, sequence, count")
}

seqtab_samples <- unique(dt$sample)
cat("[INFO] Sequence table contains", length(seqtab_samples), "samples\n")

# ---------------------------------------------------------------------------
# Read metadata file (auto-detect TSV vs CSV by extension)
# ---------------------------------------------------------------------------
ext <- tolower(tools::file_ext(meta_path))

if (ext %in% c("tsv", "txt", "tab")) {
    cat("[INFO] Reading metadata as TSV:", meta_path, "\n")
    meta <- fread(meta_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
} else if (ext %in% c("csv")) {
    cat("[INFO] Reading metadata as CSV:", meta_path, "\n")
    meta <- fread(meta_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
} else {
    cat("[INFO] Unknown extension '", ext, "', trying TSV first\n", sep = "")
    meta <- tryCatch(
        fread(meta_path, sep = "\t", header = TRUE, stringsAsFactors = FALSE),
        error = function(e) {
            cat("[INFO] TSV parse failed, trying CSV\n")
            fread(meta_path, sep = ",", header = TRUE, stringsAsFactors = FALSE)
        }
    )
}

cat("[INFO] Metadata has", nrow(meta), "rows and", ncol(meta), "columns\n")

# ---------------------------------------------------------------------------
# Check for MIMARKS fields (informational only)
# ---------------------------------------------------------------------------
mimarks_fields <- c("sample_name", "collection_date", "geo_loc_name",
                     "lat_lon", "depth", "env_broad_scale",
                     "env_local_scale", "env_medium")

present <- mimarks_fields[mimarks_fields %in% colnames(meta)]
missing <- mimarks_fields[!mimarks_fields %in% colnames(meta)]

if (length(present) > 0) {
    cat("[INFO] MIMARKS fields found:", paste(present, collapse = ", "), "\n")
}
if (length(missing) > 0) {
    cat("[INFO] MIMARKS fields absent:", paste(missing, collapse = ", "), "\n")
}

# ---------------------------------------------------------------------------
# Match metadata rows to seqtab samples
# ---------------------------------------------------------------------------
if (!id_column %in% colnames(meta)) {
    stop("[ERROR] Sample ID column '", id_column,
         "' not found in metadata. Available columns: ",
         paste(colnames(meta), collapse = ", "))
}

meta_ids <- as.character(meta[[id_column]])
cat("[INFO] Using column '", id_column, "' for sample matching\n", sep = "")

matched      <- intersect(seqtab_samples, meta_ids)
in_seq_only  <- setdiff(seqtab_samples, meta_ids)
in_meta_only <- setdiff(meta_ids, seqtab_samples)

cat("[INFO] Matched:", length(matched), "samples\n")
cat("[INFO] In seqtab but not metadata:", length(in_seq_only), "samples\n")
cat("[INFO] In metadata but not seqtab:", length(in_meta_only), "samples\n")

if (length(in_seq_only) > 0 && length(in_seq_only) <= 20) {
    cat("[INFO] Missing metadata for:", paste(in_seq_only, collapse = ", "), "\n")
}

if (length(matched) == 0) {
    stop("[ERROR] No samples matched between metadata and seqtab. ",
         "Check that the '", id_column, "' column contains the correct sample IDs.")
}

# ---------------------------------------------------------------------------
# Build output metadata data.frame
# ---------------------------------------------------------------------------
meta_df <- as.data.frame(meta)

# Set rownames to sample IDs for matched rows
row_idx <- match(matched, meta_df[[id_column]])
meta_matched <- meta_df[row_idx, , drop = FALSE]
rownames(meta_matched) <- matched

cat("[INFO] Output metadata:", nrow(meta_matched), "samples x",
    ncol(meta_matched), "columns\n")

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
saveRDS(meta_matched, "metadata.rds")

stats <- data.frame(
    category    = c("matched", "seqtab_only", "metadata_only", "total_seqtab", "total_metadata"),
    count       = c(length(matched), length(in_seq_only), length(in_meta_only),
                    length(seqtab_samples), length(meta_ids)),
    stringsAsFactors = FALSE
)
write.table(stats, "match_stats.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("[INFO] Metadata loading complete\n")
