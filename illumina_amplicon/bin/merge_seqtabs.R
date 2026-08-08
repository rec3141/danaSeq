#!/usr/bin/env Rscript
#
# merge_seqtabs.R — Merge per-plate sequence tables using long-format approach
#
# Problem: dada2::mergeSequenceTables() builds a dense matrix with columns =
# union of all unique ASV sequences. With 4K+ samples and 100K+ ASVs, this
# matrix can exceed available memory.
#
# Solution: melt each plate's small seqtab to long format (sample, sequence,
# count), rbind them all (trivial — no column alignment), and aggregate any
# duplicate (sample, sequence) pairs. The output stays in long format as an
# RDS-serialized data.table, avoiding the dense matrix entirely.
#
# Downstream scripts (chimera removal, filtering, taxonomy) operate on the
# long format directly and only cast to wide when truly needed.
#
# Usage:
#   merge_seqtabs.R
#
# Inputs (discovered automatically):
#   *.seqtab.rds   One or more per-plate sequence tables in the working directory
#
# Outputs:
#   seqtab_merged.rds   Long-format data.table (sample, sequence, count)
#   merge_stats.tsv     Summary statistics

library(data.table)

# ---------------------------------------------------------------------------
# Discover per-plate sequence tables
# ---------------------------------------------------------------------------
rds_files <- sort(list.files(".", pattern = "\\.seqtab\\.rds$",
                             full.names = TRUE))

if (length(rds_files) == 0) {
    stop("[ERROR] No .seqtab.rds files found in working directory")
}

cat("[INFO] Found", length(rds_files), "sequence tables to merge\n")

# ---------------------------------------------------------------------------
# Melt each plate to long format and stack
#
# Each seqtab is a dense matrix (samples x ASVs). Most entries are zero —
# a typical plate with 96 samples and 2K ASVs might have only ~50K non-zero
# cells out of 192K. By extracting only non-zero entries, memory use during
# merge is proportional to actual data, not to the product of all samples
# and all unique sequences.
# ---------------------------------------------------------------------------
long_list <- vector("list", length(rds_files))

for (i in seq_along(rds_files)) {
    st <- readRDS(rds_files[i])
    cat("  ", basename(rds_files[i]), "->", nrow(st), "samples,",
        ncol(st), "ASVs,", sum(st), "reads\n")

    # Melt to long format, keeping only non-zero entries
    sample_names <- rownames(st)
    seq_names    <- colnames(st)
    plate_parts  <- vector("list", nrow(st))

    for (j in seq_len(nrow(st))) {
        nonzero <- st[j, ] > 0
        if (any(nonzero)) {
            plate_parts[[j]] <- data.table(
                sample   = sample_names[j],
                sequence = seq_names[nonzero],
                count    = as.integer(st[j, nonzero])
            )
        }
    }

    long_list[[i]] <- rbindlist(plate_parts)
    rm(st, plate_parts)
    gc()
}

cat("[INFO] Concatenating long tables...\n")
dt <- rbindlist(long_list)
rm(long_list)
gc()

# ---------------------------------------------------------------------------
# Aggregate: sum counts for (sample, sequence) pairs that appear on multiple
# plates (e.g., same sample re-sequenced). In most cases there are no
# duplicates and this is a no-op, but it's necessary for correctness.
# ---------------------------------------------------------------------------
n_before <- nrow(dt)
dt <- dt[, .(count = sum(count)), by = .(sample, sequence)]
n_dupes <- n_before - nrow(dt)
if (n_dupes > 0) {
    cat("[INFO] Aggregated", n_dupes, "duplicate (sample, sequence) entries\n")
}

n_samples <- uniqueN(dt$sample)
n_asvs    <- uniqueN(dt$sequence)
n_reads   <- sum(dt$count)

cat("[INFO] Merged total:", n_samples, "samples,",
    n_asvs, "unique ASVs,", n_reads, "reads,",
    nrow(dt), "non-zero entries\n")

# ---------------------------------------------------------------------------
# Save as long-format data.table (no dense matrix)
# ---------------------------------------------------------------------------
saveRDS(dt, "seqtab_merged.rds")

stats <- data.frame(
    step           = "merged",
    samples        = n_samples,
    sequences      = n_asvs,
    total_reads    = n_reads,
    nonzero_cells  = nrow(dt)
)
write.table(stats, "merge_stats.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
