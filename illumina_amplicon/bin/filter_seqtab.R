#!/usr/bin/env Rscript
#
# filter_seqtab.R — Quality-control filtering on long-format sequence data
#
# Applies a cascade of filters to remove noise and low-confidence observations.
# Operates entirely in long format (data.table) — never builds the full dense
# matrix, so memory use is proportional to non-zero entries.
#
# Filters (in order):
#   1. Length      — ASVs shorter than a threshold (can't assign taxonomy)
#   2. Prevalence  — ASVs in fewer than N samples (likely artefacts)
#   3. Abundance   — ASVs with fewer than N total reads (singletons/doubletons)
#   4. Depth       — samples with too few total reads (insufficient coverage)
#
# Usage:
#   filter_seqtab.R <seqtab_long.rds> <min_length> <min_samples> <min_seqs> <min_reads>
#
# Outputs:
#   seqtab_final.rds          Filtered long-format data.table
#   seqtab_final_wide.rds     Filtered dense matrix (for tools that need it)
#   seqtab_orphans.rds        ASVs removed by prevalence filter
#   seqtab_small.rds          Samples removed by depth filter
#   filter_stats.tsv          Per-step filtering statistics
#   sequence_summaries.pdf    Diagnostic rank-abundance plots

library(data.table)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: filter_seqtab.R <seqtab_long.rds> ",
         "<min_length> <min_samples> <min_seqs> <min_reads>")
}

input_path  <- args[1]
min_length  <- as.integer(args[2])
min_samples <- as.integer(args[3])
min_seqs    <- as.integer(args[4])
min_reads   <- as.integer(args[5])

# ---------------------------------------------------------------------------
# Load input (long-format data.table)
# ---------------------------------------------------------------------------
dt <- readRDS(input_path)

if (!is.data.table(dt)) {
    stop("[ERROR] Expected long-format data.table with columns: sample, sequence, count")
}

total_input  <- sum(dt$count)
n_input_asvs <- uniqueN(dt$sequence)

cat("[INFO] Input:", uniqueN(dt$sample), "samples,",
    n_input_asvs, "ASVs,", total_input, "reads\n")

# ---------------------------------------------------------------------------
# 1. Remove short sequences
#    ASVs shorter than ~50 bp cannot be reliably assigned taxonomy.
#    Note: some ITS amplicons can be legitimately short — adjust if needed.
# ---------------------------------------------------------------------------
seq_lengths   <- dt[, .(len = nchar(sequence[1])), by = sequence]
short_seqs    <- seq_lengths[len < min_length, sequence]
n_short       <- length(short_seqs)

if (n_short > 0) {
    cat("[INFO] Removing", n_short, "ASVs shorter than", min_length, "bp\n")
    dt <- dt[!sequence %in% short_seqs]
}

# ---------------------------------------------------------------------------
# 2. Remove low-prevalence ASVs (orphans)
#    ASVs appearing in very few samples are likely artefacts (tag-switching,
#    index bleed, or PCR chimeras missed by earlier steps).
# ---------------------------------------------------------------------------
seq_prevalence <- dt[, .(n_samples = uniqueN(sample)), by = sequence]
orphan_seqs    <- seq_prevalence[n_samples < min_samples, sequence]
n_orphans      <- length(orphan_seqs)

# Save orphans for review before removing
dt_orphans <- dt[sequence %in% orphan_seqs]
dt          <- dt[!sequence %in% orphan_seqs]

cat("[INFO] Removed", n_orphans, "orphan ASVs (present in <",
    min_samples, "samples)\n")

# ---------------------------------------------------------------------------
# 3. Remove low-abundance ASVs
#    ASVs with very few total reads are unreliable — they could be sequencing
#    errors that survived denoising, or real but under-sampled variants.
# ---------------------------------------------------------------------------
seq_abundance <- dt[, .(total = sum(count)), by = sequence]
rare_seqs     <- seq_abundance[total < min_seqs, sequence]
n_rare        <- length(rare_seqs)

dt <- dt[!sequence %in% rare_seqs]

cat("[INFO] Removed", n_rare, "rare ASVs (<", min_seqs, "total reads)\n")

# ---------------------------------------------------------------------------
# 4. Remove shallow samples
#    Samples with too few reads lack the sequencing depth needed for reliable
#    community profiling.
# ---------------------------------------------------------------------------
sample_depth  <- dt[, .(total = sum(count)), by = sample]
small_samples <- sample_depth[total < min_reads, sample]
n_small       <- length(small_samples)

# Save small samples for review before removing
dt_small <- dt[sample %in% small_samples]
dt        <- dt[!sample %in% small_samples]

cat("[INFO] Removed", n_small, "samples (<", min_reads, "reads)\n")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
n_final_samples <- uniqueN(dt$sample)
n_final_asvs    <- uniqueN(dt$sequence)
n_final_reads   <- sum(dt$count)
pct_retained    <- round(n_final_reads / max(total_input, 1) * 100, 1)

cat("[INFO] Final:", n_final_samples, "samples,", n_final_asvs, "ASVs,",
    n_final_reads, "reads (", pct_retained, "% of input)\n")

# ---------------------------------------------------------------------------
# Diagnostic plots — four-panel rank-abundance overview
# ---------------------------------------------------------------------------
pdf("sequence_summaries.pdf", width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))

if (nrow(dt) > 0) {
    reads_per_asv    <- dt[, .(total = sum(count)), by = sequence][order(total)]
    reads_per_sample <- dt[, .(total = sum(count)), by = sample][order(total)]
    samps_per_asv    <- dt[, .(n = uniqueN(sample)), by = sequence][order(n)]
    asvs_per_sample  <- dt[, .(n = uniqueN(sequence)), by = sample][order(n)]

    plot(log10(reads_per_asv$total),
         main = "Reads per ASV", ylab = "log10(reads)", xlab = "ASV rank",
         pch = 16, cex = 0.5, col = "steelblue")

    plot(log10(reads_per_sample$total),
         main = "Reads per Sample", ylab = "log10(reads)", xlab = "Sample rank",
         pch = 16, cex = 0.5, col = "steelblue")

    plot(log10(samps_per_asv$n),
         main = "Samples per ASV", ylab = "log10(samples)", xlab = "ASV rank",
         pch = 16, cex = 0.5, col = "darkgreen")

    plot(log10(asvs_per_sample$n),
         main = "ASVs per Sample", ylab = "log10(ASVs)", xlab = "Sample rank",
         pch = 16, cex = 0.5, col = "darkgreen")
} else {
    plot.new()
    text(0.5, 0.5, "No data remaining after filtering", cex = 1.5)
}
dev.off()

# ---------------------------------------------------------------------------
# Save outputs
#
# Primary output is long format. Also produce a wide matrix for any
# downstream tools that require it (e.g., SparCC), but this is secondary.
# ---------------------------------------------------------------------------
saveRDS(dt,         "seqtab_final.rds")
saveRDS(dt_orphans, "seqtab_orphans.rds")
saveRDS(dt_small,   "seqtab_small.rds")

# Wide matrix for compatibility
cat("[INFO] Casting final table to wide matrix...\n")
dt_wide    <- dcast(dt, sample ~ sequence, value.var = "count",
                    fill = 0L, fun.aggregate = sum)
sample_col <- dt_wide$sample
dt_wide[, sample := NULL]
seqtab_mat <- as.matrix(dt_wide)
rownames(seqtab_mat) <- sample_col
saveRDS(seqtab_mat, "seqtab_final_wide.rds")

rm(dt_wide, seqtab_mat)
gc()

# Stats
stats <- data.frame(
    step               = c("length", "prevalence", "abundance", "depth", "final"),
    asvs_removed       = c(n_short, n_orphans, n_rare, NA, NA),
    samples_removed    = c(NA, NA, NA, n_small, NA),
    remaining_samples  = c(NA, NA, NA, NA, n_final_samples),
    remaining_asvs     = c(NA, NA, NA, NA, n_final_asvs),
    pct_reads_retained = c(NA, NA, NA, NA, pct_retained)
)
write.table(stats, "filter_stats.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
