#!/usr/bin/env Rscript
#
# network_sparcc.R — SparCC correlation network from normalized count data
#
# Runs SparCC (Sparse Correlations for Compositional data) to infer
# microbial co-occurrence networks. SparCC is designed specifically for
# compositional data and avoids the spurious correlations that plague
# standard correlation methods applied to relative abundances.
#
# Low-prevalence ASVs are filtered out before analysis to reduce noise
# and computation time — rare ASVs lack the statistical power for
# meaningful correlation estimates.
#
# Usage:
#   network_sparcc.R <renorm_merged.rds> <min_prevalence> <cpus>
#
# Arguments:
#   renorm_merged.rds   Long-format data.table (group, sample, sequence, count, proportion)
#   min_prevalence      Minimum number of samples an ASV must appear in
#   cpus                Threads for parallel computation
#
# Outputs:
#   sparcc_correlations.rds   Melted edge list (node1, node2, correlation, weight, color)
#   sparcc_stats.tsv          Summary statistics

library(data.table)
library(SpiecEasi)
library(Matrix)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: network_sparcc.R <renorm_merged.rds> <min_prevalence> <cpus>")
}

renorm_path    <- args[1]
min_prevalence <- as.integer(args[2])
n_cpus         <- as.integer(args[3])

# ---------------------------------------------------------------------------
# Load input
# ---------------------------------------------------------------------------
dt <- readRDS(renorm_path)
if (!is.data.table(dt)) {
    stop("[ERROR] Expected long-format data.table with columns: sample, sequence, count")
}

n_total_asvs <- uniqueN(dt$sequence)
n_samples    <- uniqueN(dt$sample)
cat("[INFO] Input:", n_samples, "samples,", n_total_asvs, "ASVs\n")

# ---------------------------------------------------------------------------
# Filter ASVs by prevalence
# ---------------------------------------------------------------------------
asv_prev <- dt[count > 0, .(prevalence = uniqueN(sample)), by = sequence]
keep_seqs <- asv_prev[prevalence >= min_prevalence, sequence]
n_kept <- length(keep_seqs)

cat("[INFO] Prevalence filter (>=", min_prevalence, "samples):",
    n_kept, "of", n_total_asvs, "ASVs retained\n")

if (n_kept < 3) {
    stop("[ERROR] Too few ASVs (", n_kept,
         ") pass the prevalence filter. Lower min_prevalence.")
}

dt_sub <- dt[sequence %in% keep_seqs]

# ---------------------------------------------------------------------------
# Cast to wide integer count matrix (SparCC needs counts, not proportions)
# ---------------------------------------------------------------------------
dt_wide    <- dcast(dt_sub, sample ~ sequence, value.var = "count",
                    fill = 0L, fun.aggregate = sum)
sample_ids <- dt_wide$sample
dt_wide[, sample := NULL]
count_mat <- as.matrix(dt_wide)
rownames(count_mat) <- sample_ids
storage.mode(count_mat) <- "integer"

cat("[INFO] Count matrix:", nrow(count_mat), "samples x",
    ncol(count_mat), "ASVs\n")

rm(dt_wide)
gc()

# ---------------------------------------------------------------------------
# Run SparCC
# ---------------------------------------------------------------------------
cat("[INFO] Running SparCC...\n")
sparcc_res <- sparcc(count_mat)

cor_mat <- sparcc_res$Cor
colnames(cor_mat) <- colnames(count_mat)
rownames(cor_mat) <- colnames(count_mat)

cat("[INFO] SparCC complete. Correlation matrix:",
    nrow(cor_mat), "x", ncol(cor_mat), "\n")

# ---------------------------------------------------------------------------
# Melt to edge list and filter by |correlation| > 0.1
# ---------------------------------------------------------------------------
# Use data.table melt on upper triangle only (avoid duplicates)
cor_dt <- as.data.table(cor_mat, keep.rownames = "node1")
cor_long <- melt(cor_dt, id.vars = "node1",
                 variable.name = "node2", value.name = "correlation")
cor_long[, node2 := as.character(node2)]

# Keep only upper triangle (node1 < node2 lexicographically)
cor_long <- cor_long[node1 < node2]

# Filter by absolute correlation
cor_long <- cor_long[abs(correlation) > 0.1]

cat("[INFO] Edges with |correlation| > 0.1:", nrow(cor_long), "\n")

# ---------------------------------------------------------------------------
# Compute edge weights and colors
# ---------------------------------------------------------------------------
# Weight: cube of correlation (preserves sign emphasis on strong edges)
cor_long[, weight := abs(correlation)^3]

# Normalize weights to [0, 1]
max_weight <- max(cor_long$weight, na.rm = TRUE)
if (max_weight > 0) {
    cor_long[, weight := weight / max_weight]
}

# Color: blue for positive, red for negative
cor_long[, color := ifelse(correlation > 0, "blue", "red")]

# ---------------------------------------------------------------------------
# Summary statistics
# ---------------------------------------------------------------------------
n_edges     <- nrow(cor_long)
n_positive  <- sum(cor_long$correlation > 0)
n_negative  <- sum(cor_long$correlation < 0)
mean_abs    <- round(mean(abs(cor_long$correlation)), 4)
max_cor     <- round(max(cor_long$correlation), 4)
min_cor     <- round(min(cor_long$correlation), 4)

cat("[INFO] Network summary:\n")
cat("[INFO]   Edges:", n_edges, "(", n_positive, "positive,",
    n_negative, "negative )\n")
cat("[INFO]   Mean |correlation|:", mean_abs, "\n")
cat("[INFO]   Range:", min_cor, "to", max_cor, "\n")

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
saveRDS(cor_long, "sparcc_correlations.rds")

stats <- data.frame(
    metric = c("input_asvs", "filtered_asvs", "samples", "min_prevalence",
               "edges", "positive_edges", "negative_edges",
               "mean_abs_correlation", "max_correlation", "min_correlation"),
    value  = c(n_total_asvs, n_kept, n_samples, min_prevalence,
               n_edges, n_positive, n_negative,
               mean_abs, max_cor, min_cor),
    stringsAsFactors = FALSE
)
write.table(stats, "sparcc_stats.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("[INFO] SparCC network analysis complete\n")
