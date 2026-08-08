#!/usr/bin/env Rscript
#
# cluster_tsne.R — t-SNE ordination of samples and ASVs via Bray-Curtis
#
# Computes Bray-Curtis distance matrices for both samples and ASVs,
# then reduces dimensionality with PCA followed by t-SNE. The PCA step
# handles the duplicate-coordinates issue that Rtsne cannot tolerate,
# and also speeds up the embedding for large datasets.
#
# Usage:
#   cluster_tsne.R <seqtab.rds> <cpus>
#
# Arguments:
#   seqtab.rds   Long-format data.table (sample, sequence, count)
#   cpus         Number of threads for parallel distance computation
#
# Outputs:
#   sample_bray_tsne.rds   t-SNE coordinates for samples (data.frame)
#   seq_bray_tsne.rds      t-SNE coordinates for ASVs (data.frame)
#   sample_bray_dist.rds   Bray-Curtis distance matrix for samples
#   seq_bray_dist.rds      Bray-Curtis distance matrix for ASVs

library(data.table)
library(parallelDist)
library(Rtsne)
library(gmodels)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: cluster_tsne.R <seqtab.rds> <cpus>")
}

seqtab_path <- args[1]
n_cpus      <- as.integer(args[2])

# ---------------------------------------------------------------------------
# Load input and cast to wide proportional matrix
# ---------------------------------------------------------------------------
dt <- readRDS(seqtab_path)
if (!is.data.table(dt)) {
    stop("[ERROR] Expected long-format data.table with columns: sample, sequence, count")
}

cat("[INFO] Input:", uniqueN(dt$sample), "samples,",
    uniqueN(dt$sequence), "ASVs\n")

# Cast to wide matrix (samples x sequences)
dt_wide    <- dcast(dt, sample ~ sequence, value.var = "count",
                    fill = 0L, fun.aggregate = sum)
sample_ids <- dt_wide$sample
dt_wide[, sample := NULL]
mat <- as.matrix(dt_wide)
rownames(mat) <- sample_ids

rm(dt_wide)
gc()

# Convert to proportions (row-wise)
row_totals <- rowSums(mat)
row_totals[row_totals == 0] <- 1  # avoid division by zero
prop_mat <- mat / row_totals

cat("[INFO] Proportional matrix:", nrow(prop_mat), "x", ncol(prop_mat), "\n")

# ---------------------------------------------------------------------------
# Helper: PCA + Rtsne on a distance matrix
#
# Rtsne requires unique input rows. We first project via PCA, then
# deduplicate, run t-SNE, and map coordinates back to all rows.
# ---------------------------------------------------------------------------
run_tsne <- function(dist_obj, labels, what = "items") {
    cat("[INFO] Running t-SNE for", what, "...\n")

    # Convert distance matrix to regular matrix for PCA
    dist_mat <- as.matrix(dist_obj)

    # PCA on distance matrix (classical MDS-like approach via PCA on distances)
    n_pcs <- min(50, nrow(dist_mat) - 1)
    pca_res <- fast.prcomp(dist_mat, center = TRUE, scale. = FALSE)
    pca_coords <- pca_res$x[, seq_len(n_pcs), drop = FALSE]

    # Deduplicate for Rtsne
    unique_coords <- unique(pca_coords)
    cat("[INFO]", what, ":", nrow(pca_coords), "total,",
        nrow(unique_coords), "unique PCA coordinates\n")

    # Perplexity must be < nrow/3
    perp <- min(30, floor((nrow(unique_coords) - 1) / 3))
    if (perp < 1) perp <- 1

    tsne_res <- Rtsne(unique_coords, dims = 2, perplexity = perp,
                      theta = 0.5, pca = FALSE, verbose = FALSE,
                      check_duplicates = FALSE)

    # Map back: find which unique row each original row corresponds to
    unique_key <- apply(unique_coords, 1, paste, collapse = "_")
    orig_key   <- apply(pca_coords, 1, paste, collapse = "_")
    map_idx    <- match(orig_key, unique_key)

    tsne_df <- data.frame(
        label = labels,
        tSNE1 = tsne_res$Y[map_idx, 1],
        tSNE2 = tsne_res$Y[map_idx, 2],
        stringsAsFactors = FALSE
    )

    cat("[INFO] t-SNE complete for", what, "\n")
    return(tsne_df)
}

# ---------------------------------------------------------------------------
# Sample ordination
# ---------------------------------------------------------------------------
cat("[INFO] Computing Bray-Curtis distances for samples (", n_cpus, " threads)...\n", sep = "")
sample_dist <- parDist(prop_mat, method = "bray", threads = n_cpus)

sample_tsne <- run_tsne(sample_dist, labels = rownames(prop_mat), what = "samples")

# ---------------------------------------------------------------------------
# ASV ordination (transpose — ASVs as rows)
# ---------------------------------------------------------------------------
cat("[INFO] Computing Bray-Curtis distances for ASVs (", n_cpus, " threads)...\n", sep = "")
asv_mat  <- t(prop_mat)
asv_dist <- parDist(asv_mat, method = "bray", threads = n_cpus)

asv_tsne <- run_tsne(asv_dist, labels = rownames(asv_mat), what = "ASVs")

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
saveRDS(sample_tsne, "sample_bray_tsne.rds")
saveRDS(asv_tsne,    "seq_bray_tsne.rds")
saveRDS(sample_dist, "sample_bray_dist.rds")
saveRDS(asv_dist,    "seq_bray_dist.rds")

cat("[INFO] Clustering complete\n")
