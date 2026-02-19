#!/usr/bin/env Rscript
# =============================================================================
# MAG Assembly Result Visualizations
# =============================================================================
#
# Usage:
#   Rscript plot_mag_results.R <results_dir>
#
# Example:
#   Rscript plot_mag_results.R nextflow/results_lorbin_test
#
# Outputs figures to <results_dir>/analysis/figures/
#
# Figures produced:
#   1. jaccard_heatmap.png       Jaccard dissimilarity among binners
#   2. tsne_gc.png               t-SNE of contigs colored by GC content
#   3. tsne_binners.png          t-SNE per binner colored by bin ID
#   4. coverage_heatmap.png      DAS Tool bin coverage per sample
# =============================================================================

suppressPackageStartupMessages({
  library(Rtsne)
  library(viridis)
  library(scales)
  library(pheatmap)
})
source("~/apps/FastPCA.R")

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat("Usage: Rscript plot_mag_results.R <results_dir>\n")
  quit(status = 1)
}

results_dir <- args[1]
binning_dir <- file.path(results_dir, "binning")
mapping_dir <- file.path(results_dir, "mapping")
assembly_dir <- file.path(results_dir, "assembly")
fig_dir     <- file.path(results_dir, "analysis", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

binners <- c("semibin", "metabat", "maxbin", "lorbin", "comebin")
labels  <- c("SemiBin2", "MetaBAT2", "MaxBin2", "LorBin", "COMEBin")

# =============================================================================
# Load shared data
# =============================================================================

# ---- Depth matrix -----------------------------------------------------------
depths_file <- file.path(mapping_dir, "depths.txt")
if (!file.exists(depths_file)) stop("depths.txt not found: ", depths_file)

depths_raw <- read.delim(depths_file, header = TRUE, check.names = FALSE)
contig_names <- depths_raw[[1]]
contig_lens  <- depths_raw[[2]]

# Extract depth columns (every other column starting at 4, skip variance)
depth_cols <- seq(4, ncol(depths_raw), by = 2)
depth_mat  <- as.matrix(depths_raw[, depth_cols])
rownames(depth_mat) <- contig_names

# Clean up sample names: strip .sorted.bam suffix
sample_names <- gsub("\\.sorted\\.bam$", "", colnames(depth_mat))
colnames(depth_mat) <- sample_names

cat("[INFO] Loaded depth matrix:", nrow(depth_mat), "contigs x",
    ncol(depth_mat), "samples\n")

# ---- GC content (approximated from TNF after loading) ----------------------
# Deferred until TNF matrix is loaded; see below

# ---- Bin assignments per binner ---------------------------------------------
bin_assignments <- list()
bin_keep <- c()
for (i in seq_along(binners)) {
  path <- file.path(binning_dir, binners[i], "contig_bins.tsv")
  if (!file.exists(path) || file.size(path) == 0) next
  df <- read.delim(path, header = FALSE, col.names = c("contig", "bin"))
  bin_assignments[[labels[i]]] <- df
  bin_keep <- c(bin_keep, labels[i])
}

# ---- DAS Tool consensus assignments ----------------------------------------
dastool_file <- file.path(binning_dir, "dastool", "contig2bin.tsv")
dastool_bins <- NULL
if (file.exists(dastool_file) && file.size(dastool_file) > 0) {
  dastool_bins <- read.delim(dastool_file, header = FALSE,
                             col.names = c("contig", "bin"))
  # Add DAS Tool to binner list for t-SNE panels
  bin_assignments[["DAS Tool"]] <- dastool_bins
  bin_keep <- c(bin_keep, "DAS Tool")
}

# ---- TNF matrix ------------------------------------------------------------
tnf_file <- file.path(assembly_dir, "tnf.tsv")
if (!file.exists(tnf_file)) stop("tnf.tsv not found: ", tnf_file)

tnf_raw <- read.delim(tnf_file, header = FALSE)
tnf_names <- tnf_raw[[1]]
tnf_mat   <- as.matrix(tnf_raw[, -1])
rownames(tnf_mat) <- tnf_names
cat("[INFO] Loaded TNF matrix:", nrow(tnf_mat), "contigs x",
    ncol(tnf_mat), "features\n")

# Align TNF to depth matrix contig order
tnf_idx <- match(contig_names, tnf_names)
tnf_aligned <- tnf_mat[tnf_idx, ]
rownames(tnf_aligned) <- contig_names

# ---- GC content from TNF (dot product) ------------------------------------
# Each canonical tetramer has a fixed GC fraction; weight by it
# Generate the 136 canonical tetranucleotides in the same order as tetramer_freqs.py
bases <- c("A", "C", "G", "T")
comp  <- c(A = "T", C = "G", G = "C", T = "A")
all_4mers <- apply(expand.grid(bases, bases, bases, bases), 1,
                   function(x) paste(x, collapse = ""))
rc_4mers  <- sapply(all_4mers, function(k)
                    paste(comp[rev(strsplit(k, "")[[1]])], collapse = ""))
seen <- character(0)
canonical <- character(0)
for (k in all_4mers) {
  rc <- rc_4mers[k]
  if (!(rc %in% seen)) {
    canonical <- c(canonical, k)
    seen <- c(seen, k)
  }
}
gc_weights <- sapply(canonical, function(k)
                     sum(strsplit(k, "")[[1]] %in% c("G", "C")) / 4)
gc_pct <- as.numeric(tnf_aligned %*% gc_weights)
names(gc_pct) <- contig_names
cat("[INFO] GC content approximated from TNF for", sum(!is.na(gc_pct)), "contigs\n")

# ---- t-SNE of TNF via FastPCA ----------------------------------------------
tsne_cache <- file.path(results_dir, "analysis", "tsne.rds")
if (file.exists(tsne_cache)) {
  cat("[INFO] Loading cached t-SNE from", tsne_cache, "\n")
  tsne_res <- readRDS(tsne_cache)
} else {
  cat("[INFO] FastPCA on TNF + length matrix...\n")
  feat <- cbind(tnf_aligned^.25, log10(contig_lens))
  tnf_pca <- FastPCA(feat, 50)
  cat("[INFO] Running t-SNE on TNF PCA (pca=F)...\n")
  set.seed(42)
  tsne_res <- Rtsne(tnf_pca$x, dims = 2, pca = FALSE,
                    perplexity = min(30, nrow(tnf_pca$x) / 4),
                    max_iter = 2500, check_duplicates = FALSE, num_threads = 8)
  saveRDS(tsne_res, tsne_cache)
  cat("[INFO] t-SNE complete, cached to", tsne_cache, "\n")
}
tsne_xy <- tsne_res$Y

# =============================================================================
# Figure 1: Jaccard dissimilarity heatmap
# =============================================================================

contig_sets <- lapply(bin_assignments, function(df) unique(df$contig))

# Jaccard compares individual binners only (not the DAS Tool consensus)
jac_keep <- setdiff(bin_keep, "DAS Tool")
if (length(jac_keep) >= 2) {
  n  <- length(jac_keep)
  jd <- matrix(0, n, n, dimnames = list(jac_keep, jac_keep))
  for (i in 1:n) for (j in 1:n) {
    a <- contig_sets[[jac_keep[i]]]
    b <- contig_sets[[jac_keep[j]]]
    inter   <- length(intersect(a, b))
    uni     <- length(union(a, b))
    jd[i,j] <- 1 - inter / uni
  }

  outfile <- file.path(fig_dir, "jaccard_heatmap.png")
  png(outfile, width = 700, height = 600, res = 150)
  par(mar = c(5, 5, 3, 6))
  cols <- colorRampPalette(c("#FFFFCC", "#FD8D3C", "#BD0026"))(100)
  image(1:n, 1:n, t(jd[n:1, ]), col = cols, zlim = c(0, 1),
        axes = FALSE, xlab = "", ylab = "")
  axis(1, 1:n, jac_keep, las = 2, cex.axis = 0.9)
  axis(2, 1:n, rev(jac_keep), las = 1, cex.axis = 0.9)
  for (i in 1:n) for (j in 1:n) {
    val <- jd[j, i]
    text(i, n - j + 1, sprintf("%.3f", val), cex = 0.85,
         col = ifelse(val > 0.5, "white", "black"))
  }
  title("Jaccard Dissimilarity of Contig Assignments", cex.main = 1.1)
  par(xpd = TRUE)
  legend_y <- seq(1, n, length.out = 5)
  legend_v <- sprintf("%.1f", seq(0, 1, 0.25))
  mtext("Dissimilarity", side = 4, line = 4.5, cex = 0.8)
  rasterImage(as.raster(matrix(rev(cols), ncol = 1)),
              n + 0.4, 1, n + 0.7, n)
  for (k in 1:5) text(n + 0.9, legend_y[k], legend_v[k], cex = 0.65)
  dev.off()
  cat("[OK]  ", outfile, "\n")
} else {
  cat("[SKIP] Jaccard heatmap: need >= 2 individual binners with output\n")
}

# =============================================================================
# Figure 2: t-SNE colored by GC content (dark background)
# =============================================================================

if (!is.null(gc_pct)) {
  outfile <- file.path(fig_dir, "tsne_gc.png")
  png(outfile, width = 1080, height = 1080, units = "px", bg = "black", res = 150)
  par(bg = "black", mar = c(3, 3, 3, 5), fg = "white",
      col.axis = "white", col.lab = "white", col.main = "white")

  gc_cols <- viridis(256)
  gc_idx  <- pmin(256, pmax(1, round(rescale(gc_pct, to = c(1, 256)))))
  gc_idx[is.na(gc_idx)] <- 1
  pt_col  <- alpha(gc_cols[gc_idx], 0.7)

  # Size by contig length
  pt_cex <- rescale(sqrt(contig_lens), to = c(0.3, 2.5))

  plot(tsne_xy, pch = 16, cex = pt_cex, col = pt_col,
       axes = FALSE, xlab = "", ylab = "",
       main = "t-SNE of Contigs by GC Content")
  axis(1, col = "gray40", col.ticks = "gray40", labels = FALSE)
  axis(2, col = "gray40", col.ticks = "gray40", labels = FALSE)

  # Color bar
  par(xpd = TRUE)
  usr <- par("usr")
  bar_x0 <- usr[2] + (usr[2] - usr[1]) * 0.02
  bar_x1 <- usr[2] + (usr[2] - usr[1]) * 0.05
  bar_y  <- seq(usr[3], usr[4], length.out = 257)
  for (k in 1:256) {
    rect(bar_x0, bar_y[k], bar_x1, bar_y[k + 1],
         col = gc_cols[k], border = NA)
  }
  gc_range <- range(gc_pct, na.rm = TRUE)
  tick_vals <- pretty(gc_range, n = 5)
  tick_y <- rescale(tick_vals, from = gc_range, to = c(usr[3], usr[4]))
  text(bar_x1 + (usr[2] - usr[1]) * 0.02, tick_y,
       sprintf("%.0f%%", tick_vals * 100), col = "white", cex = 0.6, adj = 0)
  text(bar_x0 + (bar_x1 - bar_x0) / 2, usr[4] + (usr[4] - usr[3]) * 0.02,
       "GC%", col = "white", cex = 0.7)

  dev.off()
  cat("[OK]  ", outfile, "\n")
} else {
  cat("[SKIP] t-SNE GC plot: no GC data\n")
}

# =============================================================================
# Figure 3: t-SNE per binner colored by bin ID (dark background)
# =============================================================================

if (length(bin_keep) >= 1) {
  n_binners <- length(bin_keep)
  ncols <- min(n_binners, 3)
  nrows <- ceiling(n_binners / ncols)
  cm = 20

  outfile <- file.path(fig_dir, "tsne_binners.png")
  png(outfile, width = ncols * cm, height = nrows * cm,
      units = "cm", bg = "black", res=300)
  par(mfrow = c(nrows, ncols), bg = "black",
      mar = c(1, 1, 2.5, 1), oma = c(0, 0, 2, 0),
      fg = "white", col.main = "white")

  for (bname in bin_keep) {
    df <- bin_assignments[[bname]]
    # Map contigs to bin IDs
    bin_map <- setNames(df$bin, df$contig)
    assigned <- contig_names %in% df$contig
    bin_ids  <- bin_map[contig_names]

    # Assign colors per bin
    unique_bins <- sort(unique(bin_ids[!is.na(bin_ids)]))
    n_bins <- length(unique_bins)
    if (n_bins > 0) {
      bin_pal <- setNames(
        alpha(sample(hue_pal()(min(n_bins, 256)), n_bins, replace = n_bins > 256), 0.8),
        unique_bins
      )
    }

    # Background (unassigned) contigs
    plot(tsne_xy, pch = 16, cex = 0.15, col = alpha("white", 0.5),
         axes = FALSE, xlab = "", ylab = "",
         main = paste0(bname, " (", n_bins, " bins)"))

    # Overlay assigned contigs
    if (n_bins > 0) {
      idx <- which(assigned)
      points(tsne_xy[idx, , drop = FALSE], pch = 16, cex = 0.4,
             col = bin_pal[bin_ids[idx]])
    }
  }

  # Fill empty panels if grid is not full
  remaining <- ncols * nrows - n_binners
  if (remaining > 0) for (r in seq_len(remaining)) plot.new()

  mtext("t-SNE of Contigs by Binner Assignment", outer = TRUE,
        col = "white", cex = 1.1, line = 0.3)
  dev.off()
  cat("[OK]  ", outfile, "\n")
} else {
  cat("[SKIP] t-SNE binner panels: no binner output\n")
}

# =============================================================================
# Figure 4: Coverage heatmap (DAS Tool bins x samples)
# =============================================================================

if (!is.null(dastool_bins) && nrow(dastool_bins) > 0) {
  # Map contig depths to DAS Tool bins
  bin_map <- setNames(dastool_bins$bin, dastool_bins$contig)
  matched <- contig_names %in% dastool_bins$contig
  bin_ids <- bin_map[contig_names[matched]]
  depth_matched <- depth_mat[matched, , drop = FALSE]

  # Aggregate: sum depth per bin per sample
  bin_depth <- rowsum(depth_matched, group = bin_ids)

  # Normalize to relative abundance per sample (column proportions)
  bin_prop <- sweep(bin_depth, 2, colSums(bin_depth), "/")

  # Filter bins with some signal
  keep_bins <- rowSums(bin_depth) > 0
  bin_prop <- bin_prop[keep_bins, , drop = FALSE]

  if (nrow(bin_prop) >= 2 && ncol(bin_prop) >= 2) {
    outfile <- file.path(fig_dir, "coverage_heatmap.png")
    # Scale height with number of bins
    fig_h <- max(600, nrow(bin_prop) * 14 + 200)
    png(outfile, width = max(600, ncol(bin_prop) * 120 + 300),
        height = fig_h, res = 150)
    pheatmap(
      log1p(bin_depth[keep_bins, , drop = FALSE]),
      color = viridis(100),
      clustering_distance_rows = "euclidean",
      clustering_distance_cols = "euclidean",
      clustering_method = "ward.D2",
      fontsize_row = max(4, min(8, 400 / nrow(bin_prop))),
      fontsize_col = 9,
      main = "DAS Tool Bin Coverage (log depth)",
      border_color = NA
    )
    dev.off()
    cat("[OK]  ", outfile, "\n")
  } else {
    cat("[SKIP] Coverage heatmap: not enough bins or samples\n")
  }
} else {
  cat("[SKIP] Coverage heatmap: no DAS Tool contig2bin.tsv\n")
}

cat("\n[DONE] Figures saved to:", fig_dir, "\n")
