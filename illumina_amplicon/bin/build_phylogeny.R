#!/usr/bin/env Rscript
#
# build_phylogeny.R — Multiple sequence alignment and phylogenetic tree
#
# Aligns unique ASV sequences using DECIPHER::AlignSeqs() and builds a
# neighbor-joining tree. This provides phylogenetic context for diversity
# metrics (UniFrac) and can identify ASVs that cluster unexpectedly
# (potential contaminants or chimeras missed by earlier steps).
#
# For large datasets (>10K ASVs), alignment is the bottleneck. The script
# handles this by:
#   1. Using DECIPHER (profile-to-profile, faster than MUSCLE for many seqs)
#   2. Trimming alignment ends where coverage drops below 60%
#   3. Building NJ tree from raw + indel distance (handles gaps properly)
#
# Usage:
#   build_phylogeny.R <seqtab.rds> <cpus>
#
# Arguments:
#   seqtab.rds   Long-format data.table or character vector of sequences
#   cpus         Number of threads for alignment
#
# Outputs:
#   phylo_tree.rds         Phylogenetic tree (ape::phylo object)
#   phylo_alignment.fasta  Multiple sequence alignment
#   phylo_distances.rds    Distance matrix

library(data.table)
library(DECIPHER)
library(ape)
library(Biostrings)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: build_phylogeny.R <seqtab.rds> <cpus>")
}

input_path <- args[1]
cpus       <- as.integer(args[2])

# ---------------------------------------------------------------------------
# Extract unique sequences
# ---------------------------------------------------------------------------
input <- readRDS(input_path)

if (is.data.table(input)) {
    seqs <- unique(input$sequence)
} else if (is.character(input)) {
    seqs <- input
} else if (is.matrix(input)) {
    seqs <- colnames(input)
} else {
    stop("[ERROR] Unrecognized input format")
}

cat("[INFO] Building phylogeny for", length(seqs), "sequences\n")

if (length(seqs) < 3) {
    stop("[ERROR] Need at least 3 sequences to build a tree")
}

# Create short IDs for tree tips (full sequences are unwieldy as labels)
seq_ids <- paste0("ASV_", seq_along(seqs))
names(seqs) <- seq_ids

# ---------------------------------------------------------------------------
# Multiple sequence alignment with DECIPHER
#
# DECIPHER uses a profile-to-profile approach that scales better than
# progressive alignment (MUSCLE/ClustalW) for large numbers of sequences.
# It's also more accurate for divergent sequences typical of mixed-amplicon
# datasets (16S + 18S + ITS).
# ---------------------------------------------------------------------------
cat("[INFO] Aligning sequences...\n")
dna <- DNAStringSet(seqs)
names(dna) <- seq_ids

alignment <- AlignSeqs(dna, processors = cpus, verbose = TRUE)

# Write alignment to FASTA
writeXStringSet(alignment, "phylo_alignment.fasta")
cat("[INFO] Alignment written:", width(alignment)[1], "positions\n")

# ---------------------------------------------------------------------------
# Trim alignment ends where most sequences have gaps
#
# Ragged ends from variable-length amplicons add noise to distance
# calculations. Trim positions where <60% of sequences have data.
# ---------------------------------------------------------------------------
aln_matrix <- as.matrix(alignment)
coverage <- apply(aln_matrix, 2, function(x) mean(x != "-"))
keep_cols <- coverage >= 0.6

if (sum(keep_cols) < 50) {
    cat("[WARNING] Very few positions pass coverage filter, using full alignment\n")
    aln_trimmed <- as.DNAbin(alignment)
} else {
    cat("[INFO] Trimming alignment:", sum(!keep_cols), "low-coverage positions removed,",
        sum(keep_cols), "retained\n")
    aln_trimmed <- as.DNAbin(DNAStringSet(apply(aln_matrix[, keep_cols], 1,
                                                 paste, collapse = "")))
}

# ---------------------------------------------------------------------------
# Distance matrix and neighbor-joining tree
#
# Use raw distance (proportion of differing sites) + indel distance to
# properly account for insertions/deletions. Pairwise deletion handles
# remaining gaps without discarding entire sequences.
# ---------------------------------------------------------------------------
cat("[INFO] Computing distance matrix...\n")
dist_raw   <- dist.dna(aln_trimmed, model = "raw", pairwise.deletion = TRUE)
dist_indel <- dist.dna(aln_trimmed, model = "indel")

# Combined distance: substitutions + indels
dist_combined <- dist_raw + dist_indel

# Replace any NaN distances (from zero-overlap pairs) with max distance
dist_mat <- as.matrix(dist_combined)
dist_mat[is.nan(dist_mat)] <- max(dist_mat[!is.nan(dist_mat)], na.rm = TRUE)
dist_combined <- as.dist(dist_mat)

cat("[INFO] Building neighbor-joining tree...\n")
tree <- nj(dist_combined)

# Root at midpoint for consistent orientation
tree <- midpoint(tree)

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------
saveRDS(tree, "phylo_tree.rds")
saveRDS(dist_combined, "phylo_distances.rds")

# Also save the sequence-to-ID mapping for downstream use
seq_map <- data.table(asv_id = seq_ids, sequence = unname(seqs))
saveRDS(seq_map, "phylo_seq_map.rds")

cat("[INFO] Tree built:", Ntip(tree), "tips,",
    Nnode(tree), "internal nodes\n")
