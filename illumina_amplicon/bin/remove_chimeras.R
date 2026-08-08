#!/usr/bin/env Rscript
#
# remove_chimeras.R — Sparse consensus chimera removal
#
# A memory-efficient reimplementation of dada2's removeBimeraDenovo() that
# operates on long-format count data instead of a dense sample-by-ASV matrix.
#
# Why: The standard dada2 approach materializes the full dense matrix
# (samples x ASVs), which can exceed available memory on large datasets
# (e.g., 4K samples x 180K ASVs = ~3 GB). This implementation uses
# data.table for sparse lookups and only touches non-zero entries.
#
# Algorithm (identical to dada2 "consensus" method):
#   For each query ASV (in order of increasing total abundance):
#     1. Find all samples where the query is present
#     2. In each sample, find candidate parents: ASVs with abundance >=
#        minFoldParentOverAbundance * query_abundance AND >= minParentAbundance
#     3. Check if the query is a bimera of any two parents using dada2::isBimera()
#        (calls the same C-level Needleman-Wunsch alignment as the original)
#     4. Tally: if flagged in >= minSampleFraction of samples → chimera
#
# The per-ASV checks are independent and parallelized with mclapply.
#
# Usage:
#   remove_chimeras.R <seqtab.rds> <cpus>
#
# Outputs:
#   seqtab_nochim.rds    Chimera-free sequence table (long-format data.table)
#   chimera_stats.tsv    Summary statistics

library(dada2)
library(data.table)
library(parallel)

# ---------------------------------------------------------------------------
# Sparse consensus chimera detection
#
# This function mirrors dada2::isBimeraDenovoTable() but works on a
# data.table in long format (sample, sequence, count) instead of a dense
# matrix. It reuses dada2's C-level isBimera() for the actual alignment
# and chimera check, so the biological results are identical.
#
# Arguments:
#   dt           data.table with columns: sample, sequence, count (non-zero only)
#   minFoldParentOverAbundance  Parent must be this fold more abundant than
#                               query in the same sample [default: 1.5]
#   minParentAbundance  Minimum absolute count for a parent [default: 2]
#   allowOneOff  Also flag sequences one mismatch from an exact bimera [default: FALSE]
#   minOneOffParentDistance  Min hamming distance for one-off parents [default: 4]
#   maxShift     Max alignment shift [default: 16]
#   minSampleFraction  Fraction of samples that must flag a sequence as chimeric
#                      for it to be removed [default: 0.9]
#   ignoreNNegatives  Allow this many non-flagged samples before applying
#                     minSampleFraction [default: 1]
#   cpus         Number of parallel workers [default: 1]
#
# Returns:
#   Character vector of chimeric sequences
# ---------------------------------------------------------------------------
find_chimeras_sparse <- function(dt,
                                  minFoldParentOverAbundance = 1.5,
                                  minParentAbundance = 2,
                                  allowOneOff = FALSE,
                                  minOneOffParentDistance = 4,
                                  maxShift = 16,
                                  minSampleFraction = 0.9,
                                  ignoreNNegatives = 1,
                                  cpus = 1) {

    # -----------------------------------------------------------------------
    # Build sparse indices
    #
    # sample_index: for a given sample, quickly look up all (sequence, count)
    # seq_index:    for a given sequence, quickly look up all (sample, count)
    # -----------------------------------------------------------------------
    setkey(dt, sample)
    sample_index <- dt[, .(sequence = list(sequence), count = list(count)),
                       by = sample]
    setkey(sample_index, sample)

    setkey(dt, sequence)
    seq_index <- dt[, .(samples = list(sample), counts = list(count),
                        total = sum(count)),
                    by = sequence]

    # Sort by total abundance ascending — least abundant ASVs are most likely
    # to be chimeras, and checking them first against more-abundant parents
    # is the standard UCHIME/dada2 approach
    setorder(seq_index, total)
    all_seqs <- seq_index$sequence

    cat("[INFO] Checking", length(all_seqs), "ASVs for chimeras using",
        cpus, "cores\n")

    # -----------------------------------------------------------------------
    # Per-ASV chimera check (parallelizable)
    # -----------------------------------------------------------------------
    check_one_asv <- function(idx) {
        query_seq     <- seq_index$sequence[idx]
        query_samples <- seq_index$samples[[idx]]
        query_counts  <- seq_index$counts[[idx]]
        nsam  <- length(query_samples)
        nflag <- 0L

        for (s in seq_along(query_samples)) {
            sam   <- query_samples[s]
            qcount <- query_counts[s]

            # Look up all ASVs in this sample
            sam_row <- sample_index[.(sam)]
            sam_seqs   <- sam_row$sequence[[1]]
            sam_counts <- sam_row$count[[1]]

            # Find qualifying parents: more abundant in THIS sample
            min_count <- max(minParentAbundance,
                             minFoldParentOverAbundance * qcount)
            parent_mask <- sam_counts >= min_count & sam_seqs != query_seq
            parents <- sam_seqs[parent_mask]

            if (length(parents) < 2) next  # need at least 2 parents for bimera

            # Use dada2's C-level bimera check
            # isBimera aligns query against all parents and checks if any
            # two parents can reconstruct the query as a chimera
            is_chim <- tryCatch(
                isBimera(query_seq, parents,
                         allowOneOff = allowOneOff,
                         minOneOffParentDistance = minOneOffParentDistance,
                         maxShift = maxShift),
                error = function(e) FALSE
            )

            if (is_chim) nflag <- nflag + 1L
        }

        # Consensus decision: same logic as dada2::isBimeraDenovoTable
        # Flag if chimeric in enough samples (allowing ignoreNNegatives misses)
        is_chimera <- (nflag > 0) &&
            (nflag >= nsam ||
             nflag >= (nsam - ignoreNNegatives) * minSampleFraction)

        return(is_chimera)
    }

    # Run in parallel over all ASVs
    if (cpus > 1) {
        results <- mclapply(seq_along(all_seqs), check_one_asv,
                            mc.cores = cpus)
        is_chimera <- unlist(results)
    } else {
        is_chimera <- vapply(seq_along(all_seqs), check_one_asv, logical(1))
    }

    chimeric_seqs <- all_seqs[is_chimera]

    cat("[INFO] Flagged", length(chimeric_seqs), "chimeras out of",
        length(all_seqs), "ASVs\n")

    return(chimeric_seqs)
}


# ===========================================================================
# Main script
# ===========================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: remove_chimeras.R <seqtab_long.rds> <cpus>")
}

input_path <- args[1]
cpus       <- as.integer(args[2])

# ---------------------------------------------------------------------------
# Load input (long-format data.table from merge_seqtabs.R)
# ---------------------------------------------------------------------------
dt <- readRDS(input_path)

if (!is.data.table(dt)) {
    stop("[ERROR] Expected long-format data.table with columns: sample, sequence, count")
}

n_input_asvs  <- uniqueN(dt$sequence)
n_input_reads <- sum(dt$count)
n_samples     <- uniqueN(dt$sample)

cat("[INFO] Input:", n_samples, "samples,",
    n_input_asvs, "ASVs,", n_input_reads, "reads,",
    nrow(dt), "non-zero entries\n")

# ---------------------------------------------------------------------------
# Pre-filter: remove obvious noise before the expensive chimera check
#
# Singletons can't be meaningfully assessed and ultra-short sequences
# aren't real amplicons. Filtering in long format is trivial.
# ---------------------------------------------------------------------------
seq_totals <- dt[, .(total = sum(count), len = nchar(sequence[1])),
                 by = sequence]

pre_singleton <- seq_totals[total <= 1, sequence]
pre_short     <- seq_totals[len < 20, sequence]
pre_remove    <- union(pre_singleton, pre_short)
n_pre_removed <- length(pre_remove)

if (n_pre_removed > 0) {
    cat("[INFO] Pre-filter: removing", length(pre_singleton), "singletons and",
        length(setdiff(pre_short, pre_singleton)), "ultra-short ASVs\n")
    dt <- dt[!sequence %in% pre_remove]
}

cat("[INFO] Chimera input:", uniqueN(dt$sequence), "ASVs,",
    nrow(dt), "non-zero entries\n")

# ---------------------------------------------------------------------------
# Run sparse chimera detection
# ---------------------------------------------------------------------------
chimeric_seqs <- find_chimeras_sparse(dt, cpus = cpus)

# Remove chimeric sequences
dt_clean <- dt[!sequence %in% chimeric_seqs]

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
chimeras_removed <- length(chimeric_seqs)
total_removed    <- n_pre_removed + chimeras_removed
n_output_asvs    <- uniqueN(dt_clean$sequence)
n_output_reads   <- sum(dt_clean$count)
pct_retained     <- round(n_output_reads / max(n_input_reads, 1) * 100, 1)

cat("[INFO] Removed", chimeras_removed, "chimeras +",
    n_pre_removed, "pre-filtered =", total_removed, "total\n")
cat("[INFO] Retained", pct_retained, "% of reads (",
    n_output_asvs, "ASVs across", uniqueN(dt_clean$sample), "samples)\n")

# ---------------------------------------------------------------------------
# Save output — stays in long format
# ---------------------------------------------------------------------------
saveRDS(dt_clean, "seqtab_nochim.rds")

stats <- data.frame(
    step           = "chimera_removal",
    input_asvs     = n_input_asvs,
    pre_filtered   = n_pre_removed,
    chimeras       = chimeras_removed,
    output_asvs    = n_output_asvs,
    samples        = uniqueN(dt_clean$sample),
    total_reads    = n_output_reads,
    pct_retained   = pct_retained
)
write.table(stats, "chimera_stats.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)
