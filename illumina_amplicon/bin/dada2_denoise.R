#!/usr/bin/env Rscript
#
# dada2_denoise.R — Denoise, merge pairs, and build a sequence table (per plate)
#
# For each sample on a plate:
#   1. Dereplicate forward and reverse filtered reads
#   2. Apply the DADA2 denoising algorithm using pre-learned error models
#   3. Merge denoised forward/reverse pairs (trimming read-through overhangs)
#   4. Combine all merged samples into a sequence table (samples x ASVs)
#   5. Per-plate chimera removal with removeBimeraDenovo()
#
# The sequence table is the primary output — a matrix of read counts where
# rows are samples and columns are unique amplicon sequence variants (ASVs).
#
# Running chimera removal per-plate (step 5) is a key performance optimization.
# Most chimeras are plate-specific (formed during that plate's PCR), and
# removeBimeraDenovo() scales with ASV count. Catching chimeras here, while
# plates are processed in parallel, reduces the post-merge chimera pass from
# ~180K ASVs to ~25K, cutting its runtime by an order of magnitude.
#
# Usage:
#   dada2_denoise.R <plate_id> <errF.rds> <errR.rds> <min_overlap> <cpus>
#
# Inputs (discovered automatically):
#   *_R1.filt.fastq.gz   Forward filtered reads
#   *_R2.filt.fastq.gz   Reverse filtered reads
#
# Outputs:
#   <plate_id>.seqtab.rds   Sequence table (RDS format)
#   <plate_id>.seqtab.tsv   Sequence table (tab-delimited, for inspection)

library(dada2)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
    stop("Usage: dada2_denoise.R <plate_id> <errF.rds> <errR.rds> ",
         "<min_overlap> <cpus>")
}

plate_id    <- args[1]
errF_path   <- args[2]
errR_path   <- args[3]
min_overlap <- as.integer(args[4])
cpus        <- as.integer(args[5])

# ---------------------------------------------------------------------------
# Load pre-learned error models
# ---------------------------------------------------------------------------
errF <- readRDS(errF_path)
errR <- readRDS(errR_path)

# ---------------------------------------------------------------------------
# Discover filtered FASTQ files
# ---------------------------------------------------------------------------
fwd_files <- sort(list.files(".", pattern = "_R1\\.filt\\.fastq\\.gz$",
                             full.names = TRUE))
rev_files <- sort(list.files(".", pattern = "_R2\\.filt\\.fastq\\.gz$",
                             full.names = TRUE))

# Remove empty/missing files
fwd_sizes <- file.size(fwd_files)
rev_sizes <- file.size(rev_files)
valid <- !is.na(fwd_sizes) & fwd_sizes > 20 &
         !is.na(rev_sizes) & rev_sizes > 20
fwd_files <- fwd_files[valid]
rev_files <- rev_files[valid]

# Derive sample names by stripping the _R1.filt.fastq.gz suffix
sample_names <- sub("_R1\\.filt\\.fastq\\.gz$", "", basename(fwd_files))

if (length(fwd_files) == 0) {
    stop("[ERROR] No valid filtered files found for denoising")
}

cat("[INFO] Plate", plate_id, ": denoising", length(fwd_files), "samples\n")

# ---------------------------------------------------------------------------
# Per-sample: dereplicate → denoise → merge pairs
# ---------------------------------------------------------------------------
mergers <- vector("list", length(sample_names))
names(mergers) <- sample_names

for (i in seq_along(sample_names)) {
    cat("[INFO] Processing:", sample_names[i], "\n")

    # Dereplicate: collapse identical reads, track abundances
    derepF <- derepFastq(fwd_files[i], verbose = TRUE)
    derepR <- derepFastq(rev_files[i], verbose = TRUE)

    # Denoise: apply error model to distinguish real variants from errors
    ddF <- dada(derepF, err = errF, multithread = cpus)
    ddR <- dada(derepR, err = errR, multithread = cpus)

    # Merge paired ends
    # trimOverhang=TRUE clips read-through into the opposite primer region,
    # which happens when the amplicon is shorter than 2x read length
    merger <- mergePairs(ddF, derepF, ddR, derepR,
                         trimOverhang = TRUE, verbose = TRUE,
                         minOverlap = min_overlap)

    if (nrow(merger) > 0) {
        mergers[[sample_names[i]]] <- merger
    }
}

# Drop samples that produced no merged reads
mergers <- mergers[lengths(mergers) > 0]

if (length(mergers) == 0) {
    stop("[ERROR] No merged reads produced for plate ", plate_id)
}

# ---------------------------------------------------------------------------
# Build sequence table
# ---------------------------------------------------------------------------
seqtab <- makeSequenceTable(mergers)

cat("[INFO] Plate", plate_id, ": raw table has", nrow(seqtab), "samples,",
    ncol(seqtab), "unique ASVs,", sum(seqtab), "total reads\n")

# ---------------------------------------------------------------------------
# Per-plate chimera removal
#
# Running removeBimeraDenovo() per-plate is much faster than on the merged
# table because it scales with ASV count. Most chimeras are plate-specific
# (formed during that plate's PCR), so catching them here is effective.
# A lighter second pass on the merged table catches any cross-plate chimeras.
# ---------------------------------------------------------------------------
seqtab_nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                     multithread = cpus, verbose = TRUE)

n_chimeras <- ncol(seqtab) - ncol(seqtab_nochim)
pct <- round(sum(seqtab_nochim) / max(sum(seqtab), 1) * 100, 1)
cat("[INFO] Plate", plate_id, ": removed", n_chimeras, "chimeras,",
    "retained", pct, "% of reads (", ncol(seqtab_nochim), "ASVs)\n")

# Save in both formats: RDS for downstream R steps, TSV for inspection
saveRDS(seqtab_nochim, paste0(plate_id, ".seqtab.rds"))
write.table(seqtab_nochim, paste0(plate_id, ".seqtab.tsv"),
            sep = "\t", quote = FALSE)
