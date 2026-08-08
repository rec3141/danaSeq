#!/usr/bin/env Rscript
#
# dada2_learn_errors.R — Learn error rates from filtered reads (per plate)
#
# DADA2's error model is learned per plate because samples sharing a PCR
# amplification history share error profiles. This script discovers all
# filtered FASTQ files in the working directory, learns separate forward
# and reverse error models, saves them as RDS, and generates diagnostic
# error-rate plots.
#
# Usage:
#   dada2_learn_errors.R <plate_id> <cpus>
#
# Inputs (discovered automatically):
#   *_R1.filt.fastq.gz   Forward filtered reads (one per sample)
#   *_R2.filt.fastq.gz   Reverse filtered reads (one per sample)
#
# Outputs:
#   <plate_id>_errF.rds          Forward error model
#   <plate_id>_errR.rds          Reverse error model
#   <plate_id>_error_rates.pdf   Diagnostic plots (observed vs fitted rates)

library(dada2)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
    stop("Usage: dada2_learn_errors.R <plate_id> <cpus>")
}

plate_id <- args[1]
cpus     <- as.integer(args[2])

# ---------------------------------------------------------------------------
# Discover filtered FASTQ files in the working directory
# ---------------------------------------------------------------------------
fwd_files <- sort(list.files(".", pattern = "_R1\\.filt\\.fastq\\.gz$",
                             full.names = TRUE))
rev_files <- sort(list.files(".", pattern = "_R2\\.filt\\.fastq\\.gz$",
                             full.names = TRUE))

# Drop any empty or missing files (can happen when all reads are filtered out)
fwd_sizes <- file.size(fwd_files)
rev_sizes <- file.size(rev_files)
valid <- !is.na(fwd_sizes) & fwd_sizes > 20 &
         !is.na(rev_sizes) & rev_sizes > 20
fwd_files <- fwd_files[valid]
rev_files <- rev_files[valid]

if (length(fwd_files) == 0) {
    stop("[ERROR] No valid filtered FASTQ files found for plate ", plate_id)
}

cat("[INFO] Plate", plate_id, ": learning error rates from",
    length(fwd_files), "samples\n")

# ---------------------------------------------------------------------------
# Learn error models
#
# learnErrors pools quality-score information across all input files and fits
# a parametric error model. This is the most compute-intensive DADA2 step
# because it iterates until convergence.
# ---------------------------------------------------------------------------
errF <- learnErrors(fwd_files, multithread = cpus)
errR <- learnErrors(rev_files, multithread = cpus)

saveRDS(errF, paste0(plate_id, "_errF.rds"))
saveRDS(errR, paste0(plate_id, "_errR.rds"))

# ---------------------------------------------------------------------------
# Diagnostic plots — observed error rates vs fitted model
#
# Points should track the fitted line. Large deviations suggest the error
# model is a poor fit (e.g., mixed sequencing chemistries on one plate).
# ---------------------------------------------------------------------------
pdf(paste0(plate_id, "_error_rates.pdf"), width = 10, height = 10)
print(plotErrors(errF, nominalQ = TRUE))
print(plotErrors(errR, nominalQ = TRUE))
dev.off()

cat("[INFO] Plate", plate_id, ": error models saved\n")
