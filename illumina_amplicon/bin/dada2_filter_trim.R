#!/usr/bin/env Rscript
#
# dada2_filter_trim.R — Per-sample quality filtering of paired-end reads
#
# Runs DADA2's filterAndTrim() on a single sample pair. Removes low-quality
# reads, PhiX spike-in, and optionally truncates to fixed lengths. Produces
# compressed filtered FASTQ files and a one-row stats TSV.
#
# Usage:
#   dada2_filter_trim.R <sample_id> <fwd.fastq.gz> <rev.fastq.gz> \
#       <maxEE> <truncQ> <maxN> <truncLen_fwd> <truncLen_rev> <cpus>
#
# Outputs:
#   <sample_id>_R1.filt.fastq.gz   Filtered forward reads
#   <sample_id>_R2.filt.fastq.gz   Filtered reverse reads
#   <sample_id>_filt_stats.tsv      Filtering statistics

library(dada2)

# ---------------------------------------------------------------------------
# Parse command-line arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 9) {
    stop("Usage: dada2_filter_trim.R <sample_id> <r1> <r2> ",
         "<maxEE> <truncQ> <maxN> <truncLen_fwd> <truncLen_rev> <cpus>")
}

sample_id    <- args[1]
r1_in        <- args[2]
r2_in        <- args[3]
max_ee       <- as.numeric(args[4])
trunc_q      <- as.numeric(args[5])
max_n        <- as.numeric(args[6])
trunc_fwd    <- as.integer(args[7])
trunc_rev    <- as.integer(args[8])
cpus         <- as.integer(args[9])

# Output filenames
r1_out <- paste0(sample_id, "_R1.filt.fastq.gz")
r2_out <- paste0(sample_id, "_R2.filt.fastq.gz")

# Build truncLen vector: c(0,0) means no truncation
trunc_len <- if (trunc_fwd > 0 && trunc_rev > 0) {
    c(trunc_fwd, trunc_rev)
} else {
    c(0, 0)
}

# ---------------------------------------------------------------------------
# Filter and trim
# ---------------------------------------------------------------------------
cat("[INFO]", sample_id, ": filtering with maxEE=", max_ee,
    "truncQ=", trunc_q, "truncLen=", paste(trunc_len, collapse = ","), "\n")

fto <- withCallingHandlers(
    filterAndTrim(
        fwd      = r1_in,  filt     = r1_out,
        rev      = r2_in,  filt.rev = r2_out,
        maxEE    = max_ee,
        truncQ   = trunc_q,
        maxN     = max_n,
        truncLen = trunc_len,
        rm.phix  = TRUE,
        compress = TRUE,
        verbose  = TRUE,
        multithread = cpus
    ),
    warning = function(w) {
        cat("[WARNING]", sample_id, ":", conditionMessage(w), "\n")
        invokeRestart("muffleWarning")
    }
)

# Handle samples where all reads were filtered out
reads_in  <- if (nrow(fto) > 0) fto[1, 1] else 0
reads_out <- if (nrow(fto) > 0) fto[1, 2] else 0
pct       <- round(reads_out / max(reads_in, 1) * 100, 1)

# Create empty output files if filtering removed everything
# (downstream steps filter these out by file size)
if (reads_out == 0) {
    cat("[WARNING]", sample_id, ": no reads passed filter, creating empty outputs\n")
    file.create(r1_out)
    file.create(r2_out)
}

stats <- data.frame(
    sample       = sample_id,
    reads_in     = reads_in,
    reads_out    = reads_out,
    pct_retained = pct
)
write.table(stats, paste0(sample_id, "_filt_stats.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("[INFO]", sample_id, ":", reads_out, "/", reads_in,
    "reads passed filter (", pct, "%)\n")
