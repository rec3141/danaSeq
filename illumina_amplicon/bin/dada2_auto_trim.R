#!/usr/bin/env Rscript
#
# dada2_auto_trim.R — Quality profiling and truncation length selection (R engine)
#
# The R port of `microscape auto-trim`. AUTO_TRIM used to shell out to that
# Python CLI whatever --lang was set to, so an R run still needed the Python
# microscape package installed — and envs/r.yml does not carry it. Inside the
# container both conda envs are on PATH so it happened to work; under conda
# with --lang R it could not.
#
# This must agree with microscape/quality.py exactly, or the two engines
# truncate differently on the same reads. The places that matter:
#   - numpy's default percentile ('linear') is R's quantile type 7
#   - int(np.percentile(...)) truncates toward zero, like as.integer()
#   - np.convolve(..., mode="valid") is a windowed mean over positions
#     1..(n-window+1), computed here with cumsum so it is exact
#   - the quality position is floored at min_length and THEN capped by the
#     length percentile (see quality.py: the floor exists so a quality dip
#     cannot truncate below overlap, and must never ask for bases the reads
#     do not have)
#
# Usage:
#   dada2_auto_trim.R <input_dir> <min_quality> <min_length> <output.tsv> [window] [n_reads] [n_files]
#
# Output: a TSV of key<TAB>value pairs, same keys and order as the Python CLI.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) {
    stop("Usage: dada2_auto_trim.R <input_dir> <min_quality> <min_length> ",
         "<output.tsv> [window] [n_reads] [n_files]")
}

input_dir   <- args[1]
min_quality <- as.numeric(args[2])
min_length  <- as.integer(args[3])
out_path    <- args[4]
window      <- if (length(args) >= 5) as.integer(args[5]) else 10L
n_reads     <- if (length(args) >= 6) as.integer(args[6]) else 10000L
n_files     <- if (length(args) >= 7) as.integer(args[7]) else 20L

# ---------------------------------------------------------------------------
# Find paired FASTQs — same glob patterns, dedup and sort as the Python side
# ---------------------------------------------------------------------------
.glob <- function(dir, pats) {
    hits <- unlist(lapply(pats, function(p) Sys.glob(file.path(dir, p))))
    sort(unique(hits))
}

fwd_files <- .glob(input_dir, c("*_1.fastq.gz", "*_R1*.fastq.gz", "*_R1.*.fastq.gz"))
rev_files <- .glob(input_dir, c("*_2.fastq.gz", "*_R2*.fastq.gz", "*_R2.*.fastq.gz"))

if (!length(fwd_files) || !length(rev_files)) {
    stop("No paired FASTQ files found in ", input_dir)
}

# Sample from a subset of files: files[::step][:n_sample]
.subsample <- function(files, n_files) {
    n_sample <- min(n_files, length(files))
    step <- max(1L, length(files) %/% n_sample)
    strided <- files[seq(1L, length(files), by = step)]
    strided[seq_len(min(n_sample, length(strided)))]
}

fwd_sample <- .subsample(fwd_files, n_files)
rev_sample <- .subsample(rev_files, n_files)

# ---------------------------------------------------------------------------
# Read quality scores, stopping at n_reads across files (as _read_quals does)
# ---------------------------------------------------------------------------
read_quals <- function(files, n_reads) {
    out <- vector("list", n_reads)
    count <- 0L
    for (fpath in files) {
        con <- if (grepl("\\.gz$", fpath)) gzfile(fpath, "rt") else file(fpath, "rt")
        repeat {
            if (count >= n_reads) break
            rec <- readLines(con, n = 4L, warn = FALSE)
            if (length(rec) < 4L) break            # truncated/at EOF
            qual_str <- sub("\\s+$", "", rec[4])   # Python rstrip()
            if (!nzchar(qual_str)) break
            count <- count + 1L
            out[[count]] <- utf8ToInt(qual_str) - 33L
        }
        close(con)
        if (count >= n_reads) break
    }
    if (count == 0L) list() else out[seq_len(count)]
}

# ---------------------------------------------------------------------------
# Drop fragment reads before measuring the length distribution (_amplicon_lengths)
# ---------------------------------------------------------------------------
amplicon_lengths <- function(lengths, floor_frac = 0.5) {
    if (!length(lengths)) return(lengths)
    med <- as.numeric(stats::median(lengths))
    usable <- lengths[lengths >= med * floor_frac]
    # Nearly all fragments means a broken library, not a bimodal one; measuring
    # the remainder would be less honest than measuring all of it.
    if (length(usable) >= max(10, 0.05 * length(lengths))) usable else lengths
}

# ---------------------------------------------------------------------------
# Truncation position: quality cliff, floored at min_length, capped by the
# 10th-percentile amplicon read length
# ---------------------------------------------------------------------------
find_trunc_pos <- function(quals_list, min_q, window, min_length = 0L) {
    if (!length(quals_list)) return(0L)
    lengths <- vapply(quals_list, length, integer(1))
    max_len <- max(lengths)

    p10 <- as.integer(stats::quantile(amplicon_lengths(lengths), 0.10,
                                      type = 7, names = FALSE))
    len_cap <- min(p10, max_len)

    # Per-position median over ragged rows (np.nanmedian on a NaN-padded matrix)
    mat <- matrix(NA_integer_, nrow = length(quals_list), ncol = max_len)
    for (i in seq_along(quals_list)) {
        q <- quals_list[[i]]
        mat[i, seq_along(q)] <- q
    }
    medians <- apply(mat, 2L, stats::median, na.rm = TRUE)

    n <- length(medians)
    if (n < window) {
        return(min(max(as.integer(n), min_length), len_cap))
    }

    # Windowed mean, equivalent to np.convolve(..., mode="valid")
    cs <- cumsum(c(0, medians))
    rolling <- (cs[(window + 1L):(n + 1L)] - cs[1L:(n - window + 1L)]) / window

    quality_pos <- as.integer(n)
    below <- which(rolling < min_q)
    if (length(below)) {
        # Python's 0-based i becomes below[1] - 1 here
        quality_pos <- as.integer((below[1L] - 1L) + window %/% 2L)
    }

    min(max(quality_pos, min_length), len_cap)
}

fwd_quals <- read_quals(fwd_sample, n_reads)
rev_quals <- read_quals(rev_sample, n_reads)

trunc_fwd <- find_trunc_pos(fwd_quals, min_quality, window, min_length)
trunc_rev <- find_trunc_pos(rev_quals, min_quality, window, min_length)

fwd_len <- if (length(fwd_quals)) max(vapply(fwd_quals, length, integer(1))) else 0L
rev_len <- if (length(rev_quals)) max(vapply(rev_quals, length, integer(1))) else 0L

# ---------------------------------------------------------------------------
# Emit the same keys, order and formatting as the Python CLI
# ---------------------------------------------------------------------------
# Python writes floats with a trailing .0 (f"{25.0}" -> "25.0"); match it so the
# published quality_check TSV is identical whichever engine produced it.
.pyfloat <- function(x) {
    if (x == as.integer(x)) sprintf("%.1f", x) else format(x, trim = TRUE)
}

rows <- c(
    sprintf("trunc_len_fwd\t%d",   as.integer(trunc_fwd)),
    sprintf("trunc_len_rev\t%d",   as.integer(trunc_rev)),
    sprintf("fwd_read_len\t%d",    as.integer(fwd_len)),
    sprintf("rev_read_len\t%d",    as.integer(rev_len)),
    sprintf("n_reads_sampled\t%d", length(fwd_quals)),
    sprintf("min_quality\t%s",     .pyfloat(min_quality)),
    sprintf("min_length\t%d",      as.integer(min_length))
)
writeLines(rows, out_path)

# Same stdout summary the Python CLI prints, for the run log.
cat(sprintf("trunc_len_fwd=%d\n", as.integer(trunc_fwd)))
cat(sprintf("trunc_len_rev=%d\n", as.integer(trunc_rev)))
cat(sprintf("fwd_read_len=%d\n",  as.integer(fwd_len)))
cat(sprintf("rev_read_len=%d\n",  as.integer(rev_len)))
cat(sprintf("n_reads_sampled=%d\n", length(fwd_quals)))
