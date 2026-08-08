#!/usr/bin/env Rscript
#
# build_shiny.R — Prepare pre-computed data for the Microscape Shiny app
#
# Reads all pipeline outputs and bundles them into a single app-ready data
# package (app_data.rds). This keeps the Shiny app thin — it loads one file
# and renders, with no heavy computation at runtime.
#
# Usage:
#   build_shiny.R <seqtab.rds> <renorm.rds> <taxonomy_dir> <metadata.rds> \
#                 <sample_tsne.rds> <seq_tsne.rds> <network.rds>
#
# Arguments:
#   seqtab.rds      Long-format data.table (sample, sequence, count)
#   renorm.rds      Renormalized merged data.table (sample, sequence, count,
#                   proportion, group)
#   taxonomy_dir    Directory containing *_taxonomy.rds files
#   metadata.rds    Sample metadata data.frame
#   sample_tsne.rds Sample t-SNE coordinates (matrix or data.frame, 2 cols)
#   seq_tsne.rds    ASV t-SNE coordinates (matrix or data.frame, 2 cols)
#   network.rds     Network edge list (node1, node2, corr, weight, color)
#
# Outputs (written to current directory):
#   app_data.rds    Bundled data for app.R
#   app.R           Copied from the same bin directory as this script

library(data.table)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 7) {
    stop("Usage: build_shiny.R <seqtab.rds> <renorm.rds> <taxonomy_dir> ",
         "<metadata.rds> <sample_tsne.rds> <seq_tsne.rds> <network.rds>")
}

seqtab_path     <- args[1]
renorm_path     <- args[2]
taxonomy_dir    <- args[3]
metadata_path   <- args[4]
sample_tsne_path <- args[5]
seq_tsne_path   <- args[6]
network_path    <- args[7]

cat("[INFO] Loading pipeline outputs...\n")

# ---------------------------------------------------------------------------
# 1. Load all pipeline outputs
# ---------------------------------------------------------------------------

# Long-format count data (filtered)
dt_raw <- readRDS(seqtab_path)
if (!is.data.table(dt_raw)) dt_raw <- as.data.table(dt_raw)
cat("[INFO] Seqtab:", uniqueN(dt_raw$sample), "samples,",
    uniqueN(dt_raw$sequence), "ASVs\n")

# Renormalized data with group assignments and proportions
dt_renorm <- readRDS(renorm_path)
if (!is.data.table(dt_renorm)) dt_renorm <- as.data.table(dt_renorm)
cat("[INFO] Renorm:", nrow(dt_renorm), "rows\n")

# Taxonomy: find *_taxonomy.rds files
# taxonomy_dir can be a directory OR a space-separated list of staged files
# (Nextflow stages collected files flat in the work directory)
if (file.info(taxonomy_dir)$isdir %in% TRUE) {
    tax_files <- list.files(taxonomy_dir, pattern = "_taxonomy\\.rds$",
                            full.names = TRUE)
} else {
    # Files staged flat — search the working directory
    tax_files <- list.files(".", pattern = "_taxonomy\\.rds$",
                            full.names = TRUE)
}
if (length(tax_files) == 0) {
    stop("[ERROR] No taxonomy files found")
}
taxonomy <- list()
for (tf in tax_files) {
    db_name <- sub("_taxonomy\\.rds$", "", basename(tf))
    taxonomy[[db_name]] <- readRDS(tf)
    cat("[INFO] Taxonomy db '", db_name, "':", nrow(taxonomy[[db_name]]),
        "sequences,", ncol(taxonomy[[db_name]]), "ranks\n")
}

# Metadata (optional — may be NO_METADATA placeholder)
if (file.exists(metadata_path) && metadata_path != "NO_METADATA") {
    metadata <- readRDS(metadata_path)
    if (!is.data.frame(metadata)) metadata <- as.data.frame(metadata)
    cat("[INFO] Metadata:", nrow(metadata), "rows,", ncol(metadata), "columns\n")
} else {
    cat("[INFO] No metadata provided, creating minimal metadata\n")
    samples <- unique(dt_raw$sample)
    metadata <- data.frame(
        sample_name = samples,
        row.names = samples,
        stringsAsFactors = FALSE
    )
}

# Sample t-SNE coordinates
sample_tsne_raw <- readRDS(sample_tsne_path)
if (is.matrix(sample_tsne_raw)) {
    sample_tsne <- data.frame(
        x = sample_tsne_raw[, 1],
        y = sample_tsne_raw[, 2],
        sample = if (!is.null(rownames(sample_tsne_raw))) {
            rownames(sample_tsne_raw)
        } else {
            paste0("S", seq_len(nrow(sample_tsne_raw)))
        },
        stringsAsFactors = FALSE
    )
} else {
    sample_tsne <- as.data.frame(sample_tsne_raw)
    # Map known column names to x, y, sample
    if ("tSNE1" %in% names(sample_tsne)) {
        sample_tsne$x <- sample_tsne$tSNE1
        sample_tsne$y <- sample_tsne$tSNE2
    }
    if ("label" %in% names(sample_tsne) && !"sample" %in% names(sample_tsne)) {
        sample_tsne$sample <- sample_tsne$label
    }
    if (!"sample" %in% names(sample_tsne)) {
        if (!is.null(rownames(sample_tsne))) {
            sample_tsne$sample <- rownames(sample_tsne)
        } else {
            sample_tsne$sample <- paste0("S", seq_len(nrow(sample_tsne)))
        }
    }
    if (!"x" %in% names(sample_tsne)) {
        names(sample_tsne)[1:2] <- c("x", "y")
    }
    # Keep only the columns we need
    sample_tsne <- sample_tsne[, c("x", "y", "sample")]
}
cat("[INFO] Sample t-SNE:", nrow(sample_tsne), "points\n")

# ASV t-SNE coordinates
seq_tsne_raw <- readRDS(seq_tsne_path)
if (is.matrix(seq_tsne_raw)) {
    seq_tsne <- data.frame(
        x = seq_tsne_raw[, 1],
        y = seq_tsne_raw[, 2],
        sequence = if (!is.null(rownames(seq_tsne_raw))) {
            rownames(seq_tsne_raw)
        } else {
            paste0("ESV_", seq_len(nrow(seq_tsne_raw)))
        },
        stringsAsFactors = FALSE
    )
} else {
    seq_tsne <- as.data.frame(seq_tsne_raw)
    # Map known column names
    if ("tSNE1" %in% names(seq_tsne)) {
        seq_tsne$x <- seq_tsne$tSNE1
        seq_tsne$y <- seq_tsne$tSNE2
    }
    if ("label" %in% names(seq_tsne) && !"sequence" %in% names(seq_tsne)) {
        seq_tsne$sequence <- seq_tsne$label
    }
    if (!"sequence" %in% names(seq_tsne)) {
        if (!is.null(rownames(seq_tsne))) {
            seq_tsne$sequence <- rownames(seq_tsne)
        } else {
            seq_tsne$sequence <- paste0("ESV_", seq_len(nrow(seq_tsne)))
        }
    }
    if (!"x" %in% names(seq_tsne)) {
        names(seq_tsne)[1:2] <- c("x", "y")
    }
    seq_tsne <- seq_tsne[, c("x", "y", "sequence")]
}
cat("[INFO] Sequence t-SNE:", nrow(seq_tsne), "points\n")

# Network edges
network_edges <- readRDS(network_path)
if (!is.data.frame(network_edges)) {
    network_edges <- as.data.frame(network_edges)
}
# Normalize column names to: node1, node2, correlation, weight, color
expected_cols <- c("node1", "node2", "correlation", "weight", "color")
if ("corr" %in% names(network_edges) && !"correlation" %in% names(network_edges)) {
    names(network_edges)[names(network_edges) == "corr"] <- "correlation"
}
cat("[INFO] Network:", nrow(network_edges), "edges\n")

# ---------------------------------------------------------------------------
# 2. Build lookup tables
# ---------------------------------------------------------------------------
cat("[INFO] Building lookup tables...\n")

# Use the first taxonomy database as primary for naming
primary_db <- names(taxonomy)[1]
primary_tax <- taxonomy[[primary_db]]

# Build ASV taxonomy name strings (concatenated taxonomy)
# e.g. "Bacteria;Proteobacteria;Alphaproteobacteria;Rhodobacterales;..."
all_seqs <- unique(dt_renorm$sequence)
tax_strings <- character(length(all_seqs))
names(tax_strings) <- all_seqs

for (i in seq_along(all_seqs)) {
    sq <- all_seqs[i]
    if (sq %in% rownames(primary_tax)) {
        ranks <- primary_tax[sq, ]
        ranks[is.na(ranks)] <- ""
        tax_strings[i] <- paste(ranks[ranks != ""], collapse = ";")
    } else {
        tax_strings[i] <- "unclassified"
    }
}

# ESV IDs — numeric index mapped to sequences
esv_ids <- paste0("ESV_", seq_along(all_seqs))
names(esv_ids) <- all_seqs

# Group color assignments
group_colors <- c(
    prokaryote   = "#427df4",
    eukaryote    = "#f44265",
    chloroplast  = "#0fe047",
    mitochondria = "#42e5f4",
    unknown      = "#999999"
)

# ---------------------------------------------------------------------------
# 3. Pre-compute derived data
# ---------------------------------------------------------------------------
cat("[INFO] Pre-computing derived summaries...\n")

# --- Per-ASV summaries ---
asv_summary <- dt_renorm[, .(
    total_reads = sum(count),
    n_samples   = uniqueN(sample)
), by = sequence]

# Add group (take the first group assignment per sequence)
asv_group <- dt_renorm[, .(group = group[1]), by = sequence]
asv_summary <- merge(asv_summary, asv_group, by = "sequence", all.x = TRUE)

# Add taxonomy string
asv_summary[, taxonomy_string := tax_strings[sequence]]
asv_summary[, esv_id := esv_ids[sequence]]

# Add group color
asv_summary[, group_color := group_colors[group]]
asv_summary[is.na(group_color), group_color := "#999999"]

cat("[INFO] ASV info:", nrow(asv_summary), "ASVs\n")

# --- Per-sample summaries ---
sample_summary <- dt_renorm[, .(
    total_reads = sum(count),
    n_asvs      = uniqueN(sequence)
), by = sample]

# Count primer/target groups per sample
sample_groups <- dt_renorm[, .(
    groups = paste(sort(unique(group)), collapse = ",")
), by = sample]
sample_summary <- merge(sample_summary, sample_groups, by = "sample", all.x = TRUE)

# Merge with metadata if possible
if ("cellid" %in% names(metadata)) {
    # Match samples to metadata by cellid
    meta_dt <- as.data.table(metadata)
    sample_summary <- merge(sample_summary, meta_dt,
                            by.x = "sample", by.y = "cellid", all.x = TRUE)
} else if ("sample" %in% names(metadata)) {
    meta_dt <- as.data.table(metadata)
    sample_summary <- merge(sample_summary, meta_dt,
                            by = "sample", all.x = TRUE)
} else {
    # Try matching by rownames
    meta_dt <- as.data.table(metadata, keep.rownames = "sample")
    sample_summary <- merge(sample_summary, meta_dt,
                            by = "sample", all.x = TRUE)
}

cat("[INFO] Sample info:", nrow(sample_summary), "samples\n")

# --- Prepare the long-format count data for the app ---
# Ensure dt_renorm has all needed columns
dt_counts <- dt_renorm[, .(sample, sequence, count, group, proportion)]

# Add ESV IDs to seq_tsne for network matching
seq_tsne$esv_id <- esv_ids[seq_tsne$sequence]

# Map network node IDs to t-SNE coordinates for edge drawing
# Network nodes might use ESV_N format or raw sequences
if (nrow(network_edges) > 0) {
    # Build coordinate lookup
    tsne_lookup <- seq_tsne[, c("sequence", "x", "y")]
    names(tsne_lookup) <- c("id", "x", "y")

    # Also try ESV_N format
    esv_lookup <- seq_tsne[, c("esv_id", "x", "y")]
    names(esv_lookup) <- c("id", "x", "y")
    coord_lookup <- rbind(tsne_lookup, esv_lookup)
    coord_lookup <- coord_lookup[!is.na(coord_lookup$id), ]
    rownames(coord_lookup) <- coord_lookup$id

    # Add x,y coordinates for both endpoints
    n1_match <- match(network_edges$node1, coord_lookup$id)
    n2_match <- match(network_edges$node2, coord_lookup$id)

    valid_edges <- !is.na(n1_match) & !is.na(n2_match)
    network_edges <- network_edges[valid_edges, ]
    n1_match <- n1_match[valid_edges]
    n2_match <- n2_match[valid_edges]

    network_edges$x1 <- coord_lookup$x[n1_match]
    network_edges$y1 <- coord_lookup$y[n1_match]
    network_edges$x2 <- coord_lookup$x[n2_match]
    network_edges$y2 <- coord_lookup$y[n2_match]

    cat("[INFO] Network edges with valid coordinates:", nrow(network_edges), "\n")
}

# ---------------------------------------------------------------------------
# 4. Bundle everything into app_data.rds
# ---------------------------------------------------------------------------
cat("[INFO] Bundling app data...\n")

app_data <- list(
    dt_counts     = dt_counts,
    sample_tsne   = sample_tsne,
    seq_tsne      = seq_tsne,
    taxonomy      = taxonomy,
    metadata      = as.data.frame(sample_summary),
    network_edges = network_edges,
    asv_info      = as.data.frame(asv_summary),
    sample_info   = as.data.frame(sample_summary),
    group_colors  = group_colors,
    tax_strings   = tax_strings,
    esv_ids       = esv_ids,
    primary_db    = primary_db,
    db_names      = names(taxonomy)
)

saveRDS(app_data, "app_data.rds")
cat("[INFO] Saved app_data.rds (",
    round(file.info("app_data.rds")$size / 1e6, 1), " MB)\n")

# ---------------------------------------------------------------------------
# 5. Copy app.R to the output directory
# ---------------------------------------------------------------------------
# app.R lives in the same directory as this script
script_dir <- dirname(sub("--file=", "",
                          grep("--file=", commandArgs(FALSE), value = TRUE)))
if (length(script_dir) == 0 || script_dir == "") {
    script_dir <- "."
}
app_source <- file.path(script_dir, "app.R")
if (file.exists(app_source)) {
    file.copy(app_source, "app.R", overwrite = TRUE)
    cat("[INFO] Copied app.R to output directory\n")
} else {
    cat("[WARNING] Could not find app.R at:", app_source, "\n")
    cat("[WARNING] You will need to place app.R alongside app_data.rds manually\n")
}

cat("[INFO] build_shiny.R complete\n")
