#!/usr/bin/env Rscript
#
# renormalize.R — Separate ASVs by taxonomy and normalize within groups
#
# Amplicon datasets mix biologically distinct groups (prokaryotes, organelles,
# eukaryotes) that should not be compared directly. A chloroplast ASV
# dominating 16S reads doesn't mean chloroplasts are more abundant than
# bacteria — it means the 16S primers also amplified chloroplast 16S rRNA.
#
# This script:
#   1. Classifies each ASV into a taxonomic group based on its taxonomy
#      (prokaryote, eukaryote, chloroplast, mitochondria, unknown)
#   2. Normalizes counts within each group to proportions (per-sample)
#      so that downstream analyses compare like with like
#   3. Outputs both per-group and merged normalized tables
#
# All operations work in long format (data.table) for memory efficiency.
#
# Usage:
#   renormalize.R <seqtab.rds> <taxonomy.rds> <tax_db_name>
#
# Arguments:
#   seqtab.rds      Long-format data.table (sample, sequence, count)
#   taxonomy.rds    Taxonomy matrix from assign_taxonomy.R
#   tax_db_name     Name of the taxonomy database used (for logging)
#
# Outputs:
#   renorm_by_group.rds     Named list of long-format data.tables per group
#   renorm_merged.rds       Single long-format data.table with group column
#   renorm_table_list.rds   ASV-to-group mapping
#   renorm_stats.tsv        Per-group summary statistics

library(data.table)

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
    stop("Usage: renormalize.R <seqtab.rds> <taxonomy.rds> <tax_db_name>")
}

seqtab_path  <- args[1]
tax_path     <- args[2]
tax_db_name  <- args[3]

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
dt <- readRDS(seqtab_path)
if (!is.data.table(dt)) {
    stop("[ERROR] Expected long-format data.table for seqtab")
}

tax <- readRDS(tax_path)
if (!is.matrix(tax)) {
    stop("[ERROR] Expected taxonomy matrix")
}

cat("[INFO] Input:", uniqueN(dt$sample), "samples,",
    uniqueN(dt$sequence), "ASVs,", sum(dt$count), "reads\n")
cat("[INFO] Taxonomy from:", tax_db_name, "\n")

# ---------------------------------------------------------------------------
# Classify ASVs into taxonomic groups
#
# The grouping depends on which taxonomy columns are present. SILVA-like
# databases have Domain/Phylum/Class/Order/Family/Genus. We identify:
#   - Chloroplasts: Order == "Chloroplast" (within Cyanobacteria in SILVA)
#   - Mitochondria: Family == "Mitochondria" (within Alphaproteobacteria)
#   - Eukaryotes: Domain == "Eukaryota"
#   - Prokaryotes: Domain == "Bacteria" or "Archaea" (excluding above)
#   - Unknown: everything else
# ---------------------------------------------------------------------------
cat("[INFO] Classifying ASVs into taxonomic groups...\n")

all_seqs <- unique(dt$sequence)

# Build a lookup: sequence → group
# Use the taxonomy matrix (rows = sequences, cols = ranks)
group_map <- data.table(sequence = rownames(tax), group = "unknown")

# Helper to find sequences matching a value at a given rank
find_at_rank <- function(rank_name, values) {
    if (!rank_name %in% colnames(tax)) return(character(0))
    rownames(tax)[tax[, rank_name] %in% values & !is.na(tax[, rank_name])]
}

# Identify groups (order matters — chloroplast/mitochondria before prokaryote)
chloroplast_seqs  <- find_at_rank("Order", "Chloroplast")
mitochondria_seqs <- find_at_rank("Family", "Mitochondria")
eukaryote_seqs    <- find_at_rank("Domain", "Eukaryota")
prokaryote_seqs   <- find_at_rank("Domain", c("Bacteria", "Archaea"))

# Remove organelles from prokaryote set (they're classified under Bacteria
# in SILVA but biologically belong with their eukaryotic hosts)
prokaryote_seqs <- setdiff(prokaryote_seqs,
                           c(chloroplast_seqs, mitochondria_seqs))

# Assign groups
group_map[sequence %in% chloroplast_seqs,  group := "chloroplast"]
group_map[sequence %in% mitochondria_seqs, group := "mitochondria"]
group_map[sequence %in% eukaryote_seqs,    group := "eukaryote"]
group_map[sequence %in% prokaryote_seqs,   group := "prokaryote"]

# Any sequences in dt but not in taxonomy stay "unknown"
missing_seqs <- setdiff(all_seqs, group_map$sequence)
if (length(missing_seqs) > 0) {
    group_map <- rbind(group_map,
                       data.table(sequence = missing_seqs, group = "unknown"))
}

# Summary
group_counts <- group_map[, .N, by = group][order(-N)]
cat("[INFO] ASV groups:\n")
for (i in seq_len(nrow(group_counts))) {
    cat("  ", group_counts$group[i], ":", group_counts$N[i], "ASVs\n")
}

# ---------------------------------------------------------------------------
# Add group column to count data
# ---------------------------------------------------------------------------
dt_grouped <- merge(dt, group_map, by = "sequence", all.x = TRUE)
dt_grouped[is.na(group), group := "unknown"]

# ---------------------------------------------------------------------------
# Normalize within each group
#
# For each group, compute per-sample proportions:
#   proportion = count / sum(counts for this group in this sample)
#
# This ensures that, e.g., prokaryote proportions sum to 1.0 within each
# sample, regardless of how many chloroplast or mitochondrial reads were
# also captured. This is essential for comparing relative abundances
# across samples with different organelle contamination levels.
# ---------------------------------------------------------------------------
cat("[INFO] Normalizing within groups...\n")

dt_grouped[, total_in_group := sum(count), by = .(sample, group)]
dt_grouped[, proportion := count / ifelse(total_in_group > 0, total_in_group, 1)]
dt_grouped[is.nan(proportion), proportion := 0]

# ---------------------------------------------------------------------------
# Build per-group outputs
# ---------------------------------------------------------------------------
groups <- unique(dt_grouped$group)
group_tables <- list()
stats_list <- list()

for (g in groups) {
    dt_g <- dt_grouped[group == g, .(sample, sequence, count, proportion)]
    group_tables[[g]] <- dt_g

    n_samples <- uniqueN(dt_g$sample)
    n_asvs    <- uniqueN(dt_g$sequence)
    n_reads   <- sum(dt_g$count)
    stats_list[[g]] <- data.frame(
        group = g, samples = n_samples, asvs = n_asvs, reads = n_reads
    )
    cat("[INFO]", g, ":", n_asvs, "ASVs,", n_reads, "reads across",
        n_samples, "samples\n")
}

# ---------------------------------------------------------------------------
# Save outputs
# ---------------------------------------------------------------------------

# Per-group long-format tables
saveRDS(group_tables, "renorm_by_group.rds")

# Merged table with group column (for downstream queries)
dt_out <- dt_grouped[, .(sample, sequence, count, proportion, group)]
saveRDS(dt_out, "renorm_merged.rds")

# ASV-to-group mapping
saveRDS(group_map, "renorm_table_list.rds")

# Stats
stats <- do.call(rbind, stats_list)
write.table(stats, "renorm_stats.tsv",
            sep = "\t", row.names = FALSE, quote = FALSE)

cat("[INFO] Renormalization complete\n")
