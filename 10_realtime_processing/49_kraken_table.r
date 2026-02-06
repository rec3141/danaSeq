# ABUNDANCE MATRIX GENERATOR MODULE
# Add this section to your existing eDNA mapping pipeline

# ENHANCED QUALITY FILTERING MODULE FOR KRAKEN2 RESULTS
# Add this to your eDNA mapping pipeline

# UPDATED CONFIGURATION WITH QUALITY FILTERING PARAMETERS
abundance_thresholds <- list(
  # Basic abundance thresholds
  rpm_threshold = 10,              # minimum reads per million
  read_threshold = 500,            # minimum raw read count (increased from 5)
  prevalence_threshold = 0.1,     # minimum fraction of samples (0.1 = 10%)
  
  # Quality filtering parameters
  confidence_threshold = 0.1,     # kraken2 confidence score (if available)
  min_hit_groups = 2,             # minimum hit groups for classification
  bracken_threshold = 10,         # bracken minimum read threshold
  
  # Advanced filtering
  min_unique_kmers = 50,           # minimum unique k-mers (if using KrakenUniq)
  max_contaminant_ratio = 0.95,   # filter if >95% reads are from likely contaminants
  min_taxonomic_depth = 3         # require assignment at least to family level
)

# FUNCTION TO APPLY QUALITY FILTERS TO KRAKEN DATA
apply_quality_filters <- function(kraken_df, filter_type = "standard") {
  
  # Start with original data
  filtered_df <- kraken_df
  
  # Track filtering steps
  filter_log <- tibble(
    step = character(),
    taxa_before = integer(),
    taxa_after = integer(),
    reads_before = integer(),
    reads_after = integer()
  )
  
  # Initial counts
  initial_taxa <- nrow(filtered_df)
  initial_reads <- sum(filtered_df$reads, na.rm = TRUE)
  
  # STEP 1: Basic read count filtering
  if (filter_type %in% c("standard", "strict")) {
    threshold <- ifelse(filter_type == "strict", 20, abundance_thresholds$read_threshold)
    
    before_taxa <- nrow(filtered_df)
    before_reads <- sum(filtered_df$reads, na.rm = TRUE)
    
    filtered_df <- filtered_df %>%
      filter(reads >= threshold)
    
    filter_log <- filter_log %>%
      add_row(
        step = paste0("Min ", threshold, " reads"),
        taxa_before = before_taxa,
        taxa_after = nrow(filtered_df),
        reads_before = before_reads,
        reads_after = sum(filtered_df$reads, na.rm = TRUE)
      )
  }
  
  # STEP 2: Remove likely contaminants (human, unclassified at high levels)
  if (filter_type %in% c("standard", "strict")) {
    before_taxa <- nrow(filtered_df)
    before_reads <- sum(filtered_df$reads, na.rm = TRUE)
    
    # Filter out obvious contaminants
    filtered_df <- filtered_df %>%
      filter(
        !str_detect(name, "Homo|Human|unclassified", ignore.case = TRUE),
        !str_detect(taxonomy, "root; unclassified", ignore.case = TRUE)
      )
    
    filter_log <- filter_log %>%
      add_row(
        step = "Remove contaminants",
        taxa_before = before_taxa,
        taxa_after = nrow(filtered_df),
        reads_before = before_reads,
        reads_after = sum(filtered_df$reads, na.rm = TRUE)
      )
  }
  
  # STEP 3: Taxonomic depth filtering
  if (filter_type == "strict") {
    before_taxa <- nrow(filtered_df)
    before_reads <- sum(filtered_df$reads, na.rm = TRUE)
    
    # Count taxonomic levels (root; domain; phylum; class; order; family; genus; species)
    filtered_df <- filtered_df %>%
      mutate(tax_depth = str_count(taxonomy, ";")) %>%
      filter(tax_depth >= abundance_thresholds$min_taxonomic_depth) %>%
      select(-tax_depth)
    
    filter_log <- filter_log %>%
      add_row(
        step = "Min taxonomic depth",
        taxa_before = before_taxa,
        taxa_after = nrow(filtered_df),
        reads_before = before_reads,
        reads_after = sum(filtered_df$reads, na.rm = TRUE)
      )
  }
  
  # STEP 4: Statistical outlier removal (extreme low abundance across samples)
  if (filter_type == "strict") {
    before_taxa <- nrow(filtered_df)
    before_reads <- sum(filtered_df$reads, na.rm = TRUE)
    
    # Calculate taxa statistics across samples
    taxa_stats <- filtered_df %>%
      group_by(name) %>%
      summarise(
        median_rpm = median(rpm, na.rm = TRUE),
        max_rpm = max(rpm, na.rm = TRUE),
        presence_count = sum(rpm > 0),
        .groups = 'drop'
      )
    
    # Filter out taxa that are consistently very low abundance
    taxa_to_keep <- taxa_stats %>%
      filter(
        max_rpm >= abundance_thresholds$rpm_threshold * 2 |  # At least 2x threshold in one sample
          presence_count >= 2  # Present in at least 2 samples
      ) %>%
      pull(name)
    
    filtered_df <- filtered_df %>%
      filter(name %in% taxa_to_keep)
    
    filter_log <- filter_log %>%
      add_row(
        step = "Remove statistical outliers",
        taxa_before = before_taxa,
        taxa_after = nrow(filtered_df),
        reads_before = before_reads,
        reads_after = sum(filtered_df$reads, na.rm = TRUE)
      )
  }
  
  # Print filtering summary
  message("Quality Filtering Summary (", filter_type, " mode):")
  print(filter_log)
  
  final_retention <- nrow(filtered_df) / initial_taxa * 100
  message("Overall retention: ", round(final_retention, 1), "% of taxa")
  
  return(list(
    filtered_data = filtered_df,
    filter_log = filter_log
  ))
}

# FUNCTION TO CREATE QUALITY-FILTERED ABUNDANCE MATRICES
create_quality_filtered_matrix <- function(df, tax_level, filter_type = "standard", 
                                           threshold_type = "rpm", output_file = NULL) {
  
  # Apply quality filtering first
  filtered_result <- apply_quality_filters(df, filter_type = filter_type)
  df_filtered <- filtered_result$filtered_data
  
  if (nrow(df_filtered) == 0) {
    warning("No taxa remaining after quality filtering!")
    return(NULL)
  }
  
  # Create sample identifier
  df_filtered <- df_filtered %>%
    mutate(sample_id = paste(flowcell, barcode, sep = "_"))
  
  # Apply threshold filtering
  threshold_value <- abundance_thresholds[[paste0(threshold_type, "_threshold")]]
  
  if (threshold_type == "rpm") {
    df_final <- df_filtered %>% filter(rpm >= threshold_value)
  } else if (threshold_type == "reads") {
    df_final <- df_filtered %>% filter(reads >= threshold_value)
  }
  
  # Create abundance matrix
  abundance_matrix <- df_final %>%
    select(sample_id, name, rpm) %>%
    pivot_wider(names_from = name, values_from = rpm, values_fill = 0) %>%
    column_to_rownames("sample_id")
  
  # Add metadata
  sample_metadata <- df_final %>%
    select(sample_id, flowcell, barcode, location, latitude, longitude) %>%
    distinct() %>%
    column_to_rownames("sample_id")
  
  # Combine metadata with abundance data
  full_matrix <- sample_metadata %>%
    rownames_to_column("sample_id") %>%
    left_join(
      abundance_matrix %>% rownames_to_column("sample_id"),
      by = "sample_id"
    ) %>%
    column_to_rownames("sample_id")
  
  # Replace NAs with 0s
  abundance_cols <- !names(full_matrix) %in% c("flowcell", "barcode", "location", "latitude", "longitude")
  full_matrix[abundance_cols][is.na(full_matrix[abundance_cols])] <- 0
  
  # Save to file if specified
  if (!is.null(output_file)) {
    write_tsv(full_matrix %>% rownames_to_column("sample_id"), output_file)
    
    # Also save the filtering log
    log_file <- str_replace(output_file, "\\.tsv$", "_filter_log.tsv")
    write_tsv(filtered_result$filter_log, log_file)
    
    message("Quality-filtered matrix saved to: ", output_file)
    message("Filter log saved to: ", log_file)
  }
  
  return(list(
    matrix = full_matrix,
    filter_log = filtered_result$filter_log
  ))
}

# UPDATED FUNCTION TO GENERATE ALL MATRICES WITH QUALITY FILTERING
generate_quality_filtered_matrices <- function() {
  
  # Create output directory
  matrix_folder <- file.path(output_folder, "abundance_matrices")
  dir.create(matrix_folder, showWarnings = FALSE)
  
  message("Generating quality-filtered abundance matrices...")
  
  # PHYLA LEVEL - Standard and Strict filtering
  phyla_standard <- create_quality_filtered_matrix(
    phyla_df,
    tax_level = "phyla",
    filter_type = "standard",
    threshold_type = "rpm",
    output_file = file.path(matrix_folder, "phyla_abundance_matrix_qc_standard.tsv")
  )
  
  phyla_strict <- create_quality_filtered_matrix(
    phyla_df,
    tax_level = "phyla", 
    filter_type = "strict",
    threshold_type = "rpm",
    output_file = file.path(matrix_folder, "phyla_abundance_matrix_qc_strict.tsv")
  )
  
  # GENUS LEVEL - All genera with quality filtering
  genus_all_standard <- create_quality_filtered_matrix(
    genus_all,
    tax_level = "genus_all",
    filter_type = "standard", 
    threshold_type = "rpm",
    output_file = file.path(matrix_folder, "genus_abundance_matrix_all_qc_standard.tsv")
  )
  
  genus_all_strict <- create_quality_filtered_matrix(
    genus_all,
    tax_level = "genus_all",
    filter_type = "strict",
    threshold_type = "rpm", 
    output_file = file.path(matrix_folder, "genus_abundance_matrix_all_qc_strict.tsv")
  )
  
  # TOP GENERA - Quality filtered
  genus_top_qc <- create_quality_filtered_matrix(
    genus_df,
    tax_level = "genus_top",
    filter_type = "standard",
    threshold_type = "rpm",
    output_file = file.path(matrix_folder, "genus_abundance_matrix_top_qc.tsv")
  )
  
  # HARMFUL TAXA - Quality filtered
  harmful_qc <- create_quality_filtered_matrix(
    harmful_df,
    tax_level = "harmful",
    filter_type = "strict",  # Use strict for potentially harmful organisms
    threshold_type = "rpm",
    output_file = file.path(matrix_folder, "harmful_abundance_matrix_qc.tsv")
  )
  
  # CYANOBACTERIA - Quality filtered
  cyano_qc <- create_quality_filtered_matrix(
    cyanobacteria_df,
    tax_level = "cyanobacteria",
    filter_type = "standard",
    threshold_type = "rpm",
    output_file = file.path(matrix_folder, "cyanobacteria_abundance_matrix_qc.tsv")
  )
  
  message("Quality-filtered abundance matrices generated in: ", matrix_folder)
  message("Filtering logs saved alongside each matrix file.")
}
# FUNCTION TO CREATE ABUNDANCE MATRIX
create_abundance_matrix <- function(df, tax_level, threshold_type = "rpm", 
                                    threshold_value = 1, output_file = NULL) {
  
  # Create sample identifier
  df <- df %>%
    mutate(sample_id = paste(flowcell, barcode, sep = "_"))
  
  # Filter based on threshold
  if (threshold_type == "rpm") {
    df_filtered <- df %>% filter(rpm >= threshold_value)
  } else if (threshold_type == "reads") {
    df_filtered <- df %>% filter(reads >= threshold_value)
  } else if (threshold_type == "pct") {
    df_filtered <- df %>% filter(pct >= threshold_value)
  }
  
  # Create abundance matrix (samples as rows, taxa as columns)
  abundance_matrix <- df_filtered %>%
    select(sample_id, name, rpm) %>%  # using rpm as abundance measure
    pivot_wider(names_from = name, values_from = rpm, values_fill = 0) %>%
    column_to_rownames("sample_id")
  
  # Add metadata columns
  sample_metadata <- df %>%
    select(sample_id, flowcell, barcode, location, latitude, longitude) %>%
    distinct() %>%
    column_to_rownames("sample_id")
  
  # Combine metadata with abundance data
  full_matrix <- sample_metadata %>%
    rownames_to_column("sample_id") %>%
    left_join(
      abundance_matrix %>% rownames_to_column("sample_id"),
      by = "sample_id"
    ) %>%
    column_to_rownames("sample_id")
  
  # Replace NAs with 0s in abundance columns
  abundance_cols <- !names(full_matrix) %in% c("flowcell", "barcode", "location", "latitude", "longitude")
  full_matrix[abundance_cols][is.na(full_matrix[abundance_cols])] <- 0
  
  # Save to file if specified
  if (!is.null(output_file)) {
    write_tsv(full_matrix %>% rownames_to_column("sample_id"), output_file)
    message("Abundance matrix saved to: ", output_file)
  }
  
  return(full_matrix)
}

# FUNCTION TO APPLY PREVALENCE FILTERING
filter_by_prevalence <- function(abundance_matrix, min_prevalence = 0.1) {
  # Identify metadata columns
  metadata_cols <- c("flowcell", "barcode", "location", "latitude", "longitude")
  abundance_cols <- !names(abundance_matrix) %in% metadata_cols
  
  # Calculate prevalence for each taxon
  prevalence <- abundance_matrix[abundance_cols] %>%
    summarise_all(~sum(. > 0) / nrow(.)) %>%
    pivot_longer(everything(), names_to = "taxon", values_to = "prevalence")
  
  # Filter taxa by prevalence
  taxa_to_keep <- prevalence %>%
    filter(prevalence >= min_prevalence) %>%
    pull(taxon)
  
  # Return filtered matrix
  filtered_matrix <- abundance_matrix %>%
    select(all_of(metadata_cols), all_of(taxa_to_keep))
  
  message("Filtered from ", sum(abundance_cols), " to ", length(taxa_to_keep), 
          " taxa (prevalence >= ", min_prevalence, ")")
  
  return(filtered_matrix)
}

# GENERATE ABUNDANCE MATRICES FOR DIFFERENT TAXONOMIC LEVELS
generate_all_abundance_matrices <- function() {
  
  # Create output directory
  matrix_folder <- file.path(output_folder, "abundance_matrices")
  dir.create(matrix_folder, showWarnings = FALSE)
  
  # PHYLA LEVEL MATRIX
  phyla_matrix <- create_abundance_matrix(
    phyla_df, 
    tax_level = "phyla",
    threshold_type = "rpm",
    threshold_value = abundance_thresholds$rpm_threshold,
    output_file = file.path(matrix_folder, "phyla_abundance_matrix.tsv")
  )
  
  # Apply prevalence filtering
  phyla_matrix_filtered <- filter_by_prevalence(
    phyla_matrix, 
    min_prevalence = abundance_thresholds$prevalence_threshold
  )
  write_tsv(
    phyla_matrix_filtered %>% rownames_to_column("sample_id"),
    file.path(matrix_folder, "phyla_abundance_matrix_filtered.tsv")
  )
  
  # GENUS LEVEL MATRIX (TOP GENERA ONLY)
  genus_matrix <- create_abundance_matrix(
    genus_df,
    tax_level = "genus", 
    threshold_type = "rpm",
    threshold_value = abundance_thresholds$rpm_threshold,
    output_file = file.path(matrix_folder, "genus_abundance_matrix_top.tsv")
  )
  
  # Apply prevalence filtering
  genus_matrix_filtered <- filter_by_prevalence(
    genus_matrix,
    min_prevalence = abundance_thresholds$prevalence_threshold
  )
  write_tsv(
    genus_matrix_filtered %>% rownames_to_column("sample_id"),
    file.path(matrix_folder, "genus_abundance_matrix_top_filtered.tsv")
  )
  
  # ALL GENUS LEVEL MATRIX (COMPREHENSIVE)
  all_genus_matrix <- create_abundance_matrix(
    genus_all,  # uses the unfiltered genus_all dataset
    tax_level = "genus_all", 
    threshold_type = "rpm",
    threshold_value = abundance_thresholds$rpm_threshold,
    output_file = file.path(matrix_folder, "genus_abundance_matrix_all.tsv")
  )
  
  # Apply prevalence filtering to comprehensive genus matrix
  all_genus_matrix_filtered <- filter_by_prevalence(
    all_genus_matrix,
    min_prevalence = abundance_thresholds$prevalence_threshold
  )
  write_tsv(
    all_genus_matrix_filtered %>% rownames_to_column("sample_id"),
    file.path(matrix_folder, "genus_abundance_matrix_all_filtered.tsv")
  )
  
  # HARMFUL TAXA MATRIX  
  harmful_matrix <- create_abundance_matrix(
    harmful_df,
    tax_level = "harmful",
    threshold_type = "rpm", 
    threshold_value = abundance_thresholds$rpm_threshold,
    output_file = file.path(matrix_folder, "harmful_abundance_matrix.tsv")
  )
  
  # CYANOBACTERIA MATRIX
  cyano_matrix <- create_abundance_matrix(
    cyanobacteria_df,
    tax_level = "cyanobacteria",
    threshold_type = "rpm",
    threshold_value = abundance_thresholds$rpm_threshold, 
    output_file = file.path(matrix_folder, "cyanobacteria_abundance_matrix.tsv")
  )
  
  # SUMMARY STATISTICS
  generate_matrix_summary(matrix_folder)
  
  message("All abundance matrices generated in: ", matrix_folder)
}

# FUNCTION TO GENERATE SUMMARY STATISTICS
generate_matrix_summary <- function(matrix_folder) {
  
  summary_data <- tibble(
    matrix_type = character(),
    n_samples = integer(),
    n_taxa = integer(),
    n_taxa_filtered = integer(),
    total_abundance = numeric(),
    mean_richness = numeric()
  )
  
  # List all TSV files
  tsv_files <- list.files(matrix_folder, pattern = "\\.tsv$", full.names = TRUE)
  
  for (file in tsv_files) {
    matrix_name <- tools::file_path_sans_ext(basename(file))
    
    # Skip if it's the summary file
    if (grepl("summary", matrix_name)) next
    
    # Read matrix
    mat <- read_tsv(file, show_col_types = FALSE) %>%
      column_to_rownames("sample_id")
    
    # Identify abundance columns (exclude metadata)
    metadata_cols <- c("flowcell", "barcode", "location", "latitude", "longitude")
    abundance_cols <- !names(mat) %in% metadata_cols
    
    # Calculate statistics
    n_samples <- nrow(mat)
    n_taxa <- sum(abundance_cols)
    total_abundance <- sum(mat[abundance_cols], na.rm = TRUE)
    mean_richness <- mean(rowSums(mat[abundance_cols] > 0, na.rm = TRUE))
    
    # Add to summary
    summary_data <- summary_data %>%
      add_row(
        matrix_type = matrix_name,
        n_samples = n_samples,
        n_taxa = n_taxa,
        total_abundance = round(total_abundance, 2),
        mean_richness = round(mean_richness, 2)
      )
  }
  
  # Save summary
  write_tsv(summary_data, file.path(matrix_folder, "matrix_summary.tsv"))
  
  # Print summary to console
  print(summary_data)
}

# FUNCTION TO CREATE DOWNLOAD SECTION FOR HTML DASHBOARD
create_download_section <- function() {
  
  download_links <- tagList(
    tags$h3("Download Data Matrices", style = "margin-top: 30px; border-top: 2px solid #ccc; padding-top: 20px;"),
    
    # Quality-filtered matrices section
    tags$div(
      style = "background: #e8f5e8; padding: 15px; border-radius: 8px; margin: 10px 0; border-left: 4px solid #28a745;",
      tags$h4("🔬 Quality-Filtered Matrices (Recommended)", style = "margin-top: 0; color: #28a745;"),
      tags$div(
        style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(280px, 1fr)); gap: 15px;",
        
        # QC Phyla
        tags$div(
          tags$h5("Phyla Level (QC)", style = "color: #2c3e50; margin-bottom: 8px;"),
          tags$div(
            tags$a("✅ Phyla Standard QC", 
                   href = "abundance_matrices/phyla_abundance_matrix_qc_standard.tsv", 
                   download = "phyla_abundance_matrix_qc_standard.tsv",
                   style = "display: block; padding: 6px; margin: 3px 0; background: #d4edda; text-decoration: none; border-radius: 4px; font-size: 0.9em;"),
            tags$a("✅ Phyla Strict QC", 
                   href = "abundance_matrices/phyla_abundance_matrix_qc_strict.tsv", 
                   download = "phyla_abundance_matrix_qc_strict.tsv",
                   style = "display: block; padding: 6px; margin: 3px 0; background: #c3e6cb; text-decoration: none; border-radius: 4px; font-size: 0.9em;")
          )
        ),
        
        # QC Genus
        tags$div(
          tags$h5("Genus Level (QC)", style = "color: #2c3e50; margin-bottom: 8px;"),
          tags$div(
            tags$a("✅ ALL Genera Standard QC", 
                   href = "abundance_matrices/genus_abundance_matrix_all_qc_standard.tsv", 
                   download = "genus_abundance_matrix_all_qc_standard.tsv",
                   style = "display: block; padding: 6px; margin: 3px 0; background: #d4edda; text-decoration: none; border-radius: 4px; font-size: 0.9em; font-weight: bold;"),
            tags$a("✅ ALL Genera Strict QC", 
                   href = "abundance_matrices/genus_abundance_matrix_all_qc_strict.tsv", 
                   download = "genus_abundance_matrix_all_qc_strict.tsv",
                   style = "display: block; padding: 6px; margin: 3px 0; background: #c3e6cb; text-decoration: none; border-radius: 4px; font-size: 0.9em; font-weight: bold;"),
            tags$a("✅ Top 20 Genera QC", 
                   href = "abundance_matrices/genus_abundance_matrix_top_qc.tsv", 
                   download = "genus_abundance_matrix_top_qc.tsv",
                   style = "display: block; padding: 6px; margin: 3px 0; background: #d4edda; text-decoration: none; border-radius: 4px; font-size: 0.9em;")
          )
        ),
        
        # QC Specialized
        tags$div(
          tags$h5("Specialized (QC)", style = "color: #2c3e50; margin-bottom: 8px;"),
          tags$div(
            tags$a("✅ Harmful Taxa QC", 
                   href = "abundance_matrices/harmful_abundance_matrix_qc.tsv", 
                   download = "harmful_abundance_matrix_qc.tsv",
                   style = "display: block; padding: 6px; margin: 3px 0; background: #fff3cd; text-decoration: none; border-radius: 4px; font-size: 0.9em;"),
            tags$a("✅ Cyanobacteria QC", 
                   href = "abundance_matrices/cyanobacteria_abundance_matrix_qc.tsv", 
                   download = "cyanobacteria_abundance_matrix_qc.tsv",
                   style = "display: block; padding: 6px; margin: 3px 0; background: #d1ecf1; text-decoration: none; border-radius: 4px; font-size: 0.9em;")
          )
        )
      )
    ),
    
    # Original matrices section
    tags$div(
      style = "display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; margin: 20px 0;",
      
      # Phyla matrices
      tags$div(
        tags$h4("Phyla Level (Original)", style = "color: #6c757d; margin-bottom: 10px;"),
        tags$div(
          tags$a("📊 Top Phyla Matrix", 
                 href = "abundance_matrices/phyla_abundance_matrix.tsv", 
                 download = "phyla_abundance_matrix.tsv",
                 style = "display: block; padding: 8px; margin: 5px 0; background: #ecf0f1; text-decoration: none; border-radius: 4px;"),
          tags$a("📊 Top Phyla (Filtered)", 
                 href = "abundance_matrices/phyla_abundance_matrix_filtered.tsv", 
                 download = "phyla_abundance_matrix_filtered.tsv",
                 style = "display: block; padding: 8px; margin: 5px 0; background: #ecf0f1; text-decoration: none; border-radius: 4px;")
        )
      ),
      
      # Genus matrices  
      tags$div(
        tags$h4("Genus Level (Original)", style = "color: #6c757d; margin-bottom: 10px;"),
        tags$div(
          tags$a("📊 Top 20 Genera Matrix", 
                 href = "abundance_matrices/genus_abundance_matrix_top.tsv", 
                 download = "genus_abundance_matrix_top.tsv",
                 style = "display: block; padding: 8px; margin: 5px 0; background: #ecf0f1; text-decoration: none; border-radius: 4px;"),
          tags$a("📊 Top 20 Genera (Filtered)", 
                 href = "abundance_matrices/genus_abundance_matrix_top_filtered.tsv", 
                 download = "genus_abundance_matrix_top_filtered.tsv",
                 style = "display: block; padding: 8px; margin: 5px 0; background: #ecf0f1; text-decoration: none; border-radius: 4px;"),
          tags$a("📊 ALL Genera Matrix", 
                 href = "abundance_matrices/genus_abundance_matrix_all.tsv", 
                 download = "genus_abundance_matrix_all.tsv",
                 style = "display: block; padding: 8px; margin: 5px 0; background: #d5e8d4; text-decoration: none; border-radius: 4px; font-weight: bold;"),
          tags$a("📊 ALL Genera (Filtered)", 
                 href = "abundance_matrices/genus_abundance_matrix_all_filtered.tsv", 
                 download = "genus_abundance_matrix_all_filtered.tsv",
                 style = "display: block; padding: 8px; margin: 5px 0; background: #d5e8d4; text-decoration: none; border-radius: 4px; font-weight: bold;")
        )
      ),
      
      # Specialized matrices
      tags$div(
        tags$h4("Specialized (Original)", style = "color: #6c757d; margin-bottom: 10px;"),
        tags$div(
          tags$a("⚠️ Harmful Taxa Matrix", 
                 href = "abundance_matrices/harmful_abundance_matrix.tsv", 
                 download = "harmful_abundance_matrix.tsv",
                 style = "display: block; padding: 8px; margin: 5px 0; background: #ffeaa7; text-decoration: none; border-radius: 4px;"),
          tags$a("🦠 Cyanobacteria Matrix", 
                 href = "abundance_matrices/cyanobacteria_abundance_matrix.tsv", 
                 download = "cyanobacteria_abundance_matrix.tsv",
                 style = "display: block; padding: 8px; margin: 5px 0; background: #81ecec; text-decoration: none; border-radius: 4px;")
        )
      ),
      
      # Summary and metadata
      tags$div(
        tags$h4("Summary & Info", style = "color: #2c3e50; margin-bottom: 10px;"),
        tags$div(
          tags$a("📋 Matrix Summary", 
                 href = "abundance_matrices/matrix_summary.tsv", 
                 download = "matrix_summary.tsv",
                 style = "display: block; padding: 8px; margin: 5px 0; background: #dda0dd; text-decoration: none; border-radius: 4px;"),
          tags$a("🗺️ Sample Metadata", 
                 href = "../metadata.csv", 
                 download = "metadata.csv",
                 style = "display: block; padding: 8px; margin: 5px 0; background: #dda0dd; text-decoration: none; border-radius: 4px;")
        )
      )
    ),
    
    # Info boxes
    tags$div(
      style = "background: #e8f5e8; padding: 15px; border-radius: 8px; margin: 20px 0; border-left: 4px solid #28a745;",
      tags$h4("Quality Control Information", style = "margin-top: 0; color: #28a745;"),
      tags$p("• ", tags$strong("Standard QC:"), " Minimum ", abundance_thresholds$read_threshold, " reads, removes contaminants"),
      tags$p("• ", tags$strong("Strict QC:"), " Minimum 20 reads, taxonomic depth ≥3 levels, statistical outlier removal"),
      tags$p("• ", tags$strong("Confidence filtering:"), " Based on k-mer agreement (recommended: 0.1+)"),
      tags$p("• ", tags$strong("Filter logs:"), " Detailed filtering statistics saved with each QC matrix")
    ),
    
    tags$div(
      style = "background: #f8f9fa; padding: 15px; border-radius: 8px; margin: 20px 0; border-left: 4px solid #007bff;",
      tags$h4("Matrix Information", style = "margin-top: 0; color: #007bff;"),
      tags$p("• ", tags$strong("Abundance values:"), " Reads per million (RPM)"),
      tags$p("• ", tags$strong("Quality threshold:"), " Minimum ", abundance_thresholds$rpm_threshold, " RPM"),
      tags$p("• ", tags$strong("Filtered matrices:"), " Taxa present in ≥", abundance_thresholds$prevalence_threshold * 100, "% of samples"),
      tags$p("• ", tags$strong("Matrix format:"), " Samples as rows, taxa as columns"),
      tags$p("• ", tags$strong("ALL Genera:"), " Comprehensive genus-level data (all taxa above threshold)")
    )
  )
  
  return(download_links)
}