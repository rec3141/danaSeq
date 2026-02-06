# eDNA Mapping Pipeline in R
# ---------------------------
# This script processes Kraken reports and metadata to generate interactive maps using leaflet and ggplot2.
# Output includes phyla, genera, harmful categories, and cyanobacteria data.

library(tidyverse)
library(leaflet)
library(ggplot2)
library(png)
library(htmltools)
library(htmlwidgets)
library(scales)
library(DBI)
library(duckdb)
library(crosstalk)
library(readxl)

# CONFIGURATION
# INPUT 1) BARCODE DIRECTORY
#args <- commandArgs(trailingOnly = TRUE)
#try({setwd(args[1])})

prjdir="/data/project_QEI2025/"
setwd(prjdir)
print(getwd())

#"." for all; "LWW"         "GENICE-SERF" "GENICE-CMO"  "GENICE-CI"
project = "QEI"

# CATEGORY MAPPING
category_map <- list(
  "Homo" = "Human Contamination",
  "Escherichia" = "Fecal Indicators",
  "Enterococcus" = "Fecal Indicators",
  "Shigella" = "Fecal Indicators",
  "Salmonella" = "Fecal Indicators",
  "Clostridium" = "Fecal Indicators",
  "Pseudomonas" = "Opportunistic Pathogens",
  "Acinetobacter" = "Opportunistic Pathogens",
  "Burkholderia" = "Opportunistic Pathogens",
  "Staphylococcus" = "Opportunistic Pathogens",
  "Streptococcus" = "Opportunistic Pathogens",
  "Vibrio" = "Waterborne Pathogens",
  "Campylobacter" = "Waterborne Pathogens",
  "Aeromonas" = "Waterborne Pathogens",
  "Legionella" = "Waterborne Pathogens",
  "Leptospira" = "Waterborne Pathogens",
  "Francisella" = "Waterborne Pathogens",
  "Brucella" = "Zoonotic Bacteria",
  "Pasteurella" = "Zoonotic Bacteria",
  "Mycobacterium" = "Zoonotic Bacteria",
  "Listeria" = "Zoonotic Bacteria",
  "Yersinia" = "Zoonotic Bacteria",
  "Corynebacterium" = "Zoonotic Bacteria",
  "Neisseria" = "Respiratory/Inhaled",
  "Haemophilus" = "Respiratory/Inhaled",
  "Bordetella" = "Respiratory/Inhaled",
  "Chlamydia" = "Respiratory/Inhaled",
  "Treponema" = "Respiratory/Inhaled"
)

source_tracking <- list(
  
  # Food/Agricultural crops - clear laboratory contamination
  food_agricultural = c("Capsicum", "Arachis", "Asparagus", "Benincasa", "Brassica", 
                        "Cajanus", "Camelina", "Cannabis","Carya", "Coffea", "Cucumis", "Cynara", 
                        "Dioscorea", "Elaeis", "Glycine", "Gossypium", "Helianthus", 
                        "Hevea", "Hordeum", "Humulus", "Ipomoea", "Juglans", "Lactuca", 
                        "Lolium", "Lupinus", "Lycium", "Macadamia", "Magnolia", "Malania", 
                        "Mangifera", "Miscanthus", "Musa", "Nicotiana", "Olea", "Oryza", 
                        "Panicum", "Papaver", "Phragmites", "Physcomitrium", "Pisum", 
                        "Populus", "Prunus", "Punica", "Pyrus", "Quercus", "Rhododendron", 
                        "Ricinus", "Rosa", "Rutidosis", "Salvia", "Solanum", "Sorghum", 
                        "Telopea", "Tripterygium", "Triticum", "Vicia", "Vigna", "Vitis", 
                        "Zea", "Zingiber"),
  
  # Human-associated microbes - lab worker contamination
  human_associated = c("Homo", "Staphylococcus", "Streptococcus", "Cutibacterium", 
                       "Corynebacterium", "Enterococcus", "Mammaliicoccus", "Schaalia",
                       "Leptotrichia", "Neisseria", "Moraxella", "Bacteroides"),
  
  # Clinical/pathogenic bacteria - lab contamination from medical samples
  clinical_pathogenic = c("Escherichia", "Salmonella", "Yersinia", "Citrobacter", 
                          "Morganella", "Providencia", "Klebsiella", "Campylobacter", 
                          "Helicobacter", "Legionella", "Bordetella", "Francisella", 
                          "Pasteurella", "Bartonella", "Mycobacterium", "Clostridium", 
                          "Fusobacterium"),
  
  # Soil/terrestrial bacteria - field contamination or lab cross-contamination
  soil_terrestrial = c("Bacillus", "Arthrobacter", "Microbacterium", "Nocardioides", 
                       "Streptomyces", "Paenibacillus", "Peribacillus", "Parageobacillus", 
                       "Glutamicibacter", "Agromyces", "Promicromonospora", "Nonomuraea", 
                       "Amycolatopsis", "Blastococcus", "Massilia", "Variovorax", 
                       "Hydrogenophaga", "Polaromonas", "Methylibium", "Methylobacterium", 
                       "Microvirga", "Bradyrhizobium", "Rhizobium", "Mesorhizobium", 
                       "Ensifer", "Herbaspirillum", "Azospirillum", "Azospira"),
  
  # Plant-associated bacteria/symbionts - terrestrial plant contamination
  plant_associated = c("Azoarcus", "Burkholderia", "Paraburkholderia", "Ralstonia", 
                       "Cupriavidus", "Pandoraea", "Methylobacter", "Xanthomonas", 
                       "Dickeya", "Erwinia", "Candidatus Phytoplasma", "Blattabacterium", 
                       "Buchnera", "Eremothecium", "Kwoniella"),
  
  # Reagent/kit contaminants - DNA extraction or PCR reagents
  reagent_kit = c("Chryseobacterium", "Sphingobacterium", "Pedobacter", 
                  "Mucilaginibacter", "Myroides", "Sphingomonas", 
                  "Sphingobium", "Sphingopyxis", "Brevundimonas", "Caulobacter", 
                  "Undibacterium", "Methylotenera", "Methylophilus"),
  
  # Legitimate marine bacteria - should NOT be considered contaminants
  marine = c("Candidatus Pelagibacter", "Candidatus Nitrosomarinus", 
                        "Candidatus Nitrosopelagicus", "Candidatus Thioglobus", 
                        "Candidatus Arcticimaribacter", "Alteromonas", "Pseudoalteromonas", 
                        "Vibrio", "Photobacterium", "Flavobacterium", "Shewanella", "Marinobacter", 
                        "Alcanivorax", "Halomonas", "Oceanospirillum", "Colwellia", 
                        "Psychromonas", "Psychrobacter", "Moritella", "Pseudomonas", 
                        "Acinetobacter", "Roseobacter", "Ruegeria", "Sulfitobacter", 
                        "Roseovarius", "Phaeobacter", "Loktanella", "Erythrobacter", 
                        "Paracoccus", "Planktomarina", "Maribacter", "Polaribacter", 
                        "Formosa", "Flaviramulus", "Winogradskyella", "Zobellia", 
                        "Cellulophaga", "Tenacibaculum", "Polaribacter", "Algibacter", 
                        "Dokdonia", "Kordia", "Mariniflexile", "Nonlabens", "Olleya"),
  
  # Environmental/freshwater - possible freshwater contamination
  freshwater_environmental = c("Polynucleobacter", "Limnohabitans", "Aeromonas", 
                               "Arcobacter", "Comamonas", "Hydrogenophaga", "Methylomonas", 
                               "Methylococcus", "Nitrosomonas", "Nitrosospira", 
                               "Thauera", "Aromatoleum", "Azoarcus")
)

# Function to check which category a genus belongs to
categorize_genus <- function(genus_name) {
  for (category in names(source_tracking)) {
    if (genus_name %in% source_tracking[[category]]) {
      return(category)
    }
  }
  return("unclassified")
}

# Example usage:
# categorize_genus("Capsicum")  # Returns "food_agricultural"
# categorize_genus("Homo")      # Returns "human_associated"


# KRAKEN IMPORT
read_kraken <- function(con) {
  # DB CONNECTION
  con <- dbConnect(
    duckdb::duckdb(),
    dbdir = file.path(rundir, barcode, "dana.duckdb"),
    read_only = TRUE
  )
  kraken = dbGetQuery(con, "SELECT * FROM kraken") %>% column_to_rownames("seqid") %>% mutate(genus = word(sub(
    pattern = "unclassified ",
    replacement = "",
    sub(pattern = "Candidatus ", replacement = "", taxa_name)
  ), 1))
  
  dbDisconnect(con, shutdown = TRUE)
  return(kraken)
}

#LOAD AND COMBINE KRAKEN REPORTS WITH INDENT-BASED TAXONOMY
read_kraken_report <- function(con) {
  # DB CONNECTION
  con <- dbConnect(
    duckdb::duckdb(),
    dbdir = file.path(rundir, barcode, "dana.duckdb"),
    read_only = TRUE
  )
  
  kraken_report = dbGetQuery(con, "SELECT * FROM krakenreport")
  
  dbDisconnect(con, shutdown = TRUE)
  return(kraken_report)
}


#LOAD PROKKA DATA
read_prokka <- function(con, gene) {
  # DB CONNECTION
  con <- dbConnect(
    duckdb::duckdb(),
    dbdir = file.path(rundir, barcode, "dana.duckdb"),
    read_only = TRUE
  )
  
  prokka_result <- dbGetQuery(
    con,
    "
SELECT
  k.taxa_name,
  p.product,
  COUNT(*) AS product_count
FROM
  prokka_annotations p
JOIN
  locus_index l ON p.locus_tag = l.locus_tag
JOIN
  kraken k ON l.seqid = k.seqid
WHERE
  p.product IS NOT NULL AND p.product <> ''
GROUP BY
  k.taxa_name, p.product
ORDER BY
  k.taxa_name, product_count DESC;
"
  )
  
  
  dbDisconnect(con, shutdown = TRUE)
  return(prokka_result)
}


# FILTER AND SUMMARIZE
summarize_abundance <- function(df, tax_filter, rank_col) {
  classified_df <- df %>%
    filter(rank != "U") %>%
    group_by(flowcell, barcode) %>%
    summarise(classified = sum(direct), .groups = "drop")
  
  df %>%
    filter(rank == rank_col,
           str_detect(taxonomy, tax_filter, negate = FALSE)) %>%
    group_by(name, flowcell, barcode, Category) %>%
    summarise(
      totalreads = mean(as.numeric(totalreads)),
      reads = sum(as.numeric(reads), na.rm = TRUE),
      .groups = 'drop'
    ) %>%
    left_join(meta, by = c("flowcell", "barcode")) %>%
    mutate(rpm = reads / totalreads * 1e6,
           pct = reads / totalreads * 100) %>%
    left_join(classified_df, by = c("flowcell", "barcode")) %>%
    mutate(
      rpm_classified = reads / classified * 1e6,
      pct_classified = reads / classified * 100
    )
}

read_count <- function(con) {
  con <- dbConnect(
    duckdb::duckdb(),
    dbdir = file.path(rundir, barcode, "dana.duckdb"),
    read_only = TRUE
  )
  #  reads = dbGetQuery(con, "SELECT COUNT(*) AS n FROM stats") %>% as.numeric()
  reads = dbGetQuery(
    con,
    "SELECT
  SPLIT_PART(fileid, '_', 1) AS flowcell,
  SPLIT_PART(fileid, '_', 3) AS barcode FROM sequence_index"
  )
  
  dbDisconnect(con, shutdown = TRUE)
  return(reads)
  
}


# LOAD METADATA
#metadata_file <- "~/Desktop/Amundsen-Collins/QEI_metadata.csv"
metadata_file <- "~/Desktop/Amundsen-Collins/QEI_SampleLog_2025_repaired.xlsx"
meta <- read_excel(metadata_file, sheet = "nanopore log")

#meta <- read_csv(metadata_file)
meta <- meta[!is.na(meta$flowcell) & !is.na(meta$barcode), ]
if (!(project == "."))  meta <- meta %>% filter(project == !!project)

meta$longitude = meta$longitude + sample(1:1000, nrow(meta)) / 1e6 #add jitter so we can zoom and see differences
meta$latitude = meta$latitude + sample(1:1000, nrow(meta)) / 1e6 #add jitter so we can zoom and see differences
meta$latitude[is.na(meta$latitude)] = 90
fakelong = -100 + as.numeric(as.factor(meta$flowcell)) + as.numeric(as.factor(meta$barcode)) /
  max(as.numeric(as.factor(meta$barcode))) / 3
meta$longitude[is.na(meta$longitude)] = fakelong[is.na(meta$longitude)]

output_folder <- file.path("mapping", project)
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

kraken_data <- tibble()
prokka_data <- tibble()

#dirlist = list.dirs()

host <- "concentration"; 
root <- paste0(prjdir,"out_dana/") # must have trailing slash

#cmd <- sprintf("/usr/bin/rsync -rni --progress --out-format='%%f' --contimeout=20 --timeout=120 -e 'ssh -o BatchMode=yes -o StrictHostKeyChecking=accept-new ' '%s:%s' '.'", host, root)
cmd <- sprintf("/usr/bin/rsync -rni --out-format='%%f' --timeout=10 %s:%s .", host, root)
print("fetching file list")
dirlist <- system(cmd, intern = TRUE)

# LOOP EACH BARCODE IN EACH RUN IN PROJECT
for (flowcell in unique(meta$flowcell)) {
  rundir = file.path(root,dirlist[grepl(paste0(flowcell, "$"), dirlist)])
  print(rundir)
  if (length(rundir) < 1)
    next
  for (barcode in unique(meta$barcode[meta$flowcell == flowcell])) {
    print(c(flowcell, barcode))
    if (!file.exists(file.path(rundir, barcode, "dana.duckdb"))) {
      print(paste0((file.path(rundir, barcode, "dana.duckdb")), " not found, skipping"))
      next
    }
    kraken_folder <- file.path(rundir, barcode, "kraken")
    
    kraken_genus = tryCatch({
      read_kraken()
    }, error = function(e) {
      cat("Skipping due to error:", e$message, "\n")
      return(NULL)  # Skip to next iteration instead of breaking
    })
    
    if (is.null(kraken_genus)) next
    
    if (nrow(kraken_genus) > 0) {
      totalreads_barcode = tryCatch({
        read_count() %>% filter(flowcell == flowcell &
                                  barcode == barcode) %>% nrow();
      }, error = function(e) {
        cat("Skipping due to error:", e$message, "\n")
        return(NULL)  # Skip to next iteration instead of breaking
      })
    }
    if (is.null(totalreads_barcode)) next

    print(totalreads_barcode)
    kraken_barcode <- tryCatch({
      read_kraken_report()
    }, error = function(e) {
      cat("Skipping due to error:", e$message, "\n")
      return(NULL)  # Skip to next iteration instead of breaking
    })
    
    if (is.null(kraken_barcode)) next

        if (nrow(kraken_barcode) > 0) {
      kraken_barcode$totalreads = totalreads_barcode
      kraken_barcode$Category = NA
      kraken_data <- bind_rows(kraken_data, kraken_barcode)
    }
    
    prokka_barcode = tryCatch({
      read_prokka()
    }, error = function(e) {
      cat("Skipping due to error:", e$message, "\n")
      return(NULL)  # Skip to next iteration instead of breaking
    })
    
    if (is.null(prokka_barcode)) next
    
    if (nrow(prokka_barcode) > 0) {
      prokka_barcode <- prokka_barcode %>% filter(str_detect(product, "hypothetical protein", negate =
                                                               TRUE))
      if (nrow(prokka_barcode) > 0) {
        prokka_barcode$flowcell = flowcell
        prokka_barcode$barcode = barcode
        prokka_data <- bind_rows(prokka_data, prokka_barcode)
      }
      
    }
  }
  
}

kraken_data <- kraken_data %>% filter(totalreads > 100)

phyla_df <- summarize_abundance(kraken_data, ".*", "P")

phyla_top <- phyla_df %>%
  group_by(name) %>%
  summarise(total_reads = sum(reads, na.rm = TRUE)) %>%
  arrange(desc(total_reads)) %>%
  slice_head(n = 8) %>%
  pull(name)

phyla_df <- phyla_df %>%
  filter(name %in% phyla_top) %>%
  mutate(Category = name)

genus_all <- summarize_abundance(kraken_data, ".*", "G")

genus_p1 <- genus_all %>% filter(rpm_classified > 10000)
genus_p01 <- genus_all %>% filter(rpm_classified > 1000)
genus_p001 <- genus_all %>% filter(rpm_classified > 100)

genus_allowed <- genus_p01 %>% mutate(source = sapply(name, categorize_genus)) %>% filter(! source %in% c("human_associated", "food_agricultural","clinical_pathogenic","reagent_kit")) %>% pull(name) %>% unique()

genus_top <- genus_all %>%
  group_by(name) %>%
  summarise(total_reads = sum(reads, na.rm = TRUE)) %>%
  arrange(desc(total_reads)) %>%
  filter(!name == "Homo") %>%
  slice_head(n = 20) %>%
  pull(name)

genus_df <- genus_all %>%
  filter(name %in% genus_top) %>%
  mutate(Category = name)


# Cyanobacteria identification using reconstructed taxonomy string
cyanobacteria_df <- summarize_abundance(kraken_data, "Cyanobacteriota", "F")

cyanobacteria_top <- cyanobacteria_df %>%
  group_by(name) %>%
  summarise(total_reads = sum(reads, na.rm = TRUE)) %>%
  arrange(desc(total_reads)) %>%
  slice_head(n = 8) %>%
  pull(name)

cyanobacteria_df <- cyanobacteria_df %>%
  filter(name %in% cyanobacteria_top) %>%
  mutate(Category = name)

# Potentially harmful taxa identification using reconstructed taxonomy string

harmful_regex <- paste0(names(category_map), "$", collapse = "|")
harmful_df <- summarize_abundance(kraken_data %>% mutate(Category = recode(name, !!!category_map, .default = "Other")), harmful_regex, "G")

# SCALING FUNCTION
rescale_radius <- function(x, min_r = 0, max_r = 20) {
  rng <- range(x, na.rm = TRUE)
  rng[1] = 0
  scales::rescale(x, to = c(min_r, max_r), from = rng)
}

rescale_radius_genus <- function(x, min_r = 0, max_r = 50) {
  rng <- range(x^.5, na.rm = TRUE)
  rng[1] = 0
  scales::rescale(x^.5, to = c(min_r, max_r), from = rng)
}


# MAP FUNCTION
make_map <- function(df,
                     value_col,
                     popup_label,
                     file_name,
                     color = "darkred") {
  df <- df %>%
    mutate(radius = rescale_radius(.data[[value_col]]))
  
  leaflet(df) %>%
    addProviderTiles("CartoDB.Positron") %>%
    addCircleMarkers(
      ~ longitude,
      ~ latitude,
      radius = ~ radius,
      fillColor = color,
      color = color,
      fillOpacity = 0.3,
      popup = paste(df$location, "<br>", paste0(popup_label, ": ", round(df[[value_col]], 0)))
    ) %>%
    saveWidget(file = file.path(output_folder, file_name),
               selfcontained = FALSE)
}



# GENERATE CIRCLE MAPS
make_map(
  cyanobacteria_df %>% group_by(flowcell, barcode, latitude, longitude, location) %>% summarise(rpm = sum(rpm), .groups = 'drop'),
  "rpm",
  "Cyanobacteria (per million)",
  "cyanobacteria_map.html",
  color = "teal"
)
make_map(
  harmful_df %>% group_by(flowcell, barcode, latitude, longitude, location) %>% summarise(rpm = sum(rpm), .groups = 'drop'),
  "rpm",
  "Concerning (per million)",
  "harmful_map.html",
  color = "red"
)
# make_map(phyla_df %>% group_by(barcode, latitude, longitude) %>% summarise(rpm = sum(rpm), .groups = 'drop'),
#          "rpm", "Phyla RPM", "phyla_map.html", color = "navy")
make_map(
  phyla_df %>% group_by(flowcell, barcode, latitude, longitude, location) %>% summarise(totalreads = mean(totalreads), .groups = 'drop'),
  "totalreads",
  "Data points (reads)",
  "phyla_map.html",
  color = "navy"
)

make_map(
  genus_df %>% group_by(flowcell, barcode, latitude, longitude, location) %>% summarise(totalreads = mean(totalreads), .groups = 'drop'),
  "totalreads",
  "Data points (reads)",
  "genus_map.html",
  color = "navy"
)

# STATIC PIE GENERATOR
save_pie_images <- function(df, label_col, value_col, folder) {
  dir.create(folder, showWarnings = FALSE)
  df %>% mutate(fcbc = paste0(flowcell, "-", barcode)) %>%
    split(.$fcbc) %>% purrr::walk(function(d) {
      fname <- file.path(folder, paste0(d$flowcell[1], "_", d$barcode[1], ".png"))
      ggplot(d, aes(
        x = "",
        y = !!sym(value_col),
        fill = !!sym(label_col)
      )) +
        geom_bar(stat = "identity", width = 1) +
        coord_polar("y") +
        theme_void() +
        theme(legend.position = "none")
      ggsave(
        fname,
        width = 1,
        height = 1,
        dpi = 100,
        create.dir = TRUE
      )
    })
}

# PIE MAP FUNCTION (STATIC ICONS)
make_pie_map <- function(df,
                         value_col,
                         label_col,
                         file_name,
                         pie_dir) {
  save_pie_images(df, label_col, value_col, pie_dir)
  points <- df %>% group_by(flowcell, barcode, latitude, longitude, location) %>% slice(1) %>% ungroup()
  leaflet(points) %>%
    addProviderTiles("CartoDB.Positron") %>%
    addMarkers(
      ~ longitude,
      ~ latitude,
      popup = ~ location,
      icon = ~ icons(
        iconUrl = file.path(pie_dir, paste0(flowcell, "_", barcode, ".png")),
        iconWidth = 60,
        iconHeight = 60
      )
    ) %>%
    # addCircleMarkers(~longitude, ~latitude, radius = 30,
    #                  fillColor = "white", color = "white", fillOpacity = 0,
    #                  popup = df$location) %>%
    saveWidget(file = file.path(output_folder, file_name),
               selfcontained = FALSE)
}

# GENERATE PIE MAPS
make_pie_map(
  harmful_df,
  "rpm",
  "Category",
  "harmful_pie_map.html",
  pie_dir = file.path(output_folder, "pies_harmful")
)
make_pie_map(
  phyla_df,
  "reads",
  "Category",
  "phyla_pie_map.html",
  pie_dir = file.path(output_folder, "pies_phyla")
)
make_pie_map(
  cyanobacteria_df,
  "rpm",
  "Category",
  "cyanobacteria_pie_map.html",
  pie_dir = file.path(output_folder, "pies_cyano")
)
make_pie_map(
  genus_df,
  "reads",
  "Category",
  "genus_pie_map.html",
  pie_dir = file.path(output_folder, "pies_genus")
)


# LEGEND PIES
legend_pie <- function(labels, file, colors = NULL) {
  df <- tibble(label = labels, value = 1)
  if (is.null(colors))
    colors <- scales::hue_pal()(length(labels))
  p <- ggplot(df, aes(x = "", y = value, fill = label)) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar("y") +
    scale_fill_manual(values = colors) +
    theme_void() +
    theme(legend.position = "right")
  ggsave(
    file,
    plot = p,
    width = 4,
    height = 6,
    dpi = 150
  )
}

legend_pie(unique(phyla_df$Category),
           file.path(output_folder, "legend_phyla.png"))
legend_pie(unique(genus_df$Category),
           file.path(output_folder, "legend_genus.png"))
legend_pie(unique(harmful_df$Category),
           file.path(output_folder, "legend_harmful.png"))
legend_pie(unique(cyanobacteria_df$Category),
           file.path(output_folder, "legend_cyano.png"))

# genus picker
name_order <- genus_all %>%
  filter(name %in% genus_allowed) %>%
  group_by(name) %>%
  summarise(total = sum(pct_classified), .groups = "drop") %>%
  arrange(desc(total)) %>%  # or arrange(total) for ascending
  pull(name)

gd <- genus_all %>% # one row = one genus-per-sample
  filter(name %in% genus_allowed) %>%
  mutate(radius = rescale_radius_genus(pct_classified))           # will drive circle size

# Convert to factor with levels ordered by the numeric value
gd$name <- factor(gd$name, levels = name_order)



sd_genus <- SharedData$new(gd,
                           key = ~ paste(flowcell, barcode, name),
                           group = "genus_filter")

map_genus <- leaflet(sd_genus) %>% addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(
    ~ longitude,
    ~ latitude,
    radius = ~ radius,
    fillColor = "darkblue",
    color = "darkblue",
    fillOpacity = 0.35,
    popup = ~ paste0(location, "<br><b>", name, "</b>", "<br>%: ", round(pct_classified, 3))
  )

# can't save to iframe because of picker?
# %>%
#   saveWidget(file = file.path(output_folder, "map_genus.html"),
#              selfcontained = FALSE)

genus_picker <- filter_select(
  id          = "pick_genus",
  label       = "Choose a genus to display",
  sharedData  = sd_genus,
  group       = ~ name,
  multiple    = FALSE,
  # one genus at a time
  allLevels  = TRUE            # include levels even if currently filtered out
  #  selectize   = TRUE             # Selectize gives you a type-ahead search box
)

# taxa plus product picker
# ── collapse to one row per (sample, taxon, product) ─────────────
prokka_df <- prokka_data %>%                     # your tibble shown above
  group_by(flowcell, barcode, taxa_name, product) %>%
  summarise(product_count = sum(product_count), .groups = "drop")

# ── bring in lat/lon (reuse whatever you joined for the genus map) ──
prokka_df <- prokka_df %>%
  left_join(
    meta,
    # must contain flowcell, barcode, latitude, longitude, location
    by = c("flowcell", "barcode")
  )

# ── translate counts → marker radius (reuse your helper) ──────────
prokka_df <- prokka_df %>%
  mutate(radius = rescale_radius(log1p(product_count)))      #  e.g. 0‒20 px

sd_prokka <- SharedData$new(
  prokka_df %>% head(),
  key   = ~ paste(flowcell, barcode, taxa_name, product),
  # unique row id
  group = "prokka_filter"                                    # anything you like
)

map_prokka <- leaflet(sd_prokka) %>%
  addProviderTiles("CartoDB.Positron") %>%
  addCircleMarkers(
    ~ longitude,
    ~ latitude,
    radius      = ~ radius,
    fillColor   = "steelblue",
    # pick any palette
    color       = "steelblue",
    fillOpacity = 0.35,
    popup = ~ sprintf(
      "<b>%s</b><br>%s<br>Count: %d<br>%s / %s",
      taxa_name,
      product,
      product_count,
      flowcell,
      barcode
    )
  )

picker_prokka_taxa <- filter_select(
  id         = "pick_taxa",
  label      = "Taxon",
  sharedData = sd_prokka,
  group      = ~ taxa_name,
  multiple   = FALSE
)

picker_prokka_prod <- filter_select(
  id         = "pick_prod",
  label      = "Product",
  sharedData = sd_prokka,
  group      = ~ product,
  multiple   = FALSE
)

htmltools::save_html(
  browsable(tagList(
    tags$html(
    tags$title("eDNA Mapping Dashboard"),
    tags$head(
      tags$style(HTML("
    .selectize-input {
      min-width: 480px !important;
      width: 480px !important;
    }
    .selectize-dropdown {
      min-width: 480px !important;
      width: 480px !important;
    }
  "))),
    tags$h3("Metagenomic product counts by sample"),
    tags$div(style = "display:flex; gap:0.5rem", picker_prokka_taxa, picker_prokka_prod),
    map_prokka
  ))),
  file = file.path(output_folder, "prokka_dynamic_map.html"),
  background = "white"
)

#
# # CREATE DASHBOARD INDEX PAGE (table layout with header, iframe, and legend rows)
# index_html <- file.path(output_folder, "index.html")
# htmltools::save_html(
#   htmltools::tagList(
#     tags$html(
#       tags$head(
#         tags$title("eDNA Mapping Dashboard"),
#         tags$style("body { font-family: sans-serif; }
#                     table { width: 100%; border-collapse: collapse; }
#                     td { padding: 10px; vertical-align: top; width: 25%; text-align: center; }
#                     iframe { width: 100%; height: 400px; border: 1px solid #ccc; }
#                     img { max-width: 100%; height: auto; border: 1px solid #ccc; }")
#       ),
#       tags$body(
#         tags$h2("eDNA Mapping Dashboard"),
#         tags$h3(paste0("last updated: ", date())),
#         tags$table(
#           tags$tr(
#             tags$td(
#               tags$h3("Overall Diversity"),
#               tags$p("Most abundant phyla")
#             ),
#             tags$td(
#               tags$h3("Overall Diversity"),
# #              tags$p("Choose a genus to display"),
#               genus_picker
#             ),
#             tags$td(
#               tags$h3("Of Potential Concern"),
#               tags$p("Grouped by public health relevance")
#             ),
#             tags$td(
#               tags$h3("Cyanobacteria Diversity"),
#               tags$p("Most abundant families")
#             )
#           ),
#
#           tags$tr(
#             tags$td(tags$iframe(src = "phyla_map.html")),
#             tags$td(map_genus),
#             tags$td(tags$iframe(src = "harmful_map.html")),
#             tags$td(tags$iframe(src = "cyanobacteria_map.html"))
#           ),
#           tags$tr(
#             tags$td(tags$iframe(src = "phyla_pie_map.html")),
#             tags$td(tags$iframe(src = "genus_pie_map.html")),
#             tags$td(tags$iframe(src = "harmful_pie_map.html")),
#             tags$td(tags$iframe(src = "cyanobacteria_pie_map.html"))
#           ),
#           tags$tr(
#             tags$td(tags$img(src = "legend_phyla.png")),
#             tags$td(tags$img(src = "legend_genus.png")),
#             tags$td(tags$img(src = "legend_harmful.png")),
#             tags$td(tags$img(src = "legend_cyano.png"))
#           )
#         )
#       )
#     )
#   ),
#   file = index_html
# )
  
# UPDATED CREATE DASHBOARD INDEX PAGE
# Replace your existing index.html creation with this version
source("/work/apps/dana/kraken-table.r")

# First generate both original and quality-filtered matrices
generate_all_abundance_matrices()           # Original matrices

# failing
# Generating quality-filtered abundance matrices...
# Error in `filter()`:
#   ℹ In argument: `!str_detect(name, "Homo|Human|unclassified", ignore.case = TRUE)`.
# Caused by error in `str_detect()`:
#   ! unused argument (ignore.case = TRUE)
# Run `rlang::last_trace()` to see where the error occurred.
#generate_quality_filtered_matrices()       # Quality-controlled matrices

# # Summary of what gets generated:
# message("Matrix generation complete!")
# message("Original matrices: basic filtering only")
# message("QC matrices: enhanced quality filtering to reduce false positives")
# message("- Standard QC: removes contaminants, applies read thresholds")
# message("- Strict QC: adds taxonomic depth filter, statistical outlier removal")
# message("Filter logs saved with each QC matrix for transparency")


# Create the download section
download_section <- create_download_section()

# CREATE DASHBOARD INDEX PAGE (with download section)
index_html <- file.path(output_folder, "index.html")
htmltools::save_html(htmltools::tagList(tags$html(
  tags$head(
    tags$title("eDNA Mapping Dashboard"),
    tags$style(
      "body { font-family: sans-serif; margin: 20px; }
                    table { width: 100%; border-collapse: collapse; }
                    td { padding: 10px; vertical-align: top; width: 25%; text-align: center; }
                    iframe { width: 100%; height: 800px; border: 1px solid #ccc; }
                    img { max-width: 100%; height: auto; border: 1px solid #ccc; }
                    .download-section { max-width: 1200px; margin: 0 auto; }
                    a { color: #2c3e50; }
                    a:hover { background-color: #bdc3c7 !important; }"
    )
  ),
  
  tags$head(
    tags$style(HTML("
    .selectize-input {
      min-width: 480px !important;
      width: 480px !important;
    }
    .selectize-dropdown {
      min-width: 480px !important;
      width: 480px !important;
    }
  "))
  ),
  tags$body(
    tags$h2("eDNA Mapping Dashboard"),
    tags$h3(paste0("Last updated: ", date())),
    
    # Main visualization table
    tags$table(
      tags$tr(
        tags$td(tags$h3("Overall Diversity"), tags$p("Most abundant phyla")),
        tags$td(tags$h3("Overall Diversity"), genus_picker),
        # tags$td(
        #   tags$h3("Of Potential Concern"),
        #   tags$p("Grouped by public health relevance")
        # ),
        tags$td(
          tags$h3("Cyanobacteria Diversity"),
          tags$p("Most abundant families")
        )
      ),
      
      tags$tr(
        tags$td(tags$iframe(src = "phyla_map.html")),
        tags$td(tags$div(map_genus)),
#        tags$td(tags$iframe(src = "harmful_map.html")),
        tags$td(tags$iframe(src = "cyanobacteria_map.html"))
      ),
      tags$tr(
        tags$td(tags$iframe(src = "phyla_pie_map.html")),
        tags$td(tags$iframe(src = "genus_pie_map.html")),
#        tags$td(tags$iframe(src = "harmful_pie_map.html")),
        tags$td(tags$iframe(src = "cyanobacteria_pie_map.html"))
      ),
      tags$tr(
        tags$td(tags$img(src = "legend_phyla.png")),
        tags$td(tags$img(src = "legend_genus.png")),
 #       tags$td(tags$img(src = "legend_harmful.png")),
        tags$td(tags$img(src = "legend_cyano.png"))
      )
    ),
    
    # Add download section
    tags$div(class = "download-section", download_section)
  )
)), file = index_html)

htmlfile = file.path(getwd(), index_html)
message("Dashboard index written to: ", htmlfile)

# system(
#   paste0(
#     "rsync -avz ",
#     output_folder,
#     "/* dreamhost:/home/reric/reric.org/CBM/",
#     project
#   ),
#   wait = FALSE
# )

dir.create("~/Desktop/Amundsen-Collins/edna_dashboard/", recursive = TRUE,showWarnings = FALSE)
system(
  paste0(
    "rsync -av ",
    output_folder,
    "/* ~/Desktop/Amundsen-Collins/edna_dashboard/"
  ),
  wait = FALSE
)


#gene_data = prokka_data %>% filter(str_detect(product,"peptide synthetase",negate=FALSE))



#Sys.sleep(120)
#source("~/Desktop/edna_dev/edna_mapping.r")
