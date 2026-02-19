# Load required packages
library(httr)
library(jsonlite)
library(readr)
library(dplyr)

# Your OpenAI API key (set in environment variable)
api_key <- Sys.getenv("OPENAI_API_KEY")
if (api_key == "") {
  stop("OPENAI_API_KEY environment variable not set")
}

setwd("/data/scratch/refdbs/FOAM")
# Set input and output file paths
input_file <- "FOAM-L1-L3.tsv"
output_file <- "FOAM-ESHMM-categories.tsv"

# Read your FOAM entries
foam_data <- read_tsv(input_file)

# Prepare output file with headers
output <- tibble(
  Category = character(),
  Sub_category = character(),
  Sub_subcategory = character(),
  Rationale = character()
)

# Function to call OpenAI API
classify_function <- function(context) {
  url <- "https://api.openai.com/v1/chat/completions"
  
  body <- list(
    model = "gpt-4.1",
    messages = list(
      list(role = "system", content = "You are a microbiologist classifying microbial functions into ecosystem service categories as shown below.
           🌊 Ecosystem Service Categories
① Provisioning Services

Microbial contributions to direct products obtained from aquatic ecosystems.

    Microbial biomass production

        Protein synthesis and amino acid biosynthesis pathways

    Secondary metabolite production

        Antibiotic biosynthesis (e.g., polyketide synthases [PKS], nonribosomal peptide synthetases [NRPS])

        Vitamins and cofactor synthesis (B12, biotin, riboflavin, etc.)

    Energy sources (Biofuels)

        Fatty acid biosynthesis (lipid metabolism pathways)

        Hydrogen production (hydrogenases)

② Cultural Services

Microbial functions contributing indirectly to cultural or educational value.

    Indicator species/metabolic markers

        Phylogenetic marker genes (16S rRNA, 18S rRNA, ITS regions)

        Functional gene markers linked to water quality (e.g., toxin synthesis, pathogens)

    Bioindicators of ecosystem health

        Toxin-producing cyanobacteria (microcystins, saxitoxins synthesis pathways)

        Pathogenic bacteria markers (e.g., virulence factors)

    Educational/scientific research utility

        Genes associated with extremophile metabolism (e.g., heat-shock proteins, osmolyte synthesis pathways)

③ Habitat Supporting Services

Microbial activities critical to maintaining habitat structure and function.

    Biofilm formation and stabilization

        Extracellular polymeric substances (EPS) biosynthesis pathways

        Quorum sensing genes (luxI/luxR, autoinducer synthases)

    Primary productivity

        Photosynthesis (psbA, rbcL, chlorophyll biosynthesis pathways)

        Chemosynthesis (e.g., sulfur oxidation, sulfide oxidation, ammonia oxidation)

    Microbial symbiosis and mutualism

        Symbiotic nitrogen fixation genes (nifH, nifD, nifK)

        Ectosymbiont/endosymbiont-specific markers (symbiosis-related genes)

④ Regulatory Services

Microbial activities regulating environmental quality and ecosystem stability.

    Nutrient Cycling

        Nitrogen cycle

            Nitrogen fixation (nifH, nifD, nifK)

            Nitrification (amoA, hao, nxrA/B)

            Denitrification (nirS, nirK, norB, nosZ)

            Ammonification (urease, ammonium transporter genes)

            Anammox (hydrazine synthase, hydrazine oxidoreductase)

        Phosphorus cycling

            Phosphate solubilization (phoD, phoX)

            Phosphorus uptake and metabolism (pst genes)

        Sulfur cycling

            Sulfate reduction (dsrAB)

            Sulfur oxidation (sox genes, sulfide:quinone oxidoreductase [sqr])

    Carbon Cycling and Sequestration

        Carbon fixation pathways

            Calvin-Benson-Bassham (RuBisCO genes rbcL, rbcS)

            Reverse tricarboxylic acid cycle (rTCA)

            Wood-Ljungdahl pathway (carbon monoxide dehydrogenase, CODH)

        Methane metabolism

            Methanogenesis (mcrA, mcrB)

            Methane oxidation (pmoA, mmoX)

        Organic matter degradation

            Complex polysaccharide breakdown (CAZymes, cellulases, hemicellulases)

            Aromatic compound degradation (ring cleavage dioxygenases)

    Pollutant Detoxification and Bioremediation

        Hydrocarbon degradation (alkane hydroxylases, aromatic oxygenases)

        Xenobiotic metabolism (cytochrome P450 monooxygenases, glutathione transferases)

        Metal/metalloid detoxification (mercury reductase, arsenate reductase)
           
           Follow a strict output format: Category<TAB>Sub-category<TAB>Sub-subcategory<TAB>Rationale"),
      list(role = "user", content = paste0("Context: ", context))
    ),
    temperature = 0
  )
  
  body$messages[[1]]$content = paste(rep(x = body$messages[[1]]$content,10), collapse="\n")
  
  response <- POST(
    url,
    add_headers(`Authorization` = paste("Bearer", api_key),
                `Content-Type` = "application/json"),
    body = toJSON(body, auto_unbox = TRUE)
  )
  
  if (status_code(response) == 200) {
    content_text <- content(response, as = "text", encoding = "UTF-8")
    message_content <- fromJSON(content_text)$choices$message$content
    return(strsplit(message_content, "\t")[[1]])
  } else {
    print(content(response, as = "text"))
    stop("Request failed.")
  }
}

# Loop through each row, classify, and store the results
for (i in 1:nrow(foam_data)) {
  context <- paste(foam_data$L1[i], foam_data$L2[i], foam_data$L3[i], sep = " > ")
  cat(paste0("Processing row ", i, ": ", context, "\n"))
  
  result <- tryCatch(
    {
      classify_function(context)
    },
    error = function(e) {
      cat("Error at row ", i, ": ", e$message, "\n")
      return(rep(NA, 4))
    }
  )
  
  output <- output %>% add_row(
    Category = result[1],
    Sub_category = result[2],
    Sub_subcategory = result[3],
    Rationale = result[4]
  )
  
  Sys.sleep(20)  # Be respectful of API rate limits
}

# Write output to TSV
write_tsv(output, output_file)

cat("Classification complete. Output written to:", output_file, "\n")

