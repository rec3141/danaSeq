#!/usr/bin/env Rscript
#
# app.R — Microscape Shiny Visualization App
#
# High-performance interactive explorer for amplicon sequencing data.
# Designed for large datasets (4K+ samples, 100K+ ASVs) using:
#   - plotly with WebGL (scattergl) for client-side rendering
#   - data.table for all filtering operations
#   - Pre-computed data from build_shiny.R (no heavy computation at runtime)
#   - Debounced inputs to avoid redundant re-renders
#
# Usage:
#   Place app.R alongside app_data.rds, then:
#     shiny::runApp(".")
#
# Dependencies: shiny, plotly, DT, data.table, htmlwidgets

library(shiny)
library(plotly)
library(DT)
library(data.table)

# ---------------------------------------------------------------------------
# Load pre-computed data
# ---------------------------------------------------------------------------
app_dir <- getwd()
data_path <- file.path(app_dir, "app_data.rds")
if (!file.exists(data_path)) {
    stop("app_data.rds not found in ", app_dir,
         "\nRun build_shiny.R first to generate it.")
}

D <- readRDS(data_path)

# Convert to data.tables for fast filtering
dt_counts   <- as.data.table(D$dt_counts)
dt_sample   <- as.data.table(D$sample_tsne)
dt_seq      <- as.data.table(D$seq_tsne)
dt_asv_info <- as.data.table(D$asv_info)
dt_sample_info <- as.data.table(D$sample_info)
network_edges  <- as.data.table(D$network_edges)

group_colors <- D$group_colors
tax_strings  <- D$tax_strings
esv_ids      <- D$esv_ids
db_names     <- D$db_names
primary_db   <- D$primary_db
taxonomy     <- D$taxonomy

# Pre-compute sample-level stats for filters
sample_reads <- dt_counts[, .(total_reads = sum(count),
                               n_asvs = uniqueN(sequence)), by = sample]
setkey(sample_reads, sample)

# Pre-compute ASV prevalence for filters
asv_prevalence <- dt_counts[, .(n_samples = uniqueN(sample),
                                 total_reads = sum(count)), by = sequence]
setkey(asv_prevalence, sequence)

# Merge sample stats into sample t-SNE data
dt_sample <- merge(dt_sample, sample_reads, by = "sample", all.x = TRUE)
dt_sample[is.na(total_reads), total_reads := 0]
dt_sample[is.na(n_asvs), n_asvs := 0]

# Build dt_seq with all needed columns by merging t-SNE coords with ASV info
# Use a fresh merge to avoid column name collisions
dt_seq_base <- dt_seq[, .(x, y, sequence)]
dt_seq <- merge(dt_seq_base, dt_asv_info, by = "sequence", all.x = TRUE)
dt_seq[is.na(total_reads), total_reads := 0]
dt_seq[is.na(n_samples), n_samples := 0]
dt_seq[is.na(group), group := "unknown"]
dt_seq[is.na(group_color), group_color := "#999999"]
dt_seq[is.na(taxonomy_string), taxonomy_string := "unclassified"]
if (!"esv_id" %in% names(dt_seq)) {
    dt_seq[, esv_id := paste0("ASV_", .I)]
}

# Identify metadata columns suitable for coloring the sample plot
meta_cols <- names(dt_sample_info)
meta_cols <- meta_cols[!meta_cols %in% c("sample", "total_reads", "n_asvs")]
# Keep columns with reasonable cardinality for coloring (< 200 unique values)
meta_color_cols <- c("total_reads", "n_asvs")
for (mc in meta_cols) {
    vals <- dt_sample_info[[mc]]
    if (!is.null(vals) && length(unique(vals[!is.na(vals)])) < 200) {
        meta_color_cols <- c(meta_color_cols, mc)
    }
}

# Merge metadata into sample t-SNE for hover and coloring
if (nrow(dt_sample_info) > 0) {
    # Only merge columns not already in dt_sample
    extra_cols <- setdiff(names(dt_sample_info), names(dt_sample))
    if (length(extra_cols) > 0) {
        merge_cols <- c("sample", extra_cols)
        dt_sample <- merge(dt_sample, dt_sample_info[, ..merge_cols],
                           by = "sample", all.x = TRUE)
    }
}

# Identify available taxonomic groups for checkbox filter
available_groups <- sort(unique(dt_seq$group))

# Slider ranges
max_log_reads <- max(1, ceiling(log2(max(dt_sample$total_reads, na.rm = TRUE))))
max_log_richness <- max(1, ceiling(log2(max(dt_sample$n_asvs, na.rm = TRUE))))
max_prevalence <- max(1, ceiling(log2(max(dt_seq$n_samples, na.rm = TRUE))))

# ===========================================================================
# UI
# ===========================================================================
ui <- fluidPage(
    tags$head(
        tags$style(HTML("
            .sidebar-panel { padding: 10px; }
            .nav-tabs > li > a { padding: 8px 12px; }
            body { font-size: 13px; }
            .dataTables_wrapper { font-size: 12px; }
        "))
    ),

    titlePanel("Microscape Explorer"),

    sidebarLayout(
        sidebarPanel(
            width = 3,

            # --- Shared controls ---
            selectInput("tax_db", "Taxonomy Database",
                        choices = db_names, selected = primary_db),

            checkboxGroupInput("groups", "Taxonomic Groups",
                               choices = available_groups,
                               selected = available_groups,
                               inline = TRUE),

            hr(),

            # --- Selected items display ---
            h5("Selected Sample"),
            verbatimTextOutput("selected_sample_display", placeholder = TRUE),

            h5("Selected ASV"),
            verbatimTextOutput("selected_asv_display", placeholder = TRUE),

            hr(),

            # --- Tab-specific controls rendered server-side ---
            uiOutput("tab_controls")
        ),

        mainPanel(
            width = 9,
            tabsetPanel(
                id = "main_tabs",

                # Tab 1: Sample Explorer
                tabPanel(
                    "Sample Explorer",
                    value = "samples",
                    plotlyOutput("sample_plot", height = "650px"),
                    hr(),
                    h5("Top taxa in selected sample"),
                    DT::dataTableOutput("sample_taxa_table")
                ),

                # Tab 2: ASV Network
                tabPanel(
                    "ASV Network",
                    value = "network",
                    plotlyOutput("network_plot", height = "650px")
                ),

                # Tab 3: Summary Tables
                tabPanel(
                    "Summary Tables",
                    value = "tables",
                    h4("ASV Summary"),
                    downloadButton("download_asv", "Download ASV table (CSV)"),
                    DT::dataTableOutput("asv_table"),
                    hr(),
                    h4("Sample Summary"),
                    downloadButton("download_sample", "Download Sample table (CSV)"),
                    DT::dataTableOutput("sample_table")
                )
            )
        )
    )
)

# ===========================================================================
# Server
# ===========================================================================
server <- function(input, output, session) {

    # -----------------------------------------------------------------------
    # Reactive state
    # -----------------------------------------------------------------------
    selected_sample <- reactiveVal(NULL)
    selected_asv    <- reactiveVal(NULL)

    # Debounced filter inputs (avoid re-render while slider is being dragged)
    min_reads_d     <- debounce(reactive(input$min_reads), 300)
    min_richness_d  <- debounce(reactive(input$min_richness), 300)
    min_prevalence_d <- debounce(reactive(input$min_prevalence), 300)
    min_corr_d      <- debounce(reactive(input$min_corr), 300)

    # -----------------------------------------------------------------------
    # Tab-specific controls
    # -----------------------------------------------------------------------
    output$tab_controls <- renderUI({
        tab <- input$main_tabs

        if (tab == "samples") {
            tagList(
                h5("Sample Filters"),
                sliderInput("min_reads", "Min reads (log2)",
                            min = 0, max = max_log_reads,
                            value = 0, step = 0.5),
                sliderInput("min_richness", "Min richness (log2)",
                            min = 0, max = max_log_richness,
                            value = 0, step = 0.5),
                textInput("sample_regex", "Sample name filter (regex)", ""),
                hr(),
                h5("Taxa Overlay"),
                textInput("taxa_overlay_regex", "Taxonomy filter (regex)", ""),
                sliderInput("taxa_overlay_min_hits", "Min samples per ASV (log2)",
                            min = 0, max = max_prevalence,
                            value = 0, step = 0.5)
            )
        } else if (tab == "network") {
            tagList(
                h5("ASV Filters"),
                sliderInput("min_prevalence", "Min samples per ASV (log2)",
                            min = 0, max = max_prevalence,
                            value = 1, step = 0.5),
                sliderInput("min_corr", "Min SparCC |correlation|",
                            min = 0, max = 1, value = 0.3, step = 0.01),
                textInput("taxa_regex", "Taxonomy filter (regex)", "")
            )
        } else {
            tagList()
        }
    })

    # -----------------------------------------------------------------------
    # Display for selected items
    # -----------------------------------------------------------------------
    output$selected_sample_display <- renderText({
        s <- selected_sample()
        if (is.null(s)) "(none)" else s
    })

    output$selected_asv_display <- renderText({
        a <- selected_asv()
        if (is.null(a)) {
            "(none)"
        } else {
            ts <- tax_strings[a]
            if (is.na(ts)) ts <- "unclassified"
            eid <- esv_ids[a]
            if (is.na(eid)) eid <- ""
            paste0(eid, "\n", ts)
        }
    })

    # -----------------------------------------------------------------------
    # Tab 1: Sample Explorer
    # -----------------------------------------------------------------------

    # Filtered sample data
    filtered_samples <- reactive({
        mr <- min_reads_d()
        mrich <- min_richness_d()
        regex <- input$sample_regex

        dt <- copy(dt_sample)

        if (!is.null(mr) && mr > 0) {
            dt <- dt[total_reads >= 2^mr]
        }
        if (!is.null(mrich) && mrich > 0) {
            dt <- dt[n_asvs >= 2^mrich]
        }
        if (!is.null(regex) && nzchar(regex)) {
            tryCatch({
                dt <- dt[grepl(regex, sample, ignore.case = TRUE)]
            }, error = function(e) {
                # Invalid regex — skip filter
            })
        }

        dt
    })

    output$sample_plot <- renderPlotly({
        dt_samp <- filtered_samples()
        if (nrow(dt_samp) == 0) {
            return(plotly_empty(type = "scatter", mode = "markers") %>%
                       layout(title = "No samples match filters"))
        }

        # Get selected taxonomic groups
        sel_groups <- input$groups

        # ---- Build taxa overlay data ----
        # For each (sample, ASV) pair, place a point at the sample's t-SNE
        # position, sized by that ASV's proportion in that sample, colored
        # by taxonomic group. This is the "microscape" visualization.

        # Filter counts to selected samples and groups
        visible_samples <- dt_samp$sample
        overlay <- dt_counts[sample %in% visible_samples]
        if (!is.null(sel_groups) && length(sel_groups) > 0) {
            overlay <- overlay[group %in% sel_groups]
        }

        # Filter by taxonomy regex if provided
        taxa_regex <- input$taxa_overlay_regex
        if (!is.null(taxa_regex) && nzchar(taxa_regex)) {
            tryCatch({
                matching_seqs <- names(tax_strings)[grepl(taxa_regex, tax_strings,
                                                          ignore.case = TRUE)]
                overlay <- overlay[sequence %in% matching_seqs]
            }, error = function(e) NULL)
        }

        # Filter by min prevalence
        taxa_min <- input$taxa_overlay_min_hits
        if (!is.null(taxa_min) && taxa_min > 0) {
            prev_seqs <- asv_prevalence[n_samples >= 2^taxa_min, sequence]
            overlay <- overlay[sequence %in% prev_seqs]
        }

        # Only keep ASVs with nonzero counts
        overlay <- overlay[count > 0]

        if (nrow(overlay) > 0) {
            # Join sample t-SNE coordinates
            overlay <- merge(overlay, dt_samp[, .(sample, x, y)],
                             by = "sample", all.x = TRUE)

            # Compute point sizes from proportion (4th-root scaling like OG app)
            overlay[, size := pmax(1, 15 * sqrt(sqrt(proportion)))]

            # Assign colors by group
            group_color_map <- c(
                prokaryote   = "#427df4",
                eukaryote    = "#f44265",
                chloroplast  = "#0fe047",
                mitochondria = "#42e5f4",
                unknown      = "#999999"
            )
            overlay[, color := group_color_map[group]]
            overlay[is.na(color), color := "#999999"]

            # Build hover text
            overlay[, hover := {
                tax <- tax_strings[sequence]
                tax[is.na(tax)] <- "unclassified"
                paste0(tax, "\n", round(proportion * 1000, 1), " permille")
            }]
        }

        # ---- Base plot: grey sample outlines ----
        dt_samp[, base_size := pmax(3, log10(total_reads + 1) * 3)]
        dt_samp[, hover_base := paste0(
            "Sample: ", sample,
            "\nReads: ", format(total_reads, big.mark = ","),
            "\nASVs: ", n_asvs
        )]

        p <- plot_ly(source = "sample_plot") %>%
            layout(
                dragmode = "pan",
                xaxis = list(title = "", zeroline = FALSE, showgrid = FALSE),
                yaxis = list(title = "", zeroline = FALSE, showgrid = FALSE,
                             scaleanchor = "x"),
                title = paste0("Microscape (", nrow(dt_samp), " samples)")
            )

        # Grey base layer — sample positions
        p <- p %>% add_trace(
            data = dt_samp, x = ~x, y = ~y,
            type = "scattergl", mode = "markers",
            marker = list(size = ~base_size, color = "lightgrey",
                          opacity = 0.5,
                          line = list(width = 0.5, color = "grey")),
            text = ~hover_base, hoverinfo = "text",
            key = ~sample, name = "samples",
            showlegend = FALSE
        )

        # ---- Taxa overlay: one trace per group ----
        if (nrow(overlay) > 0) {
            for (grp in unique(overlay$group)) {
                grp_data <- overlay[group == grp]
                grp_color <- group_color_map[grp]
                if (is.na(grp_color)) grp_color <- "#999999"

                p <- p %>% add_trace(
                    data = grp_data, x = ~x, y = ~y,
                    type = "scattergl", mode = "markers",
                    marker = list(size = ~size, color = grp_color,
                                  opacity = 0.6,
                                  line = list(width = 0)),
                    text = ~hover, hoverinfo = "text",
                    key = ~sample, name = grp
                )
            }
        }

        p %>% config(scrollZoom = TRUE)
    })

    # Capture click events from sample plot
    observe({
        click_data <- event_data("plotly_click", source = "sample_plot")
        if (!is.null(click_data) && !is.null(click_data$key)) {
            selected_sample(click_data$key[[1]])
        }
    })

    # Top taxa table for selected sample
    output$sample_taxa_table <- DT::renderDataTable({
        s <- selected_sample()
        if (is.null(s)) {
            return(DT::datatable(
                data.frame(Message = "Click a sample to view its top taxa"),
                options = list(dom = "t"), rownames = FALSE
            ))
        }

        sample_data <- dt_counts[sample == s]
        if (nrow(sample_data) == 0) {
            return(DT::datatable(
                data.frame(Message = "No data for selected sample"),
                options = list(dom = "t"), rownames = FALSE
            ))
        }

        sample_data <- sample_data[order(-count)]
        sample_data <- sample_data[1:min(50, nrow(sample_data))]
        sample_data[, taxonomy := tax_strings[sequence]]
        sample_data[, esv := esv_ids[sequence]]

        out <- sample_data[, .(ESV = esv, Taxonomy = taxonomy,
                               Reads = count, Group = group,
                               Proportion = round(proportion, 4))]

        DT::datatable(out, rownames = FALSE,
                      options = list(pageLength = 15, dom = "ftp",
                                     scrollX = TRUE),
                      selection = "single")
    })

    # -----------------------------------------------------------------------
    # Tab 2: ASV Network
    # -----------------------------------------------------------------------

    # Filtered ASV data
    filtered_asvs <- reactive({
        mp <- min_prevalence_d()
        regex <- input$taxa_regex
        groups <- input$groups

        dt <- copy(dt_seq)

        if (!is.null(mp) && mp > 0) {
            dt <- dt[n_samples >= 2^mp]
        }
        if (!is.null(regex) && nzchar(regex)) {
            tryCatch({
                dt <- dt[grepl(regex, taxonomy_string, ignore.case = TRUE)]
            }, error = function(e) {
                # Invalid regex — skip filter
            })
        }
        if (!is.null(groups)) {
            dt <- dt[group %in% groups]
        }

        dt
    })

    # Filtered edges
    filtered_edges <- reactive({
        mc <- min_corr_d()
        if (is.null(mc)) mc <- 0.3

        dt_asvs <- filtered_asvs()
        if (nrow(dt_asvs) == 0 || nrow(network_edges) == 0) {
            return(data.table())
        }

        edges <- network_edges[abs(correlation) >= mc]
        if (nrow(edges) == 0) return(data.table())

        # Filter to only include edges where both nodes are visible
        visible_seqs <- dt_asvs$sequence
        visible_esvs <- dt_asvs$esv_id

        visible_ids <- c(visible_seqs, visible_esvs)
        visible_ids <- visible_ids[!is.na(visible_ids)]

        edges <- edges[node1 %in% visible_ids & node2 %in% visible_ids]
        edges
    })

    output$network_plot <- renderPlotly({
        # Return a simple placeholder plot — full network plot needs debugging
        dt <- copy(filtered_asvs())
        if (TRUE) {  # TODO: fix network plot data.table column issues
        return(plot_ly(dt, x = ~x, y = ~y,
                       type = "scattergl", mode = "markers",
                       marker = list(size = 5, color = ~group_color, opacity = 0.7),
                       text = ~paste0(sequence),
                       hoverinfo = "text") %>%
            layout(title = paste0("ASV t-SNE (", nrow(dt), " ASVs)"),
                   xaxis = list(title = "", zeroline = FALSE),
                   yaxis = list(title = "", zeroline = FALSE, scaleanchor = "x")) %>%
            config(scrollZoom = TRUE))
        }
        tryCatch({
        if (nrow(dt) == 0) {
            return(plotly_empty(type = "scatter", mode = "markers") %>%
                       layout(title = "No ASVs match filters"))
        }

        # Base point size proportional to log(total abundance)
        dt[, base_size := pmax(2, log10(total_reads + 1) * 2)]

        # If a sample is selected, adjust sizes to show that sample's composition
        sel_sample <- selected_sample()
        if (!is.null(sel_sample)) {
            # Look up per-sequence counts for this sample (no merge needed)
            sample_counts <- dt_counts[sample == sel_sample, .(sequence, count)]
            setkey(sample_counts, sequence)
            dt[, seq_count := sample_counts[sequence, sum(count, na.rm = TRUE),
                                             by = .EACHI]$V1]
            dt[is.na(seq_count), seq_count := 0L]
            dt[, size := ifelse(seq_count > 0,
                                pmax(4, log10(seq_count + 1) * 5),
                                base_size * 0.5)]
            dt[, opacity := ifelse(seq_count > 0, 0.9, 0.3)]
        } else {
            dt[, size := base_size]
            dt[, opacity := 0.7]
            dt[, seq_count := 0L]
        }

        # Hover text
        dt[, hover := paste0(
            esv_id, "\n", taxonomy_string,
            "\nTotal reads: ", format(total_reads, big.mark = ","),
            "\nSamples: ", n_samples,
            "\nGroup: ", group
        )]
        if (!is.null(sel_sample)) {
            dt[seq_count > 0,
               hover := paste0(hover, "\nIn ", sel_sample, ": ",
                                format(seq_count, big.mark = ","), " reads")]
        }

        # Start plot with ASV points
        p <- plot_ly() %>%
            layout(
                dragmode = "pan",
                xaxis = list(title = "", zeroline = FALSE, showgrid = FALSE),
                yaxis = list(title = "", zeroline = FALSE, showgrid = FALSE,
                             scaleanchor = "x"),
                title = paste0("ASV Network (", nrow(dt), " ASVs)")
            )

        # Draw edges first (behind points) using add_segments
        edges <- filtered_edges()
        if (nrow(edges) > 0) {
            # Limit edges for performance — draw at most 50K edges
            if (nrow(edges) > 50000) {
                edges <- edges[order(-abs(correlation))][1:50000]
            }
            p <- p %>% add_segments(
                data = edges,
                x = ~x1, xend = ~x2,
                y = ~y1, yend = ~y2,
                line = list(color = ~color, width = 0.5),
                hoverinfo = "none",
                showlegend = FALSE,
                inherit = FALSE
            )
        }

        # Draw ASV points, one trace per group for legend
        for (grp in unique(dt$group)) {
            dt_grp <- dt[group == grp]
            if (nrow(dt_grp) == 0) next

            grp_color <- group_colors[grp]
            if (is.na(grp_color)) grp_color <- "#999999"

            p <- p %>% add_trace(
                data = dt_grp,
                x = ~x, y = ~y,
                type = "scattergl",
                mode = "markers",
                marker = list(
                    size = ~size,
                    color = grp_color,
                    opacity = ~opacity,
                    line = list(width = 0.2, color = "white")
                ),
                text = ~hover,
                hoverinfo = "text",
                key = ~sequence,
                name = grp,
                legendgroup = grp,
                source = "network_plot",
                inherit = FALSE
            )
        }

        p %>% config(scrollZoom = TRUE)
        }, error = function(e) {
            message("[WARNING] Network plot error: ", e$message)
            plotly_empty(type = "scatter", mode = "markers") %>%
                layout(title = paste("Network plot error:", e$message))
        })
    })

    # Capture click events from network plot
    observe({
        click_data <- event_data("plotly_click", source = "network_plot")
        if (!is.null(click_data) && !is.null(click_data$key)) {
            selected_asv(click_data$key[[1]])
        }
    })

    # -----------------------------------------------------------------------
    # Tab 3: Summary Tables
    # -----------------------------------------------------------------------

    # ASV table
    asv_table_data <- reactive({
        dt <- dt_asv_info[, .(esv_id, taxonomy_string, total_reads,
                               n_samples, group)]
        setnames(dt, c("ESV", "Taxonomy", "Total Reads", "Samples", "Group"))
        dt[order(-`Total Reads`)]
    })

    output$asv_table <- DT::renderDataTable({
        DT::datatable(
            asv_table_data(),
            rownames = FALSE,
            filter = "top",
            options = list(
                pageLength = 25,
                scrollX = TRUE,
                dom = "lftipr",
                search = list(regex = TRUE, caseInsensitive = TRUE)
            ),
            selection = "single"
        )
    })

    # Sample table
    sample_table_data <- reactive({
        dt <- copy(dt_sample_info)
        dt[order(-total_reads)]
    })

    output$sample_table <- DT::renderDataTable({
        DT::datatable(
            sample_table_data(),
            rownames = FALSE,
            filter = "top",
            options = list(
                pageLength = 25,
                scrollX = TRUE,
                dom = "lftipr",
                search = list(regex = TRUE, caseInsensitive = TRUE)
            ),
            selection = "single"
        )
    })

    # Download handlers
    output$download_asv <- downloadHandler(
        filename = function() {
            paste0("microscape_asv_summary_", Sys.Date(), ".csv")
        },
        content = function(file) {
            # Get filtered rows from DT
            rows <- input$asv_table_rows_all
            dt <- asv_table_data()
            if (!is.null(rows)) {
                dt <- dt[rows, ]
            }
            fwrite(dt, file)
        }
    )

    output$download_sample <- downloadHandler(
        filename = function() {
            paste0("microscape_sample_summary_", Sys.Date(), ".csv")
        },
        content = function(file) {
            rows <- input$sample_table_rows_all
            dt <- sample_table_data()
            if (!is.null(rows)) {
                dt <- dt[rows, ]
            }
            fwrite(dt, file)
        }
    )
}

# ---------------------------------------------------------------------------
# Launch
# ---------------------------------------------------------------------------
shinyApp(ui = ui, server = server)
