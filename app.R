# app.R — Genetic Association Heatmap
#
# Interactive Shiny application for visualising beta values from genetic
# association studies.  Data sources:
#   • GWAS Catalog REST API v2  (https://www.ebi.ac.uk/gwas/rest/api/v2)
#   • HuGeAMP BioIndex REST API (https://bioindex.hugeamp.org/api/bio)
#
# Usage (from RStudio):  shiny::runApp()
# Usage (from R console): shiny::runApp(".")

# ---- Dependencies ------------------------------------------------------------
library(shiny)
library(httr)
library(jsonlite)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(circlize)
library(DT)

# Load helper modules
source("R/api_utils.R")
source("R/heatmap_utils.R")

# ---- UI ----------------------------------------------------------------------

ui <- fluidPage(
  title = "Genetic Association Heatmap",

  tags$head(
    tags$style(HTML("
      body { font-family: 'Helvetica Neue', Arial, sans-serif; }
      .app-subtitle { color: #555; font-size: 13px; margin-top: -6px; margin-bottom: 14px; }
      .section-title { font-size: 14px; font-weight: bold; margin-bottom: 6px;
                        border-bottom: 1px solid #ddd; padding-bottom: 4px; }
      .status-box { padding: 6px 10px; border-radius: 4px; margin: 6px 0;
                    font-size: 13px; }
      .status-ok   { background:#dff0d8; color:#3c763d; border:1px solid #c3e6cb; }
      .status-err  { background:#f8d7da; color:#721c24; border:1px solid #f5c6cb; }
      .status-info { background:#d1ecf1; color:#0c5460; border:1px solid #bee5eb; }
      .btn-block   { width: 100%; margin-bottom: 4px; }
      #heatmap_container { min-height: 400px; }
      .chip { display:inline-block; background:#e9ecef; border-radius:12px;
              padding:2px 8px; font-size:12px; margin:2px; }
      textarea.bulk-input { font-family: monospace; font-size: 12px;
                            resize: vertical; }
      .bulk-note { color:#666; font-size:11px; margin-top:2px; }
    "))
  ),

  # ---- App header ----
  fluidRow(
    column(12,
      h2("Genetic Association Heatmap"),
      p(class = "app-subtitle",
        "Beta values for genetic variant–study associations sourced from ",
        a("GWAS Catalog v2", href = "https://www.ebi.ac.uk/gwas/",
          target = "_blank"),
        " or ",
        a("HuGeAMP BioIndex", href = "https://hugeamp.org/",
          target = "_blank"), "."
      ),
      hr(style = "margin-top:4px;")
    )
  ),

  sidebarLayout(

    # ---- Sidebar ----
    sidebarPanel(
      width = 4,

      # Data source selector
      div(class = "section-title", "Data Source"),
      radioButtons(
        "data_source", label = NULL,
        choices = c(
          "GWAS Catalog v2"  = "gwas",
          "HuGeAMP BioIndex" = "hugeamp"
        ),
        selected = "gwas"
      ),

      hr(),

      # Search panel
      div(class = "section-title", icon("magnifying-glass"), " Search Variant"),
      textInput("search_variant_id",
                label   = NULL,
                value   = "",
                placeholder = "Variant ID, e.g. rs7903146"),
      actionButton("btn_search", "Search",
                   class = "btn-primary btn-sm",
                   icon  = icon("magnifying-glass"),
                   width = "100%"),
      br(),
      uiOutput("search_status_ui"),

      # Search results
      conditionalPanel(
        condition = "output.has_search_results",
        br(),
        div(class = "section-title", "Associations found"),
        DT::dataTableOutput("search_results_table"),
        br(),
        actionButton("btn_add_all_heatmap_studies",
                     "Add variant to all studies in heatmap",
                     class = "btn-info btn-sm btn-block",
                     icon  = icon("layer-group")),
        actionButton("btn_add_all_found",
                     "Add variant to ALL found studies",
                     class = "btn-success btn-sm btn-block",
                     icon  = icon("plus-circle")),
        actionButton("btn_add_selected",
                     "Add variant to selected rows",
                     class = "btn-secondary btn-sm btn-block",
                     icon  = icon("check")),
        br()
      ),

      hr(),

      # Manual entry
      div(class = "section-title", icon("keyboard"), " Manual Entry"),
      textInput("manual_variant", label = NULL,
                placeholder = "Variant ID, e.g. rs7903146"),
      textInput("manual_study",   label = NULL,
                placeholder = "Study ID, e.g. GCST001234"),
      actionButton("btn_add_manual", "Add to Heatmap",
                   class = "btn-warning btn-sm",
                   icon  = icon("plus"),
                   width = "100%"),

      hr(),

      # ---- Bulk Variant Input ----
      div(class = "section-title", icon("list"), " Bulk Variant Input"),
      p(class = "bulk-note",
        "Paste one rsID per line to add multiple variants at once."),
      tags$textarea(
        id          = "bulk_variants",
        class       = "form-control bulk-input",
        rows        = "6",
        placeholder = "rs7903146\nrs1801282\nrs10811661\n..."
      ),
      br(),
      checkboxInput(
        "bulk_heatmap_studies_only",
        label = "Limit to studies already in heatmap",
        value = FALSE
      ),
      actionButton("btn_add_bulk", "Add All Variants",
                   class = "btn-primary btn-sm",
                   icon  = icon("layer-group"),
                   width = "100%"),
      br(),
      uiOutput("add_status_ui"),

      hr(),

      # Current heatmap summary
      div(class = "section-title", icon("table-cells"), " Current Heatmap"),
      uiOutput("heatmap_summary_ui"),
      br(),
      actionButton("btn_clear", "Clear All",
                   class = "btn-danger btn-sm",
                   icon  = icon("trash"),
                   width = "100%")
    ),

    # ---- Main panel ----
    mainPanel(
      width = 8,
      tabsetPanel(
        id = "main_tabs",

        # Heatmap tab
        tabPanel(
          title = tagList(icon("fire"), " Heatmap"),
          value = "tab_heatmap",
          br(),
          conditionalPanel(
            condition = "output.has_heatmap_data",
            id = "heatmap_container",
            InteractiveComplexHeatmapOutput(
              "ht_main",
              height           = 520,
              width            = "100%",
              response         = c("click", "brush"),
              brush_opt        = list(stroke = "#1A73E8", opacity = 0.6),
              output_ui        = htmlOutput("ht_info")
            ),
            br(),
            downloadButton("btn_download", "Download PNG",
                           class = "btn-sm btn-secondary")
          ),
          conditionalPanel(
            condition = "!output.has_heatmap_data",
            div(
              style = "text-align:center; padding:80px 40px; color:#888;",
              tags$i(class = "fa fa-chart-area fa-4x",
                     style = "color:#ccc; display:block; margin-bottom:16px;"),
              h4("No data yet"),
              p("Search for a variant and add it to the heatmap,",
                br(),
                "or enter a variant ID and study ID manually.")
            )
          )
        ),

        # Data table tab
        tabPanel(
          title = tagList(icon("table"), " Data Table"),
          value = "tab_data",
          br(),
          DT::dataTableOutput("data_table")
        ),

        # About tab
        tabPanel(
          title = tagList(icon("circle-info"), " Help"),
          value = "tab_help",
          br(),
          uiOutput("help_ui")
        )
      )
    )
  )
)

# ---- Server ------------------------------------------------------------------

server <- function(input, output, session) {

  # ---- Reactive state -------------------------------------------------------

  # Main data store: one row per (variant, study) combination
  beta_data <- reactiveVal(
    data.frame(
      variant_id = character(0),
      study_id   = character(0),
      beta       = numeric(0),
      source     = character(0),
      stringsAsFactors = FALSE
    )
  )

  search_results  <- reactiveVal(NULL)
  search_msg      <- reactiveVal(NULL)
  add_msg         <- reactiveVal(NULL)

  # Track which source belongs to each study already in the heatmap
  study_sources <- reactiveVal(setNames(character(0), character(0)))

  # ---- Helpers ---------------------------------------------------------------

  #' Add a vector of (variant, study, beta, source) rows to beta_data,
  #' skipping any that already exist.
  add_rows <- function(rows_df) {
    existing <- beta_data()
    if (nrow(rows_df) == 0) return(invisible(NULL))

    # Deduplicate against existing entries
    existing_keys <- paste(existing$variant_id, existing$study_id, sep = "|||")
    new_keys      <- paste(rows_df$variant_id,  rows_df$study_id,  sep = "|||")
    rows_df <- rows_df[!new_keys %in% existing_keys, , drop = FALSE]

    if (nrow(rows_df) > 0) {
      beta_data(rbind(existing, rows_df))
      # Update study → source mapping
      ss <- study_sources()
      for (i in seq_len(nrow(rows_df))) {
        sid <- rows_df$study_id[i]
        if (!sid %in% names(ss)) ss[[sid]] <- rows_df$source[i]
      }
      study_sources(ss)
    }
    invisible(nrow(rows_df))
  }

  #' Ensure beta_data contains cells for every (variant × study) combination
  #' currently represented in the heatmap, fetching missing values from the
  #' appropriate source.
  fill_matrix <- function() {
    data   <- beta_data()
    if (nrow(data) == 0) return(invisible(NULL))

    variants <- unique(data$variant_id)
    studies  <- unique(data$study_id)
    ss       <- study_sources()

    new_rows <- list()
    existing_keys <- paste(data$variant_id, data$study_id, sep = "|||")

    for (v in variants) {
      for (s in studies) {
        key <- paste(v, s, sep = "|||")
        if (!key %in% existing_keys) {
          src  <- ss[[s]] %||% input$data_source
          beta <- get_beta(v, s, src)
          new_rows[[length(new_rows) + 1]] <- data.frame(
            variant_id = v,
            study_id   = s,
            beta       = beta,
            source     = src,
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (length(new_rows) > 0) {
      add_rows(do.call(rbind, new_rows))
    }
  }

  # ---- Search ---------------------------------------------------------------

  observeEvent(input$btn_search, {
    vid <- trimws(input$search_variant_id)
    req(nzchar(vid))

    search_msg(list(type = "info",
                    msg  = paste("Searching for", vid, "…")))
    search_results(NULL)

    withProgress(message = paste("Querying", input$data_source, "…"), {
      res <- search_variant(vid, input$data_source)
    })

    if (is.null(res) || nrow(res) == 0) {
      src_name <- if (input$data_source == "gwas") "GWAS Catalog" else "HuGeAMP"
      search_msg(list(type = "err",
                      msg  = paste0("No associations found for '", vid,
                                    "' in ", src_name, ".")))
    } else {
      search_msg(list(type = "ok",
                      msg  = paste(nrow(res), "associations found for", vid)))
      search_results(res)
    }
  })

  output$search_status_ui <- renderUI({
    m <- search_msg()
    if (is.null(m)) return(NULL)
    css <- paste0("status-box status-", m$type)
    div(class = css, m$msg)
  })

  output$has_search_results <- reactive({
    res <- search_results()
    !is.null(res) && nrow(res) > 0
  })
  outputOptions(output, "has_search_results", suspendWhenHidden = FALSE)

  output$search_results_table <- DT::renderDataTable({
    res <- search_results()
    req(res)
    show_cols <- intersect(c("study_id", "beta", "pvalue", "trait"),
                           names(res))
    DT::datatable(
      res[, show_cols, drop = FALSE],
      selection  = "multiple",
      rownames   = FALSE,
      options    = list(
        pageLength = 5,
        scrollX    = TRUE,
        dom        = "tip"
      )
    )
  })

  # ---- Add: variant to all studies already in the heatmap ------------------

  observeEvent(input$btn_add_all_heatmap_studies, {
    res <- search_results()
    req(res)
    vid  <- trimws(input$search_variant_id)
    data <- beta_data()

    heatmap_studies <- unique(data$study_id)
    if (length(heatmap_studies) == 0) {
      add_msg(list(type = "err",
                   msg  = "No studies in the heatmap yet. Use 'Add variant to ALL found studies' first."))
      return()
    }

    to_add <- res[res$study_id %in% heatmap_studies, , drop = FALSE]
    if (nrow(to_add) == 0) {
      # Fall back: fetch betas individually for existing studies
      new_rows <- lapply(heatmap_studies, function(s) {
        src  <- (study_sources()[[s]]) %||% input$data_source
        beta <- get_beta(vid, s, src)
        data.frame(variant_id = vid, study_id = s,
                   beta = beta, source = src, stringsAsFactors = FALSE)
      })
      to_add <- do.call(rbind, new_rows)
    }

    n <- add_rows(to_add)
    # Fill any remaining gaps
    fill_matrix()

    add_msg(list(type = "ok",
                 msg  = paste("Added", vid,
                              "for", length(heatmap_studies),
                              "heatmap study/studies.")))
  })

  # ---- Add: variant to ALL found studies ------------------------------------

  observeEvent(input$btn_add_all_found, {
    res <- search_results()
    req(res)
    vid <- trimws(input$search_variant_id)

    withProgress(message = "Adding data …", {
      n <- add_rows(res[, c("variant_id","study_id","beta","source")])
      fill_matrix()
    })

    add_msg(list(type = "ok",
                 msg  = paste("Added", vid, "for", nrow(res), "studies.")))
  })

  # ---- Add: variant to selected rows in table --------------------------------

  observeEvent(input$btn_add_selected, {
    res      <- search_results()
    selected <- input$search_results_table_rows_selected
    if (length(selected) == 0) {
      add_msg(list(type = "err", msg = "No rows selected in the table."))
      return()
    }
    to_add <- res[selected, c("variant_id","study_id","beta","source"),
                  drop = FALSE]
    n <- add_rows(to_add)
    fill_matrix()
    add_msg(list(type = "ok",
                 msg  = paste("Added", nrow(to_add), "selected association(s).")))
  })

  # ---- Add: manual entry ----------------------------------------------------

  observeEvent(input$btn_add_manual, {
    vid <- trimws(input$manual_variant)
    sid <- trimws(input$manual_study)

    if (!nzchar(vid) || !nzchar(sid)) {
      add_msg(list(type = "err",
                   msg  = "Please enter both a Variant ID and a Study ID."))
      return()
    }

    add_msg(list(type = "info",
                 msg  = paste("Fetching beta for", vid, "×", sid, "…")))

    src  <- input$data_source
    beta <- withProgress(message = paste("Querying", src, "…"), {
      get_beta(vid, sid, src)
    })

    new_row <- data.frame(variant_id = vid, study_id = sid,
                          beta = beta, source = src,
                          stringsAsFactors = FALSE)
    add_rows(new_row)

    # Fill missing cells created by this new (variant, study) combination
    fill_matrix()

    msg <- if (is.na(beta)) {
      paste0("Added ", vid, " × ", sid,
             " (beta not available from ", src, ").")
    } else {
      paste0("Added ", vid, " × ", sid, "  —  beta = ", round(beta, 5))
    }
    add_msg(list(type = "ok", msg = msg))
  })

  # ---- Add: bulk pasted variant list ----------------------------------------

  observeEvent(input$btn_add_bulk, {
    raw_text <- input$bulk_variants
    if (is.null(raw_text) || !nzchar(trimws(raw_text))) {
      add_msg(list(type = "err", msg = "Please paste at least one variant ID."))
      return()
    }

    # Parse: handle \r\n (Windows), \n (Unix), \r (old Mac); strip comments
    vids <- readLines(textConnection(raw_text), warn = FALSE)
    vids <- trimws(vids)
    vids <- vids[nzchar(vids) & !grepl("^#", vids)]
    vids <- unique(vids)

    if (length(vids) == 0) {
      add_msg(list(type = "err", msg = "No valid variant IDs found in the input."))
      return()
    }

    heatmap_only    <- isTRUE(input$bulk_heatmap_studies_only)
    heatmap_studies <- unique(beta_data()$study_id)
    src             <- input$data_source

    if (heatmap_only && length(heatmap_studies) == 0) {
      add_msg(list(
        type = "err",
        msg  = paste0("No studies in the heatmap yet. ",
                      "Uncheck 'Limit to studies already in heatmap' to add ",
                      "all found associations.")
      ))
      return()
    }

    add_msg(list(type = "info",
                 msg  = paste("Processing", length(vids), "variant(s)…")))

    n_added    <- 0L
    not_found  <- character(0)

    withProgress(message = "Adding variants…", value = 0, {
      for (i in seq_along(vids)) {
        vid <- vids[[i]]
        setProgress(
          value  = i / length(vids),
          detail = sprintf("(%d/%d) %s", i, length(vids), vid)
        )

        res <- search_variant(vid, src)

        if (is.null(res) || nrow(res) == 0) {
          not_found <- c(not_found, vid)
          next
        }

        if (heatmap_only) {
          res <- res[res$study_id %in% heatmap_studies, , drop = FALSE]
          if (nrow(res) == 0) {
            not_found <- c(not_found, vid)
            next
          }
        }

        n_new   <- add_rows(res[, c("variant_id", "study_id", "beta", "source"),
                                drop = FALSE])
        n_added <- n_added + (if (is.null(n_new)) 0L else n_new)
      }

      # Fill any gaps in the complete (variant × study) matrix
      fill_matrix()
    })

    # Build human-readable summary
    n_found <- length(vids) - length(not_found)
    parts   <- paste(n_found, "of", length(vids), "variant(s) found")
    if (n_added > 0)
      parts <- paste(parts, "—", n_added, "association(s) added")
    if (length(not_found) > 0)
      parts <- paste(parts,
                     "— not found:", paste(not_found, collapse = ", "))

    add_msg(list(
      type = if (n_added > 0) "ok" else "err",
      msg  = parts
    ))
  })

  output$add_status_ui <- renderUI({
    m <- add_msg()
    if (is.null(m)) return(NULL)
    css <- paste0("status-box status-", m$type)
    div(class = css, m$msg)
  })

  # ---- Clear ----------------------------------------------------------------

  observeEvent(input$btn_clear, {
    beta_data(data.frame(
      variant_id = character(0), study_id = character(0),
      beta       = numeric(0),   source   = character(0),
      stringsAsFactors = FALSE
    ))
    study_sources(setNames(character(0), character(0)))
    search_results(NULL)
    add_msg(list(type = "ok", msg = "Heatmap cleared."))
  })

  # ---- Heatmap summary sidebar widget ----------------------------------------

  output$heatmap_summary_ui <- renderUI({
    data <- beta_data()
    if (nrow(data) == 0)
      return(p(style = "color:#888; font-size:13px;", "No data yet."))

    variants <- unique(data$variant_id)
    studies  <- unique(data$study_id)
    n_na     <- sum(is.na(data$beta))

    tagList(
      p(strong(length(variants)), "variant(s),",
        strong(length(studies)),  "stud(y/ies)"),
      p(style = "color:#888; font-size:12px;",
        n_na, "of", nrow(data), "cells have no beta value"),
      tags$div(
        lapply(variants, function(v)
          tags$span(class = "chip", v))
      )
    )
  })

  # ---- Beta matrix ----------------------------------------------------------

  beta_matrix <- reactive({
    data <- beta_data()
    if (nrow(data) == 0) return(NULL)

    variants <- unique(data$variant_id)
    studies  <- unique(data$study_id)

    mat <- matrix(
      NA_real_,
      nrow     = length(variants),
      ncol     = length(studies),
      dimnames = list(variants, studies)
    )
    for (i in seq_len(nrow(data))) {
      v <- data$variant_id[i]
      s <- data$study_id[i]
      mat[v, s] <- data$beta[i]
    }
    mat
  })

  output$has_heatmap_data <- reactive({
    mat <- beta_matrix()
    !is.null(mat) && nrow(mat) >= 1 && ncol(mat) >= 1
  })
  outputOptions(output, "has_heatmap_data", suspendWhenHidden = FALSE)

  # ---- Interactive Heatmap --------------------------------------------------

  heatmap_reactive <- reactive({
    mat <- beta_matrix()
    req(mat)
    build_heatmap(mat)
  })

  observe({
    ht <- heatmap_reactive()
    req(ht)
    makeInteractiveComplexHeatmap(
      input, output, session,
      ht, "ht_main",
      click_action = function(df, output) {
        output[["ht_info"]] <- renderUI({
          if (is.null(df)) return(NULL)
          div(
            style = "padding: 8px; font-size: 13px;",
            tags$b("Selected cell"),
            tags$ul(
              tags$li(tags$b("Variant: "), df$row_label),
              tags$li(tags$b("Study: "),   df$column_label),
              tags$li(tags$b("Beta: "),
                      if (is.na(df$value)) "NA" else round(df$value, 6))
            )
          )
        })
      },
      brush_action = function(df, output) {
        output[["ht_info"]] <- renderUI({
          if (is.null(df) || nrow(df) == 0) return(NULL)
          n_non_na <- sum(!is.na(df$value))
          mean_val <- if (n_non_na > 0) round(mean(df$value, na.rm=TRUE), 4)
                      else NA
          div(
            style = "padding: 8px; font-size: 13px;",
            tags$b("Selected region"),
            tags$ul(
              tags$li(tags$b("Cells: "),    nrow(df)),
              tags$li(tags$b("Non-NA: "),   n_non_na),
              tags$li(tags$b("Mean beta: "), mean_val)
            )
          )
        })
      }
    )
  })

  # ---- Data table tab -------------------------------------------------------

  output$data_table <- DT::renderDataTable({
    data <- beta_data()
    if (nrow(data) == 0)
      return(DT::datatable(data.frame(Message = "No data yet.")))

    DT::datatable(
      data,
      rownames = FALSE,
      options  = list(pageLength = 15, scrollX = TRUE),
      caption  = "All variant × study beta values"
    )
  })

  # ---- Download -------------------------------------------------------------

  output$btn_download <- downloadHandler(
    filename = function() {
      paste0("genetic_association_heatmap_", Sys.Date(), ".png")
    },
    content = function(file) {
      mat <- beta_matrix()
      req(mat)
      ht  <- build_heatmap(mat)
      req(ht)
      height_px <- max(600, nrow(mat) * 40 + 280)
      width_px  <- max(900, ncol(mat) * 90 + 280)
      grDevices::png(file, width = width_px, height = height_px,
                     res = 120, type = "cairo-png")
      ComplexHeatmap::draw(ht, heatmap_legend_side = "bottom")
      grDevices::dev.off()
    }
  )

  # ---- Help tab -------------------------------------------------------------

  output$help_ui <- renderUI({
    tagList(
      h4("Quick Start"),
      tags$ol(
        tags$li("Select a ", tags$b("Data Source"), " (GWAS Catalog or HuGeAMP)."),
        tags$li("Enter a variant rsID in the ", tags$b("Search Variant"), " box and click ", tags$b("Search"), "."),
        tags$li("Review the associations found and click one of the ", tags$b("Add"), " buttons."),
        tags$li("Repeat for additional variants. The heatmap fills in beta values for all variant × study combinations automatically."),
        tags$li("Click a cell in the heatmap to see details, or brush to select a region."),
        tags$li("Download the heatmap as a PNG with the ", tags$b("Download PNG"), " button.")
      ),
      h4("Bulk Variant Input"),
      p("Paste a list of rsIDs — one per line — into the ",
        tags$b("Bulk Variant Input"), " box and click ",
        tags$b("Add All Variants"), ". Lines starting with ", code("#"),
        " are treated as comments and ignored. Duplicate IDs are silently deduplicated."),
      p("Enable ", tags$b("Limit to studies already in heatmap"),
        " to restrict results to studies that are already represented,",
        " e.g. to expand the heatmap by adding more variants across the same study set."),
      h4("Manual Entry"),
      p("You can also type a Variant ID and Study ID directly (e.g. ",
        code("rs7903146"), " and ", code("GCST001234"),
        ") and click ", tags$b("Add to Heatmap"), "."),
      h4("Colour Scale"),
      p("Blue = negative beta, white = zero, orange = positive beta. Grey cells have no available beta value."),
      h4("Data Sources"),
      tags$ul(
        tags$li(a("GWAS Catalog REST API v2",
                  href = "https://www.ebi.ac.uk/gwas/rest/api/v2/docs/reference",
                  target = "_blank")),
        tags$li(a("HuGeAMP BioIndex REST API",
                  href = "https://bioindex.hugeamp.org/docs",
                  target = "_blank")),
        tags$li(a("T2D Knowledge Portal examples",
                  href = "https://t2d.hugeamp.org/help.html?page=914",
                  target = "_blank"))
      )
    )
  })
}

# ---- Launch ------------------------------------------------------------------
shinyApp(ui = ui, server = server)
