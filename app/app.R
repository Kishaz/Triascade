# =============================================================================
# ICI scorer - MVP Shiny app. Upload a gene x sample expression table, score each
# patient against the frozen locked 98-gene reference (single-patient, selfnorm),
# and explore cohort- and patient-level views of the 3 layer scores.
# Run:  shiny::runApp("app")
# =============================================================================
library(shiny); library(bslib); library(plotly); library(DT); library(ggplot2)
source("ici_locked_reference.R"); source("score_cohort.R")
assets <- readRDS("ici_reference.rds")

LAYCOL  <- c(Sensory="#9C7A3E", Integration="#6A2C5A", Execution="#B2567E")
ARCHCOL <- c(intact="#6A2C5A", cold="#C9A66B", blocked="#8A6D30", `no-kill`="#B2567E")
SANS <- font_collection("Inter", "system-ui", "Helvetica Neue", "Arial", "sans-serif")

theme <- bs_theme(version = 5, bg = "#ffffff", fg = "#262220",
  primary = "#6A2C5A", secondary = "#9C7A3E",
  base_font = SANS, heading_font = SANS, "border-radius" = "0.65rem")

ui <- page_sidebar(
  title = tagList(span("ICI", style="font-weight:700;color:#6A2C5A"),
                  span("  ·  Immune Cascade Index scorer", style="font-weight:500")),
  theme = theme,
  sidebar = sidebar(width = 320, open = "open",
    markdown("Upload a **gene × sample** expression table (CSV; first column = gene symbol, columns = samples). Each patient is scored against the **locked 98-gene reference** — one patient at a time, platform-agnostic."),
    fileInput("file", NULL, accept = c(".csv", ".tsv", ".txt"), buttonLabel = "Browse…", placeholder = "expression .csv"),
    actionButton("example", "Load example cohort", class = "btn-outline-secondary btn-sm", width = "100%"),
    hr(),
    uiOutput("qc"),
    hr(),
    downloadButton("dl", "Download scores", class = "btn-primary btn-sm"),
    div(style="margin-top:auto;font-size:0.72rem;color:#8A7F71",
        textOutput("reflab"))
  ),
  layout_columns(
    fill = FALSE,
    value_box("Patients scored", textOutput("n_pat", inline = TRUE), theme = "primary"),
    value_box("Panel coverage", textOutput("cov", inline = TRUE), theme = value_box_theme(bg="#9C7A3E", fg="white")),
    value_box("Cascade-intact", textOutput("intact", inline = TRUE), theme = value_box_theme(bg="#C9A66B", fg="#262220"))
  ),
  layout_columns(
    col_widths = c(6, 6),
    card(full_screen = TRUE,
      card_header("Cohort distribution — the three layer scores"),
      plotlyOutput("cohort", height = "330px")),
    card(full_screen = TRUE,
      card_header("Per-patient profiles  ·  hover a point for patient detail"),
      plotlyOutput("patients", height = "330px"))
  ),
  card(full_screen = TRUE, card_header("Per-patient scores"), DTOutput("tbl"))
)

server <- function(input, output, session) {
  res <- reactiveVal(NULL)
  load_df <- function(path) {
    df <- tryCatch(read.csv(path, row.names = 1, check.names = FALSE), error = function(e) NULL)
    validate(need(!is.null(df) && ncol(df) >= 1, "Could not read a gene x sample table from that file."))
    res(score_uploaded(df, assets))
  }
  observeEvent(input$file, { req(input$file); load_df(input$file$datapath) })
  observeEvent(input$example, { load_df("example_cohort.csv") })

  long <- reactive({ req(res()); s <- res()$scores
    L <- do.call(rbind, lapply(c("Sensory","Integration","Execution"), function(l)
      data.frame(sample=s$sample, layer=l, z=s[[l]], archetype=s$archetype, floor=s$floor_layer)))
    L$layer <- factor(L$layer, c("Sensory","Integration","Execution"))
    L$tt <- sprintf("<b>%s</b><br>Sensory %.2f · Integration %.2f · Execution %.2f<br>floor: %s  ·  archetype: %s",
                    L$sample, s$Sensory[match(L$sample,s$sample)], s$Integration[match(L$sample,s$sample)],
                    s$Execution[match(L$sample,s$sample)], L$floor, L$archetype)
    L })

  output$n_pat  <- renderText({ req(res()); res()$n_samples })
  output$cov    <- renderText({ req(res()); sprintf("%d / %d genes", res()$coverage, res()$n_panel) })
  output$intact <- renderText({ req(res()); s<-res()$scores; sprintf("%d of %d (%.0f%%)",
                    sum(s$archetype=="intact"), nrow(s), 100*mean(s$archetype=="intact")) })
  output$reflab <- renderText({ sprintf("Locked reference: %s", assets$ref_label) })
  output$qc <- renderUI({ req(res())
    cov <- res()$coverage; np <- res()$n_panel; pct <- round(100*cov/np)
    col <- if (pct>=90) "#2e7d32" else if (pct>=75) "#b26a00" else "#b00020"
    tagList(
      div(style=sprintf("font-size:0.8rem;color:%s;font-weight:600", col),
          sprintf("Coverage %d%%  (%d/%d panel genes)", pct, cov, np)),
      if (res()$logged) div(style="font-size:0.72rem;color:#8A7F71", "Input looked linear — log2 applied."),
      if (pct<75) div(style="font-size:0.72rem;color:#b00020", "Low coverage — interpret scores with caution.")
    )
  })

  base_theme <- function(p) ggplotly(p, tooltip = "text") |>
    layout(font = list(family = "Inter, system-ui, Arial, sans-serif"),
           margin = list(t = 24, l = 64, b = 48, r = 16)) |>
    config(displayModeBar = FALSE)

  output$cohort <- renderPlotly({ L <- long()
    p <- ggplot(L, aes(layer, z, fill = layer)) +
      geom_hline(yintercept = 0, colour = "grey70", linewidth = 0.3) +
      geom_violin(colour = NA, alpha = 0.35, scale = "width") +
      geom_jitter(aes(text = tt), width = 0.12, height = 0, size = 0.9, alpha = 0.6,
                  colour = "#4C1F41") +
      scale_fill_manual(values = LAYCOL, guide = "none") +
      labs(x = NULL, y = "score (z)") +
      theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
    base_theme(p) })

  output$patients <- renderPlotly({ L <- long()
    p <- ggplot(L, aes(z, layer, colour = archetype, text = tt)) +
      geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.3) +
      geom_jitter(width = 0, height = 0.26, size = 1.7, alpha = 0.85) +
      scale_colour_manual(values = ARCHCOL, name = NULL) +
      labs(x = "layer score (z vs reference)", y = NULL) +
      theme_minimal(base_size = 12) + theme(panel.grid.minor = element_blank())
    base_theme(p) |> layout(legend = list(orientation = "h", y = 1.14, x = 0.5, xanchor = "center"),
                            margin = list(t = 40, l = 76, b = 44, r = 16)) })

  output$tbl <- renderDT({ req(res()); s <- res()$scores
    d <- data.frame(Patient=s$sample, Sensory=round(s$Sensory,2), Integration=round(s$Integration,2),
                    Execution=round(s$Execution,2), Floor=s$floor_layer, Archetype=s$archetype)
    datatable(d, rownames = FALSE, options = list(pageLength = 8, dom = "tp"), class = "compact stripe") |>
      formatStyle("Archetype", backgroundColor = styleEqual(names(ARCHCOL), paste0(ARCHCOL,"33")))
  })

  output$dl <- downloadHandler(
    filename = function() sprintf("ICI_scores_%s.csv", Sys.Date()),
    content = function(f) { req(res()); write.csv(res()$scores, f, row.names = FALSE) })
}

shinyApp(ui, server)
