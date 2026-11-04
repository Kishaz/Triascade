# =============================================================================
# ICI scorer - Shiny app. Score tumour RNA against a locked Immune Cascade Index
# reference (single-patient, selfnorm, cohort-independent).
#   * Patient report  : one-patient clinical-style read-out (call + 3-layer gauge
#                       + named block + biological rationale + QC + version stamp)
#   * Cohort explorer : distributions + per-patient beeswarm for trial/biomarker use
# Future-proof: references live in app/references/*.rds, each self-describing its
# panel version / gene count / reference cohort. Drop a new reference (e.g. an
# expanded panel, or one built on a newer platform) into that folder and it
# appears in the selector automatically - no code change. Input: RNA-seq (TPM /
# log-TPM / counts), gene x sample CSV (first col = symbol). Run: shiny::runApp("app")
# =============================================================================
library(shiny); library(bslib); library(plotly); library(DT); library(ggplot2)
source("ici_locked_reference.R"); source("score_cohort.R")
`%||%` <- function(a, b) if (is.null(a)) b else a

# ---- reference registry (self-describing, versioned, swappable) -------------
ref_files <- list.files("references", pattern="\\.rds$", full.names=TRUE)
stopifnot(length(ref_files) >= 1)
ref_choices <- setNames(ref_files, vapply(ref_files, function(f){ m<-readRDS(f)
  sprintf("%s  ·  %d genes", m$version %||% basename(f), m$n_genes %||% length(m$panel_genes)) }, character(1)))

LAYCOL  <- c(Sensory="#9C7A3E", Integration="#6A2C5A", Execution="#B2567E")
ARCHCOL <- c(intact="#2e7d32", cold="#9C7A3E", blocked="#6A2C5A", `no-kill`="#B2567E")
SANS <- font_collection("Inter", "system-ui", "Helvetica Neue", "Arial", "sans-serif")
theme <- bs_theme(version = 5, bg = "#ffffff", fg = "#262220", primary = "#6A2C5A",
  secondary = "#9C7A3E", base_font = SANS, heading_font = SANS, "border-radius" = "0.65rem")

CALL <- list(
  intact  = list(t="Immune-engaged — all three layers intact",
    s="All three cascade layers are intact relative to the reference, a profile consistent with checkpoint-blockade benefit. If this tumour does not respond, resistance likely lies outside the measured immune cascade."),
  cold    = list(t="Cold — Sensory is the limiting step",
    s="Antigen presentation / IFN signalling is the rate-limiting layer: the tumour is immunologically cold. Checkpoint blockade alone is unlikely to engage an absent T-cell response."),
  blocked = list(t="Blocked — Integration is the limiting step",
    s="T-cell trafficking / infiltration is the limiting layer, a pattern consistent with stromal / TGF-β-mediated exclusion. Checkpoint blockade alone is unlikely to overcome exclusion."),
  `no-kill`= list(t="No-kill — Execution is the limiting step",
    s="A T-cell compartment is present but cytotoxic effector output is the limiting layer (effector-deficient / exhausted)."))
BIO <- list(
  intact  = "No single limiting layer. If non-responsive, consider non-cascade escape (antigen loss, β2M/MHC defects, non-immune resistance).",
  cold    = "Mechanistic rationale (investigational): priming strategies — antigen-presentation / IFN induction, innate agonists, or radiation / vaccine combinations.",
  blocked = "Mechanistic rationale (investigational): exclusion-targeting — anti-TGF-β, anti-angiogenic, or stroma-directed combinations.",
  `no-kill`= "Mechanistic rationale (investigational): effector support / exhaustion reversal — e.g. LAG3 / TIM-3 / TIGIT combinations or cytokine support.")
RUO <- "Research use only. Investigational biomarker — not a validated clinical test. Does not replace approved biomarkers (PD-L1, TMB, MSI) or clinical judgement. Mechanistic rationale is hypothesis-generating, not a treatment recommendation."

ui <- page_navbar(
  title = tagList(span("ICI", style="font-weight:700;color:#6A2C5A"),
                  span("  ·  Immune Cascade Index", style="font-weight:500")),
  theme = theme, fillable = FALSE,
  sidebar = sidebar(width = 330, open = "open",
    selectInput("refver", "Panel / reference", choices = ref_choices, width = "100%"),
    markdown("Score tumour **RNA-seq** (TPM, log-TPM or counts) against the selected locked reference. Upload a **gene × sample** CSV (first column = gene symbol); a single patient is just one sample."),
    fileInput("file", NULL, accept = c(".csv", ".tsv", ".txt"), buttonLabel = "Browse…", placeholder = "expression .csv"),
    actionButton("example", "Load example cohort", class = "btn-outline-secondary btn-sm", width = "100%"),
    hr(), uiOutput("qc"), hr(),
    downloadButton("dl", "Download scores", class = "btn-primary btn-sm"),
    div(style="margin-top:auto;font-size:0.72rem;color:#8A7F71", textOutput("reflab"))),

  nav_panel("Patient report",
    layout_columns(col_widths = c(5,7), fill = FALSE,
      card(selectInput("patient", "Patient", choices = NULL, width = "100%")),
      value_box("Overall read", textOutput("vb_call", inline = TRUE), theme = "primary")),
    card(card_header(textOutput("rep_hdr")),
      uiOutput("rep_call"),
      plotOutput("gauge", height = "230px"),
      uiOutput("rep_bio"),
      uiOutput("rep_stamp"),
      div(class="text-muted", style="font-size:0.72rem;border-top:1px solid #eee;padding-top:8px;margin-top:6px", RUO))),

  nav_panel("Cohort explorer",
    layout_columns(fill = FALSE,
      value_box("Patients scored", textOutput("n_pat", inline = TRUE), theme = "primary"),
      value_box("Panel coverage", textOutput("cov", inline = TRUE), theme = value_box_theme(bg="#9C7A3E", fg="white")),
      value_box("Cascade-intact", textOutput("intact", inline = TRUE), theme = value_box_theme(bg="#C9A66B", fg="#262220"))),
    layout_columns(col_widths = c(6,6),
      card(full_screen=TRUE, card_header("Cohort distribution — the three layer scores"), plotlyOutput("cohort", height="320px")),
      card(full_screen=TRUE, card_header("Per-patient profiles  ·  hover for detail"), plotlyOutput("patients", height="320px"))),
    card(full_screen=TRUE, card_header("Per-patient scores"), DTOutput("tbl")))
)

server <- function(input, output, session) {
  assets <- reactive({ req(input$refver); readRDS(input$refver) })
  raw <- reactiveVal(NULL)
  read_expr <- function(path) {
    df <- tryCatch(read.csv(path, row.names = 1, check.names = FALSE), error = function(e) NULL)
    validate(need(!is.null(df) && ncol(df) >= 1, "Could not read a gene x sample table from that file."))
    raw(df)
  }
  observeEvent(input$file, { req(input$file); read_expr(input$file$datapath) })
  observeEvent(input$example, { read_expr("example_cohort.csv") })
  res <- reactive({ req(raw(), assets()); score_uploaded(raw(), assets()) })
  observeEvent(res(), { updateSelectInput(session, "patient", choices = res()$scores$sample) })
  patient <- reactive({ req(res(), input$patient); s<-res()$scores; s[s$sample==input$patient, ][1, ] })

  output$reflab <- renderText({ a<-assets(); sprintf("Reference cohort: %s", a$ref_label) })
  output$qc <- renderUI({ req(res()); cov<-res()$coverage; np<-res()$n_panel; pct<-round(100*cov/np)
    col <- if (pct>=90) "#2e7d32" else if (pct>=75) "#b26a00" else "#b00020"
    tagList(div(style=sprintf("font-size:0.8rem;color:%s;font-weight:600",col), sprintf("Coverage %d%%  (%d/%d panel genes)",pct,cov,np)),
      if (res()$logged) div(style="font-size:0.72rem;color:#8A7F71","Input looked linear — log2 applied."),
      if (pct<75) div(style="font-size:0.72rem;color:#b00020","Low coverage — interpret with caution.")) })

  # ---- Patient report ----
  output$vb_call <- renderText({ req(patient()); CALL[[patient()$archetype]]$t })
  output$rep_hdr <- renderText({ req(patient()); sprintf("ICI report — %s", patient()$sample) })
  output$rep_call <- renderUI({ req(patient()); a<-patient()$archetype
    div(div(style=sprintf("font-size:1.15rem;font-weight:700;color:%s",ARCHCOL[a]), CALL[[a]]$t),
        div(style="font-size:0.92rem;margin-top:3px", CALL[[a]]$s)) })
  output$rep_bio <- renderUI({ req(patient()); div(style="margin-top:8px;font-size:0.86rem",
      tags$b("Limiting-layer biology.  "), BIO[[patient()$archetype]]) })
  output$rep_stamp <- renderUI({ req(res()); a<-assets()
    div(style="margin-top:8px;font-size:0.72rem;color:#6b6258",
        sprintf("Scored with %s (%d genes)  ·  reference %s  ·  coverage %d/%d  ·  %s",
                a$version, a$n_genes, a$ref_label, res()$coverage, res()$n_panel, as.character(a$built))) })
  output$gauge <- renderPlot({ req(patient()); p<-patient()
    g <- data.frame(layer=factor(c("Sensory","Integration","Execution"), levels=c("Execution","Integration","Sensory")),
                    z=c(p$Sensory,p$Integration,p$Execution))
    g$pct <- round(pnorm(g$z)*100); g$floor <- (as.character(g$layer)==p$floor_layer)
    ggplot(g, aes(z, layer)) +
      annotate("rect", xmin=-1, xmax=1, ymin=0.4, ymax=3.6, fill="#F2ECE2") +
      geom_vline(xintercept=0, colour="grey65", linewidth=0.4) +
      geom_segment(aes(x=0, xend=z, yend=layer, colour=floor), linewidth=2.4) +
      geom_point(aes(colour=floor), size=5) +
      geom_text(aes(x=z, label=sprintf("%+.2f", z)), vjust=-1.1, size=3.6, colour="#262220", fontface=2) +
      geom_text(aes(x=-3.1, label=paste0(pct, "th pct")), hjust=0, size=3.2, colour="#6b6258") +
      { if(any(g$floor)) annotate("text", x=g$z[g$floor], y=as.numeric(g$layer[g$floor])-0.34,
          label="← limiting step", colour=ARCHCOL[[p$archetype]], fontface=2, size=3.3) } +
      scale_colour_manual(values=c(`TRUE`=unname(ARCHCOL[[p$archetype]]), `FALSE`="#C9BEB0"), guide="none") +
      scale_x_continuous("layer score (z vs reference population)", limits=c(-3.3,3.3), breaks=seq(-3,3,1)) +
      scale_y_discrete(NULL) + theme_minimal(base_size=13) +
      theme(panel.grid.major.y=element_blank(), panel.grid.minor=element_blank(),
            axis.text.y=element_text(face=2, colour="#262220")) }, res=108)

  # ---- Cohort explorer ----
  long <- reactive({ req(res()); s<-res()$scores
    L <- do.call(rbind, lapply(c("Sensory","Integration","Execution"), function(l)
      data.frame(sample=s$sample, layer=l, z=s[[l]], archetype=s$archetype, floor=s$floor_layer)))
    L$layer <- factor(L$layer, c("Sensory","Integration","Execution"))
    L$tt <- sprintf("<b>%s</b><br>S %.2f · I %.2f · E %.2f<br>floor: %s · %s",
      L$sample, s$Sensory[match(L$sample,s$sample)], s$Integration[match(L$sample,s$sample)],
      s$Execution[match(L$sample,s$sample)], L$floor, L$archetype); L })
  output$n_pat  <- renderText({ req(res()); res()$n_samples })
  output$cov    <- renderText({ req(res()); sprintf("%d / %d genes", res()$coverage, res()$n_panel) })
  output$intact <- renderText({ req(res()); s<-res()$scores; sprintf("%d of %d (%.0f%%)", sum(s$archetype=="intact"), nrow(s), 100*mean(s$archetype=="intact")) })
  base_theme <- function(p) ggplotly(p, tooltip="text") |>
    layout(font=list(family="Inter, system-ui, Arial, sans-serif"), margin=list(t=24,l=64,b=48,r=16)) |> config(displayModeBar=FALSE)
  output$cohort <- renderPlotly({ L<-long()
    p <- ggplot(L, aes(layer, z, fill=layer)) + geom_hline(yintercept=0, colour="grey70", linewidth=0.3) +
      geom_violin(colour=NA, alpha=0.35, scale="width") +
      geom_jitter(aes(text=tt), width=0.12, height=0, size=0.9, alpha=0.6, colour="#4C1F41") +
      scale_fill_manual(values=LAYCOL, guide="none") + labs(x=NULL, y="score (z)") +
      theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank()); base_theme(p) })
  output$patients <- renderPlotly({ L<-long()
    p <- ggplot(L, aes(z, layer, colour=archetype, text=tt)) + geom_vline(xintercept=0, colour="grey70", linewidth=0.3) +
      geom_jitter(width=0, height=0.26, size=1.7, alpha=0.85) + scale_colour_manual(values=ARCHCOL, name=NULL) +
      labs(x="score (z vs reference)", y=NULL) + theme_minimal(base_size=12) + theme(panel.grid.minor=element_blank())
    base_theme(p) |> layout(legend=list(orientation="h", y=1.14, x=0.5, xanchor="center"), margin=list(t=40,l=76,b=44,r=16)) })
  output$tbl <- renderDT({ req(res()); s<-res()$scores
    d <- data.frame(Patient=s$sample, Sensory=round(s$Sensory,2), Integration=round(s$Integration,2),
      Execution=round(s$Execution,2), Floor=s$floor_layer, Archetype=s$archetype)
    datatable(d, rownames=FALSE, options=list(pageLength=8, dom="tp"), class="compact stripe") |>
      formatStyle("Archetype", backgroundColor=styleEqual(names(ARCHCOL), paste0(ARCHCOL,"33"))) })
  output$dl <- downloadHandler(filename=function() sprintf("ICI_scores_%s.csv", Sys.Date()),
    content=function(f){ req(res()); write.csv(res()$scores, f, row.names=FALSE) })
}

shinyApp(ui, server)
