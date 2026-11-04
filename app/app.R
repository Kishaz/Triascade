# =============================================================================
# ICI scorer - Shiny app. Score tumour RNA against a locked Immune Cascade Index
# reference (single-patient, selfnorm, cohort-independent).
#   * Patient report  : one-patient clinical read-out + PDF export
#   * Cohort explorer : archetype proportions, distributions, beeswarm + PDF report
# Input: RNA-seq (TPM / log-TPM / counts) as CSV, TSV or Excel (gene x sample;
# first column = gene symbol). Reference registry: app/references/*.rds (versioned,
# swappable). User-selectable colour palette (colour-blind-safe options).
# Run: shiny::runApp("app")
# =============================================================================
library(shiny); library(bslib); library(plotly); library(DT); library(ggplot2); library(patchwork); library(readxl)
source("ici_locked_reference.R"); source("score_cohort.R")
`%||%` <- function(a, b) if (is.null(a)) b else a

ref_files <- list.files("references", pattern = "[.]rds$", full.names = TRUE)
stopifnot(length(ref_files) >= 1)
ref_choices <- setNames(ref_files, vapply(ref_files, function(f){ m<-readRDS(f)
  sprintf("%s  ·  %d genes", m$version %||% basename(f), m$n_genes %||% length(m$panel_genes)) }, character(1)))

# ---- colour palettes (archetype + layer); user-selectable --------------------
PALETTES <- list(
  "Okabe-Ito (colour-blind safe)" = list(
    arch = c(intact="#009E73", cold="#0072B2", blocked="#D55E00", `no-kill`="#CC79A7"),
    lay  = c(Sensory="#E69F00", Integration="#0072B2", Execution="#CC79A7")),
  "House (plum / sand)" = list(
    arch = c(intact="#2e7d32", cold="#9C7A3E", blocked="#6A2C5A", `no-kill`="#B2567E"),
    lay  = c(Sensory="#9C7A3E", Integration="#6A2C5A", Execution="#B2567E")),
  "High contrast" = list(
    arch = c(intact="#1b9e77", cold="#1f78b4", blocked="#e6550d", `no-kill`="#6a51a3"),
    lay  = c(Sensory="#d95f02", Integration="#1f78b4", Execution="#6a51a3")),
  "Viridis" = list(
    arch = c(intact="#440154", cold="#31688e", blocked="#35b779", `no-kill`="#fde725"),
    lay  = c(Sensory="#443983", Integration="#21918c", Execution="#7ad151")))
DEF <- PALETTES[[1]]

SANS <- font_collection("Inter", "system-ui", "Helvetica Neue", "Arial", "sans-serif")
app_theme <- bs_theme(version = 5, bg = "#ffffff", fg = "#000000", primary = "#6A2C5A",
  secondary = "#9C7A3E", base_font = SANS, heading_font = SANS, "border-radius" = "0.65rem")

CALL <- list(
  intact  = list(t="Immune-engaged: all three layers intact",
    s="All three cascade layers are intact relative to the reference, a profile consistent with checkpoint-blockade benefit. If this tumour does not respond, resistance likely lies outside the measured immune cascade."),
  cold    = list(t="Cold: Sensory is the limiting step",
    s="Antigen presentation and IFN signalling are the rate-limiting layer: the tumour is immunologically cold. Checkpoint blockade alone is unlikely to engage an absent T-cell response."),
  blocked = list(t="Blocked: Integration is the limiting step",
    s="T-cell trafficking and infiltration are the limiting layer, a pattern consistent with stromal or TGF-beta mediated exclusion. Checkpoint blockade alone is unlikely to overcome exclusion."),
  `no-kill`= list(t="No-kill: Execution is the limiting step",
    s="A T-cell compartment is present but cytotoxic effector output is the limiting layer (effector-deficient or exhausted)."))
BIO <- list(
  intact  = "No single limiting layer. If non-responsive, consider non-cascade escape (antigen loss, B2M or MHC defects, non-immune resistance).",
  cold    = "Mechanistic rationale (investigational): priming strategies such as antigen-presentation or IFN induction, innate agonists, or radiation and vaccine combinations.",
  blocked = "Mechanistic rationale (investigational): exclusion-targeting such as anti-TGF-beta, anti-angiogenic, or stroma-directed combinations.",
  `no-kill`= "Mechanistic rationale (investigational): effector support or exhaustion reversal, for example LAG3, TIM-3 or TIGIT combinations, or cytokine support.")
GAUGE_CAP <- "Each bar is the patient's score for that layer as a percentile of the reference population (50 = median; shaded band = typical middle 50%). The lowest layer is the rate-limiting step."
RUO <- "Research use only. Investigational biomarker; not a validated clinical test. Does not replace approved biomarkers (PD-L1, TMB, MSI) or clinical judgement. Mechanistic rationale is hypothesis-generating, not a treatment recommendation."

# ---- multi-format reader (CSV / TSV / Excel) --------------------------------
read_expr_file <- function(path, name){
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xlsx","xls")) {
    df <- as.data.frame(readxl::read_excel(path))
    rn <- make.unique(as.character(df[[1]])); df[[1]] <- NULL; rownames(df) <- rn
  } else {
    l1 <- readLines(path, n = 1, warn = FALSE)
    sep <- if (grepl("\t", l1)) "\t" else if (grepl(";", l1) && !grepl(",", l1)) ";" else ","
    df <- read.delim(path, sep = sep, row.names = 1, check.names = FALSE)
  }
  df
}

# ---- plot builders (take colour vectors) ------------------------------------
build_gauge <- function(p, lay = DEF$lay){
  L <- c("Sensory","Integration","Execution"); z <- c(p$Sensory,p$Integration,p$Execution)
  g <- data.frame(layer=factor(L, levels=rev(L)), pct=round(pnorm(z)*100))
  g$floor <- as.character(g$layer)==p$floor_layer; g$fill <- unname(lay[as.character(g$layer)])
  fy <- which(rev(L)==p$floor_layer)
  ggplot(g, aes(pct, layer)) +
    annotate("rect", xmin=25, xmax=75, ymin=0.4, ymax=3.6, fill="#F1EADF") +
    geom_vline(xintercept=50, colour="black", linewidth=0.4) +
    geom_segment(aes(x=0, xend=pct, yend=layer), colour=g$fill, linewidth=3) +
    geom_point(aes(x=pct), colour=g$fill, size=7) +
    geom_point(data=g[g$floor,], aes(x=pct), shape=21, size=8.4, stroke=1.3, colour="black", fill=NA) +
    geom_text(aes(x=pct, label=pct), colour="white", fontface=2, size=3.3) +
    { if(length(fy)) annotate("text", x=50, y=fy-0.42, label="rate-limiting step", colour="black", fontface=2, size=3.5) } +
    scale_x_continuous("Percentile vs reference population  (50 = median)", limits=c(0,100), breaks=c(0,25,50,75,100)) +
    scale_y_discrete(NULL) + theme_minimal(base_size=13) +
    theme(panel.grid=element_blank(), axis.text=element_text(colour="black"),
          axis.text.y=element_text(face=2), axis.title.x=element_text(colour="black"),
          axis.line.x=element_line(colour="black", linewidth=0.3))
}
build_archetype_bar <- function(scores, arch = DEF$arch){
  tab <- as.data.frame(table(factor(scores$archetype, names(arch)))); names(tab)<-c("archetype","n")
  tab <- tab[tab$n>0,]; tab$pct <- round(100*tab$n/sum(tab$n))
  tab$archetype <- factor(tab$archetype, levels=rev(names(arch)[names(arch) %in% as.character(tab$archetype)]))
  ggplot(tab, aes(n, archetype, fill=archetype)) + geom_col(width=0.7) +
    geom_text(aes(label=sprintf("%d  (%d%%)", n, pct)), hjust=-0.1, colour="black", fontface=2, size=4.2) +
    scale_fill_manual(values=arch, guide="none") +
    scale_x_continuous("patients", expand=expansion(mult=c(0,0.34))) + scale_y_discrete(NULL) +
    theme_minimal(base_size=13) + theme(panel.grid=element_blank(),
      axis.text=element_text(colour="black", face=2), axis.title=element_text(colour="black"))
}
build_layer_violin <- function(scores, lay = DEF$lay){
  L <- do.call(rbind, lapply(c("Sensory","Integration","Execution"), function(l) data.frame(layer=l, z=scores[[l]])))
  L$layer <- factor(L$layer, c("Sensory","Integration","Execution"))
  ggplot(L, aes(layer, z, fill=layer)) + geom_hline(yintercept=0, colour="black", linewidth=0.3) +
    geom_violin(colour=NA, alpha=0.45, scale="width") + geom_jitter(width=0.1, size=0.6, alpha=0.55, colour="black") +
    scale_fill_manual(values=lay, guide="none") + labs(x=NULL, y="score (z vs reference)") +
    theme_minimal(base_size=13) + theme(panel.grid=element_blank(),
      axis.text=element_text(colour="black"), axis.title=element_text(colour="black"))
}
build_report <- function(p, a, cov, np, arch = DEF$arch, lay = DEF$lay){
  wrap <- function(s,w=104) paste(strwrap(s,w), collapse="\n"); ar <- p$archetype
  tp <- function() ggplot()+xlim(0,1)+ylim(0,1)+theme_void()
  top <- tp() +
    annotate("text",x=0,y=0.97,hjust=0,vjust=1,label="Immune Cascade Index: patient report",fontface=2,size=6.3,colour="#6A2C5A")+
    annotate("text",x=0,y=0.72,hjust=0,vjust=1,label=sprintf("Patient: %s",p$sample),size=4.3,colour="black")+
    annotate("text",x=0,y=0.5,hjust=0,vjust=1,label=CALL[[ar]]$t,fontface=2,size=5.3,colour=unname(arch[ar]))+
    annotate("text",x=0,y=0.32,hjust=0,vjust=1,label=wrap(CALL[[ar]]$s),size=3.8,colour="black",lineheight=0.95)
  bottom <- tp() +
    annotate("text",x=0,y=0.97,hjust=0,vjust=1,label=wrap(paste0("Limiting-layer biology: ",BIO[[ar]])),size=3.7,colour="black",lineheight=0.95)+
    annotate("text",x=0,y=0.63,hjust=0,vjust=1,label=wrap(GAUGE_CAP),size=3.3,colour="black",lineheight=0.95)+
    annotate("text",x=0,y=0.33,hjust=0,vjust=1,label=sprintf("Scored with %s (%d genes)  |  reference %s  |  coverage %d/%d  |  %s",
             a$version,a$n_genes,a$ref_label,cov,np,as.character(a$built)),size=3,colour="black")+
    annotate("text",x=0,y=0.2,hjust=0,vjust=1,label=wrap(RUO,110),size=2.9,colour="black",lineheight=0.95)
  top / build_gauge(p, lay) / bottom + plot_layout(heights=c(1.05,1.0,1.2))
}
build_cohort_report <- function(scores, a, cov, np, arch = DEF$arch, lay = DEF$lay){
  wrap <- function(s,w=120) paste(strwrap(s,w), collapse="\n")
  n<-nrow(scores); tab<-table(factor(scores$archetype, names(arch))); pct<-round(100*as.numeric(tab)/n)
  medz <- sapply(c("Sensory","Integration","Execution"), function(l) median(scores[[l]]))
  arch_line <- paste(sprintf("%s %d (%d%%)", names(tab), as.numeric(tab), pct)[as.numeric(tab)>0], collapse="     ")
  top <- ggplot()+xlim(0,1)+ylim(0,1)+theme_void()+
    annotate("text",x=0,y=0.95,hjust=0,vjust=1,label="Immune Cascade Index: cohort report",fontface=2,size=6.3,colour="#6A2C5A")+
    annotate("text",x=0,y=0.6,hjust=0,vjust=1,size=3.7,colour="black",
      label=sprintf("%d patients   |   coverage %d/%d panel genes   |   %s (%d genes)   |   reference %s", n, cov, np, a$version, a$n_genes, a$ref_label))+
    annotate("text",x=0,y=0.34,hjust=0,vjust=1,size=3.7,colour="black", label=paste("Archetype breakdown:   ", arch_line))+
    annotate("text",x=0,y=0.1,hjust=0,vjust=1,size=3.7,colour="black",
      label=sprintf("Median layer scores (z):   Sensory %.2f    Integration %.2f    Execution %.2f", medz[1], medz[2], medz[3]))
  bottom <- ggplot()+xlim(0,1)+ylim(0,1)+theme_void()+
    annotate("text",x=0,y=0.85,hjust=0,vjust=1,size=2.9,colour="black", label=sprintf("%s  |  %s", a$version, as.character(a$built)))+
    annotate("text",x=0,y=0.5,hjust=0,vjust=1,size=2.9,colour="black",lineheight=0.95,label=wrap(RUO))
  top / (build_archetype_bar(scores, arch) | build_layer_violin(scores, lay)) / bottom + plot_layout(heights=c(0.9,1.3,0.5))
}

ui <- page_navbar(
  title = tagList(span("ICI", style="font-weight:700;color:#6A2C5A"),
                  span("  ·  Immune Cascade Index", style="font-weight:500")),
  theme = app_theme, fillable = FALSE,
  sidebar = sidebar(width = 330, open = "open",
    selectInput("refver", "Panel / reference", choices = ref_choices, width = "100%"),
    markdown("Score tumour **RNA-seq** (TPM, log-TPM or counts) against the selected locked reference. Upload a **gene × sample** table (**CSV, TSV or Excel**; first column = gene symbol); a single patient is just one sample."),
    fileInput("file", NULL, accept = c(".csv",".tsv",".txt",".xlsx",".xls"), buttonLabel = "Browse…", placeholder = "CSV / TSV / Excel"),
    div(style="margin-top:-8px;margin-bottom:8px;font-size:0.74rem;color:#000000",
        "Unsure of the format? ", downloadLink("dl_template", "Download an input template")),
    actionButton("example", "Load example cohort", class = "btn-outline-secondary btn-sm", width = "100%"),
    selectInput("palette", "Colour palette", choices = names(PALETTES), width = "100%"),
    hr(), uiOutput("qc"), hr(),
    selectInput("dl_format", "Scores download format", choices = c("CSV","TSV","Excel"), width = "100%"),
    downloadButton("dl", "Download scores", class = "btn-primary btn-sm"),
    div(style="margin-top:auto;font-size:0.72rem;color:#000000", textOutput("reflab"))),

  nav_panel("Patient report",
    layout_columns(col_widths = c(5,7), fill = FALSE,
      card(selectInput("patient", "Patient", choices = NULL, width = "100%"),
           downloadButton("dl_pdf", "Download report (PDF)", class = "btn-primary btn-sm")),
      value_box("Overall read", textOutput("vb_call", inline = TRUE), theme = "primary")),
    card(card_header(textOutput("rep_hdr")),
      uiOutput("rep_call"),
      plotOutput("gauge", height = "220px"),
      div(style="font-size:0.78rem;color:#000000;margin-top:2px", GAUGE_CAP),
      uiOutput("rep_bio"), uiOutput("rep_stamp"),
      div(style="font-size:0.72rem;color:#000000;border-top:1px solid #000000;padding-top:8px;margin-top:6px", RUO))),

  nav_panel("Cohort explorer",
    layout_columns(fill = FALSE,
      value_box("Patients scored", textOutput("n_pat", inline = TRUE), theme = "primary"),
      value_box("Panel coverage", textOutput("cov", inline = TRUE), theme = value_box_theme(bg="#9C7A3E", fg="white")),
      value_box("Cascade-intact", textOutput("intact", inline = TRUE), theme = value_box_theme(bg="#C9A66B", fg="#000000"))),
    div(style="margin:2px 0 8px 0", downloadButton("dl_cohort", "Download cohort report (PDF)", class = "btn-primary btn-sm")),
    layout_columns(col_widths = c(5,7),
      card(full_screen=TRUE, card_header("Archetype proportions"), plotlyOutput("donut", height="300px")),
      card(full_screen=TRUE, card_header("Cohort distribution: the three layer scores"), plotlyOutput("cohort", height="300px"))),
    card(full_screen=TRUE, card_header("Per-patient profiles (hover for detail)"), plotlyOutput("patients", height="300px")),
    card(full_screen=TRUE, card_header("Per-patient scores"), DTOutput("tbl")))
)

server <- function(input, output, session) {
  assets <- reactive({ req(input$refver); readRDS(input$refver) })
  pal <- reactive({ PALETTES[[input$palette %||% names(PALETTES)[1]]] })
  raw <- reactiveVal(NULL)
  load_file <- function(path, name) {
    df <- tryCatch(read_expr_file(path, name), error = function(e) NULL)
    validate(need(!is.null(df) && ncol(df) >= 1, "Could not read a gene x sample table (CSV, TSV or Excel)."))
    raw(df)
  }
  observeEvent(input$file, { req(input$file); load_file(input$file$datapath, input$file$name) })
  observeEvent(input$example, { load_file("example_cohort.csv", "example_cohort.csv") })
  res <- reactive({ req(raw(), assets()); score_uploaded(raw(), assets()) })
  observeEvent(res(), { updateSelectInput(session, "patient", choices = res()$scores$sample) })
  patient <- reactive({ req(res(), input$patient); s<-res()$scores; s[s$sample==input$patient, ][1, ] })

  output$dl_template <- downloadHandler(
    filename = function() "ICI_input_template.csv",
    content = function(f){ pg <- assets()$panel_genes
      ex <- tryCatch(read.csv("example_cohort.csv", row.names=1, check.names=FALSE), error=function(e) NULL)
      if (!is.null(ex)) { rn<-toupper(rownames(ex)); g<-intersect(pg, rn)
        m <- ex[match(g, rn), seq_len(min(3,ncol(ex))), drop=FALSE]
        out <- data.frame(gene=g, round(m,3), check.names=FALSE); names(out)[-1] <- paste0("Sample", seq_len(ncol(out)-1))
      } else out <- data.frame(gene=pg, Sample1=NA, Sample2=NA, Sample3=NA)
      write.csv(out, f, row.names=FALSE) })
  output$reflab <- renderText({ sprintf("Reference cohort: %s", assets()$ref_label) })
  output$qc <- renderUI({ req(res()); cov<-res()$coverage; np<-res()$n_panel; pct<-round(100*cov/np)
    col <- if (pct>=90) "#1b5e20" else if (pct>=75) "#8a4b00" else "#8a0017"
    tagList(div(style=sprintf("font-size:0.8rem;color:%s;font-weight:600",col), sprintf("Coverage %d%%  (%d/%d panel genes)",pct,cov,np)),
      if (res()$logged) div(style="font-size:0.72rem;color:#000000","Input looked linear; log2 applied."),
      if (pct<75) div(style="font-size:0.72rem;color:#8a0017","Low coverage; interpret with caution.")) })

  output$vb_call <- renderText({ req(patient()); CALL[[patient()$archetype]]$t })
  output$rep_hdr <- renderText({ req(patient()); sprintf("ICI report: %s", patient()$sample) })
  output$rep_call <- renderUI({ req(patient()); a<-patient()$archetype
    div(div(style=sprintf("font-size:1.15rem;font-weight:700;color:%s",pal()$arch[a]), CALL[[a]]$t),
        div(style="font-size:0.92rem;margin-top:3px;color:#000000", CALL[[a]]$s)) })
  output$rep_bio <- renderUI({ req(patient()); div(style="margin-top:8px;font-size:0.86rem;color:#000000",
      tags$b("Limiting-layer biology.  "), BIO[[patient()$archetype]]) })
  output$rep_stamp <- renderUI({ req(res()); a<-assets()
    div(style="margin-top:8px;font-size:0.72rem;color:#000000",
        sprintf("Scored with %s (%d genes)  ·  reference %s  ·  coverage %d/%d  ·  %s",
                a$version, a$n_genes, a$ref_label, res()$coverage, res()$n_panel, as.character(a$built))) })
  output$gauge <- renderPlot({ req(patient()); build_gauge(patient(), pal()$lay) }, res=108)
  output$dl_pdf <- downloadHandler(
    filename = function() sprintf("ICI_report_%s.pdf", gsub("[^A-Za-z0-9]", "", patient()$sample)),
    content = function(f){ req(patient()); ggsave(f, build_report(patient(), assets(), res()$coverage, res()$n_panel, pal()$arch, pal()$lay),
                                                  width = 8.27, height = 8.8, device = "pdf") })

  output$donut <- renderPlotly({ req(res()); s<-res()$scores; arch<-pal()$arch
    tab <- as.data.frame(table(factor(s$archetype, names(arch)))); names(tab)<-c("a","n"); tab<-tab[tab$n>0,]
    plot_ly(tab, labels=~a, values=~n, type="pie", hole=0.58, sort=FALSE,
            marker=list(colors=unname(arch[as.character(tab$a)]), line=list(color="#ffffff", width=2)),
            textinfo="label+percent", textfont=list(color="black", size=14),
            hovertemplate="%{label}: %{value} patients (%{percent})<extra></extra>") |>
      layout(showlegend=FALSE, margin=list(t=12,b=12,l=12,r=12),
             font=list(color="black", family="Inter, system-ui, Arial, sans-serif")) |> config(displayModeBar=FALSE) })
  output$dl_cohort <- downloadHandler(
    filename = function() sprintf("ICI_cohort_report_%s.pdf", Sys.Date()),
    content = function(f){ req(res()); ggsave(f, build_cohort_report(res()$scores, assets(), res()$coverage, res()$n_panel, pal()$arch, pal()$lay),
                                              width = 10, height = 8, device = "pdf") })

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
    layout(font=list(family="Inter, system-ui, Arial, sans-serif", color="black"), margin=list(t=24,l=64,b=48,r=16)) |> config(displayModeBar=FALSE)
  output$cohort <- renderPlotly({ base_theme(build_layer_violin_jit(long(), pal()$lay)) })
  output$patients <- renderPlotly({ L<-long()
    p <- ggplot(L, aes(z, layer, colour=archetype, text=tt)) + geom_vline(xintercept=0, colour="black", linewidth=0.3) +
      geom_jitter(width=0, height=0.26, size=1.7, alpha=0.85) + scale_colour_manual(values=pal()$arch, name=NULL) +
      labs(x="score (z vs reference)", y=NULL) + theme_minimal(base_size=12) + theme(panel.grid=element_blank(), axis.text=element_text(colour="black"), axis.title=element_text(colour="black"))
    base_theme(p) |> layout(legend=list(orientation="h", y=1.14, x=0.5, xanchor="center"), margin=list(t=40,l=76,b=44,r=16)) })
  output$tbl <- renderDT({ req(res()); s<-res()$scores; arch<-pal()$arch
    d <- data.frame(Patient=s$sample, Sensory=round(s$Sensory,2), Integration=round(s$Integration,2),
      Execution=round(s$Execution,2), Floor=s$floor_layer, Archetype=s$archetype)
    datatable(d, rownames=FALSE, options=list(pageLength=8, dom="tp"), class="compact stripe") |>
      formatStyle("Archetype", backgroundColor=styleEqual(names(arch), paste0(arch,"33"))) })
  output$dl <- downloadHandler(
    filename = function(){ ext <- switch(input$dl_format %||% "CSV", CSV="csv", TSV="tsv", Excel="xlsx")
                           sprintf("ICI_scores_%s.%s", Sys.Date(), ext) },
    content = function(f){ req(res()); s<-res()$scores; fmt<-input$dl_format %||% "CSV"
      if (fmt=="Excel") openxlsx::write.xlsx(s, f)
      else if (fmt=="TSV") write.table(s, f, sep="\t", row.names=FALSE, quote=FALSE)
      else write.csv(s, f, row.names=FALSE) })
}

# cohort violin for the live (plotly) view, with hover text
build_layer_violin_jit <- function(L, lay){
  ggplot(L, aes(layer, z, fill=layer)) + geom_hline(yintercept=0, colour="black", linewidth=0.3) +
    geom_violin(colour=NA, alpha=0.35, scale="width") +
    geom_jitter(aes(text=tt), width=0.12, height=0, size=0.9, alpha=0.65, colour="black") +
    scale_fill_manual(values=lay, guide="none") + labs(x=NULL, y="score (z)") +
    theme_minimal(base_size=12) + theme(panel.grid=element_blank(), axis.text=element_text(colour="black"), axis.title=element_text(colour="black"))
}

shinyApp(ui, server)
