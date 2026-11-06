# =============================================================================
# CascadeScope - score tumor RNA against a locked Immune Cascade Index (ICI)
# reference (single-patient, selfnorm, cohort-independent). One consistent scale
# throughout: the STANDARDIZED SCORE vs the reference population (layer scores:
# 0 = reference mean, in SD units; left of the typical-range band = lowest-scoring layer).
#   * Patient report  : one-patient read-out (call, gauge, drill-down) + PDF
#   * Gene drivers     : how all 98 panel genes scored for the selected patient
#   * Cohort explorer  : archetype proportions, distributions, click-through beeswarm
#   * About            : what the tool is, the cascade model, the scale, RUO
# Input: RNA-seq (TPM / log-TPM / counts) as CSV, TSV or Excel (gene x sample;
# first column = gene symbol). Reference registry: app/references/*.rds.
# Run: shiny::runApp("app")
# =============================================================================
library(shiny); library(bslib); library(plotly); library(DT); library(ggplot2); library(patchwork)
# readxl (Excel input) and openxlsx (Excel download) are used namespaced and are
# optional: CSV / TSV work without them, so the app still loads where they are
# unavailable (e.g. the in-browser webR build).
source("ici_locked_reference.R"); source("score_cohort.R")
`%||%` <- function(a, b) if (is.null(a)) b else a


ref_files <- list.files("references", pattern = "[.]rds$", full.names = TRUE)
stopifnot(length(ref_files) >= 1)
ref_choices <- setNames(ref_files, vapply(ref_files, function(f){ m<-readRDS(f)
  sprintf("%s  ·  %d genes", m$version %||% basename(f), m$n_genes %||% length(m$panel_genes)) }, character(1)))

# ---- colour palettes (archetype + layer); user-selectable --------------------
PALETTES <- list(
  "Okabe-Ito (color-blind safe)" = list(
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
BAND <- 0.674                                                # typical range = middle 50% of the reference (SD units)

# ---- contrast helpers (WCAG relative luminance) -----------------------------
# ink(bg): the higher-contrast of black/white for text on a colour.
# onwhite(hex): darken a hue just enough to read sharply on a white background.
.lum1 <- function(v){ f <- function(c) ifelse(c <= 0.03928, c/12.92, ((c+0.055)/1.055)^2.4)
  w <- f(v); 0.2126*w[1] + 0.7152*w[2] + 0.0722*w[3] }
ink <- function(bg){ L <- .lum1(grDevices::col2rgb(bg)[,1]/255)
  if ((L+0.05)/0.05 >= 1.05/(L+0.05)) "#000000" else "#ffffff" }
onwhite <- function(hex, target = 4.5){ v <- grDevices::col2rgb(hex)[,1]/255
  for (f in seq(1, 0, by = -0.02)){ cc <- v*f
    if (1.05/(.lum1(cc)+0.05) >= target) return(grDevices::rgb(cc[1],cc[2],cc[3])) }
  "#000000" }

SANS <- font_collection("Inter", "system-ui", "Helvetica Neue", "Arial", "sans-serif")
app_theme <- bs_theme(version = 5, bg = "#ffffff", fg = "#000000", primary = "#6A2C5A",
  secondary = "#9C7A3E", base_font = SANS, heading_font = SANS, "border-radius" = "0.65rem")

# ---- interpretation text -----------------------------------------------------
GLOSS <- c(intact="inflamed / hot", cold="immune-desert", blocked="immune-excluded", `no-kill`="inflamed-dysfunctional")
CALL <- list(
  intact  = list(short="Immune-intact", t="Immune-intact (inflamed / hot): no substantially reduced layer",
    s="All three cascade layers are within or above the typical range of the reference. If this tumor is non-responsive, resistance likely lies outside the measured immune cascade. The biological and therapeutic implications are investigational."),
  cold    = list(short="Cold", t="Cold (immune-desert): Sensory is the lowest-scoring layer",
    s="Antigen presentation and IFN signaling are the lowest-scoring layer relative to the reference, a profile consistent with an immunologically cold tumor. The biological and therapeutic implications are investigational."),
  blocked = list(short="Blocked", t="Blocked (immune-excluded): Integration is the lowest-scoring layer",
    s="T-cell trafficking and infiltration are the lowest-scoring layer relative to the reference, a pattern consistent with stromal or TGF-beta mediated exclusion. The biological and therapeutic implications are investigational."),
  `no-kill`= list(short="No-kill", t="No-kill (inflamed-dysfunctional): Execution is the lowest-scoring layer",
    s="Immune cells have reached the tumor, but cytotoxic effector output (the shared killing program of CD8 T cells and NK cells) is the lowest-scoring layer (effector-deficient or exhausted)."))
BIO_LEAD <- list(
  intact  = "No layer stands out as low; all three are within range of the reference. If the tumor is non-responsive, look outside the cascade.",
  cold    = "Priming is the gap: the immune response is not being started.",
  blocked = "Entry is the gap: immune cells are held at the margin and kept out of the tumor.",
  `no-kill`= "Killing is the gap: cytotoxic lymphocytes (CD8 T and NK cells) are present but their killing output is low.")
BIO_CHIPS <- list(
  intact  = c("Antigen loss / B2M / MHC defects", "Non-immune resistance", "Confirm with PD-L1 / TMB / MSI"),
  cold    = c("Antigen presentation / IFN signaling", "Innate sensing (STING / TLR)"),
  blocked = c("TGF-β signaling", "Stromal exclusion", "Vascular barriers", "Impaired trafficking"),
  `no-kill`= c("Effector exhaustion", "Inhibitory checkpoints (LAG-3 / TIM-3 / TIGIT)", "Low cytotoxic output"))
GAUGE_CAP <- "Each bar is the patient's standardized layer score versus the reference population (0 = reference mean, in SD units; shaded band = typical range, the middle 50% of the reference). The starred layer is the patient's relative weakest layer within the analyzed cohort (for a single uploaded sample, versus the reference)."
RUO <- "Research use only. Investigational biomarker; not a validated clinical test. Does not replace approved biomarkers (PD-L1, TMB, MSI) or clinical judgment. Mechanistic rationale is hypothesis-generating, not a treatment recommendation."

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

# ---- sub-module metadata ----------------------------------------------------
SM_ORD   <- c("SENSORY","INTPOS","INTNEG","EXECPOS","EXECNEG")
SM_LAB   <- c(SENSORY="Antigen presentation / IFN", INTPOS="T-cell trafficking",
              INTNEG="Stroma / Treg suppression", EXECPOS="Cytotoxic effector", EXECNEG="Exhaustion / checkpoint")
SM_LAYER <- c(SENSORY="Sensory", INTPOS="Integration", INTNEG="Integration", EXECPOS="Execution", EXECNEG="Execution")
SM_SIGN  <- c(SENSORY="+", INTPOS="+", INTNEG="−", EXECPOS="+", EXECNEG="−")
LANE_LAB <- c(SENSORY="Sensory", INTPOS="Integration (+)", INTNEG="Integration (−)",
              EXECPOS="Execution (+)", EXECNEG="Execution (−)")
SM_SHORT <- c(SENSORY="Sensory +", INTPOS="Trafficking +", INTNEG="Stroma/Treg −",
              EXECPOS="Effector +", EXECNEG="Exhaustion −")
CASC_LOW <- "#9C7A3E"; CASC_MID <- "#F2ECE0"; CASC_HIGH <- "#6A2C5A"   # gold -> cream -> plum

# ---- gauge: three layer standardized scores (diverging) ---------------------
build_gauge <- function(p, lay = DEF$lay){
  L <- c("Sensory","Integration","Execution"); z <- c(p$Sensory,p$Integration,p$Execution)
  g <- data.frame(layer=factor(L, levels=rev(L)), z=z)
  g$floor <- as.character(g$layer)==p$floor_layer; g$fill <- unname(lay[as.character(g$layer)])
  g$bar <- ifelse(g$z >= 0, "#6A2C5A", "#9C7A3E")   # above reference = plum, below = gold (matches gene colours)
  lim <- max(2.2, max(abs(z),na.rm=TRUE)+0.85); fl <- g[g$floor,,drop=FALSE]
  g$lx <- g$z + ifelse(g$z>=0, 0.15, -0.15); g$hj <- ifelse(g$z>=0, 0, 1)  # labels sit just outside each dot
  pl <- ggplot(g, aes(z, layer)) +
    annotate("rect", xmin=-BAND, xmax=BAND, ymin=0.4, ymax=3.6, fill="#ECECEC") +
    geom_vline(xintercept=0, colour="black", linewidth=0.4) +
    geom_segment(aes(x=0, xend=z, yend=layer), colour=g$bar, linewidth=2.6, lineend="round") +
    geom_point(aes(x=z), fill=g$bar, colour="white", size=6.8, shape=21, stroke=1.1) +
    geom_text(aes(x=lx, label=sprintf("%+.1f", z), hjust=hj), colour="black", fontface=2, size=3.6)
  if (nrow(fl)) {                                      # floor marker: a true 5-point star drawn as a polygon
    fx <- fl$z[1]; fy <- as.integer(fl$layer)[1]       # (font-independent, so it renders in webR too)
    aa <- (90 + seq(0, 324, by=36)) * pi/180; rr <- rep(c(1, 0.4), 5)
    pl <- pl + annotate("polygon", x = fx + 0.052*rr*cos(aa), y = fy + 0.095*rr*sin(aa),
                        fill="black", colour="white", linewidth=0.35)
  }
  pl +
    scale_x_continuous("Standardized score vs reference  (0 = mean)",
                       limits=c(-lim,lim), breaks=seq(-floor(lim), floor(lim), 1)) +
    scale_y_discrete(NULL, labels=c(Sensory="Sensory (See)",
                                    Integration="Integration (Reach)",
                                    Execution="Execution (Kill)")) + theme_minimal(base_size=13) +
    theme(panel.grid=element_blank(), axis.text=element_text(colour="black"),
          axis.text.y=element_text(face=2, size=10.5), axis.title.x=element_text(colour="black"),
          axis.line.x=element_line(colour="black", linewidth=0.3),
          legend.position="bottom", legend.margin=margin(t=-2,b=0),
          legend.text=element_text(colour="black", size=12.5),
          legend.key.size=unit(1.2, "lines"))
}

# ---- per-patient cascade drill-down (sub-module decomposition) --------------
drill_text <- function(smz, floor_layer){
  if (is.null(smz) || floor_layer=="none") return("All three layers are within the typical range; no sub-module stands out.")
  subs <- SM_ORD[SM_LAYER[SM_ORD]==floor_layer]
  z <- as.numeric(smz[subs]); names(z) <- subs
  unfav <- ifelse(SM_SIGN[subs]=="+", -z, z); w <- subs[which.max(unfav)]
  sprintf("In the %s layer, %s %s is furthest from the reference (standardized score %+.2f).", floor_layer,
          if (SM_SIGN[w]=="+") "low" else "high", SM_LAB[w], z[w])
}
build_drilldown <- function(smz, floor_layer, lay = DEF$lay){
  d <- data.frame(sm=SM_ORD, z=as.numeric(smz[SM_ORD]))
  d$layer <- SM_LAYER[d$sm]; d$sign <- SM_SIGN[d$sm]
  d$lab <- sprintf("%s (%s)", SM_LAB[d$sm], d$sign); d$fill <- unname(lay[d$layer])
  d$unfav <- ifelse(d$sign=="+", -d$z, d$z)
  d$driver <- d$layer==floor_layer & floor_layer!="none" & d$unfav==ave(d$unfav, d$layer, FUN=max) & d$unfav>0
  d$lab <- factor(d$lab, levels=rev(d$lab)); lim <- max(1.5, max(abs(d$z),na.rm=TRUE)+0.3)
  ggplot(d, aes(z, lab)) +
    annotate("rect", xmin=-BAND, xmax=BAND, ymin=0.4, ymax=5.6, fill="#ECECEC") +
    geom_vline(xintercept=0, colour="black", linewidth=0.4) +
    geom_segment(aes(x=0, xend=z, yend=lab), colour=d$fill, linewidth=3) +
    geom_point(aes(x=z), colour=d$fill, size=6.5) +
    geom_point(data=d[d$driver,,drop=FALSE], aes(x=z), shape=21, size=8, stroke=1.3, colour="black", fill=NA) +
    scale_x_continuous("Standardized score vs reference  (0 = mean)", limits=c(-lim,lim), breaks=seq(-2,2,1)) +
    scale_y_discrete(NULL) + theme_minimal(base_size=12) +
    theme(panel.grid=element_blank(), axis.text=element_text(colour="black"),
          axis.text.y=element_text(face=2), axis.title.x=element_text(colour="black"),
          axis.line.x=element_line(colour="black", linewidth=0.3))
}

# ---- gene drivers: barcode (all layers) / labeled heatmap (one layer) -------
build_gene_heat <- function(df){
  df <- df[order(df$module_ord, df$dev),]
  df$gene <- factor(df$gene, levels=unique(df$gene))
  df$col  <- factor(SM_SHORT[df$module], levels=SM_SHORT[SM_ORD[SM_ORD %in% df$module]])
  lim <- max(0.6, max(abs(df$dev), na.rm=TRUE))
  ggplot(df, aes(col, gene)) +
    geom_tile(aes(fill=dev, colour=factor(on)), width=0.9, height=0.92, linewidth=0.5) +
    scale_fill_gradient2(low=CASC_LOW, mid=CASC_MID, high=CASC_HIGH, midpoint=0, limits=c(-lim,lim),
                         name="z vs reference") +
    scale_colour_manual(values=c(`0`="#D9D9D9", `1`="black"), guide="none") +
    scale_x_discrete(NULL, position="top") + scale_y_discrete(NULL) +
    guides(fill=guide_colourbar(barwidth=unit(7,"lines"), barheight=unit(0.7,"lines"), title.position="top")) +
    theme_minimal(base_size=12) +
    theme(panel.grid=element_blank(), axis.text=element_text(colour="black"),
          axis.text.x=element_text(face=2, size=12.5), axis.text.y=element_text(size=9.5),
          panel.border=element_rect(colour="black", fill=NA, linewidth=1),
          plot.margin=margin(6,14,6,14),
          legend.position="bottom", legend.title=element_text(colour="black", size=11.5, face=2),
          legend.text=element_text(colour="black", size=10.5))
}

# ---- per-patient cascade tree (top-down): gene panel -> sub-modules -> layers -
build_cascade_tree <- function(df, p, smz, lay = DEF$lay){
  if (is.null(smz)) smz <- setNames(rep(NA_real_, length(SM_ORD)), SM_ORD)
  mods <- SM_ORD[SM_ORD %in% df$module]
  df <- df[order(match(df$module, SM_ORD)), ]
  gap <- 4; xc <- 0; centers <- numeric(0); df$x <- NA_real_
  for (m in mods){ rows <- which(df$module==m); n <- length(rows)
    df$x[rows] <- xc + seq_len(n); centers[m] <- xc + (n+1)/2; xc <- xc + n + gap }
  df$cx <- unname(centers[df$module])
  lyr_x <- c(Sensory=unname(centers["SENSORY"]),
             Integration=mean(centers[c("INTPOS","INTNEG")]),
             Execution=mean(centers[c("EXECPOS","EXECNEG")]))
  gy <- 10; sy <- 5.5; ly <- 1          # gene row (top) / sub-module row / layer row (bottom)
  snode <- data.frame(x=unname(centers[mods]), y=sy, z=as.numeric(smz[mods]),
                      sign=unname(SM_SIGN[mods]), layer=unname(SM_LAYER[mods]),
                      lab=unname(SM_SHORT[mods]), stringsAsFactors=FALSE)
  lnode <- data.frame(layer=c("Sensory","Integration","Execution"),
                      x=unname(lyr_x[c("Sensory","Integration","Execution")]), y=ly,
                      z=c(p$Sensory,p$Integration,p$Execution), letter=c("S","I","E"),
                      stringsAsFactors=FALSE)
  lnode$lab <- sprintf("%s  %+.2f", lnode$layer, lnode$z)
  e2 <- data.frame(x=snode$x, xend=unname(lyr_x[snode$layer]), y=sy, yend=ly,
                   col=ifelse(snode$sign=="+", CASC_HIGH, CASC_LOW), lwd=1.4 + 1.8*pmin(1.2, abs(snode$z)))
  e1 <- data.frame(x=df$x, xend=df$cx, y=gy, yend=sy, dev=df$dev, on=df$on)
  lim <- max(1, max(abs(c(df$dev, snode$z, lnode$z)), na.rm=TRUE)); xr <- range(df$x)
  lfill <- scales::gradient_n_pal(c(CASC_LOW, CASC_MID, CASC_HIGH))(pmax(0, pmin(1, (lnode$z + lim)/(2*lim))))
  lletter <- vapply(lfill, ink, character(1))   # S/I/E letter colour contrasts its node
  ggplot() +
    geom_segment(data=e2, aes(x=x,xend=xend,y=y,yend=yend), colour=e2$col, linewidth=e2$lwd, lineend="round") +
    geom_segment(data=e1[e1$on==0,,drop=FALSE], aes(x=x,xend=xend,y=y,yend=yend, colour=dev), linewidth=0.25, alpha=0.30) +
    geom_segment(data=e1[e1$on==1,,drop=FALSE], aes(x=x,xend=xend,y=y,yend=yend, colour=dev), linewidth=0.35, alpha=0.85) +
    geom_point(data=df, aes(x=x, y=gy, colour=dev, shape=factor(on)), size=3.3, stroke=1.05) +
    geom_text(data=df, aes(x=x, y=gy+0.5, label=gene), angle=90, hjust=0, vjust=0.5, size=2.9, colour="black") +
    geom_point(data=snode, aes(x=x, y=y, colour=z), size=9) +
    geom_label(data=snode, aes(x=x, y=y+1.3, label=lab), fontface=2, size=3.5, colour="black",
               fill="white", linewidth=0, label.padding=unit(0.1,"lines")) +
    geom_point(data=lnode, aes(x=x, y=y, colour=z), size=17) +
    geom_text(data=lnode, aes(x=x, y=y, label=letter), colour=lletter, fontface=2, size=5.4) +
    geom_label(data=lnode, aes(x=x, y=y-1.5, label=lab), fontface=2, size=4.3, colour="black",
               fill="white", linewidth=0, label.padding=unit(0.12,"lines")) +
    scale_colour_gradient2(low=CASC_LOW, mid=CASC_MID, high=CASC_HIGH, midpoint=0, limits=c(-lim,lim),
                           na.value="#BBBBBB", name="Standardized score (z) vs reference") +
    scale_shape_manual(values=c(`0`=0, `1`=15),
                       labels=c(`0`="off  (< ref median)", `1`="on  (≥ ref median)"), name=NULL) +
    scale_x_continuous(limits=c(xr[1]-1.5, xr[2]+1.5)) +
    scale_y_continuous(limits=c(ly-2.8, gy+8.4)) +
    coord_cartesian(clip="off") +
    guides(colour=guide_colourbar(order=1, barwidth=unit(13,"lines"), barheight=unit(0.85,"lines"), title.position="top"),
           shape=guide_legend(order=2, override.aes=list(size=5.5))) +
    theme_void(base_size=13) +
    theme(legend.position="bottom", legend.box="horizontal", legend.spacing.x=unit(0.6,"lines"),
          plot.margin=margin(16,12,8,12),
          legend.title=element_text(colour="black", size=13, face=2),
          legend.text=element_text(colour="black", size=12.5))
}
# ---- three layer heatmaps side by side, each on its OWN scale -----------------
build_layer_heatmaps <- function(df){
  ps <- lapply(c("Sensory","Integration","Execution"), function(L){
    build_gene_heat(df[df$layer==L,,drop=FALSE]) +
      ggtitle(L) + theme(plot.title=element_text(face=2, size=15, hjust=0.5, colour="black")) })
  ps[[1]] | ps[[2]] | ps[[3]]
}

# ---- cohort plot builders ---------------------------------------------------
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
    scale_fill_manual(values=lay, guide="none") + labs(x=NULL, y="standardized score vs reference") +
    theme_minimal(base_size=13) + theme(panel.grid=element_blank(),
      axis.text=element_text(colour="black"), axis.text.x=element_text(face=2, size=13),
      axis.title=element_text(colour="black"))
}
# cohort violin for the live (plotly) view, with hover text
build_layer_violin_jit <- function(L, lay){
  ggplot(L, aes(layer, z, fill=layer)) + geom_hline(yintercept=0, colour="black", linewidth=0.3) +
    geom_violin(colour=NA, alpha=0.35, scale="width") +
    geom_jitter(aes(text=tt), width=0.12, height=0, size=0.9, alpha=0.65, colour="black") +
    scale_fill_manual(values=lay, guide="none") + labs(x=NULL, y="standardized score") +
    theme_minimal(base_size=12) + theme(panel.grid=element_blank(), axis.text=element_text(colour="black"),
      axis.text.x=element_text(face=2, size=13), axis.title=element_text(colour="black"))
}

# ---- PDF reports ------------------------------------------------------------
build_report <- function(p, a, cov, np, arch = DEF$arch, lay = DEF$lay, smz = NULL, study = NULL){
  wrap <- function(s,w=104) paste(strwrap(s,w), collapse="\n"); ar <- p$archetype
  tp <- function() ggplot()+xlim(0,1)+ylim(0,1)+theme_void()
  readout <- sprintf("Sensory %+.2f      ·      Integration %+.2f      ·      Execution %+.2f        (standardized vs reference; 0 = mean)",
                     p$Sensory, p$Integration, p$Execution)
  hdr <- if (!is.null(study) && nzchar(study)) sprintf("Study: %s        Generated: %s", study, Sys.Date()) else sprintf("Generated: %s", Sys.Date())
  drv <- if (!is.null(smz)) drill_text(smz, p$floor_layer) else NULL
  top <- tp() +
    annotate("text",x=0,y=0.98,hjust=0,vjust=1,label="CascadeScope: patient report",fontface=2,size=5.9,colour="#6A2C5A")+
    annotate("text",x=0,y=0.82,hjust=0,vjust=1,label=hdr,size=3.2,colour="black")+
    annotate("text",x=0,y=0.70,hjust=0,vjust=1,label=sprintf("Patient: %s",p$sample),fontface=2,size=3.9,colour="black")+
    annotate("text",x=0,y=0.55,hjust=0,vjust=1,label="RESULT",fontface=2,size=3.0,colour="#6A2C5A")+
    annotate("text",x=0,y=0.47,hjust=0,vjust=1,label=CALL[[ar]]$t,fontface=2,size=4.2,colour=onwhite(unname(arch[ar])))+
    annotate("text",x=0,y=0.33,hjust=0,vjust=1,label=wrap(paste0("What this means: ",CALL[[ar]]$s),88),size=3.2,colour="black",lineheight=0.95)+
    annotate("text",x=0,y=0.08,hjust=0,vjust=1,label=wrap(paste0("Layer scores:  ",readout),92),size=2.95,colour="black",lineheight=0.95)
  bottom <- tp() +
    (if (!is.null(drv)) annotate("text",x=0,y=0.98,hjust=0,vjust=1,label=wrap(drv,84),fontface=2,size=3.3,colour="black")) +
    annotate("text",x=0,y=0.80,hjust=0,vjust=1,label=wrap(paste0("Lowest-layer biology (investigational): ",BIO_LEAD[[ar]],"  Programs for investigation: ",paste(BIO_CHIPS[[ar]],collapse=", "),"."),92),size=3.0,colour="black",lineheight=0.95)+
    annotate("text",x=0,y=0.52,hjust=0,vjust=1,label=wrap(GAUGE_CAP,100),size=2.9,colour="black",lineheight=0.95)+
    annotate("text",x=0,y=0.30,hjust=0,vjust=1,label=wrap(sprintf("Method: scored with %s (%d genes); reference %s; within-sample standardization; panel coverage %d/%d; reference built %s.",
             a$version,a$n_genes,a$ref_label,cov,np,as.character(a$built)),104),size=2.85,colour="black",lineheight=0.95)+
    annotate("text",x=0,y=0.13,hjust=0,vjust=1,label=wrap(RUO,104),size=2.8,colour="black",lineheight=0.95)
  if (!is.null(smz)) {
    dd <- build_drilldown(smz, p$floor_layer, lay) +
          ggtitle("Cascade drill-down: sub-modules of the lowest-scoring layer") +
          theme(plot.title = element_text(face=2, size=11, colour="black"))
    return(top / build_gauge(p, lay) / dd / bottom + plot_layout(heights=c(1.5,1.0,1.2,1.55)))
  }
  top / build_gauge(p, lay) / bottom + plot_layout(heights=c(1.5,1.0,1.5))
}
build_cohort_report <- function(scores, a, cov, np, arch = DEF$arch, lay = DEF$lay, study = NULL){
  wrap <- function(s,w=120) paste(strwrap(s,w), collapse="\n")
  n<-nrow(scores); tab<-table(factor(scores$archetype, names(arch))); pct<-round(100*as.numeric(tab)/n)
  medz <- sapply(c("Sensory","Integration","Execution"), function(l) median(scores[[l]]))
  arch_line <- paste(sprintf("%s %d (%d%%)", names(tab), as.numeric(tab), pct)[as.numeric(tab)>0], collapse="     ")
  hdr <- if (!is.null(study) && nzchar(study)) sprintf("Study: %s        Generated: %s", study, Sys.Date()) else sprintf("Generated: %s", Sys.Date())
  top <- ggplot()+xlim(0,1)+ylim(0,1)+theme_void()+
    annotate("text",x=0,y=0.96,hjust=0,vjust=1,label="CascadeScope: cohort report",fontface=2,size=6.0,colour="#6A2C5A")+
    annotate("text",x=0,y=0.76,hjust=0,vjust=1,size=3.5,colour="black",label=hdr)+
    annotate("text",x=0,y=0.56,hjust=0,vjust=1,size=3.5,colour="black",
      label=sprintf("%d patients   |   coverage %d/%d panel genes   |   %s (%d genes)   |   reference %s", n, cov, np, a$version, a$n_genes, a$ref_label))+
    annotate("text",x=0,y=0.36,hjust=0,vjust=1,size=3.5,colour="black", label=paste("Archetype breakdown:   ", arch_line))+
    annotate("text",x=0,y=0.16,hjust=0,vjust=1,size=3.5,colour="black",
      label=sprintf("Median layer scores (standardized vs reference):   Sensory %+.2f    Integration %+.2f    Execution %+.2f", medz[1], medz[2], medz[3]))
  bottom <- ggplot()+xlim(0,1)+ylim(0,1)+theme_void()+
    annotate("text",x=0,y=0.85,hjust=0,vjust=1,size=2.9,colour="black", label=sprintf("%s  |  %s", a$version, as.character(a$built)))+
    annotate("text",x=0,y=0.5,hjust=0,vjust=1,size=2.9,colour="black",lineheight=0.95,label=wrap(RUO))
  top / (build_archetype_bar(scores, arch) | build_layer_violin(scores, lay)) / bottom + plot_layout(heights=c(0.9,1.3,0.5))
}

# ---- small UI helpers -------------------------------------------------------
chips <- function(items, col){
  div(style="display:flex;flex-wrap:wrap;gap:6px;margin-top:8px",
      lapply(items, function(x) span(x, style=sprintf(
        "display:inline-block;padding:3px 11px;border:1px solid %s;border-radius:999px;font-size:0.82rem;color:#000;background:%s22", col, col))))
}

about_tab <- nav_panel(title = tagList(icon("circle-info"), " About"), value = "about",
  layout_columns(col_widths = c(8,4),
    card(card_header("CascadeScope and the Immune Cascade Index"),
      markdown(paste(
        "**CascadeScope** is the software; the **Immune Cascade Index (ICI)** is the method it implements.",
        "",
        "- **Immune Cascade Index (ICI)** is the *science*: the 98-gene, three-layer framework for reading the antitumor immune response, developed and validated in the manuscript.",
        "- **CascadeScope** is the *software*: it computes and visualizes the ICI for one patient at a time. The method can evolve independently; the tool is the interface built around it.",
        "",
        "CascadeScope scores a tumor's RNA against a locked ICI reference and shows, for one patient, **which stage of the antitumor immune response reads weakest**.",
        "",
        "The ICI reads the response as a three-stage **cascade**, from sensing to killing, and scores each stage against the reference:",
        "",
        "1. **Sensory** - antigen presentation and IFN signaling (is the tumor being *seen*?)",
        "2. **Integration** - T-cell trafficking and infiltration vs stromal / Treg exclusion (do immune cells *get in*?)",
        "3. **Execution** - cytotoxic effector output vs exhaustion / checkpoints (do the cells that got in *kill*?)",
        "",
        "**Execution** measures the shared cytotoxic machinery of **CD8 T cells and NK cells** (granzymes, perforin, granulysin, NKG7): it reflects their combined killing output and does not separate the two lineages.",
        "",
        "The patient's **relative weakest layer** is flagged for attention - the stage that reads lowest within the analyzed cohort (or versus the reference for a single uploaded sample). The archetype is a relative, population-dependent label, not a proven causal bottleneck; the continuous layer scores above are population-independent.",
        sep="\n")),
      card_header("One scale throughout"),
      markdown(paste(
        "Every score - layers, sub-modules and individual genes - is shown as a **standardized score versus the reference population**:",
        "",
        "- **0** = the reference mean (for per-gene calls, the reference median).",
        "- Units are **standard deviations** of the reference.",
        "- The **shaded band** on each gauge is the *typical range* (middle 50% of the reference).",
        "- A layer **left of the band** sits below the reference typical range; the starred layer is the patient's relative weakest layer within the analyzed cohort.",
        "",
        "Genes are colored by the same idea: <span style='color:#6A2C5A'>**plum = above the reference**</span>,",
        "<span style='color:#9C7A3E'>**gold = below**</span>. A filled square (black outline) marks a gene the score counts as \"on\".",
        sep="\n"))),
    card(card_header("The four archetypes"),
      tags$table(class="table table-sm", style="font-size:0.85rem",
        tags$thead(tags$tr(tags$th("Call"), tags$th("Phenotype"), tags$th("Lowest-scoring layer"))),
        tags$tbody(
          tags$tr(tags$td(tags$b("Immune-intact")), tags$td("inflamed / hot"), tags$td("none")),
          tags$tr(tags$td(tags$b("Cold")), tags$td("immune-desert"), tags$td("Sensory")),
          tags$tr(tags$td(tags$b("Blocked")), tags$td("immune-excluded"), tags$td("Integration")),
          tags$tr(tags$td(tags$b("No-kill")), tags$td("inflamed-dysfunctional"), tags$td("Execution")))),
      card_header("Panel & input"),
      markdown(paste(
        "- **98-gene panel**, locked reference (versioned; swappable).",
        "- Input: **RNA-seq** (TPM, log-TPM or counts), gene × sample, as CSV / TSV / Excel.",
        "- A single patient is just a one-sample file; scoring is **cohort-independent**.",
        sep="\n")),
      div(style="font-size:0.74rem;color:#000;border-top:1px solid #000;padding-top:8px;margin-top:4px", RUO))))

# ---- first-load landing (shown until a file or the example is loaded) --------
landing_card <- card(class = "border-0 shadow-none",
  div(style="max-width:760px;margin:26px auto 8px auto;padding:8px 16px;text-align:center",
    HTML('<svg width="54" height="54" viewBox="0 0 120 120" style="display:block;margin:0 auto 10px auto">
      <rect width="120" height="120" rx="26" fill="#6A2C5A"/>
      <rect x="28" y="32" width="64" height="14" rx="7" fill="#ffffff"/>
      <rect x="28" y="54" width="46" height="14" rx="7" fill="#ffffff" fill-opacity="0.82"/>
      <rect x="28" y="76" width="30" height="14" rx="7" fill="#E8B54B"/></svg>'),
    h3(style="font-weight:800;color:#6A2C5A;margin-bottom:4px", "Profile antitumor immune state from tumor RNA"),
    div(style="font-size:0.95rem;color:#333;max-width:620px;margin:0 auto 14px auto",
        "CascadeScope applies the 98-gene Immune Cascade Index to quantify three functional states against a frozen reference. Scores are computed for a single sample or a whole cohort."),
    div(style="font-weight:800;letter-spacing:.3px;color:#6A2C5A;font-size:1.0rem;margin-bottom:16px",
        "Recognition  →  Access  →  Effector function"),
    div(style="display:flex;gap:14px;justify-content:center;align-items:flex-start;flex-wrap:wrap;text-align:left",
        actionButton("example_land", "Try example cohort", class = "btn-primary", icon = icon("play")),
        div(style="min-width:280px",
            fileInput("file_land", NULL, accept = c(".csv",".tsv",".txt",".xlsx",".xls"),
                      buttonLabel = "Upload RNA-seq…", placeholder = "CSV / TSV / Excel"),
            div(style="margin-top:-10px;font-size:0.72rem;color:#777;text-align:center",
                "gene × sample table (first column = gene symbol)"))),
    div(style="margin-top:16px;font-size:0.72rem;color:#000000", "Research use only. Investigational biomarker; not a validated clinical test.")))

ui <- page_navbar(
  id = "nav",
  window_title = "CascadeScope",
  title = div(style="display:flex;align-items:center;gap:9px",
    HTML('<svg width="38" height="38" viewBox="0 0 120 120" style="display:block;flex:none">
      <rect width="120" height="120" rx="26" fill="#6A2C5A"/>
      <rect x="28" y="32" width="64" height="14" rx="7" fill="#ffffff"/>
      <rect x="28" y="54" width="46" height="14" rx="7" fill="#ffffff" fill-opacity="0.82"/>
      <rect x="28" y="76" width="30" height="14" rx="7" fill="#E8B54B"/>
    </svg>'),
    div(style="width:2px;align-self:stretch;margin:4px 0;background:#6A2C5A;border-radius:1px;opacity:.55"),
    div(style="line-height:1.12",
      div(style="font-weight:800;font-size:1.15rem;letter-spacing:.2px",
          HTML('<span style="color:#6A2C5A">Cascade</span><span style="color:#9C7A3E">Scope</span>')),
      div("Immune Cascade Index", style="font-size:0.74rem;font-weight:600;color:#777777;letter-spacing:.3px"))),
  theme = app_theme, fillable = FALSE,
  header = tagList(
    # Prime the Plotly.js load at app start (see output$`.plotly_primer`): the main
    # plotly library is attached only when a plot first renders and loads async, so
    # without this the cohort plots - inserted together when their lazily-rendered
    # tab opens - can run before it loads ("Plotly is not defined", a blank widget
    # that never retries). This always-present 1px plot triggers the load up front.
    div(style="position:absolute;width:1px;height:1px;overflow:hidden;opacity:0;pointer-events:none",
        plotlyOutput(".plotly_primer", width="1px", height="1px")),
    tags$head(
      tags$link(rel="icon", type="image/svg+xml", href="icon.svg"),
      tags$style(HTML(paste(
        # highlight the active navbar tab
        ".navbar-nav .nav-link.active, .nav-link.active{font-weight:700 !important;color:#6A2C5A !important;}",
        ".navbar-nav .nav-link.active{border-bottom:3px solid #6A2C5A;}",
        sep="\n")))),
    tags$script(HTML(
      "document.addEventListener('shown.bs.collapse', function(){ setTimeout(function(){ window.dispatchEvent(new Event('resize')); }, 50); });"))),
  sidebar = sidebar(width = 335, title = "Data & settings",
                    open = list(desktop = "open", mobile = "closed"),
    selectInput("refver", "Panel / reference", choices = ref_choices, width = "100%"),
    markdown("Score tumor **RNA-seq** (TPM, log-TPM or counts) against the selected locked reference. Upload a **gene × sample** table (**CSV, TSV or Excel**; first column = gene symbol); a single patient is just one sample."),
    fileInput("file", NULL, accept = c(".csv",".tsv",".txt",".xlsx",".xls"), buttonLabel = "Browse…", placeholder = "CSV / TSV / Excel"),
    div(style="margin-top:-8px;margin-bottom:8px;font-size:0.74rem;color:#000000",
        "Unsure of the format? ", downloadLink("dl_template", "Download an input template")),
    actionButton("example", "Load example cohort", class = "btn-outline-secondary btn-sm", width = "100%"),
    selectInput("patient", "Selected patient", choices = NULL, width = "100%"),
    textInput("study", "Study / cohort label (optional)", width = "100%", placeholder = "e.g. KEYNOTE-052 baseline"),
    div(style="margin-top:-8px;margin-bottom:8px;font-size:0.72rem;color:#000000",
        "Tags every report and download so you can tell runs apart."),
    selectInput("palette", "Colour palette", choices = names(PALETTES), width = "100%"),
    hr(), uiOutput("qc"), hr(),
    selectInput("dl_format", "Scores download format", choices = c("CSV","TSV","Excel"), width = "100%"),
    downloadButton("dl", "Download scores", class = "btn-primary btn-sm"),
    div(style="margin-top:auto;font-size:0.72rem;color:#000000", textOutput("reflab"))),

  nav_panel(title = tagList(icon("user-check"), " Patient report"), value = "patient",
    conditionalPanel("!output.has_data", landing_card),
    conditionalPanel("output.has_data", uiOutput("qc_report"),
    card(
      card_header(div(style="display:flex;justify-content:space-between;align-items:center;gap:12px;flex-wrap:wrap",
        div(textOutput("rep_hdr", inline = TRUE), style="font-weight:600"),
        uiOutput("rep_floor_badge"))),
      uiOutput("rep_call"),
      plotOutput("gauge", height = "230px"),
      uiOutput("rep_readout"),
      div(style="display:flex;justify-content:space-between;align-items:flex-end;gap:12px;flex-wrap:wrap;margin-top:4px",
        div(style="font-size:0.76rem;color:#000000;max-width:72%", GAUGE_CAP),
        downloadButton("dl_pdf", "Download report (PDF)", class = "btn-primary btn-sm"))),
    accordion(id = "acc", open = c("meaning","biology"),
      accordion_panel("What this result means", value = "meaning",
        uiOutput("rep_meaning")),
      accordion_panel("Lowest-layer biology (programs for investigation)", value = "biology",
        uiOutput("rep_bio")),
      accordion_panel("Cascade drill-down: sub-modules within the lowest-scoring layer", value = "drill",
        plotOutput("drill", height = "230px"),
        uiOutput("drill_text"),
        div(style="font-size:0.74rem;color:#000000;margin-top:4px",
            "Activating programs (+): higher is better. Suppressive programs (−): lower is better. The ringed point is the sub-module furthest from the reference in the lowest-scoring layer.")),
      accordion_panel("Gene drivers: how all 98 panel genes scored", value = "genes",
        div(style="display:flex;gap:16px;align-items:center;flex-wrap:wrap;margin-bottom:4px",
          div(style="min-width:260px", selectInput("gene_layer", "Show", width="100%",
              choices = c("Cascade tree (all layers)"="tree","Layer heatmaps (compare, own scale)"="heat"))),
          div(style="font-size:0.82rem;color:#000000",
              "Genes feed sub-modules, which feed the three layers. Node & gene colour = standardized score vs reference (",
              tags$b(style=sprintf("color:%s", CASC_LOW), "gold = low"), ", ",
              tags$b(style=sprintf("color:%s", CASC_HIGH), "plum = high"),
              "); filled square = gene “on” (≥ reference median), open square = “off”.")),
        card(full_screen = TRUE,
          div(style="overflow-x:auto", plotOutput("gene_plot", width = "1850px", height = "720px")),
          div(style="font-size:0.72rem;color:#000000;margin-top:2px",
              "Scroll horizontally to see all 98 genes at full size; use the expand control (top-right) for a full-screen view.")),
        div(style="font-size:0.78rem;color:#000000;margin-top:4px",
            tags$b("Cascade tree:"), " top row = the 98-gene panel; thick edges from sub-module to layer are ",
            tags$b(style=sprintf("color:%s", CASC_HIGH), "plum = activating (+)"), ", ",
            tags$b(style=sprintf("color:%s", CASC_LOW), "gold = suppressive (−)"),
            "; the layer node shows its standardized score. ",
            tags$b("Layer heatmaps:"), " the three layers side by side, each on its own colour scale for head-to-head comparison. Use the expand control at the panel's top-right for a full-screen view.")),
      accordion_panel("Method, reference & disclaimer", value = "method",
        uiOutput("rep_stamp"),
        div(style="font-size:0.74rem;color:#000000;border-top:1px solid #000000;padding-top:8px;margin-top:6px", RUO))))),

  nav_panel(title = tagList(icon("chart-column"), " Cohort explorer"), value = "cohort",
    layout_columns(fill = FALSE,
      value_box("Patients scored", textOutput("n_pat", inline = TRUE), theme = "primary"),
      value_box("Panel coverage", textOutput("cov", inline = TRUE), theme = value_box_theme(bg="#9C7A3E", fg=ink("#9C7A3E"))),
      value_box("Immune-intact", textOutput("intact", inline = TRUE), theme = value_box_theme(bg="#C9A66B", fg="#000000"))),
    div(style="margin:2px 0 8px 0", downloadButton("dl_cohort", "Download cohort report (PDF)", class = "btn-primary btn-sm")),
    layout_columns(col_widths = c(5,7),
      card(full_screen=TRUE, card_header("Archetype proportions"), plotlyOutput("donut", height="300px"),
           div(style="font-size:0.78rem;color:#000000",
               div(tags$b("Canonical phenotypes:")),
               div(tags$b("cold"), " = immune-desert  ·  ", tags$b("blocked"), " = immune-excluded  ·  ",
                   tags$b("no-kill"), " = inflamed-dysfunctional  ·  ", tags$b("intact"), " = inflamed / hot."))),
      card(full_screen=TRUE, card_header("Cohort distribution: the three layer scores"), plotlyOutput("cohort", height="300px"))),
    card(full_screen=TRUE, card_header("Per-patient profiles"),
      div(style="font-size:0.8rem;color:#000000;margin-bottom:4px",
          tags$b("Tip:"), " click a point, or a table row below, to open that patient's report."),
      plotlyOutput("patients", height="300px")),
    card(full_screen=TRUE, card_header("Per-patient scores (click a row to open the report)"), DTOutput("tbl"))),

  about_tab
)

server <- function(input, output, session) {
  assets <- reactive({ req(input$refver); readRDS(input$refver) })
  pal <- reactive({ PALETTES[[input$palette %||% names(PALETTES)[1]]] })
  studylab <- reactive({ s <- input$study; if (is.null(s)) "" else trimws(s) })
  study_slug <- function(){ sl <- gsub("[^A-Za-z0-9]+","_", studylab()); sl <- gsub("^_+|_+$","",sl); sl }
  raw <- reactiveVal(NULL)
  load_file <- function(path, name) {
    df <- tryCatch(read_expr_file(path, name), error = function(e) NULL)
    validate(need(!is.null(df) && ncol(df) >= 1, "Could not read a gene x sample table (CSV, TSV or Excel)."))
    raw(df)
  }
  observeEvent(input$file, { req(input$file); load_file(input$file$datapath, input$file$name) })
  observeEvent(input$example, { load_file("example_cohort.csv", "example_cohort.csv") })
  observeEvent(input$example_land, { load_file("example_cohort.csv", "example_cohort.csv") })
  observeEvent(input$file_land, { req(input$file_land); load_file(input$file_land$datapath, input$file_land$name) })
  output$has_data <- reactive({ !is.null(raw()) })
  outputOptions(output, "has_data", suspendWhenHidden = FALSE)
  res <- reactive({ req(raw(), assets()); score_uploaded(raw(), assets()) })
  observeEvent(res(), { updateSelectInput(session, "patient", choices = res()$scores$sample) })
  patient <- reactive({ req(res(), input$patient); s<-res()$scores; s[s$sample==input$patient, ][1, ] })
  psub <- reactive({ req(res(), input$patient); sm<-res()$submod
                     if (is.null(sm) || !(input$patient %in% rownames(sm))) return(NULL); sm[input$patient, ] })

  # per-patient gene drivers (one row per module x present gene)
  gene_df <- reactive({ gd <- res()$gene; req(gd, input$patient); PANEL <- assets()$ref$PANEL
    s <- input$patient; if (!(s %in% colnames(gd$dev))) return(NULL)
    dv <- gd$dev[, s]; onv <- gd$on[, s]
    rows <- do.call(rbind, lapply(SM_ORD, function(m){
      g <- intersect(toupper(PANEL[[m]]), rownames(gd$dev)); if (!length(g)) return(NULL)
      data.frame(module=m, module_ord=match(m,SM_ORD), gene=g, dev=as.numeric(dv[g]),
                 on=as.integer(onv[g]), stringsAsFactors=FALSE) }))
    rows$layer <- unname(SM_LAYER[rows$module]); rows$sign <- unname(SM_SIGN[rows$module])
    rows$lane  <- unname(LANE_LAB[rows$module]); rows$module_name <- unname(SM_LAB[rows$module])
    rows$contrib <- ifelse(rows$sign=="+", ifelse(rows$on==1,"supports this layer","below reference"),
                           ifelse(rows$on==1,"suppresses this layer","not suppressing"))
    rows$tt <- sprintf("<b>%s</b><br>%s (%s)<br>vs reference: %+.2f<br>call: %s<br>%s",
        rows$gene, rows$module_name, ifelse(rows$sign=="+","activating","suppressive"),
        rows$dev, ifelse(rows$on==1,"on","off"), rows$contrib)
    rows })

  output$dl_template <- downloadHandler(
    filename = function() "ICI_input_template.csv",
    content = function(f){ pg <- assets()$panel_genes
      ex <- tryCatch(read.csv("example_cohort.csv", row.names=1, check.names=FALSE), error=function(e) NULL)
      if (!is.null(ex)) { rn<-toupper(rownames(ex)); g<-intersect(pg, rn)
        m <- ex[match(g, rn), seq_len(min(3,ncol(ex))), drop=FALSE]
        out <- data.frame(gene=g, round(m,3), check.names=FALSE); names(out)[-1] <- paste0("Sample", seq_len(ncol(out)-1))
      } else out <- data.frame(gene=pg, Sample1=NA, Sample2=NA, Sample3=NA)
      write.csv(out, f, row.names=FALSE) })
  output$reflab <- renderText({ sprintf("Reference: %s", assets()$ref_label) })
  output$qc <- renderUI({ req(res()); cov<-res()$coverage; np<-res()$n_panel; pct<-round(100*cov/np)
    col <- if (pct>=90) "#1b5e20" else if (pct>=75) "#8a4b00" else "#8a0017"
    tagList(div(style=sprintf("font-size:0.8rem;color:%s;font-weight:600",col), sprintf("Coverage %d%%  (%d/%d panel genes)",pct,cov,np)),
      if ((res()$mapped %||% 0) > 0) div(style="font-size:0.72rem;color:#000000", sprintf("%d gene(s) matched via alias / Ensembl mapping.", res()$mapped)),
      if (res()$logged) div(style="font-size:0.72rem;color:#000000","Input looked linear; log2 applied."),
      if (pct<75) div(style="font-size:0.72rem;color:#8a0017","Low coverage; interpret with caution.")) })

  # compact "input quality" strip shown above the patient result
  output$qc_report <- renderUI({ req(res()); r<-res(); cov<-r$coverage; np<-r$n_panel; pct<-round(100*cov/np)
    lc <- r$layer_cov
    pill <- function(label, val, col="#1b5e20")
      div(style="padding:4px 11px;border:1px solid #e3dccb;border-radius:14px;background:#faf7f1;display:inline-flex;gap:6px;align-items:center;font-size:0.8rem",
          tags$span(style=sprintf("color:%s;font-weight:700",col), val), tags$span(style="color:#555", label))
    layer_pill <- function(nm,key){ if (is.null(lc) || !(key %in% colnames(lc))) return(NULL)
      have<-lc["have",key]; tot<-lc["total",key]
      col <- if (have>=0.9*tot) "#1b5e20" else if (have>=0.8*tot) "#8a4b00" else "#8a0017"
      pill(nm, sprintf("%d/%d", have, tot), col) }
    covcol <- if (pct>=90) "#1b5e20" else if (pct>=75) "#8a4b00" else "#8a0017"
    card(class="mb-2", card_body(style="padding:9px 14px",
      div(style="display:flex;gap:9px;flex-wrap:wrap;align-items:center",
        tags$span(style="font-weight:700;color:#000;margin-right:2px", "Input quality"),
        pill("panel coverage", sprintf("%d/%d", cov, np), covcol),
        layer_pill("Sensory","Sensory"), layer_pill("Integration","Integration"), layer_pill("Execution","Execution"),
        if ((r$mapped %||% 0)>0) pill("identifiers mapped", sprintf("+%d", r$mapped), "#555"),
        if (isTRUE(r$logged)) pill("input", "log2 applied", "#555"),
        pill("reference", sprintf("%s locked", assets()$version), "#6A2C5A"),
        if (pct<75) tags$span(style="color:#8a0017;font-size:0.8rem;font-weight:600;margin-left:4px",
                              "Low coverage — interpret with caution.")))) })

  # ---- Patient report ----
  output$rep_floor_badge <- renderUI({ req(patient()); p<-patient(); col<-unname(pal()$arch[p$archetype])
    lab <- if (p$floor_layer=="none") "All three layers intact" else sprintf("Lowest-scoring layer: %s ★", p$floor_layer)
    span(lab, style=sprintf("display:inline-block;padding:4px 13px;border-radius:999px;font-weight:700;font-size:0.9rem;color:%s;background:%s", ink(col), col)) })
  output$rep_hdr <- renderText({ req(patient()); st<-studylab()
    if (nzchar(st)) sprintf("Patient %s  ·  %s", patient()$sample, st) else sprintf("Patient %s", patient()$sample) })
  output$rep_readout <- renderUI({ req(patient()); p<-patient()
    txt <- sprintf("Sensory %+.2f  ·  Integration %+.2f  ·  Execution %+.2f", p$Sensory, p$Integration, p$Execution)
    div(style="font-size:0.84rem;color:#000000;margin-top:4px", tags$b("Layer scores: "), txt,
        span(style="color:#000000", "  (standardized vs reference; 0 = mean)")) })
  output$rep_call <- renderUI({ req(patient()); a<-patient()$archetype
    div(div(style=sprintf("font-size:1.2rem;font-weight:800;color:%s", onwhite(unname(pal()$arch[a]))), CALL[[a]]$t)) })
  output$rep_meaning <- renderUI({ req(patient()); div(style="font-size:0.92rem;color:#000000", CALL[[patient()$archetype]]$s) })
  output$rep_bio <- renderUI({ req(patient()); a<-patient()$archetype; col<-unname(pal()$arch[a])
    tagList(
      div(style="font-size:0.92rem;color:#000000", tags$b(BIO_LEAD[[a]])),
      div(style="font-size:0.78rem;color:#000000;margin-top:6px",
          if (a=="intact") "If this tumor is non-responsive, check:" else "Biological programs for investigation (hypothesis-generating, not treatment recommendations):"),
      chips(BIO_CHIPS[[a]], col)) })
  output$rep_stamp <- renderUI({ req(res()); a<-assets(); st<-studylab()
    pre <- if (nzchar(st)) sprintf("Study: %s  ·  ", st) else ""
    div(style="font-size:0.74rem;color:#000000",
        sprintf("%sGenerated %s  ·  Scored with %s (%d genes)  ·  reference %s  ·  within-sample standardization  ·  coverage %d/%d  ·  reference built %s",
                pre, Sys.Date(), a$version, a$n_genes, a$ref_label, res()$coverage, res()$n_panel, as.character(a$built))) })
  output$gauge <- renderPlot({ req(patient()); build_gauge(patient(), pal()$lay) }, res=108)
  output$drill <- renderPlot({ req(psub()); build_drilldown(psub(), patient()$floor_layer, pal()$lay) }, res=108)
  output$drill_text <- renderUI({ req(psub()); div(style="font-size:0.88rem;color:#000000;margin-top:6px",
      drill_text(psub(), patient()$floor_layer)) })
  output$dl_pdf <- downloadHandler(
    filename = function(){ base <- gsub("[^A-Za-z0-9]", "", patient()$sample); sl <- study_slug()
      if (nzchar(sl)) sprintf("ICI_report_%s_%s.pdf", sl, base) else sprintf("ICI_report_%s.pdf", base) },
    content = function(f){ req(patient()); ggsave(f, build_report(patient(), assets(), res()$coverage, res()$n_panel, pal()$arch, pal()$lay, psub(), studylab()),
                                                  width = 8.27, height = if (is.null(psub())) 9.6 else 11.69, device = "pdf") })

  # ---- Gene drivers (in the patient report) ----
  output$gene_plot <- renderPlot({ df <- gene_df(); req(df); v <- input$gene_layer %||% "tree"
    if (v=="heat") build_layer_heatmaps(df)
    else build_cascade_tree(df, patient(), psub(), pal()$lay) }, res = 108)

  # ---- Cohort explorer ----
  # Primer: renders at startup (suspendWhenHidden=FALSE) to trigger the async
  # Plotly.js load before the user reaches the cohort tab, so its plots never race it.
  output$`.plotly_primer` <- renderPlotly({
    plotly::plot_ly(x=0, y=0, type="scatter", mode="markers") |> plotly::config(displayModeBar=FALSE)
  })
  outputOptions(output, ".plotly_primer", suspendWhenHidden=FALSE)
  output$donut <- renderPlotly({ req(res()); s<-res()$scores; arch<-pal()$arch
    tab <- as.data.frame(table(factor(s$archetype, names(arch)))); names(tab)<-c("a","n"); tab<-tab[tab$n>0,]
    plot_ly(tab, labels=~a, values=~n, type="pie", hole=0.58, sort=FALSE,
            marker=list(colors=unname(arch[as.character(tab$a)]), line=list(color="#ffffff", width=2)),
            textinfo="label+percent", textfont=list(color=vapply(unname(arch[as.character(tab$a)]), ink, character(1)), size=16),
            hovertemplate="%{label}: %{value} patients (%{percent})<extra></extra>") |>
      layout(showlegend=FALSE, margin=list(t=12,b=12,l=12,r=12),
             font=list(color="black", family="Inter, system-ui, Arial, sans-serif")) |> config(displayModeBar=FALSE) })
  output$dl_cohort <- downloadHandler(
    filename = function(){ sl <- study_slug(); if (nzchar(sl)) sprintf("ICI_cohort_report_%s_%s.pdf", sl, Sys.Date()) else sprintf("ICI_cohort_report_%s.pdf", Sys.Date()) },
    content = function(f){ req(res()); ggsave(f, build_cohort_report(res()$scores, assets(), res()$coverage, res()$n_panel, pal()$arch, pal()$lay, studylab()),
                                              width = 10, height = 8, device = "pdf") })

  long <- reactive({ req(res()); s<-res()$scores
    L <- do.call(rbind, lapply(c("Sensory","Integration","Execution"), function(l)
      data.frame(sample=s$sample, layer=l, z=s[[l]], archetype=s$archetype, floor=s$floor_layer)))
    L$layer <- factor(L$layer, c("Sensory","Integration","Execution"))
    L$tt <- sprintf("<b>%s</b><br>S %+.2f · I %+.2f · E %+.2f<br>lowest: %s · %s",
      L$sample, s$Sensory[match(L$sample,s$sample)], s$Integration[match(L$sample,s$sample)],
      s$Execution[match(L$sample,s$sample)], L$floor, L$archetype)
    set.seed(1); L$jit <- runif(nrow(L), -0.26, 0.26); L })
  output$n_pat  <- renderText({ req(res()); res()$n_samples })
  output$cov    <- renderText({ req(res()); sprintf("%d / %d genes", res()$coverage, res()$n_panel) })
  output$intact <- renderText({ req(res()); s<-res()$scores; sprintf("%d of %d (%.0f%%)", sum(s$archetype=="intact"), nrow(s), 100*mean(s$archetype=="intact")) })
  # ggplotly() measures text on the active graphics device; measure on a throwaway
  # null PDF device so text sizing is stable regardless of the display device.
  base_theme <- function(p){
    grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add=TRUE)
    ggplotly(p, tooltip="text") |>
      layout(font=list(family="Inter, system-ui, Arial, sans-serif", color="black"), margin=list(t=24,l=88,b=48,r=16)) |> config(displayModeBar=FALSE)
  }
  output$cohort <- renderPlotly({ base_theme(build_layer_violin_jit(long(), pal()$lay)) })
  output$patients <- renderPlotly({ req(res()); L<-long(); arch<-pal()$arch
    L$yj <- as.integer(L$layer) + L$jit
    plot_ly(L, x=~z, y=~yj, type="scatter", mode="markers", source="beeswarm",
            color=~archetype, colors=arch, customdata=~sample, text=~tt, hoverinfo="text",
            marker=list(size=11, line=list(width=1.4, color="#ffffff"))) |>
      layout(shapes=list(list(type="line", x0=0, x1=0, y0=0.5, y1=3.5, line=list(color="black", width=1))),
             xaxis=list(title=list(text="standardized score vs reference (0 = mean)", font=list(size=14)), zeroline=FALSE),
             yaxis=list(tickvals=c(1,2,3), ticktext=c("<b>Sensory</b>","<b>Integration</b>","<b>Execution</b>"),
                        tickfont=list(size=15, color="black"), title="", range=c(0.4,3.6)),
             legend=list(orientation="h", y=1.17, x=0.5, xanchor="center", font=list(size=15, color="black"), itemsizing="constant"),
             font=list(color="black", family="Inter, system-ui, Arial, sans-serif"),
             margin=list(t=44,l=100,b=48,r=16)) |>
      config(displayModeBar=FALSE) |> event_register("plotly_click") })
  output$tbl <- renderDT({ req(res()); s<-res()$scores; arch<-pal()$arch
    d <- data.frame(Patient=s$sample, Sensory=round(s$Sensory,2), Integration=round(s$Integration,2),
      Execution=round(s$Execution,2), Floor=s$floor_layer, Archetype=s$archetype)
    datatable(d, rownames=FALSE, selection="single", options=list(pageLength=8, dom="tp"), class="compact stripe row-border hover") |>
      formatStyle("Archetype", backgroundColor=styleEqual(names(arch), paste0(arch,"33"))) })
  output$dl <- downloadHandler(
    filename = function(){ ext <- switch(input$dl_format %||% "CSV", CSV="csv", TSV="tsv", Excel="xlsx"); sl <- study_slug()
                           if (nzchar(sl)) sprintf("ICI_scores_%s_%s.%s", sl, Sys.Date(), ext) else sprintf("ICI_scores_%s.%s", Sys.Date(), ext) },
    content = function(f){ req(res()); s<-res()$scores; fmt<-input$dl_format %||% "CSV"; st<-studylab()
      if (nzchar(st)) s <- cbind(study=st, s)
      if (fmt=="Excel") openxlsx::write.xlsx(s, f)
      else if (fmt=="TSV") write.table(s, f, sep="\t", row.names=FALSE, quote=FALSE)
      else write.csv(s, f, row.names=FALSE) })

  # ---- click-through: beeswarm point or table row -> that patient's report ----
  go_to_patient <- function(s){ req(s); updateSelectInput(session, "patient", selected = s); bslib::nav_select("nav", "patient") }
  observeEvent(event_data("plotly_click", source="beeswarm"), {
    e <- event_data("plotly_click", source="beeswarm"); if (!is.null(e$customdata)) go_to_patient(e$customdata) })
  observeEvent(input$tbl_rows_selected, { r <- input$tbl_rows_selected; req(r, res())
    go_to_patient(res()$scores$sample[r]) })
}

shinyApp(ui, server)
