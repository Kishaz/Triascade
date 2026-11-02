###############################################################################
# Figure 1 — ICI three-layer framework and translational interpretation
# Cancer Cell–style multi-panel concept figure
#
# DESIGN
#   a. Fixed a priori tumor-immune architecture (SEE → REACH → KILL)
#   b. Deterministic patient scoring pipeline (no model training)
#   c. Pan-cancer layer performance and dominant bottlenecks
#   d. Translational interpretation: select, combine, stratify
#
# OUTPUTS
#   Figure1_ICI_framework_CancerCell_V2.pdf   (vector; primary manuscript file)
#   Figure1_ICI_framework_CancerCell_V2.svg   (editable vector)
#   Figure1_ICI_framework_CancerCell_V2.tiff  (600 dpi)
#
# The figure intentionally has NO final full-width slogan/banner.
###############################################################################

SCRIPT_VERSION <- "Figure1_ICI_framework  v3.4  2026-10-28"
# v3.4: panel d p-values by likelihood-ratio test, as in the text and Figure 3
#       (I-SPY2 Integration x arm 0.018; lowest layer 0.065). Nothing else changed.
# v3.3: the two NSCLC studies shown as separate rows (13 cohort-by-arm analyses);
#       Composite read from Oct 16/source_data/percohort_AUC_13.csv.
# v3.2: Composite AUC column added to panel c, read from your per-cohort AUC table
#       (Oct 16/source_data/SuppFigS2_source_FINAL.csv); S/I/E values checked against it.
# v3.1: layout fixes only (text overflow, arrow through bar, unreadable cell text,
#       empty-box glyph) and a Unicode-safe PDF device. No values changed.
# v3.0 changes from CancerCell_V2:
#   c: I-SPY2 pembro at n=248 values; I-SPY2 paclitaxel-only and Friedrich rows added
#      (all cohort-by-arm analyses; 13 from v3.3); "Pooled by indication" removed; Riaz labelled by drug
#   b: clipped "patient score" text inside the ICI circle removed
#   d: Select/Combine/Stratify replaced by predictive / prognostic / descriptive summary;
#      failure-mode rows marked as hypotheses
#   output to Oct 16/figures
cat("\n", SCRIPT_VERSION, "\n", sep = "")
suppressPackageStartupMessages(library(grid))

# -----------------------------------------------------------------------------#
# Output settings
# -----------------------------------------------------------------------------#
OUTDIR <- if (dir.exists("/Users/dbedi")) "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript/Oct 16/figures" else file.path(getwd(), "figures")
dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

FIG_W_MM <- 180
FIG_H_MM <- 145

# Use a base-PDF-safe built-in font.
# Helvetica is supported by R's standard pdf() device on macOS.
FONT <- if (dir.exists("/Users")) "Helvetica" else "Liberation Sans"

# -----------------------------------------------------------------------------#
# Palette
# -----------------------------------------------------------------------------#
COL <- list(
  blue       = "#2E6DA4",
  blue_bg    = "#EDF3F9",
  amber      = "#D9770E",
  amber_bg   = "#FDF3E6",
  red        = "#C0392B",
  red_bg     = "#FBEDEC",
  navy       = "#1F3B57",
  salmon     = "#E8A79C",
  salmon_dk  = "#C56F61",
  tcell      = "#7FB3D9",
  tcell_dk   = "#3979A8",
  stroma     = "#C8A97E",
  dc         = "#B49BD1",
  exhausted  = "#D6DEE4",
  grey       = "#6E7781",
  grey_lt    = "#F3F4F5",
  black      = "#111111",
  white      = "#FFFFFF"
)

# -----------------------------------------------------------------------------#
# Data embedded exactly as supplied
# -----------------------------------------------------------------------------#
cohort_auc <- data.frame(
  cohort = c(
    "Urothelial — IMvigor210 (n=298)",
    "Urothelial — Rose BLCA (n=87)",
    "Breast — I-SPY2 pembro (n=69)",
    "Breast — I-SPY2 paclitaxel only (n=177)",
    "Melanoma — Gide pembro (n=32)",
    "Melanoma — Gide ipi+pembro (n=25)",
    "Melanoma — Hugo (n=26)",
    "Esophageal — EAC (n=32)",
    "Melanoma — Riaz nivolumab (n=49)",
    "NSCLC — GSE126044 nivolumab (n=16)",
    "NSCLC — GSE135222 anti-PD-(L)1 (n=27)",
    "Myeloma — Friedrich (n=9)",
    "Renal — Braun ccRCC (n=181)"
  ),
  Sensory = c(.553,.590,.754,.616,.910,.861,.642,.800,.536,.945,.618,.722,.532),
  Integration = c(.616,.692,.777,.574,.890,.649,.624,.588,.581,.655,.693,.722,.496),
  Execution = c(.568,.545,.687,.673,.635,.590,.544,.513,.645,.982,.643,.778,.557),
  Dominant = c(
    "Integration","Integration","Integration","Execution",
    "Sensory","Sensory","Sensory","Sensory",
    "Execution","Execution","Integration","Execution","None"
  ),
  stringsAsFactors = FALSE
)

# Composite AUC per cohort-by-arm analysis, read from your per-cohort AUC table
# (Supp Fig S2 source: I-SPY2 at n=248, Friedrich in). The S/I/E values above are
# checked against the same table.
S2_CSV <- Sys.getenv("ICI_S2_CSV",
  "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript/Oct 16/source_data/percohort_AUC_13.csv")
if (!file.exists(S2_CSV)) stop("Per-cohort AUC table not found: ", S2_CSV)
s2 <- read.csv(S2_CSV, stringsAsFactors = FALSE, check.names = FALSE)
s2$key <- tolower(trimws(s2$cohort))
cohort_auc$key <- c("imvigor", "rose", "ispy2 [pembro]", "ispy2 [control]", "gide [pembro]",
                    "gide [ipi+pembro]", "hugo", "eac", "riaz", "nsclc gse126044", "nsclc gse135222",
                    "friedrich", "braun")
get_auc <- function(k, L) {
  v <- s2$AUC[s2$key == k & tolower(s2$layer) == tolower(L)]
  if (length(v)) v[1] else NA_real_
}
cohort_auc$Composite <- sapply(cohort_auc$key, get_auc, L = "Composite")
chk <- data.frame(
  cohort    = cohort_auc$key,
  d_S       = round(cohort_auc$Sensory     - sapply(cohort_auc$key, get_auc, L = "Sensory"), 3),
  d_I       = round(cohort_auc$Integration - sapply(cohort_auc$key, get_auc, L = "Integration"), 3),
  d_E       = round(cohort_auc$Execution   - sapply(cohort_auc$key, get_auc, L = "Execution"), 3),
  Composite = round(cohort_auc$Composite, 3))
cat("\nComposite AUC from", basename(S2_CSV), "(d_ = figure value minus table value):\n")
print(chk, row.names = FALSE)
if (any(is.na(cohort_auc$Composite)))
  stop("Composite AUC missing for: ", paste(cohort_auc$key[is.na(cohort_auc$Composite)], collapse = ", "))
if (any(abs(unlist(chk[, c("d_S", "d_I", "d_E")])) > 0.0015, na.rm = TRUE))
  cat("NOTE: a layer value in the figure differs from the table - paste this output to Claude\n")

gene_text <- list(
  sensory = paste0(
    "21 genes: HLA-A/B/C/E/F/G, B2M, TAP1/2, TAPBP,\n",
    "PSMB8/9/10, IRF1, STAT1/2, JAK1/2, IFNG, IFNGR1/2"
  ),
  integration = paste0(
    "41 genes, incl. CXCL9/10/11/13, CCL4/5/19/21, ICAM1,\n",
    "CD80/86, CD40 (recruit) vs VEGFA/B/C, TGFB1/2/3,\n",
    "COL1A1, FAP, MMP2/9 (exclude)"
  ),
  execution = paste0(
    "38 genes, incl. GZMA/B/H/K/M, PRF1, GNLY, NKG7,\n",
    "CD8A/B (kill) vs PDCD1, CTLA4, LAG3, HAVCR2,\n",
    "TIGIT, TOX (suppress)"
  )
)

# -----------------------------------------------------------------------------#
# Low-level drawing helpers
# -----------------------------------------------------------------------------#
gp_txt <- function(size = 7, col = COL$black, face = "plain", lineheight = 1.0) {
  gpar(fontfamily = FONT, fontsize = size, col = col,
       fontface = face, lineheight = lineheight)
}

wrap_text <- function(x, width) {
  paste(strwrap(x, width = width), collapse = "\n")
}

panel_vp <- function(x, y, w, h) {
  viewport(x = x, y = y, width = w, height = h,
           just = c("left","bottom"), clip = "on")
}

panel_border <- function(col = "#C8CDD2") {
  grid.rect(gp = gpar(fill = NA, col = col, lwd = 0.8))
}

panel_label <- function(letter) {
  grid.text(letter, x = unit(0.012, "npc"), y = unit(0.982, "npc"),
            just = c("left","top"), gp = gp_txt(12, COL$black, "bold"))
}

panel_title <- function(title, x = 0.065, y = 0.976, size = 8.2) {
  grid.text(title, x = unit(x, "npc"), y = unit(y, "npc"),
            just = c("left","top"), gp = gp_txt(size, COL$black, "bold"))
}

draw_arrow <- function(x0, y0, x1, y1, col = COL$black, lwd = 1.1,
                       dashed = FALSE) {
  grid.lines(
    x = unit(c(x0,x1), "npc"),
    y = unit(c(y0,y1), "npc"),
    gp = gpar(col = col, lwd = lwd, lty = if (dashed) 2 else 1),
    arrow = arrow(length = unit(1.8, "mm"), type = "closed")
  )
}

draw_x <- function(x, y, size = 0.025, col = COL$red, lwd = 2.5) {
  grid.lines(unit(c(x-size,x+size),"npc"), unit(c(y-size,y+size),"npc"),
             gp = gpar(col = col, lwd = lwd))
  grid.lines(unit(c(x-size,x+size),"npc"), unit(c(y+size,y-size),"npc"),
             gp = gpar(col = col, lwd = lwd))
}

draw_tumor_cell <- function(x, y, r = 0.035, fill = COL$salmon,
                            outline = COL$salmon_dk, dashed = FALSE) {
  grid.circle(x = unit(x,"npc"), y = unit(y,"npc"), r = unit(r,"npc"),
              gp = gpar(fill = fill, col = outline, lwd = 1.0,
                        lty = if (dashed) 2 else 1))
  if (!is.na(fill)) {
    grid.circle(x = unit(x,"npc"), y = unit(y,"npc"), r = unit(r*0.42,"npc"),
                gp = gpar(fill = "#CD7F70", col = NA))
  }
}

draw_tcell <- function(x, y, r = 0.032, exhausted = FALSE) {
  fc <- if (exhausted) COL$exhausted else COL$tcell
  ec <- if (exhausted) COL$grey else COL$tcell_dk
  grid.circle(unit(x,"npc"), unit(y,"npc"), r = unit(r,"npc"),
              gp = gpar(fill = fc, col = ec, lwd = 1))
  grid.circle(unit(x,"npc"), unit(y,"npc"), r = unit(r*0.48,"npc"),
              gp = gpar(fill = if (exhausted) "#BFC8CE" else "#4E89B6",
                        col = NA))
}

draw_dc <- function(x, y, r = 0.028) {
  ang <- seq(0, 2*pi, length.out = 17)[-17]
  rr <- rep(c(r*1.55, r), 8)
  xx <- x + rr*cos(ang)
  yy <- y + rr*sin(ang)
  grid.polygon(unit(xx,"npc"), unit(yy,"npc"),
               gp = gpar(fill = COL$dc, col = "#765C9D", lwd = 0.9))
  grid.circle(unit(x,"npc"), unit(y,"npc"), r = unit(r*0.42,"npc"),
              gp = gpar(fill = "#8068A7", col = NA))
}

draw_mhci <- function(x, y, scale = 1) {
  grid.lines(unit(c(x,x),"npc"), unit(c(y,y+0.018*scale),"npc"),
             gp = gpar(col = COL$blue, lwd = 1.6))
  grid.lines(unit(c(x-0.009*scale,x+0.009*scale),"npc"),
             unit(c(y+0.018*scale,y+0.018*scale),"npc"),
             gp = gpar(col = COL$blue, lwd = 3.0, lineend = "round"))
}

draw_vessel <- function(x, y, w, h) {
  grid.rect(unit(x,"npc"), unit(y,"npc"), unit(w,"npc"), unit(h,"npc"),
            just = c("left","bottom"),
            gp = gpar(fill = "#F3B2A7", col = COL$salmon_dk, lwd = 1))
  grid.lines(unit(c(x+0.12*w,x+0.12*w),"npc"),
             unit(c(y,y+h),"npc"), gp = gpar(col = "#B85E53", lwd = 1))
  grid.lines(unit(c(x+0.88*w,x+0.88*w),"npc"),
             unit(c(y,y+h),"npc"), gp = gpar(col = "#B85E53", lwd = 1))
}

draw_collagen <- function(x, y, w, h, n = 9) {
  for (i in seq_len(n)) {
    xx <- x + (i-0.5)*w/n
    yy <- seq(y, y+h, length.out = 80)
    wav <- 0.012*sin(seq(0, 5*pi, length.out = 80) + i*0.7)
    grid.lines(unit(xx + wav,"npc"), unit(yy,"npc"),
               gp = gpar(col = COL$stroma, lwd = 1.05))
  }
}

draw_inhibitory_bar <- function(x0, y0, x1, y1, col = COL$red) {
  grid.lines(unit(c(x0,x1),"npc"), unit(c(y0,y1),"npc"),
             gp = gpar(col = col, lwd = 2.0))
  grid.lines(unit(c(x1,x1),"npc"), unit(c(y1-0.018,y1+0.018),"npc"),
             gp = gpar(col = col, lwd = 2.0))
}

mix_with_white <- function(hex, amount) {
  rgb0 <- grDevices::col2rgb(hex)/255
  amount <- max(0, min(1, amount))
  out <- rgb0*amount + 1*(1-amount)
  grDevices::rgb(out[1], out[2], out[3])
}

txt_col <- function(fill) {
  v <- grDevices::col2rgb(fill)[, 1] / 255
  if (sum(c(0.2126, 0.7152, 0.0722) * v) < 0.55) COL$white else COL$black
}

auc_fill <- function(v, layer_col) {
  strength <- 0.10 + 0.78*max(0, min(1, (v - 0.48)/(0.92 - 0.48)))
  mix_with_white(layer_col, strength)
}

# -----------------------------------------------------------------------------#
# Panel a — biological architecture
# -----------------------------------------------------------------------------#
draw_panel_a <- function() {
  panel_border()
  panel_label("a")
  panel_title("Fixed a priori tumor–immune response architecture", size = 7.4)

  # Flow line
  draw_arrow(0.318, 0.67, 0.352, 0.67, lwd = 1.0)
  draw_arrow(0.648, 0.67, 0.682, 0.67, lwd = 1.0)
  grid.text("then", unit(0.335,"npc"), unit(0.705,"npc"),
            gp = gp_txt(5.0, COL$grey, "italic"))
  grid.text("then", unit(0.665,"npc"), unit(0.705,"npc"),
            gp = gp_txt(5.0, COL$grey, "italic"))

  # Three layer regions, no large card boxes
  xs <- c(0.17, 0.50, 0.83)
  layer_cols <- c(COL$blue, COL$amber, COL$red)
  layer_bgs  <- c(COL$blue_bg, COL$amber_bg, COL$red_bg)
  titles <- c("SENSORY (SEE)", "INTEGRATION (REACH)", "EXECUTION (KILL)")
  questions <- c(
    "Can the immune system\nsee the tumor?",
    "Can activated T cells\nreach the tumor?",
    "Can T cells kill despite\ncheckpoint suppression?"
  )

  for (i in 1:3) {
    grid.roundrect(unit(xs[i]-0.145,"npc"), unit(0.20,"npc"),
                   unit(0.29,"npc"), unit(0.69,"npc"),
                   r = unit(1.5,"mm"),
                   just = c("left","bottom"),
                   gp = gpar(fill = layer_bgs[i], col = layer_cols[i], lwd = 0.9))
    grid.rect(unit(xs[i]-0.145,"npc"), unit(0.82,"npc"),
              unit(0.29,"npc"), unit(0.07,"npc"),
              just = c("left","bottom"),
              gp = gpar(fill = layer_cols[i], col = layer_cols[i]))
    grid.text(titles[i], unit(xs[i],"npc"), unit(0.855,"npc"),
              gp = gp_txt(7.1, COL$white, "bold"))
    grid.text(questions[i], unit(xs[i],"npc"), unit(0.785,"npc"),
              gp = gp_txt(4.85, layer_cols[i], "bold", lineheight = 0.95))
  }

  # Sensory cartoon
  draw_tumor_cell(0.125,0.63,0.034); draw_tumor_cell(0.175,0.66,0.036)
  draw_tumor_cell(0.215,0.61,0.034); draw_tumor_cell(0.165,0.57,0.034)
  draw_mhci(0.115,0.664); draw_mhci(0.17,0.697)
  draw_mhci(0.217,0.645); draw_mhci(0.158,0.603)
  draw_dc(0.075,0.48,0.022)
  draw_arrow(0.095,0.50,0.137,0.555, col = COL$black, lwd = 0.8)
  draw_tumor_cell(0.252,0.55,0.029, fill = NA, outline = COL$grey, dashed = TRUE)
  draw_x(0.252,0.55,0.018)
  grid.text("no MHC-I", unit(0.252,"npc"), unit(0.505,"npc"),
            gp = gp_txt(5.2, COL$red, "italic"))
  grid.text("antigen presentation • MHC-I • IFN-γ",
            unit(0.17,"npc"), unit(0.435,"npc"),
            gp = gp_txt(5.1, COL$black))
  grid.text(gene_text$sensory, unit(0.035,"npc"), unit(0.355,"npc"),
            just = c("left","top"), gp = gp_txt(3.3, COL$black, lineheight = 0.90))
  grid.roundrect(unit(0.035,"npc"), unit(0.215,"npc"),
                 unit(0.270,"npc"), unit(0.058,"npc"),
                 r = unit(1.2,"mm"), just = c("left","bottom"),
                 gp = gpar(fill = COL$white, col = COL$blue, lwd = 0.8))
  grid.text("Score_S = Σ 1[gene ≥ cohort median]",
            unit(0.17,"npc"), unit(0.244,"npc"),
            gp = gp_txt(4.0, COL$blue, "bold"))

  # Integration cartoon
  draw_vessel(0.385,0.47,0.055,0.22)
  draw_tcell(0.412,0.62,0.020); draw_tcell(0.412,0.54,0.020)
  draw_tcell(0.455,0.58,0.020)
  draw_arrow(0.43,0.58,0.475,0.58, lwd = 0.8)
  draw_collagen(0.49,0.46,0.075,0.25,8)
  # chemokine gradient
  gradx <- seq(0.61,0.46,length.out=12)
  grady <- 0.68 - 0.05*(seq_along(gradx)-1)/11
  for (k in seq_along(gradx)) {
    grid.circle(unit(gradx[k],"npc"), unit(grady[k],"npc"),
                r = unit(0.0045*(1-k/20),"npc"),
                gp = gpar(fill = COL$amber, col = NA))
  }
  grid.text("CXCL9/10/11, CCL5 — attract",
            unit(0.50,"npc"), unit(0.735,"npc"),
            gp = gp_txt(4.7, COL$amber, "bold"))
  draw_tcell(0.60,0.56,0.023)
  draw_x(0.60,0.56,0.017)
  grid.text("TGF-β, VEGF, FAP, COL1A1 — exclude",
            unit(0.50,"npc"), unit(0.425,"npc"),
            gp = gp_txt(4.4, "#8A6137"))
  grid.text(gene_text$integration, unit(0.365,"npc"), unit(0.355,"npc"),
            just = c("left","top"), gp = gp_txt(3.2, COL$black, lineheight = 0.88))
  grid.roundrect(unit(0.365,"npc"), unit(0.215,"npc"),
                 unit(0.270,"npc"), unit(0.058,"npc"),
                 r = unit(1.2,"mm"), just = c("left","bottom"),
                 gp = gpar(fill = COL$white, col = COL$amber, lwd = 0.8))
  grid.text("Score_I = Σ IntPos (21) − Σ IntNeg (20)",
            unit(0.50,"npc"), unit(0.244,"npc"),
            gp = gp_txt(4.0, COL$amber, "bold"))

  # Execution cartoon
  draw_tumor_cell(0.765,0.60,0.045)
  draw_tcell(0.895,0.62,0.034)
  # granules
  for (k in 1:12) {
    xx <- 0.85 - 0.006*k
    yy <- 0.62 + 0.018*sin(k*1.5)
    grid.circle(unit(xx,"npc"), unit(yy,"npc"), r = unit(0.0035,"npc"),
                gp = gpar(fill = COL$red, col = NA))
  }
  grid.text("granzymes • perforin", unit(0.84,"npc"), unit(0.695,"npc"),
            gp = gp_txt(5.2, COL$red, "bold"))
  draw_inhibitory_bar(0.805,0.54,0.875,0.54)
  grid.text("PD-1 / PD-L1", unit(0.84,"npc"), unit(0.505,"npc"),
            gp = gp_txt(5.1, COL$red))
  draw_tcell(0.89,0.44,0.026, exhausted = TRUE)
  grid.text("TOX, LAG3 — exhausted", unit(0.83,"npc"), unit(0.40,"npc"),
            gp = gp_txt(4.9, COL$grey))
  draw_x(0.93,0.44,0.016)
  grid.text(gene_text$execution, unit(0.695,"npc"), unit(0.355,"npc"),
            just = c("left","top"), gp = gp_txt(3.2, COL$black, lineheight = 0.90))
  grid.roundrect(unit(0.695,"npc"), unit(0.215,"npc"),
                 unit(0.270,"npc"), unit(0.058,"npc"),
                 r = unit(1.2,"mm"), just = c("left","bottom"),
                 gp = gpar(fill = COL$white, col = COL$red, lwd = 0.8))
  grid.text("Score_E = Σ ExecPos (19) − Σ ExecNeg (19)",
            unit(0.83,"npc"), unit(0.244,"npc"),
            gp = gp_txt(4.0, COL$red, "bold"))

  grid.text("Failure at any step can limit response; the rate-limiting layer varies by tumor context.",
            unit(0.5,"npc"), unit(0.135,"npc"),
            gp = gp_txt(6.0, COL$navy, "bold"))
  grid.text("Red × marks denote the matched failure state at each layer.",
            unit(0.5,"npc"), unit(0.092,"npc"),
            gp = gp_txt(5.3, COL$grey, "italic"))
}

# -----------------------------------------------------------------------------#
# Panel b — deterministic scoring pipeline
# -----------------------------------------------------------------------------#
draw_matrix_icon <- function(cx, cy, w = 0.11, h = 0.12, binary = FALSE) {
  nr <- 7; nc <- 5
  for (r in seq_len(nr)) for (c in seq_len(nc)) {
    xx <- cx - w/2 + (c-0.5)*w/nc
    yy <- cy - h/2 + (r-0.5)*h/nr
    if (binary) {
      val <- ((r+c) %% 3) != 0
      fc <- if (val) COL$navy else COL$white
    } else {
      pal <- c(COL$blue_bg, COL$blue, COL$amber_bg, COL$amber,
               COL$red_bg, COL$red)
      fc <- pal[((r*2+c) %% length(pal))+1]
    }
    grid.rect(unit(xx,"npc"), unit(yy,"npc"),
              unit(w/nc*0.84,"npc"), unit(h/nr*0.82,"npc"),
              gp = gpar(fill = fc, col = "#D7DADF", lwd = 0.25))
  }
}

draw_panel_b <- function() {
  panel_border()
  panel_label("b")
  panel_title("Deterministic scoring pipeline", size = 7.2)

  # Workflow nodes
  ypos <- 0.67
  xpos <- c(0.11,0.30,0.485,0.68,0.87)
  labs <- c("Input\ntranscriptome", "Median\nbinarization",
            "Layer\nscores", "Within-cohort\nz-scores", "Composite\nICI")

  draw_matrix_icon(xpos[1], ypos, 0.13, 0.18, FALSE)
  draw_matrix_icon(xpos[2], ypos, 0.13, 0.18, TRUE)

  # Three layer scores
  grid.circle(unit(xpos[3],"npc"),unit(ypos+0.07,"npc"),r=unit(0.035,"npc"),
              gp=gpar(fill=COL$blue_bg,col=COL$blue,lwd=1))
  grid.circle(unit(xpos[3],"npc"),unit(ypos,"npc"),r=unit(0.035,"npc"),
              gp=gpar(fill=COL$amber_bg,col=COL$amber,lwd=1))
  grid.circle(unit(xpos[3],"npc"),unit(ypos-0.07,"npc"),r=unit(0.035,"npc"),
              gp=gpar(fill=COL$red_bg,col=COL$red,lwd=1))
  grid.text("S",unit(xpos[3],"npc"),unit(ypos+0.07,"npc"),gp=gp_txt(7,COL$blue,"bold"))
  grid.text("I",unit(xpos[3],"npc"),unit(ypos,"npc"),gp=gp_txt(7,COL$amber,"bold"))
  grid.text("E",unit(xpos[3],"npc"),unit(ypos-0.07,"npc"),gp=gp_txt(7,COL$red,"bold"))

  # z bars
  zvals <- c(0.06,0.10,0.04)
  zcols <- c(COL$blue,COL$amber,COL$red)
  for (i in 1:3) {
    x0 <- xpos[4]-0.058+(i-1)*0.046
    grid.rect(unit(x0,"npc"), unit(ypos-0.08,"npc"),
              unit(0.030,"npc"), unit(0.16,"npc"),
              just=c("left","bottom"),
              gp=gpar(fill=COL$grey_lt,col="#D0D4D8",lwd=.4))
    grid.rect(unit(x0,"npc"), unit(ypos-0.08,"npc"),
              unit(0.030,"npc"), unit(0.06+zvals[i],"npc"),
              just=c("left","bottom"),
              gp=gpar(fill=zcols[i],col=NA))
  }

  # Composite output
  grid.circle(unit(xpos[5],"npc"),unit(ypos,"npc"),r=unit(0.060,"npc"),
              gp=gpar(fill=COL$navy,col=COL$navy))
  grid.text("ICI",unit(xpos[5],"npc"),unit(ypos,"npc"),
            gp=gp_txt(11,COL$white,"bold"))

  for (i in 1:4) draw_arrow(xpos[i]+0.072,ypos,xpos[i+1]-0.072,ypos,lwd=0.9)

  for (i in 1:5) {
    grid.text(labs[i],unit(xpos[i],"npc"),unit(0.475,"npc"),
              gp=gp_txt(4.55,COL$black,"bold", lineheight = 0.95))
  }

  # Formula and design notes
  grid.roundrect(unit(0.06,"npc"),unit(0.26,"npc"),
                 unit(0.88,"npc"),unit(0.15,"npc"),
                 r=unit(1.5,"mm"),just=c("left","bottom"),
                 gp=gpar(fill="#F7F8FA",col="#B8BEC5",lwd=.7))
  grid.text("Composite ICI = mean of within-cohort z-scores",
            unit(0.50,"npc"),unit(0.355,"npc"),
            gp=gp_txt(5.7,COL$navy,"bold"))
  grid.text("(Score_S, Score_I, Score_E)",
            unit(0.50,"npc"),unit(0.315,"npc"),
            gp=gp_txt(5.2,COL$navy,"bold"))

  grid.text("98 unique genes  •  100 layer assignments",
            unit(0.08,"npc"),unit(0.18,"npc"),just=c("left","center"),
            gp=gp_txt(5.6,COL$black,"bold"))
  grid.text("IFNG and TGFB1 are assigned to two biological layers",
            unit(0.08,"npc"),unit(0.135,"npc"),just=c("left","center"),
            gp=gp_txt(5.0,COL$grey))
  grid.text("Fixed a priori  •  no training  •  no fitted parameters",
            unit(0.08,"npc"),unit(0.085,"npc"),just=c("left","center"),
            gp=gp_txt(5.6,COL$navy,"bold"))
}

# -----------------------------------------------------------------------------#
# Panel c — pan-cancer AUC map
# -----------------------------------------------------------------------------#
draw_panel_c <- function() {
  panel_border()
  panel_label("c")
  panel_title("Layer performance in each cohort-by-arm analysis", size = 7.0)

  left <- 0.025
  rowlab_w <- 0.36
  cell_w <- 0.115
  xcols <- left + rowlab_w + c(0.5,1.5,2.5)*cell_w
  xcomp <- left + rowlab_w + 3.5*cell_w + 0.025
  xdom <- xcomp + cell_w + 0.005

  grid.text("Cohort", unit(left,"npc"),unit(0.875,"npc"),
            just=c("left","center"),gp=gp_txt(5.6,COL$black,"bold"))
  grid.text("Sensory",unit(xcols[1],"npc"),unit(0.875,"npc"),
            gp=gp_txt(5.5,COL$blue,"bold"))
  grid.text("Integration",unit(xcols[2],"npc"),unit(0.875,"npc"),
            gp=gp_txt(5.5,COL$amber,"bold"))
  grid.text("Execution",unit(xcols[3],"npc"),unit(0.875,"npc"),
            gp=gp_txt(5.5,COL$red,"bold"))
  grid.text("Composite",unit(xcomp,"npc"),unit(0.875,"npc"),
            gp=gp_txt(5.5,COL$navy,"bold"))
  grid.text("Dominant",unit(xdom,"npc"),unit(0.875,"npc"),
            gp=gp_txt(5.3,COL$black,"bold"))
  grid.lines(unit(c(left,0.98),"npc"),unit(c(0.85,0.85),"npc"),
             gp=gpar(col="#8C939A",lwd=.65))

  n <- nrow(cohort_auc)
  ytop <- 0.80
  dy <- 0.0555
  layers <- c("Sensory","Integration","Execution")
  layer_cols <- c(COL$blue,COL$amber,COL$red)

  for (r in seq_len(n)) {
    y <- ytop - (r-1)*dy
    if (r %% 2 == 0) {
      grid.rect(unit(left,"npc"),unit(y-dy*0.45,"npc"),
                unit(0.955,"npc"),unit(dy*0.9,"npc"),
                just=c("left","bottom"),
                gp=gpar(fill="#FAFAFA",col=NA))
    }
    grid.text(cohort_auc$cohort[r],unit(left,"npc"),unit(y,"npc"),
              just=c("left","center"),gp=gp_txt(4.2,COL$black))
    for (j in 1:3) {
      val <- cohort_auc[r,layers[j]]
      is_dom <- cohort_auc$Dominant[r] == layers[j]
      grid.rect(unit(xcols[j],"npc"),unit(y,"npc"),
                unit(cell_w*0.86,"npc"),unit(dy*0.72,"npc"),
                gp=gpar(fill=auc_fill(val,layer_cols[j]),
                        col=if(is_dom) layer_cols[j] else "#D3D7DB",
                        lwd=if(is_dom) 1.2 else .35))
      grid.text(sprintf("%.3f",val),unit(xcols[j],"npc"),unit(y,"npc"),
                gp=gp_txt(4.25,txt_col(auc_fill(val,layer_cols[j])),
                          if(is_dom) "bold" else "plain"))
    }
    vc <- cohort_auc$Composite[r]
    grid.rect(unit(xcomp,"npc"),unit(y,"npc"),
              unit(cell_w*0.86,"npc"),unit(dy*0.72,"npc"),
              gp=gpar(fill=auc_fill(vc,COL$navy),col="#D3D7DB",lwd=.35))
    grid.text(sprintf("%.3f",vc),unit(xcomp,"npc"),unit(y,"npc"),
              gp=gp_txt(4.25,txt_col(auc_fill(vc,COL$navy))))
    dom_col <- switch(cohort_auc$Dominant[r],
                      Sensory=COL$blue,Integration=COL$amber,
                      Execution=COL$red,None=COL$grey)
    grid.text(cohort_auc$Dominant[r],unit(xdom,"npc"),unit(y,"npc"),
              gp=gp_txt(4.2,dom_col,
                        if(cohort_auc$Dominant[r]!="None") "bold" else "plain"))
  }

  grid.text("AUC for response within each cohort-by-arm analysis; highest layer outlined. Composite = mean of the three z-scored layers. Cohorts are not pooled.",
            unit(left,"npc"),unit(0.020,"npc"),just=c("left","bottom"),
            gp=gp_txt(4.0,COL$grey,"italic"))
}

# -----------------------------------------------------------------------------#
# Panel d — translational interpretation
# -----------------------------------------------------------------------------#
draw_failure_row <- function(y, col, bg, title, tumors, strategy) {
  grid.rect(unit(0.035,"npc"),unit(y-0.075,"npc"),
            unit(0.93,"npc"),unit(0.145,"npc"),
            just=c("left","bottom"),
            gp=gpar(fill=bg,col="#D5D8DC",lwd=.55))
  grid.rect(unit(0.035,"npc"),unit(y-0.075,"npc"),
            unit(0.012,"npc"),unit(0.145,"npc"),
            just=c("left","bottom"),
            gp=gpar(fill=col,col=col))
  grid.text(title,unit(0.07,"npc"),unit(y+0.043,"npc"),
            just=c("left","center"),gp=gp_txt(5.7,col,"bold"))
  grid.text(wrap_text(tumors, 60),unit(0.07,"npc"),unit(y+0.005,"npc"),
            just=c("left","center"),gp=gp_txt(4.25,COL$black, lineheight = 0.95))
  grid.text(wrap_text(strategy, 66),unit(0.07,"npc"),unit(y-0.043,"npc"),
            just=c("left","center"),gp=gp_txt(4.05,col,"italic", lineheight = 0.92))
}

draw_panel_d <- function() {
  panel_border()
  panel_label("d")
  panel_title("Interpretation", size = 7.2)

  draw_failure_row(
    0.81, COL$blue, COL$blue_bg,
    "SENSORY-BLIND",
    "Melanoma (Gide, Hugo) • Esophageal",
    "Restore antigen presentation: IFN-γ agonism, HDAC/DNMT modulation,\nvaccines or oncolytic virus"
  )
  draw_failure_row(
    0.63, COL$amber, COL$amber_bg,
    "INTEGRATION-BLOCKED",
    "Urothelial • Breast (pembrolizumab arm) • NSCLC (GSE135222)",
    "Remodel stroma/vasculature: anti-VEGF, TGF-β blockade,\nCXCR4 or CSF1R inhibition, FAP-directed agents"
  )
  draw_failure_row(
    0.45, COL$red, COL$red_bg,
    "EXECUTION-FAILED",
    "NSCLC (GSE126044) • Melanoma (Riaz)",
    "Release or bypass the brake: PD-1/PD-L1, LAG3, TIGIT,\nbispecific T-cell engagers"
  )

  # Negative control
  grid.rect(unit(0.035,"npc"),unit(0.245,"npc"),
            unit(0.93,"npc"),unit(0.120,"npc"),
            just=c("left","bottom"),
            gp=gpar(fill=COL$grey_lt,col="#AEB4BA",lwd=.65))
  grid.text("NO CASCADE SIGNAL",unit(0.065,"npc"),unit(0.338,"npc"),
            just=c("left","center"),gp=gp_txt(5.7,COL$grey,"bold"))
  grid.text("Clear-cell renal carcinoma",unit(0.065,"npc"),unit(0.303,"npc"),
            just=c("left","center"),gp=gp_txt(4.8,COL$black))
  grid.text("Angiogenesis-driven context; all layers near-null (AUC 0.50–0.56)",
            unit(0.065,"npc"),unit(0.268,"npc"),
            just=c("left","center"),gp=gp_txt(4.45,COL$grey,"italic"))

  # What the data establish
  x0 <- c(0.18,0.50,0.82)
  heads <- c("Predictive","Prognostic","Descriptive")
  body <- c(
    "Integration × arm\nI-SPY2 p = 0.018\nNeoTRIP p = 0.040",
    "Sensory, Execution\nassociated with response\nin both arms",
    "Lowest layer did not\nstratify benefit\n(p = 0.065)"
  )
  hcols <- c(COL$amber,COL$navy,COL$grey)
  for (i in 1:3) {
    grid.text(heads[i],unit(x0[i],"npc"),unit(0.19,"npc"),
              gp=gp_txt(6.0,hcols[i],"bold"))
    grid.text(body[i],unit(x0[i],"npc"),unit(0.12,"npc"),
              gp=gp_txt(4.4,COL$black, lineheight = 0.95))
  }
  grid.lines(unit(c(0.345,0.345),"npc"),unit(c(0.07,0.21),"npc"),
             gp=gpar(col="#C6CACF",lwd=.7))
  grid.lines(unit(c(0.655,0.655),"npc"),unit(c(0.07,0.21),"npc"),
             gp=gpar(col="#C6CACF",lwd=.7))

  grid.text(
    "Layer failure modes above are hypotheses, not treatment assignments;\nnot a companion diagnostic and not a basis for withholding standard care.",
    unit(0.50,"npc"),unit(0.026,"npc"),
    gp=gp_txt(3.75,COL$grey,"italic", lineheight = 0.95)
  )
}

# -----------------------------------------------------------------------------#
# Full figure
# -----------------------------------------------------------------------------#
draw_figure <- function() {
  grid.newpage()
  grid.rect(gp = gpar(fill = COL$white, col = NA))

  # Multi-panel layout deliberately differs from a three-column slide:
  # top: biological architecture + scoring
  # bottom: pan-cancer data + translational interpretation
  pushViewport(panel_vp(0.012,0.520,0.615,0.465))
  draw_panel_a()
  popViewport()

  pushViewport(panel_vp(0.64,0.520,0.348,0.465))
  draw_panel_b()
  popViewport()

  pushViewport(panel_vp(0.012,0.018,0.615,0.485))
  draw_panel_c()
  popViewport()

  pushViewport(panel_vp(0.64,0.018,0.348,0.485))
  draw_panel_d()
  popViewport()
}

# -----------------------------------------------------------------------------#
# Export
# -----------------------------------------------------------------------------#
pdf_file <- file.path(OUTDIR, "Figure1_concept_v3_4.pdf")
if (capabilities("cairo")) {
  cairo_pdf(pdf_file, width = FIG_W_MM/25.4, height = FIG_H_MM/25.4, family = FONT)
} else {
  pdf(pdf_file, width = FIG_W_MM/25.4, height = FIG_H_MM/25.4,
      family = FONT, useDingbats = FALSE, onefile = TRUE)
}
draw_figure()
dev.off()

svg_file <- file.path(OUTDIR, "Figure1_concept_v3_4.svg")
if (requireNamespace("svglite", quietly = TRUE)) {
  svglite::svglite(svg_file,
                   width = FIG_W_MM/25.4, height = FIG_H_MM/25.4,
                   system_fonts = list(sans = FONT))
} else {
  svg(svg_file,
      width = FIG_W_MM/25.4, height = FIG_H_MM/25.4,
      family = FONT)
}
draw_figure()
dev.off()

tiff_file <- file.path(OUTDIR, "Figure1_concept_v3_4.tiff")
if (requireNamespace("ragg", quietly = TRUE)) {
  ragg::agg_tiff(tiff_file,
                 width = FIG_W_MM, height = FIG_H_MM,
                 units = "mm", res = 600, compression = "lzw",
                 background = "white")
} else {
  tiff(tiff_file,
       width = FIG_W_MM/25.4*600,
       height = FIG_H_MM/25.4*600,
       units = "px", res = 600, compression = "lzw",
       family = FONT, bg = "white")
}
draw_figure()
dev.off()

cat("\nCancer Cell–style Figure 1 written to:\n",
    "  ", pdf_file, "\n",
    "  ", svg_file, "\n",
    "  ", tiff_file, "\n\n",
    "Primary manuscript file: vector PDF\n",
    "No final full-width slogan/banner was added.\n", sep = "")
