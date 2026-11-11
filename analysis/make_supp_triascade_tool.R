# =============================================================================
# SUPPLEMENTARY: the Triascade tool as a data-flow schematic - tumor RNA in,
# per-patient cascade read-outs out. Shapes with geometric (square-corner) black
# connectors, read left -> right: INPUT -> HARMONIZE -> SCORE (against the frozen
# locked reference) -> CALL the archetype -> three OUTPUTS. House style: Arial,
# black ink, warm plum/sand/rose palette, editable SVG (text joined, no textLength).
# No results are drawn; this is what the software DOES, not what it found.
# Drawing helpers ported from the pangenome house kit (fig_design_helpers.R) and
# recoloured to the ICI theme, per the "copy, don't depend on the path" rule.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R")
.pt <- 2.845276                              # geom_text size (mm) = pt / .pt

# ---- canvas (millimetres; 1 x-unit == 1 y-unit so corners stay round) --------
W <- 174; H <- 120; YB <- 3; ASP <- 1
PTS <- 9.5; SM <- 8; TINY <- 7; HEAD <- 10; TITLE <- 13

# ---- ported drawing helpers (house kit, ICI palette) -------------------------
L <- list(); add <- function(g) L[[length(L) + 1]] <<- g
rr <- function(x0, x1, y0, y1, r = 2, n = 10) {          # rounded-rectangle outline
  ry <- r * ASP; a <- seq(0, pi/2, length.out = n)
  data.frame(x = c(x1-r + r*cos(a), x0+r - r*sin(a), x0+r - r*cos(a), x1-r + r*sin(a)),
             y = c(y1-ry + ry*sin(a), y1-ry + ry*cos(a), y0+ry - ry*sin(a), y0+ry - ry*cos(a)))
}
dbox <- function(x0, x1, y0, y1, fill = "white", col = INK, lwd = 0.5, r = 2.2)
  add(geom_polygon(data = rr(x0, x1, y0, y1, r), aes(x, y), fill = fill, colour = col,
                   linewidth = lwd, inherit.aes = FALSE))
dlab <- function(x, y, txt, size = PTS, hj = 0.5, vj = 0.5, col = INK, face = 1, lh = 1.14)
  add(annotate("text", x = x, y = y, label = txt, size = size/.pt, hjust = hj, vjust = vj,
               colour = col, family = PUB_FONT, fontface = face, lineheight = lh))
darw <- function(x0, x1, y0, y1 = y0, col = INK, lwd = 0.6, len = 1.7)              # straight arrow
  add(annotate("segment", x = x0, xend = x1, y = y0, yend = y1, colour = col, linewidth = lwd,
               arrow = arrow(type = "closed", length = unit(len, "mm")), lineend = "butt", linejoin = "mitre"))
dpath <- function(x, y, col = INK, lwd = 0.6, ends = "last", len = 1.7)             # ONE mitred path, arrow at `ends`
  add(annotate("path", x = x, y = y, colour = col, linewidth = lwd, linejoin = "mitre", lineend = "butt",
               arrow = arrow(type = "closed", length = unit(len, "mm"), ends = ends)))
dseg <- function(x0, x1, y0, y1 = y0, col = INK, lwd = 0.6)
  add(annotate("segment", x = x0, xend = x1, y = y0, yend = y1, colour = col, linewidth = lwd, lineend = "round"))
dpt <- function(x, y, size = 2, col = INK, fill = col, shape = 21, stroke = 0.3)
  add(annotate("point", x = x, y = y, size = size, colour = "white", fill = fill, shape = shape, stroke = stroke))
ell <- function(x, y, rx, ry, fill = NA, col = INK, lwd = 0.45) {
  a <- seq(0, 2*pi, length.out = 60)
  add(geom_polygon(data = data.frame(x = x + rx*cos(a), y = y + ry*sin(a)), aes(x, y),
                   fill = fill, colour = col, linewidth = lwd, inherit.aes = FALSE))
}
# a database cylinder (the classic data-store icon): top ellipse, two sides, a
# front bottom arc, and one platter line below the top.
dbcyl <- function(xc, y0, y1, rx, fill = "white", col = INK, lwd = 0.6) {
  ry <- rx * 0.34; a <- seq(0, pi, length.out = 30)
  add(geom_polygon(data = data.frame(x = c(xc-rx, xc+rx, xc+rx, xc-rx), y = c(y0, y0, y1, y1)),
                   aes(x, y), fill = fill, colour = NA, inherit.aes = FALSE))   # body
  ell(xc, y0, rx, ry, fill = fill, col = NA, lwd = 0); ell(xc, y1, rx, ry, fill = fill, col = NA, lwd = 0)
  dseg(xc-rx, xc-rx, y0, y1, col = col, lwd = lwd); dseg(xc+rx, xc+rx, y0, y1, col = col, lwd = lwd)
  ell(xc, y1, rx, ry, fill = NA, col = col, lwd = lwd)                          # top ellipse (full)
  add(geom_path(data = data.frame(x = xc+rx*cos(a), y = y0 - ry*sin(a)), aes(x, y), colour = col, linewidth = lwd, inherit.aes = FALSE))             # front bottom arc
  add(geom_path(data = data.frame(x = xc+rx*cos(a), y = (y1 - (y1-y0)*0.3) - ry*sin(a)), aes(x, y), colour = col, linewidth = lwd*0.85, inherit.aes = FALSE))  # platter line
}
# a chip: small rounded square with a white letter (used for the three layers)
chip <- function(xc, yc, s, fill, letter) { dbox(xc-s, xc+s, yc-s, yc+s, fill = fill, col = "white", lwd = 0.4, r = 0.7)
  dlab(xc, yc, letter, TINY, col = "white", face = 2) }
# a titled box: bold title (one line) with body hanging below it, top-anchored
boxlab <- function(cx, y1, title, body, tcol = INK) {
  dlab(cx, y1-4.2, title, HEAD, col = tcol, face = 2)
  dlab(cx, y1-8.6, body, SM, vj = 1)
}

# ---- palette roles -----------------------------------------------------------
BORDER <- PLUM_DEEP                                    # box outline
FILL_STEP <- "white"; FILL_REF <- "#F4ECDB"; FILL_OUT <- CREAM_LT
ARCH <- c(intact = "#2E7D32", cold = SAND_DEEP, blocked = PLUM, "no-kill" = ROSE)  # archetype = its floor layer; intact = none

# =============================================================================
# header: the Triascade logo in place of a text title (aspect preserved)
logo <- png::readPNG("../app/www/logo.png")
LGH <- 15.5; LGW <- LGH * (dim(logo)[2] / dim(logo)[1])
add(annotation_raster(logo, xmin = 3.5, xmax = 3.5 + LGW, ymin = 118.5 - LGH, ymax = 118.5, interpolate = TRUE))

# ---- the pipeline band (4 steps, left -> right; equal widths) -----------------
PY0 <- 74; PY1 <- 102
bw <- 38.5; gp <- 4; x0 <- 4
B <- list(); for (k in 0:3) B[[c("input","harm","score","call")[k+1]]] <- x0 + k*(bw+gp) + c(0, bw)
for (nm in names(B)) dbox(B[[nm]][1], B[[nm]][2], PY0, PY1, fill = FILL_STEP, col = BORDER, lwd = 0.6)
ctr <- function(b) mean(B[[b]]); PYM <- (PY0 + PY1)/2

boxlab(ctr("input"), PY1, "Tumor RNA-seq", "gene × sample matrix\nTPM, log-TPM or counts\n1 patient = 1 column")
boxlab(ctr("harm"),  PY1, "Harmonize",     "map IDs to gene symbols\nlog2 if linear-scale\n98-gene coverage")
boxlab(ctr("score"), PY1, "Score the cascade", "self-normalize vs\nthe locked reference\ngene on / off (≥ median)")
chip(ctr("score")-9, PY0+4, 2.3, SAND_DEEP, "S"); chip(ctr("score"), PY0+4, 2.3, PLUM, "I"); chip(ctr("score")+9, PY0+4, 2.3, ROSE, "E")
boxlab(ctr("call"),  PY1, "Call the archetype", "lowest layer\nnames the type")
ax <- ctr("call") + c(-13.5, -4.5, 4.5, 13.5); for (i in seq_along(ARCH)) dpt(ax[i], PY0+7.8, 3.4, fill = ARCH[i])
for (i in seq_along(ARCH)) dlab(ax[i], PY0+3.2, names(ARCH)[i], TINY-0.8)

darw(B$input[2], B$harm[1], PYM); darw(B$harm[2], B$score[1], PYM); darw(B$score[2], B$call[1], PYM)

# ---- frozen reference (data store) feeding the SCORE step --------------------
RX <- c(ctr("score")-30, ctr("score")+30); RY0 <- 44; RY1 <- 64; RXc <- mean(RX)
dbox(RX[1], RX[2], RY0, RY1, fill = FILL_REF, col = SAND_DEEP, lwd = 0.6)
dbcyl(RX[1] + 12, RY0 + 4, RY1 - 4, 5, fill = "white", col = SAND_DEEP, lwd = 0.6)   # database icon
TXc <- (RX[1] + 20 + RX[2]) / 2                                     # centre of the text, right of the icon
dlab(TXc, RY1-4.2, "Locked reference", HEAD, col = SAND_DEEP, face = 2)
dlab(TXc, RY1-8.6, "frozen ICI reference (v1)\ndeterministic scoring\ncohort-independent (n = 1)", SM, vj = 1)
darw(ctr("score"), ctr("score"), RY1, PY0, len = 1.7)              # up into SCORE

# ---- outputs band (fan out from CALL) ---------------------------------------
OY0 <- 8; OY1 <- 34
O <- list(report = c(4, 56), cohort = c(61, 113), export = c(118, 170))
for (nm in names(O)) dbox(O[[nm]][1], O[[nm]][2], OY0, OY1, fill = FILL_OUT, col = BORDER, lwd = 0.6)
octr <- function(o) mean(O[[o]])
boxlab(octr("report"), OY1, "Patient report",  "archetype call + gauge\ncascade drill-down, 98-gene tree\none-page PDF")
boxlab(octr("cohort"), OY1, "Cohort explorer", "archetype proportions\nlayer distributions, beeswarm\nclick a point for that patient")
boxlab(octr("export"), OY1, "Exports",         "scores: CSV / TSV / Excel\nPDF reports\nversion-stamped output")

# fan: trunk down from CALL, a bus, an arrow into each output (no arrow past the last step)
BUS <- 38; trunk <- ctr("call")
dseg(trunk, trunk, PY0, BUS); dseg(octr("report"), trunk, BUS, BUS)
for (o in c("report","cohort","export")) darw(octr(o), octr(o), BUS, OY1)

# =============================================================================
p <- ggplot() + L +
  coord_cartesian(xlim = c(0, W), ylim = c(YB, H), expand = FALSE, clip = "off") +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA),
                       plot.margin = margin(1,1,1,1))
save_fig(p, "ICI_SuppFig_triascade_tool", width_mm = W, height_mm = H - YB,
         dir = "../results/figures/supplementary")
cat("done - Triascade tool schematic\n")
