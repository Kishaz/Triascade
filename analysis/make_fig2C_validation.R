# =============================================================================
# Figure 2 Panel C - the ICI framework in action on a real randomized ICB trial
# (I-SPY2).  Story (user narrative 2026-11-02): a treated cohort yields observed
# R/NR; the framework (i) TELLS R / NR APART (layer-score discrimination, ROC/AUC
# in the ICB arm) and (ii) CLUSTERS patients into archetypes that EXPLAIN why the
# non-responders failed (limiting layer).  Performance block therefore shows:
#   - ROC/AUC (ICB arm, cohort scoring)                 = tells R/NR apart
#   - ICB-ARM non-responder composition by limiting layer = why they failed
# NR composition uses the ICB(pembro)-arm denominator (the framework explains ICB
# failure, not chemo failure).  The Integration x treatment interaction (ICB-
# specificity) is a one-line caption.  Numbers from fig2C_ISPY2_perf.rds and
# nr_taxonomy_legacy98_ISPY2.rds; never typed.  Legacy 98 only.
# House style via ici_theme.R; all text BLACK except layer-colored marks; no em-dashes.
#
# Builds C_layers (a list) so it draws standalone OR sourced into the combined
# Figure 2 at an x/y offset (FIG2C_XOFF / FIG2C_YOFF; FIG2C_ASLAYERS to suppress save).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R")
.pt <- 2.845276
col <- pal_layer
FSCALE <- 1.22
P   <- readRDS("fig2C_ISPY2_perf.rds")
NRT <- readRDS("nr_taxonomy_legacy98_ISPY2.rds")

CXO <- if (exists("FIG2C_XOFF")) FIG2C_XOFF else 0    # x-offset (combined layout)
CYO <- if (exists("FIG2C_YOFF")) FIG2C_YOFF else 0    # y-offset (combined layout)

txt <- function(x,y,lab,size=9,face=1,colour=INK,hjust=0.5,parse=FALSE,angle=0)
  annotate("text", x=x+CXO, y=y+CYO, label=lab, family=PUB_FONT, size=size*FSCALE/.pt, fontface=face, colour=colour, hjust=hjust, parse=parse, angle=angle)
rect <- function(xmin,xmax,ymin,ymax,fill,colour,lwd=0.6)
  annotate("rect", xmin=xmin+CXO,xmax=xmax+CXO,ymin=ymin+CYO,ymax=ymax+CYO,fill=fill,colour=colour,linewidth=lwd)
arrw <- function(x,xend,y,yend=y,colour=INK,lwd=0.7,len=1.5)
  annotate("segment", x=x+CXO,xend=xend+CXO,y=y+CYO,yend=yend+CYO,colour=colour,linewidth=lwd,arrow=arrow(length=unit(len,"mm"),type="closed"))
person <- function(cx,cy,fill,s=1) list(
  annotate("polygon", x=c(cx-1.6*s,cx-1.3*s,cx-0.55*s,cx+0.55*s,cx+1.3*s,cx+1.6*s)+CXO,
           y=c(cy-2.0*s,cy+0.4*s,cy+0.95*s,cy+0.95*s,cy+0.4*s,cy-2.0*s)+CYO, fill=fill, colour=fill, linewidth=0),
  annotate("point", x=cx+CXO, y=cy+1.8*s+CYO, size=1.7*s, colour=fill))
cluster <- function(x0,y0,nx,ny,dx,dy,fill,s=1){ L<-list()
  for(i in 0:(nx-1)) for(j in 0:(ny-1)) L<-c(L, person(x0+i*dx, y0-j*dy, fill, s)); L }

# ---- ROC square (treated arm): tells R / NR apart ----
rx0<-168; rx1<-200; ry0<-34; ry1<-70
mapx<-function(f) rx0+f*(rx1-rx0)+CXO; mapy<-function(t) ry0+t*(ry1-ry0)+CYO
roc <- P$roc; roc$px<-mapx(roc$fpr); roc$py<-mapy(roc$tpr)
roc$layer <- factor(roc$layer, c("Sensory","Integration","Execution"))
roccol <- c(Sensory=unname(col["Sense"]), Integration=unname(col["Integration"]), Execution=unname(col["Execution"]))
roc_layers <- lapply(levels(roc$layer), function(L){ d<-roc[roc$layer==L,]; geom_path(data=d, aes(px,py), colour=roccol[L], linewidth=0.9) })

# ---- ICB-arm NR composition bars (why ICB non-responders fail; legacy-98 taxonomy) ----
nrp <- NRT$nr_comp_pembro
ord <- c("Sense-low (cold)","Integration-suppressed","Execution-low","Cascade-intact")
nrp <- nrp[match(ord, as.character(nrp$class)), ]
nnrp <- sum(nrp$n_NR)
compcol <- c(unname(col["Sense"]), unname(col["Integration"]), unname(col["Execution"]), NEUTRAL)
complab <- c("cold","blocked","no-kill","none")
cbx <- c(240, 250, 260, 270); cbw <- 7; cbb <- 34; cbsc <- 36/max(nrp$n_NR)  # tallest bar matches ROC height (y34-70)
comp_layers <- c(list(
  annotate("segment", x=236+CXO, xend=236+CXO, y=cbb+CYO, yend=cbb+36+CYO, colour=INK, linewidth=0.4),
  annotate("segment", x=236+CXO, xend=273+CXO, y=cbb+CYO, yend=cbb+CYO, colour=INK, linewidth=0.4)),
  unlist(lapply(1:4, function(i) list(
    rect(cbx[i]-cbw/2, cbx[i]+cbw/2, cbb, cbb+nrp$n_NR[i]*cbsc, compcol[i], INK, 0.4),
    txt(cbx[i], cbb+nrp$n_NR[i]*cbsc+2, sprintf("%d", nrp$n_NR[i]), size=7, face=2,
        colour=if(i==4) INK else compcol[i]),
    txt(cbx[i], cbb-3, complab[i], size=6, colour=INK))), recursive=FALSE))

C_layers <- list(
  annotate("text", x=1+CXO, y=88+CYO, label="C", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0),
  txt(9, 88, sprintf("Independent-cohort test on a randomized ICB trial: I-SPY2 (n = %d)", P$n_pembro + P$n_control), size=10, face=2, hjust=0),

  # ===== COHORT (branch point) =====
  txt(15, 74, "cohort", size=9, face=2, hjust=0.5),
  cluster(9, 69, 4, 2, 4, 4, NEUTRAL, 0.95),
  # branch 1: treat -> observed outcome
  arrw(24, 48, 66), txt(36, 69, "treat with ICB", size=8.3, hjust=0.5),
  # branch 2: drop down to the transcriptome (label rotated, runs along the shaft)
  arrw(15, 15, 61, 45), txt(19, 53, "profile", size=8, hjust=0.5, face=3, angle=90),

  # ===== OBSERVED OUTCOME box (top row) =====
  rect(50, 140, 53, 79, "white", MUTED, 0.5),
  txt(95, 74, "Observed outcome", size=9.5, face=2),
  cluster(68, 69, 3, 2, 4.2, 4, unname(col["Integration"]), 0.95),
  txt(72, 60, "Responders (R)", size=8, colour=unname(col["Integration"])),
  txt(72, 56, sprintf("pCR %.0f%%", P$pcr_pembro), size=8),
  cluster(111, 69, 4, 2, 4.2, 4, MUTED, 0.95),
  txt(117, 60, "Non-responders (NR)", size=8),

  # ===== ICI PIPELINE (bottom row): scores -> archetype -> R/NR + why =====
  rect(2, 28, 33, 43, CREAM, MUTED), txt(15, 40, "bulk", size=8.3), txt(15, 36, "transcriptome", size=7.6),
  arrw(28, 38, 38),
  rect(38, 66, 32, 44, "white", MUTED, 0.5), txt(52, 46.5, "layer scores", size=8, face=3),
  txt(41.5, 41, "S", size=8.5, face=2, colour=unname(col["Sense"])),
  txt(41.5, 38, "I", size=8.5, face=2, colour=unname(col["Integration"])),
  txt(41.5, 35, "E", size=8.5, face=2, colour=unname(col["Execution"])),
  annotate("segment", x=45+CXO, xend=63+CXO, y=41+CYO, yend=41+CYO, colour=unname(col["Sense"]), linewidth=0.5),
  annotate("segment", x=45+CXO, xend=63+CXO, y=38+CYO, yend=38+CYO, colour=unname(col["Integration"]), linewidth=0.5),
  annotate("segment", x=45+CXO, xend=63+CXO, y=35+CYO, yend=35+CYO, colour=unname(col["Execution"]), linewidth=0.5),
  annotate("point", x=57+CXO, y=41+CYO, size=1.6, colour=unname(col["Sense"])),
  annotate("point", x=52+CXO, y=38+CYO, size=1.6, colour=unname(col["Integration"])),
  annotate("point", x=60+CXO, y=35+CYO, size=1.6, colour=unname(col["Execution"])),
  arrw(66, 76, 38),
  rect(76, 102, 32, 44, "white", INK, 0.5), txt(89, 46.5, "archetype", size=8, face=3),
  rect(80,83,37,40, unname(col["Sense"]), unname(col["Sense"]),0.4),
  rect(84,87,37,40, unname(col["Integration"]), unname(col["Integration"]),0.4),
  rect(88,91,37,40, "white", unname(col["Execution"]),0.4),
  txt(96, 38.5, "type", size=8),
  arrw(102, 112, 38),
  rect(112, 140, 32, 44, "white", unname(col["Integration"]), 0.7),
  txt(126, 41, "R / NR ranking", size=7.5, face=2, colour=INK),
  txt(126, 36, "+ limiting layer", size=7, face=3, colour=INK),

  # ===== observed + predicted -> one Performance block =====
  arrw(140, 162, 66), txt(151, 69, "observed", size=8, hjust=0.5),
  arrw(140, 162, 38), txt(151, 41, "predicted", size=8, hjust=0.5),

  # ---- LEFT message: tells R / NR apart (ROC/AUC, ICB arm) ----
  txt(184, 85, "Tells R / NR apart", size=9, face=2, hjust=0.5),
  annotate("segment", x=rx0+CXO, xend=rx0+CXO, y=ry0+CYO, yend=ry1+CYO, colour=INK, linewidth=0.6),
  annotate("segment", x=rx0+CXO, xend=rx1+CXO, y=ry0+CYO, yend=ry0+CYO, colour=INK, linewidth=0.6),
  annotate("segment", x=rx0+CXO, xend=rx1+CXO, y=ry0+CYO, yend=ry1+CYO, colour="grey70", linewidth=0.4, linetype="22"),
  roc_layers,
  txt((rx0+rx1)/2, ry1+3, "ROC, ICB arm", size=8.3, face=2),
  txt((rx0+rx1)/2, ry0-4, "false positive rate", size=8),
  txt(rx0-4, (ry0+ry1)/2, "true positive rate", size=8, angle=90),
  txt(202, 60, "AUC, ICB arm", size=7.2, face=2, hjust=0),
  txt(202, 55.5, "(cohort scoring)", size=6, face=3, hjust=0),
  txt(202, 50.5, sprintf("Integration  %.2f", P$auc["Integration"]), size=7.4, face=2, hjust=0, colour=unname(col["Integration"])),
  txt(202, 46, sprintf("Sensory  %.2f", P$auc["Sensory"]), size=7.2, hjust=0, colour=unname(col["Sense"])),
  txt(202, 41.5, sprintf("Execution  %.2f", P$auc["Execution"]), size=7.2, hjust=0, colour=unname(col["Execution"])),
  txt(202, 37, sprintf("Int × treatment p = %.3f", P$int_p), size=6, face=3, hjust=0),

  # ---- RIGHT message: why non-responders fail (ICB-arm limiting layer) ----
  txt(255, 85, "Why NRs fail", size=9, face=2, hjust=0.5),
  txt(255, 80.5, sprintf("ICB arm, n = %d", nnrp), size=6.5, face=3, hjust=0.5),
  comp_layers
)

if (exists("FIG2C_ASLAYERS") && isTRUE(FIG2C_ASLAYERS)) {
  # combined mode: C_layers is exposed for make_fig2_combined.R; do not save here.
  invisible(NULL)
} else {
  p <- ggplot() + C_layers +
    coord_fixed(ratio=1, xlim=c(0,276), ylim=c(24,91), expand=FALSE, clip="off") +
    theme_void() + theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(1,1,1,1,"mm"))
  FIG_W_DOUBLE <- 320
  save_fig(p, "ICI_Fig2C", width_mm=276*174/186, height_mm=67*174/186)
  cat("done\n")
}
