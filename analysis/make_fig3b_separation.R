# =============================================================================
# Figure 3b: does the single-patient Integration score separate RESPONDERS from
# NON-RESPONDERS?  Per cohort, the score (z-scored within cohort) for R vs NR,
# shown as dodged patient-level clouds + median/IQR. Grouped: pure PD-(L)1
# (validated) | CTLA-4 | non-ICB/cold (null). AUC annotated per cohort.
# The separation is visible where it works and collapses in the null block.
# Needs per-sample scores in pancancer_validation.rds (icb_eval $scores).
# House style (ici_theme.R); single canvas; editable SVG.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(grid) })
source("ici_theme.R"); .pt <- 2.845276
set.seed(7)
R <- readRDS("pancancer_validation.rds")

WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
CTLA4 <- c("VanAllen2015"); NULLB <- c("GSE165252","PRINCE","Braun")
## author/trial + year so every label identifies its study on its own
AUTHOR <- c(IMvigor210="IMvigor210", GSE78220="Hugo 2016", GSE91061="Riaz 2017",
          GSE126044="Cho 2020", GSE135222="Jung 2019", Gide2019="Gide 2019", Liu2019="Liu 2019",
          Puch="PUCH 2022", Kim2018gas="Kim 2018", Rose="Rose 2021", ISPY2="I-SPY2 trial",
          VanAllen2015="Van Allen 2015", GSE165252="PERFECT trial", PRINCE="PRINCE trial", Braun="Braun 2020")
CANC <- c(melanoma="Melanoma", urothelial="Urothelial", NSCLC="NSCLC", gastric="Gastric",
          breast="Breast", EAC="Esophageal", pancreatic="Pancreatic", ccRCC="Renal")
cdisp <- function(x){ v<-unname(CANC[x]); ifelse(is.na(v), x, v) }

## ---- order within works by single-patient Integration AUC (desc) ------------
aucs <- sapply(WORKS, function(k) R[[k]]$int_ci_single[1])
WORKS <- WORKS[order(-aucs)]
ORDER <- c(WORKS, CTLA4, NULLB)
bucket <- setNames(c(rep("works",length(WORKS)), rep("ctla4",length(CTLA4)),
                     rep("null",length(NULLB))), ORDER)

## ---- x positions with gaps between buckets ---------------------------------
xpos <- numeric(0); xc <- 1; prev <- "works"
for(k in ORDER){ if(bucket[k]!=prev){ xc <- xc + 1.0; prev <- bucket[k] }  # tighter inter-group gap
  xpos[k] <- xc; xc <- xc + 1 }

RCOL <- PLUM; NRCOL <- SAND            # responder = plum, non-responder = sand (house two-group)
ytop <- 3.1; yAUC <- 3.6; yband <- c(-4.2, 4.3); yHEAD <- yband[2]+1.25; yLEG <- yband[2]+2.1

## ---- bucket background bands + headers -------------------------------------
band <- function(ks, fill, lab){ xs <- xpos[ks]
  list(annotate("rect", xmin=min(xs)-0.45, xmax=max(xs)+0.45, ymin=yband[1], ymax=yband[2],
                fill=fill, colour=NA),
       annotate("text", x=mean(range(xs)), y=yHEAD, label=lab, lineheight=0.9,
                family=PUB_FONT, size=(PT_BASE-1)/.pt, fontface=2, colour=INK, hjust=0.5, vjust=1)) }

## =============================================================================
## build one layer's separation plot. layer = "Sensory"|"Integration"|"Execution"
## (cohort ORDER is fixed = Integration-AUC order, so panels line up for comparison)
## =============================================================================
build <- function(layer, score_col, tag, outname, outdir="../results/figures/composite"){
  set.seed(7)
  ## long per-patient frame, z-scored within cohort (this layer's score)
  DF <- do.call(rbind, lapply(ORDER, function(k){ s <- R[[k]]$scores
    sc <- s[[score_col]]; z <- (sc - mean(sc)) / sd(sc)
    resp <- factor(ifelse(s$y==1,"R","NR"), levels=c("NR","R"))
    off  <- ifelse(resp=="R", 0.205, -0.205)
    data.frame(cohort=k, x=xpos[k], resp=resp,
               px = xpos[k] + off + runif(length(z), -0.11, 0.11),
               z  = pmax(pmin(z,3.1),-3.1), bucket=bucket[k], stringsAsFactors=FALSE) }))
  DF$vx <- DF$x + ifelse(DF$resp=="R", 0.21, -0.21)
  agg <- do.call(rbind, lapply(split(DF, list(DF$cohort, DF$resp), drop=TRUE), function(d){
    data.frame(xc=d$x[1]+ifelse(d$resp[1]=="R",0.205,-0.205), resp=d$resp[1], med=median(d$z)) }))
  ## per-cohort AUC for THIS layer (single-patient selfnorm) + x-axis label
  meta <- do.call(rbind, lapply(ORDER, function(k){
    data.frame(x=xpos[k], auc=R[[k]]$auc["single_self", layer], cohort=k,
               label=sprintf("%s · %s (n=%d)", cdisp(R[[k]]$cancer), AUTHOR[[k]], R[[k]]$n)) }))

  p <- ggplot() +
    band(WORKS, "#F3EAF0", "Validated\n(anti-PD-1/PD-L1)") +
    band(CTLA4, "#F5EEDD", "CTLA-4") +
    band(NULLB, "#EFECE7", "Null\n(cold / resistant)") +
    annotate("segment", x=0.2, xend=max(xpos)+0.8, y=0, yend=0, linewidth=0.3, linetype="22", colour=RULE) +
    geom_violin(data=DF, aes(x=vx, y=z, group=interaction(cohort,resp), fill=resp, colour=resp),
                width=0.36, alpha=0.5, linewidth=0.3, trim=TRUE, scale="width", adjust=0.5) +
    geom_point(data=DF, aes(px, z, colour=resp), size=0.3, alpha=0.25, stroke=0) +
    scale_fill_manual(values=c(R=RCOL, NR=NRCOL), guide="none") +
    scale_colour_manual(values=c(R=RCOL, NR=NRCOL), guide="none") +
    geom_segment(data=agg, aes(x=xc-0.18, xend=xc+0.18, y=med, yend=med, colour=resp), linewidth=0.7)

  ## AUC per cohort (bold where it clears chance)
  p <- p + geom_text(data=meta, aes(x=x, y=yAUC, label=sprintf("%.2f", auc)),
                     family=PUB_FONT, size=(PT_BASE-2)/.pt, fontface=ifelse(meta$auc>=0.6,2,1), colour=INK) +
    annotate("text", x=0.1, y=yAUC, label="AUC", family=PUB_FONT, size=(PT_BASE-2.5)/.pt, colour=MUTED, hjust=1) +
    geom_text(data=meta, aes(x=x, y=-3.35, label=label), family=PUB_FONT,
              size=(PT_BASE-0.5)/.pt, angle=90, hjust=1, vjust=0.5, colour=INK)

  ## manual y axis + title (this layer)
  for(t in c(-3,-2,-1,0,1,2,3)) p <- p +
    annotate("segment", x=-0.06, xend=0.0, y=t, yend=t, linewidth=0.3, colour=RULE) +
    annotate("text", x=-0.18, y=t, label=t, family=PUB_FONT, size=(PT_BASE-2.5)/.pt, colour=INK, hjust=1)
  p <- p + annotate("segment", x=0.0, xend=0.0, y=-3.1, yend=3.1, linewidth=0.3, colour=RULE) +
    annotate("text", x=-0.95, y=0, label=sprintf("%s score  (z within cohort)", layer),
             family=PUB_FONT, size=PT_AXIS_TITLE/.pt, colour=INK, angle=90, hjust=0.5)

  ## legend + panel tag + frame
  cxL <- mean(range(xpos[WORKS]))
  p <- p +
    annotate("point", x=cxL-2.6, y=yLEG, colour=RCOL,  size=3.4) +
    annotate("text",  x=cxL-2.35, y=yLEG, label="responder",     family=PUB_FONT, size=PT_BASE/.pt, hjust=0, colour=INK) +
    annotate("point", x=cxL+0.5, y=yLEG, colour=NRCOL, size=3.4) +
    annotate("text",  x=cxL+0.75, y=yLEG, label="non-responder", family=PUB_FONT, size=PT_BASE/.pt, hjust=0, colour=INK) +
    annotate("text", x=-0.95, y=yLEG, label=tag, family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0, colour=INK) +
    coord_cartesian(xlim=c(-1.05, max(xpos)+0.6), ylim=c(-10.8, yLEG+0.6), clip="off") +
    theme_void() +
    theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(2,3,2,2,"mm"))

  dir.create(outdir, showWarnings=FALSE, recursive=TRUE)
  save_fig(p, outname, width_mm=174, height_mm=168, dir=outdir)
  cat(sprintf("  %s -> %s/%s\n", layer, basename(outdir), outname))
}

## Integration = the main Fig 3b (composite/). Sensory + Execution = SUPPLEMENTARY
## figures: same plot for the other two layers, showing Integration separates best.
SUPP <- "../results/figures/supplementary"
build("Integration","integ_single","B","ICI_Fig3b_separation")
build("Sensory",    "sens_single", "a","ICI_SuppFig_separation_sensory",   SUPP)
build("Execution",  "exec_single", "b","ICI_SuppFig_separation_execution", SUPP)
cat("done\n")
