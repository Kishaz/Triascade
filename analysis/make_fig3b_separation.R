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
          Puch="Puch", Kim2018gas="Kim 2018", Rose="Rose 2021", ISPY2="I-SPY2 trial",
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

## ---- long per-patient frame, z-scored within cohort ------------------------
DF <- do.call(rbind, lapply(ORDER, function(k){ s <- R[[k]]$scores
  z <- (s$integ_single - mean(s$integ_single)) / sd(s$integ_single)
  resp <- factor(ifelse(s$y==1,"R","NR"), levels=c("NR","R"))
  off  <- ifelse(resp=="R", 0.205, -0.205)
  data.frame(cohort=k, x=xpos[k], resp=resp,
             px = xpos[k] + off + runif(length(z), -0.11, 0.11),
             z  = z, bucket=bucket[k], stringsAsFactors=FALSE) }))
DF$z <- pmax(pmin(DF$z, 3.1), -3.1)                     # clip rare outliers for framing

## ---- medians / IQR per cohort x group --------------------------------------
agg <- do.call(rbind, lapply(split(DF, list(DF$cohort, DF$resp), drop=TRUE), function(d){
  data.frame(cohort=d$cohort[1], x=d$x[1], resp=d$resp[1],
             xc=d$x[1]+ifelse(d$resp[1]=="R",0.205,-0.205),
             med=median(d$z), q1=quantile(d$z,.25), q3=quantile(d$z,.75)) }))

## ---- per-cohort meta (AUC, n, label, y-top) --------------------------------
## x-axis label: CANCER TYPE first, then author/trial, then n
meta <- do.call(rbind, lapply(ORDER, function(k){
  data.frame(cohort=k, x=xpos[k], auc=R[[k]]$int_ci_single[1],
             n=R[[k]]$n, cancer=R[[k]]$cancer,
             label=sprintf("%s · %s (n=%d)", cdisp(R[[k]]$cancer), AUTHOR[[k]], R[[k]]$n),
             bucket=bucket[k]) }))

RCOL <- PLUM; NRCOL <- SAND            # responder = plum, non-responder = sand (house two-group)
ytop <- 3.1; yAUC <- 3.6; yband <- c(-4.2, 4.3); yHEAD <- yband[2]+0.5; yLEG <- yband[2]+1.3

## ---- bucket background bands + headers -------------------------------------
band <- function(ks, fill, lab){ xs <- xpos[ks]
  list(annotate("rect", xmin=min(xs)-0.45, xmax=max(xs)+0.45, ymin=yband[1], ymax=yband[2],
                fill=fill, colour=NA),
       annotate("text", x=mean(range(xs)), y=yHEAD, label=lab, lineheight=0.9,
                family=PUB_FONT, size=(PT_BASE-1)/.pt, fontface=2, colour=INK, hjust=0.5, vjust=1)) }

p <- ggplot() +
  band(WORKS, "#F3EAF0", "Validated\n(anti-PD-1/PD-L1)") +
  band(CTLA4, "#F5EEDD", "CTLA-4") +
  band(NULLB, "#EFECE7", "Null\n(cold / resistant)") +
  annotate("segment", x=0.2, xend=max(xpos)+0.8, y=0, yend=0,
           linewidth=0.3, linetype="22", colour=RULE)

## patient clouds
p <- p + geom_point(data=DF, aes(px, z, colour=resp), size=0.55, alpha=0.5, stroke=0) +
  scale_colour_manual(values=c(R=RCOL, NR=NRCOL), guide="none")

## median + IQR crossbars
p <- p +
  geom_segment(data=agg, aes(x=xc-0.17, xend=xc+0.17, y=med, yend=med, colour=resp), linewidth=0.7) +
  geom_segment(data=agg, aes(x=xc, xend=xc, y=q1, yend=q3, colour=resp), linewidth=0.4)

## AUC per cohort (bold where it clears chance meaningfully)
p <- p + geom_text(data=meta, aes(x=x, y=yAUC, label=sprintf("%.2f", auc)),
                   family=PUB_FONT, size=(PT_BASE-2)/.pt,
                   fontface=ifelse(meta$auc>=0.6,2,1), colour=INK) +
  annotate("text", x=0.1, y=yAUC, label="AUC", family=PUB_FONT,
           size=(PT_BASE-2.5)/.pt, colour=MUTED, hjust=1)

## cohort labels (cancer, author, n) under the axis — vertical, house-style font
p <- p +
  geom_text(data=meta, aes(x=x, y=-3.35, label=label), family=PUB_FONT,
            size=(PT_BASE-0.5)/.pt, angle=90, hjust=1, vjust=0.5, colour=INK)

## y axis (manual, left side) + R/NR legend swatches
yt <- c(-3,-2,-1,0,1,2,3)
p <- p + annotate("segment", x=0.0, xend=0.0, y=-3.1, yend=3.1, linewidth=0.3, colour=RULE)
for(t in yt) p <- p +
  annotate("segment", x=-0.06, xend=0.0, y=t, yend=t, linewidth=0.3, colour=RULE) +
  annotate("text", x=-0.18, y=t, label=t, family=PUB_FONT, size=(PT_BASE-2.5)/.pt, colour=INK, hjust=1)
p <- p + annotate("text", x=-0.95, y=0, label="Integration score  (z within cohort)",
                  family=PUB_FONT, size=PT_AXIS_TITLE/.pt, colour=INK, angle=90, hjust=0.5)

## legend (R / NR) — own line, centred over the validated band
cxL <- mean(range(xpos[WORKS]))
p <- p +
  annotate("point", x=cxL-1.9, y=yLEG, colour=RCOL,  size=1.5) +
  annotate("text",  x=cxL-1.7, y=yLEG, label="responder",     family=PUB_FONT, size=(PT_BASE-2)/.pt, hjust=0, colour=INK) +
  annotate("point", x=cxL+0.3, y=yLEG, colour=NRCOL, size=1.5) +
  annotate("text",  x=cxL+0.5, y=yLEG, label="non-responder", family=PUB_FONT, size=(PT_BASE-2)/.pt, hjust=0, colour=INK)

## panel tag + frame
p <- p +
  annotate("text", x=-0.95, y=yLEG, label="B", family=PUB_FONT,
           size=PT_PANEL_TAG/.pt, fontface=2, hjust=0, colour=INK) +
  coord_cartesian(xlim=c(-1.05, max(xpos)+0.6), ylim=c(-10.8, yLEG+0.6), clip="off") +
  theme_void() +
  theme(plot.background=element_rect(fill="white", colour=NA),
        plot.margin=margin(2,3,2,2,"mm"))

save_fig(p, "ICI_Fig3b_separation", width_mm=174, height_mm=168, dir="../results/figures/composite")
cat("done\n")
