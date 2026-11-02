# =============================================================================
# Figure 3b: pan-cancer ICB validation of the Integration layer, as a forest
# plot. One row per cohort = single-patient Integration AUC with DeLong 95% CI,
# grouped: pure PD-(L)1 (validated) -> pooled random-effects diamond; CTLA-4
# monotherapy; non-ICB / cold tumours (null, incl. the Braun ccRCC negative
# control). Chance line at 0.5. House style (ici_theme.R); editable SVG.
# Reads pancancer_validation.rds (run_pancancer_extra.R). Single-patient =
# the paper's thesis metric (self-normalised, one patient scorable alone).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(grid) })
source("ici_theme.R"); .pt <- 2.845276
R <- readRDS("pancancer_validation.rds")

## ---- buckets (match run_pooled_auc.R) --------------------------------------
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose")
CTLA4 <- c("VanAllen2015"); NULLB <- c("GSE165252","PRINCE","Braun")

## ---- friendly display names (no fabricated years where unsure) -------------
DISP <- c(IMvigor210="IMvigor210", GSE78220="Hugo 2016", GSE91061="Riaz 2017",
          GSE126044="Cho 2020", GSE135222="Jung 2019", Gide2019="Gide 2019",
          Liu2019="Liu 2019", Puch="Puch", Kim2018gas="Kim 2018", Rose="Rose",
          VanAllen2015="Van Allen 2015", GSE165252="EAC (PERFECT)",
          PRINCE="PRINCE", Braun="Braun 2020")

row_of <- function(k){ x<-R[[k]]; ci<-x$int_ci_single
  data.frame(key=k, name=DISP[[k]], cancer=x$cancer, n=x$n, nR=x$nR, nNR=x$nNR,
             auc=ci[1], lo=ci[2], hi=ci[3], stringsAsFactors=FALSE) }

## ---- DeLong-CI -> SE -> logit -> DerSimonian-Laird pooled (works) ----------
se_a <- function(lo,hi)(hi-lo)/(2*qnorm(0.975)); lg<-function(a)log(a/(1-a)); ilg<-function(x)1/(1+exp(-x))
dw <- do.call(rbind, lapply(WORKS, row_of))
y  <- lg(dw$auc); sy <- se_a(dw$lo,dw$hi)/(dw$auc*(1-dw$auc)); w<-1/sy^2
yf <- sum(w*y)/sum(w); Q<-sum(w*(y-yf)^2); k<-nrow(dw); C<-sum(w)-sum(w^2)/sum(w)
tau2<-max(0,(Q-(k-1))/C); wr<-1/(sy^2+tau2); yr<-sum(wr*y)/sum(wr); ser<-sqrt(1/sum(wr))
POOL <- list(auc=ilg(yr), lo=ilg(yr-1.96*ser), hi=ilg(yr+1.96*ser),
             n=sum(dw$n), nR=sum(dw$nR), nNR=sum(dw$nNR), I2=max(0,100*(Q-(k-1))/Q))

## ---- assemble rows top -> bottom, assign y ---------------------------------
dw <- dw[order(-dw$auc),]                      # best AUC on top within works
items <- list()
add <- function(type, ...) items[[length(items)+1]] <<- c(list(type=type), list(...))
add("header", text="Pure PD-(L)1  — validated")
for(i in seq_len(nrow(dw))) add("cohort", r=dw[i,], col=PLUM)
add("pooled")
add("gap"); add("header", text="CTLA-4 monotherapy")
for(k in CTLA4) add("cohort", r=row_of(k), col=SAND_DEEP)
add("gap"); add("header", text="Non-ICB / cold tumours  — null")
for(k in NULLB) add("cohort", r=row_of(k), col=MUTED)

## y-coordinates (top rows get the highest y)
yy <- length(items); ypos <- numeric(length(items))
for(i in seq_along(items)){ t<-items[[i]]$type
  if(t=="gap"){ yy <- yy-0.6; ypos[i]<-NA; next }
  if(t=="header") yy <- yy-0.4
  ypos[i]<-yy; yy<-yy-1 }

## ---- geometry / layout columns ---------------------------------------------
xL1 <- -1.18          # display name column
xL2 <- -0.50          # cancer column
xR  <-  1.03          # stats column
xlo <- 0.30; xhi <- 1.00                 # plotted AUC axis range
axY <- min(ypos,na.rm=TRUE) - 1.1        # x-axis baseline
ytop<- max(ypos,na.rm=TRUE)

fmt <- function(a,lo,hi,n) sprintf("%.2f  [%.2f–%.2f]   n=%d", a,lo,hi,n)

p <- ggplot() +
  ## chance line
  annotate("segment", x=0.5, xend=0.5, y=axY+0.4, yend=ytop+0.4,
           linetype="22", linewidth=0.3, colour=ACCENT) +
  annotate("text", x=0.5, y=ytop+1.0, label="chance", family=PUB_FONT,
           size=(PT_BASE-2.5)/.pt, colour=ACCENT, hjust=0.5)

## rows
for(i in seq_along(items)){ it<-items[[i]]; yv<-ypos[i]; if(is.na(yv)) next
  if(it$type=="header"){
    p <- p + annotate("text", x=xL1, y=yv, label=it$text, family=PUB_FONT,
                      size=(PT_BASE-1)/.pt, fontface=2, hjust=0, colour=INK)
  } else if(it$type=="cohort"){
    r<-it$r
    p <- p +
      annotate("segment", x=r$lo, xend=r$hi, y=yv, yend=yv, linewidth=0.5, colour=it$col) +
      annotate("point", x=r$auc, y=yv, size=2.2, shape=21, fill=it$col, colour="white", stroke=0.5) +
      annotate("text", x=xL1, y=yv, label=r$name, family=PUB_FONT,
               size=(PT_BASE-1.5)/.pt, hjust=0, colour=INK) +
      annotate("text", x=xL2, y=yv, label=r$cancer, family=PUB_FONT,
               size=(PT_BASE-2.5)/.pt, hjust=0, colour=MUTED) +
      annotate("text", x=xR, y=yv, label=fmt(r$auc,r$lo,r$hi,r$n), family=PUB_FONT,
               size=(PT_BASE-2.5)/.pt, hjust=0, colour=INK)
  } else if(it$type=="pooled"){
    d<-data.frame(x=c(POOL$lo,POOL$auc,POOL$hi,POOL$auc),
                  y=yv+c(0,0.32,0,-0.32))
    p <- p +
      geom_polygon(data=d, aes(x,y), fill=PLUM_DEEP, colour=PLUM_DEEP) +
      annotate("text", x=xL1, y=yv, label="Pooled (random effects)", family=PUB_FONT,
               size=(PT_BASE-1.5)/.pt, fontface=2, hjust=0, colour=PLUM_DEEP) +
      annotate("text", x=xL2, y=yv, label=sprintf("%d cohorts", nrow(dw)), family=PUB_FONT,
               size=(PT_BASE-2.5)/.pt, hjust=0, colour=MUTED) +
      annotate("text", x=xR, y=yv,
               label=sprintf("%.2f  [%.2f–%.2f]   n=%d", POOL$auc,POOL$lo,POOL$hi,POOL$n),
               family=PUB_FONT, size=(PT_BASE-2.5)/.pt, fontface=2, hjust=0, colour=PLUM_DEEP) }
}

## ---- manual x-axis (only across the plotted AUC range) ----------------------
ticks <- seq(0.3,1.0,0.1)
p <- p + annotate("segment", x=xlo, xend=xhi, y=axY, yend=axY, linewidth=0.3, colour=RULE)
for(t in ticks) p <- p +
  annotate("segment", x=t, xend=t, y=axY, yend=axY-0.18, linewidth=0.3, colour=RULE) +
  annotate("text", x=t, y=axY-0.55, label=sprintf("%.1f",t), family=PUB_FONT,
           size=(PT_BASE-2.5)/.pt, colour=INK, hjust=0.5)
p <- p + annotate("text", x=(xlo+xhi)/2, y=axY-1.5,
                  label="Single-patient Integration AUC", family=PUB_FONT,
                  size=PT_AXIS_TITLE/.pt, colour=INK, hjust=0.5)

## ---- panel tag + frame ------------------------------------------------------
p <- p +
  annotate("text", x=xL1, y=ytop+1.9, label="B", family=PUB_FONT,
           size=PT_PANEL_TAG/.pt, fontface=2, hjust=0, colour=INK) +
  coord_cartesian(xlim=c(xL1, 1.62), ylim=c(axY-2.1, ytop+2.3), clip="off") +
  theme_void() +
  theme(plot.background=element_rect(fill="white", colour=NA),
        plot.margin=margin(2,3,2,2,"mm"))

n_rows <- sum(!is.na(ypos))
save_fig(p, "ICI_Fig3b_auc", width_mm=174, height_mm=118, dir="../results/figures/composite")
cat(sprintf("done — pooled works single-patient Integration AUC %.3f [%.3f-%.3f], I2=%.0f%%\n",
            POOL$auc,POOL$lo,POOL$hi,POOL$I2))
