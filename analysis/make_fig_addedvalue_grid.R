# =============================================================================
# Added-value as 2x2 grids, one per established marker (GEP, IFN-gamma, CYT, PD-L1 mRNA).
# Each validated patient is split at the within-cohort median on the marker (columns)
# and on the ICI Integration score (rows); each cell = ICB response rate (%) + n,
# fill = rate. Arrows: plum (up the marker-high column) = ICI effect; sand (across the
# ICI-low row) = marker effect. Reading is explained in the caption / results, not on
# the figure. Reads pancancer_validation.rds. House style; black text only; editable SVG.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(grid) })
source("ici_theme.R"); .pt <- 2.845276
R <- readRDS("pancancer_validation.rds")
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
MK <- c(GEP="T-cell-inflamed GEP", IFNG="IFN-γ signature",
        CYT="Cytolytic activity", CD274="PD-L1 mRNA (CD274)")

hi <- function(v) ifelse(v > median(v), "High", "Low")
D <- do.call(rbind, lapply(WORKS, function(k){ s<-R[[k]]$scores
  data.frame(y=s$y, ICI=hi(s$integ_single),
             GEP=hi(s$GEP), IFNG=hi(s$IFNG), CYT=hi(s$CYT), CD274=hi(s$CD274)) }))
rate <- function(sub){ v<-D$y[sub]; c(pct=100*mean(v), n=length(v)) }

## four cells per marker -------------------------------------------------------
grid <- expand.grid(m=names(MK), GEP=c("Low","High"), ICI=c("Low","High"), stringsAsFactors=FALSE)
AG <- do.call(rbind, Map(function(m, ms, is){
  v <- rate(D[[m]]==ms & D$ICI==is)
  data.frame(marker=MK[m], cx=ifelse(ms=="Low",1,2), cy=ifelse(is=="Low",1,2),
             pct=v["pct"], n=v["n"], row.names=NULL)
}, grid$m, grid$GEP, grid$ICI))
AG$marker <- factor(AG$marker, levels=unname(MK))

p <- ggplot(AG, aes(cx, cy)) +
  geom_tile(aes(fill=pct), colour="white", linewidth=2.0, width=0.96, height=0.96) +
  geom_text(aes(label=sprintf("%.0f%%", pct), colour=pct>32), vjust=-0.15,
            family=PUB_FONT, size=(PT_BASE+5)/.pt, fontface=2) +
  geom_text(aes(label=sprintf("n=%d", n), colour=pct>32), vjust=3.0,
            family=PUB_FONT, size=(PT_BASE-2)/.pt) +
  ## plum arrow: up the marker-high column (ICI effect); sand arrow: across ICI-low row (marker effect)
  annotate("segment", x=2.58, xend=2.58, y=1.2, yend=1.8,
           arrow=arrow(length=unit(1.9,"mm"), type="closed"), linewidth=0.8, colour=PLUM) +
  annotate("segment", x=1.2, xend=1.8, y=0.42, yend=0.42,
           arrow=arrow(length=unit(1.7,"mm"), type="closed"), linewidth=0.7, colour=SAND_DEEP) +
  facet_wrap(~marker, nrow=2) +
  scale_fill_gradient(low="#F3EAF0", high=PLUM, limits=c(12,50), oob=scales::squish, guide="none") +
  scale_colour_manual(values=c(`TRUE`="white",`FALSE`=INK), guide="none") +
  scale_x_continuous(breaks=c(1,2), labels=c("Low","High"),
                     limits=c(0.45,2.75), expand=expansion(add=0)) +
  scale_y_continuous("ICI Integration score", breaks=c(1,2), labels=c("Low","High"),
                     limits=c(0.3,2.55), expand=expansion(add=0)) +
  coord_fixed(ratio=1, clip="off") +
  theme_publication() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=PT_AXIS_TITLE+1, colour=INK, face=2),
        axis.text=element_text(size=PT_BASE+1, colour=INK),
        axis.ticks=element_blank(), panel.border=element_blank(), panel.grid=element_blank(),
        strip.text=element_text(size=PT_BASE+2, face=2, colour=INK),
        strip.background=element_blank(),
        panel.spacing=unit(6,"mm"),
        plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(3,4,3,3,"mm"))

save_fig(p, "ICI_Fig_addedvalue_grid", width_mm=150, height_mm=150, dir="../results/figures/composite")
cat("done — 2x2 grids for", paste(names(MK), collapse=", "), "\n")
