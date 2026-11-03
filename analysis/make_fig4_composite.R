# =============================================================================
# Figure 4 composite: ICI Integration vs established RNA biomarkers.
#   a  pooled-AUC ranking (ICI leads all competitors)
#   b  per-cohort AUC heatmap (ICI leads cohort by cohort, never fails)
#   c  added value: 2x2 grids x4 markers (ICI subsumes each marker)
# Each panel is rendered as its own editable SVG (house style, no textLength) for
# Illustrator assembly; a raster composite (this file's PNG) is for LAYOUT REVIEW.
# Reads benchmark tables + pancancer_validation.rds. Numbers computed, never typed.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(grid); library(png) })
source("ici_theme.R"); .pt <- 2.845276
CDIR <- "../results/figures/composite"

## =============================================================================
## PANEL a — pooled-AUC ranking
## =============================================================================
POOL <- read.csv("../results/tables/ICI_benchmark_pooled_auc.csv", check.names=FALSE)
is_ici <- function(v) v=="ICI Integration"
POOL$biom  <- factor(POOL$label, levels=rev(POOL$label[order(POOL$auc)]))
POOL$fillc <- ifelse(is_ici(POOL$label), PLUM, MUTED)
pa <- ggplot(POOL, aes(auc, biom)) +
  geom_vline(xintercept=0.5, linetype="22", colour=RULE, linewidth=0.35) +
  annotate("text", x=0.5, y=0.5, label="chance", family=PUB_FONT, size=(PT_BASE-1)/.pt,
           colour=INK, vjust=1, hjust=0.5) +
  geom_errorbarh(aes(xmin=lo, xmax=hi, colour=fillc), height=0.32, linewidth=1.2) +
  geom_point(aes(colour=fillc), size=4.2) +
  geom_text(data=~subset(.x, !is_ici(label)), aes(x=hi+0.006, label=sprintf("%.2f", auc)),
            family=PUB_FONT, size=(PT_BASE+1.5)/.pt, colour=INK, hjust=0, fontface=1) +
  geom_text(data=~subset(.x, is_ici(label)), aes(x=hi+0.006, label=sprintf("%.2f", auc)),
            family=PUB_FONT, size=(PT_BASE+1.5)/.pt, colour=INK, hjust=0, fontface=2) +
  scale_colour_identity() +
  scale_x_continuous("Pooled AUC  (11 validated cohorts, random effects)",
                     limits=c(0.48, 0.84), breaks=seq(0.5,0.8,0.1), expand=expansion(mult=0)) +
  scale_y_discrete(expand=expansion(add=c(0.6,0.7))) +
  theme_publication(grid="x") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(size=PT_AXIS_TITLE+1, colour=INK),
        axis.text.x=element_text(size=PT_BASE+1, colour=INK),
        axis.text.y=element_text(size=PT_BASE+2.5, colour=INK),
        plot.tag=element_text(size=PT_PANEL_TAG+2, face="bold", family=PUB_FONT, colour=INK)) +
  labs(tag="a")
save_fig(pa, "ICI_Fig4_a_ranking", width_mm=174, height_mm=46, dir=CDIR)

## =============================================================================
## PANEL b — per-cohort AUC heatmap
## =============================================================================
PC <- read.csv("../results/tables/ICI_benchmark_per_cohort_auc.csv", check.names=FALSE)
AUTHOR <- c(IMvigor210="IMvigor210", GSE78220="Hugo 2016", GSE91061="Riaz 2017",
  GSE126044="Cho 2020", GSE135222="Jung 2019", Gide2019="Gide 2019", Liu2019="Liu 2019",
  Puch="PUCH 2022", Kim2018gas="Kim 2018", Rose="Rose 2021", ISPY2="I-SPY2")
CANC <- c(melanoma="Melanoma", urothelial="Urothelial", NSCLC="NSCLC", gastric="Gastric", breast="Breast")
clab <- function(k, cn) sprintf("%s · %s", ifelse(is.na(CANC[cn]), cn, CANC[cn]), AUTHOR[k])
BIOM <- c(ICI="ICI Integration", GEP="T-cell-inflamed GEP", IFNG="IFN-γ signature",
          CYT="Cytolytic activity", CD274="PD-L1 mRNA (CD274)")
BIOM_ORD <- c("ICI Integration","IFN-γ signature","T-cell-inflamed GEP",
              "Cytolytic activity","PD-L1 mRNA (CD274)")
PC$biom   <- factor(unname(BIOM[PC$biom]), levels=rev(BIOM_ORD))
ici       <- PC[PC$biom=="ICI Integration", ]
ord_k     <- as.character(ici$cohort[order(-ici$auc)])
cancer_of <- tapply(as.character(PC$cancer), as.character(PC$cohort), `[`, 1)
labmap    <- setNames(clab(ord_k, cancer_of[ord_k]), ord_k)
PC$clab   <- factor(unname(labmap[as.character(PC$cohort)]), levels=unname(labmap[ord_k]))
pb <- ggplot(PC, aes(clab, biom, fill=auc)) +
  geom_tile(colour="white", linewidth=1.1) +
  geom_text(data=~subset(.x, biom!="ICI Integration"), aes(label=sprintf("%.2f", auc), colour=auc>0.70),
            family=PUB_FONT, size=(PT_BASE+0.5)/.pt, fontface=1) +
  geom_text(data=~subset(.x, biom=="ICI Integration"), aes(label=sprintf("%.2f", auc), colour=auc>0.70),
            family=PUB_FONT, size=(PT_BASE+0.5)/.pt, fontface=2) +
  annotate("segment", x=0.4, xend=length(ord_k)+0.6, y=4.5, yend=4.5, colour=INK, linewidth=0.6) +
  scale_fill_gradient(low="#F3EAF0", high=PLUM, limits=c(0.45,1), oob=scales::squish, name="AUC",
                      breaks=c(0.5,0.75,1.0), guide=guide_colourbar(barheight=unit(26,"mm"), barwidth=unit(3.5,"mm"))) +
  scale_colour_manual(values=c(`TRUE`="white",`FALSE`=INK), guide="none") +
  scale_x_discrete(position="top", expand=expansion(add=0)) +
  scale_y_discrete(expand=expansion(add=0)) +
  coord_cartesian(clip="off") +
  theme_publication() +
  theme(axis.title=element_blank(), panel.border=element_blank(), panel.grid=element_blank(),
        axis.text.x.top=element_text(angle=45, hjust=0, vjust=0, size=PT_BASE+0.5, colour=INK),
        axis.text.y=element_text(size=PT_BASE+2, colour=INK),
        legend.title=element_text(size=PT_LEGEND+1.5, colour=INK),
        legend.text=element_text(size=PT_LEGEND+1, colour=INK),
        axis.ticks=element_blank(), legend.position="right",
        plot.tag=element_text(size=PT_PANEL_TAG+2, face="bold", family=PUB_FONT, colour=INK),
        plot.margin=margin(2,3,2,3,"mm")) +
  labs(tag="b")
save_fig(pb, "ICI_Fig4_b_heatmap", width_mm=174, height_mm=86, dir=CDIR)

## =============================================================================
## PANEL c — added-value 2x2 grids, four markers in one row
## =============================================================================
R <- readRDS("pancancer_validation.rds")
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
MK <- c(GEP="T-cell-inflamed GEP", IFNG="IFN-γ signature",
        CYT="Cytolytic activity", CD274="PD-L1 mRNA")   # (CD274) spelled out in panel b
hi <- function(v) ifelse(v > median(v), "High", "Low")
D <- do.call(rbind, lapply(WORKS, function(k){ s<-R[[k]]$scores
  data.frame(y=s$y, ICI=hi(s$integ_single),
             GEP=hi(s$GEP), IFNG=hi(s$IFNG), CYT=hi(s$CYT), CD274=hi(s$CD274)) }))
rate <- function(sub){ v<-D$y[sub]; c(pct=100*mean(v), n=length(v)) }
g4 <- expand.grid(m=names(MK), M=c("Low","High"), I=c("Low","High"), stringsAsFactors=FALSE)
AG <- do.call(rbind, Map(function(m, ms, is){
  v <- rate(D[[m]]==ms & D$ICI==is)
  data.frame(marker=MK[m], cx=ifelse(ms=="Low",1,2), cy=ifelse(is=="Low",1,2),
             pct=v["pct"], n=v["n"], row.names=NULL)
}, g4$m, g4$M, g4$I))
AG$marker <- factor(AG$marker, levels=unname(MK))
pc <- ggplot(AG, aes(cx, cy)) +
  geom_tile(aes(fill=pct), colour="white", linewidth=1.6, width=0.96, height=0.96) +
  geom_text(aes(label=sprintf("%.0f%%", pct), colour=pct>32), vjust=-0.2,
            family=PUB_FONT, size=(PT_BASE+2)/.pt, fontface=2) +
  geom_text(aes(label=sprintf("n=%d", n), colour=pct>32), vjust=3.1,
            family=PUB_FONT, size=(PT_BASE-3)/.pt) +
  annotate("segment", x=2.6, xend=2.6, y=1.2, yend=1.8,
           arrow=arrow(length=unit(1.6,"mm"), type="closed"), linewidth=0.7, colour=PLUM) +
  annotate("segment", x=1.2, xend=1.8, y=0.42, yend=0.42,
           arrow=arrow(length=unit(1.4,"mm"), type="closed"), linewidth=0.6, colour=SAND_DEEP) +
  facet_wrap(~marker, nrow=1) +
  scale_fill_gradient(low="#F3EAF0", high=PLUM, limits=c(12,50), oob=scales::squish, guide="none") +
  scale_colour_manual(values=c(`TRUE`="white",`FALSE`=INK), guide="none") +
  scale_x_continuous(breaks=c(1,2), labels=c("Low","High"), limits=c(0.45,2.75), expand=expansion(add=0)) +
  scale_y_continuous("ICI Integration", breaks=c(1,2), labels=c("Low","High"),
                     limits=c(0.3,2.55), expand=expansion(add=0)) +
  theme_publication() +
  theme(aspect.ratio=1, axis.title.x=element_blank(),
        axis.title.y=element_text(size=PT_AXIS_TITLE, colour=INK, face=2),
        axis.text=element_text(size=PT_BASE, colour=INK),
        axis.ticks=element_blank(), panel.border=element_blank(), panel.grid=element_blank(),
        strip.text=element_text(size=PT_BASE+0.5, face=2, colour=INK), strip.background=element_blank(),
        panel.spacing=unit(3,"mm"),
        plot.tag=element_text(size=PT_PANEL_TAG+2, face="bold", family=PUB_FONT, colour=INK),
        plot.margin=margin(2,3,2,3,"mm")) +
  labs(tag="c")
save_fig(pc, "ICI_Fig4_c_grids", width_mm=174, height_mm=58, dir=CDIR)

## =============================================================================
## RASTER COMPOSITE (layout review): stack a / b / c, full width
## =============================================================================
rd <- function(f){ a<-readPNG(file.path(CDIR, paste0(f,".png")))
  list(g=rasterGrob(a, interpolate=TRUE, width=unit(1,"npc"), height=unit(1,"npc")), asp=dim(a)[2]/dim(a)[1]) }
A <- rd("ICI_Fig4_a_ranking"); B <- rd("ICI_Fig4_b_heatmap"); C <- rd("ICI_Fig4_c_grids")
W <- 174; gap <- 2
hA <- W/A$asp; hB <- W/B$asp; hC <- W/C$asp; H <- hA+hB+hC+2*gap
place <- function(p, topy, h){ pushViewport(viewport(x=unit(0,"mm"), y=unit(topy,"mm"),
             width=unit(W,"mm"), height=unit(h,"mm"), just=c("left","top")))
  grid.draw(p$g); popViewport() }
MM <- 1/25.4
png_dev <- if (requireNamespace("ragg", quietly=TRUE)) ragg::agg_png else png
png_dev(file.path(CDIR,"ICI_Fig4_composite.png"), width=W*MM, height=H*MM, units="in", res=300, background="white")
grid.newpage()
place(A, H, hA); place(B, H-hA-gap, hB); place(C, H-hA-hB-2*gap, hC)
invisible(dev.off())
cat(sprintf("done — ICI_Fig4_composite.png  (%.0f x %.0f mm)\n", W, H))
