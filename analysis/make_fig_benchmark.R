# =============================================================================
# Benchmark figure: ICI Integration vs established RNA biomarkers, on the SAME
# patients. (a) pooled AUC across the 11 validated cohorts (random effects) as a
# ranked dot+CI; (b) per-cohort AUC heatmap so the reader sees ICI leads not just
# on average but cohort by cohort. Reads the tables from run_benchmark.R.
# House style (ici_theme.R); single canvas; editable SVG (no textLength).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("ici_theme.R"); .pt <- 2.845276
PC   <- read.csv("../results/tables/ICI_benchmark_per_cohort_auc.csv", check.names=FALSE)
POOL <- read.csv("../results/tables/ICI_benchmark_pooled_auc.csv", check.names=FALSE)

## cohort display labels: Cancer · Author (year) ------------------------------
AUTHOR <- c(IMvigor210="IMvigor210", GSE78220="Hugo 2016", GSE91061="Riaz 2017",
  GSE126044="Cho 2020", GSE135222="Jung 2019", Gide2019="Gide 2019", Liu2019="Liu 2019",
  Puch="PUCH 2022", Kim2018gas="Kim 2018", Rose="Rose 2021", ISPY2="I-SPY2")
CANC <- c(melanoma="Melanoma", urothelial="Urothelial", NSCLC="NSCLC", gastric="Gastric", breast="Breast")
clab <- function(k, cn) sprintf("%s · %s", ifelse(is.na(CANC[cn]), cn, CANC[cn]), AUTHOR[k])

BIOM <- c(ICI="ICI Integration", GEP="T-cell-inflamed GEP", IFNG="IFN-γ signature",
          CYT="Cytolytic activity", CD274="PD-L1 mRNA (CD274)")      # short key -> display label
BIOM_ORD <- c("ICI Integration","IFN-γ signature","T-cell-inflamed GEP",
              "Cytolytic activity","PD-L1 mRNA (CD274)")             # ICI first (highlighted)
is_ici <- function(v) v=="ICI Integration"

## =============================================================================
## (a) pooled AUC ranking — dot + 95% CI, ICI in plum on top
## =============================================================================
POOL$biom  <- factor(POOL$label, levels=rev(POOL$label[order(POOL$auc)]))   # best on top
POOL$fillc <- ifelse(is_ici(POOL$label), PLUM, MUTED)
pa <- ggplot(POOL, aes(auc, biom)) +
  geom_vline(xintercept=0.5, linetype="22", colour=RULE, linewidth=0.35) +
  annotate("text", x=0.5, y=0.5, label="chance", family=PUB_FONT, size=(PT_BASE-1)/.pt,
           colour=MUTED, vjust=1, hjust=0.5) +
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
        axis.title.x=element_text(size=PT_AXIS_TITLE+1.5, colour=INK),
        axis.text.x=element_text(size=PT_BASE+1, colour=INK),
        axis.text.y=element_text(size=PT_BASE+2.5, colour=INK),
        plot.tag=element_text(size=PT_PANEL_TAG+2, face="bold", family=PUB_FONT, colour=INK)) +
  labs(tag="a")

## =============================================================================
## (b) per-cohort AUC heatmap — biomarker x cohort, fill = AUC, value printed
## =============================================================================
PC$biom   <- factor(unname(BIOM[PC$biom]), levels=rev(BIOM_ORD))      # short key -> long label; ICI top row
ici       <- PC[PC$biom=="ICI Integration", ]
ord_k     <- as.character(ici$cohort[order(-ici$auc)])                # cohorts by ICI AUC (desc)
cancer_of <- tapply(as.character(PC$cancer), as.character(PC$cohort), `[`, 1)
labmap    <- setNames(clab(ord_k, cancer_of[ord_k]), ord_k)           # cohort -> display label
PC$clab   <- factor(unname(labmap[as.character(PC$cohort)]), levels=unname(labmap[ord_k]))

pb <- ggplot(PC, aes(clab, biom, fill=auc)) +
  geom_tile(colour="white", linewidth=1.1) +
  geom_text(data=~subset(.x, biom!="ICI Integration"), aes(label=sprintf("%.2f", auc), colour=auc>0.70),
            family=PUB_FONT, size=(PT_BASE+0.5)/.pt, fontface=1) +
  geom_text(data=~subset(.x, biom=="ICI Integration"), aes(label=sprintf("%.2f", auc), colour=auc>0.70),
            family=PUB_FONT, size=(PT_BASE+0.5)/.pt, fontface=2) +
  annotate("segment", x=0.4, xend=length(ord_k)+0.6, y=4.5, yend=4.5, colour=INK, linewidth=0.6) + # under ICI row
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

## ---- compose: a on top (short), b below (tall) -----------------------------
fig <- pa / pb + plot_layout(heights=c(1, 1.5))
dir.create("../results/figures/composite", showWarnings=FALSE, recursive=TRUE)
save_fig(fig, "ICI_Fig_benchmark", width_mm=174, height_mm=150, dir="../results/figures/composite")
cat("done — benchmark figure (a pooled ranking, b per-cohort heatmap)\n")
