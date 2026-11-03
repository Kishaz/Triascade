# =============================================================================
# Added-value figure (incremental predictive value). For each established RNA
# biomarker, two cohort-stratified logistic ORs per SD:
#   plum  = ICI Integration added ON TOP of that biomarker  (large, significant)
#   sand  = that biomarker added ON TOP of ICI Integration  (near 1 = adds nothing)
# A dumbbell joins the pair; OR=1 dashed = "no added value". The asymmetry is the
# message: ICI subsumes the inflammation markers. Reads run_added_value.R table.
# House style; single canvas; editable SVG (no textLength).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R"); .pt <- 2.845276
AV <- read.csv("../results/tables/ICI_added_value.csv", check.names=FALSE)

ICI_ON <- AV[AV$added_key=="ICI", ]                       # ICI added on top of each competitor
CMP_ON <- AV[AV$added_key!="ICI", ]                       # competitor added on top of ICI
CMP_ON$marker <- CMP_ON$added                             # the competitor name
ICI_ON$marker <- ICI_ON$base

ord <- ICI_ON$marker[order(ICI_ON$OR)]                    # weakest ICI-gain at bottom
lev <- ord
ICI_ON$marker <- factor(ICI_ON$marker, levels=lev)
CMP_ON$marker <- factor(CMP_ON$marker, levels=lev)
star <- function(p) ifelse(p<1e-3,"***",ifelse(p<1e-2,"**",ifelse(p<0.05,"*","ns")))

DB <- merge(data.frame(marker=ICI_ON$marker, x_ici=ICI_ON$OR),
            data.frame(marker=CMP_ON$marker, x_cmp=CMP_ON$OR), by="marker")
S_ICI <- "+ ICI Integration beyond marker"; S_CMP <- "+ marker beyond ICI Integration"
ICI_ON$series <- factor(S_ICI, levels=c(S_ICI, S_CMP))
CMP_ON$series <- factor(S_CMP, levels=c(S_ICI, S_CMP))

p <- ggplot() +
  geom_vline(xintercept=1, linetype="22", colour=RULE, linewidth=0.4) +
  annotate("text", x=1, y=4.72, label="no added value", family=PUB_FONT,
           size=(PT_BASE-1)/.pt, colour=MUTED, hjust=0.5, vjust=0) +
  ## dumbbell connector
  geom_segment(data=DB, aes(x=x_cmp, xend=x_ici, y=marker, yend=marker),
               colour="grey70", linewidth=0.9) +
  ## competitor-added-on-top-of-ICI (sand) then ICI-added (plum), with CI
  geom_errorbarh(data=CMP_ON, aes(xmin=lo, xmax=hi, y=marker, colour=series), height=0.16, linewidth=0.9) +
  geom_errorbarh(data=ICI_ON, aes(xmin=lo, xmax=hi, y=marker, colour=series), height=0.16, linewidth=1.1) +
  geom_point(data=CMP_ON, aes(x=OR, y=marker, colour=series), size=4.0) +
  geom_point(data=ICI_ON, aes(x=OR, y=marker, colour=series), size=4.6) +
  ## significance stars above each point
  geom_text(data=ICI_ON, aes(x=OR, y=marker, label=star(p)), vjust=-1.1,
            family=PUB_FONT, size=(PT_BASE-1)/.pt, fontface=2, colour=PLUM) +
  geom_text(data=CMP_ON, aes(x=OR, y=marker, label=star(p)), vjust=-1.1,
            family=PUB_FONT, size=(PT_BASE-1.5)/.pt, colour=SAND_DEEP) +
  scale_colour_manual(values=setNames(c(PLUM, SAND_DEEP), c(S_ICI, S_CMP)), name=NULL,
                      guide=guide_legend(override.aes=list(linewidth=0, size=4.4))) +
  scale_x_continuous("Adjusted odds ratio of ICB response per SD  (cohort-stratified, n=768)",
                     limits=c(0.8, 2.75), breaks=seq(1,2.5,0.5), expand=expansion(mult=0)) +
  scale_y_discrete(expand=expansion(add=c(0.55,0.95))) +
  coord_cartesian(clip="off") +
  theme_publication(grid="x") +
  theme(axis.title.y=element_blank(),
        axis.title.x=element_text(size=PT_AXIS_TITLE+1, colour=INK),
        axis.text.x=element_text(size=PT_BASE+1, colour=INK),
        axis.text.y=element_text(size=PT_BASE+2.5, colour=INK),
        legend.position="top", legend.justification="right", legend.box.spacing=unit(1,"mm"),
        legend.text=element_text(size=PT_BASE+0.5, colour=INK),
        plot.margin=margin(3,6,3,3,"mm"))

dir.create("../results/figures/composite", showWarnings=FALSE, recursive=TRUE)
save_fig(p, "ICI_Fig_addedvalue", width_mm=174, height_mm=100, dir="../results/figures/composite")
cat("done — added-value dumbbell figure\n")
