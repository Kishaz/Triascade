# =============================================================================
# SUPPLEMENTARY: leave-one-gene-out stability of the predictive signal.
# Each of the 41 Integration-layer genes is dropped in turn; the I-SPY2
# Integration x treatment interaction P and the pembrolizumab-arm AUC are
# recomputed. Shows no single gene drives the predictive result.
# Source: stability_leavegeneout.csv (run_stability.R). Full-panel reference
# lines are READ from ICI_predictive_interaction.csv (never typed). House style.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
setwd("/Users/smwambu1/Desktop/Projects/ICI/analysis")
source("ici_theme.R"); .pt <- 2.845276

S <- read.csv("../results/tables/stability_leavegeneout.csv", stringsAsFactors=FALSE)
ref <- read.csv("../results/tables/ICI_predictive_interaction.csv", stringsAsFactors=FALSE)
r <- ref[grepl("I-SPY2", ref$trial) & ref$scoring=="locked single-patient" & ref$layer=="Integration", ][1, ]
full_p <- r$interaction_p; full_auc <- r$AUC_icb

S <- S[order(S$interaction_p), ]
S$dropped <- factor(S$dropped, levels=S$dropped)          # shared gene order (ascending P)
ng <- nrow(S)

ytheme <- theme(axis.text.y=element_text(size=PT_BASE-3.5, colour=INK),
                axis.text.x=element_text(size=PT_BASE-2, colour=INK),
                axis.title.x=element_text(size=PT_AXIS_TITLE-1, colour=INK, margin=margin(t=2)),
                plot.title=element_text(size=PT_BASE, face=2, colour=INK, margin=margin(b=3)),
                plot.margin=margin(3,4,3,3,"mm"))

# in-plot key placed in the empty right-hand zone (clear of the lines and the title); black text
lsz <- (PT_BASE-3)/.pt; sw0 <- 0.0540; sw1 <- 0.0575; stx <- 0.0590
pA <- ggplot(S, aes(interaction_p, dropped)) +
  geom_vline(xintercept=0.05,   linetype="dashed", colour=INK, linewidth=0.4) +
  geom_vline(xintercept=full_p, linetype="dotted", colour=INK, linewidth=0.5) +
  geom_segment(aes(x=0, xend=interaction_p, yend=dropped), colour="grey82", linewidth=0.3) +
  geom_point(colour=PLUM, size=1.6) +
  annotate("segment", x=sw0, xend=sw1, y=ng-0.3, yend=ng-0.3, linetype="dashed", colour=INK, linewidth=0.4) +
  annotate("text", x=stx, y=ng-0.3, label="P = 0.05", hjust=0, vjust=0.5, family=PUB_FONT, size=lsz, colour=INK) +
  annotate("segment", x=sw0, xend=sw1, y=ng-2.1, yend=ng-2.1, linetype="dotted", colour=INK, linewidth=0.5) +
  annotate("text", x=stx, y=ng-2.1, label=sprintf("full panel (%.3f)", full_p), hjust=0, vjust=0.5, family=PUB_FONT, size=lsz, colour=INK) +
  scale_x_continuous("Integration × treatment interaction P",
                     limits=c(0, 0.078), breaks=seq(0, 0.05, 0.01), expand=expansion(mult=c(0, 0.02))) +
  scale_y_discrete(NULL) +
  coord_cartesian(clip="off") +
  ggtitle("a  Interaction significance") +
  theme_publication(grid="x") + ytheme

auclo <- 0.763
pB <- ggplot(S, aes(AUC_pembro, dropped)) +
  geom_vline(xintercept=full_auc, linetype="dotted", colour=INK, linewidth=0.5) +
  geom_segment(aes(xend=AUC_pembro, yend=dropped), x=auclo, colour="grey82", linewidth=0.3) +
  geom_point(colour=SAND_DEEP, size=1.6) +
  annotate("text", x=0.7885, y=ng-0.3, label=sprintf("full panel (%.3f)", full_auc), hjust=1, vjust=0.5,
           family=PUB_FONT, size=lsz, colour=INK) +
  annotate("segment", x=0.7895, xend=0.7925, y=ng-0.3, yend=ng-0.3, linetype="dotted", colour=INK, linewidth=0.5) +
  scale_x_continuous("Pembrolizumab-arm Integration AUC", limits=c(auclo, 0.796),
                     breaks=seq(0.76, 0.79, 0.01), expand=expansion(mult=c(0.01, 0.02))) +
  scale_y_discrete(NULL) +
  coord_cartesian(clip="off") +
  ggtitle("b  Discrimination") +
  theme_publication(grid="x") +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=PT_BASE-2, colour=INK),
        axis.title.x=element_text(size=PT_AXIS_TITLE-1, colour=INK, margin=margin(t=2)),
        plot.title=element_text(size=PT_BASE, face=2, colour=INK, margin=margin(b=3)),
        plot.margin=margin(3,4,3,2,"mm"))

p <- pA + pB + plot_layout(widths=c(1.25, 1))
save_fig(p, "ICI_SuppFig_stability", width_mm=174, height_mm=190, dir="../results/figures/supplementary")
cat(sprintf("done - ICI_SuppFig_stability: %d genes, interaction P %.3f-%.3f (all < 0.05 ? %s), AUC %.3f-%.3f, full panel P=%.3f AUC=%.3f\n",
            ng, min(S$interaction_p), max(S$interaction_p), all(S$interaction_p < 0.05),
            min(S$AUC_pembro), max(S$AUC_pembro), full_p, full_auc))
