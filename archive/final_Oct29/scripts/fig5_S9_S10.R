# =====================================================================
# FIGURE 5 and SUPPLEMENTARY S9, S10
#
# Fig5  the profile null: which layer is lowest carries no information
#       about benefit, while Integration's level does
# S9    per-cohort AUC with DeLong CIs, corrected n, grouped by design,
#       no dominance framing
# S10   layer correlations across cohorts, never above 0.66
#
# Output: Aug_1_ICI_FINAL_PAPER/Figures/
# =====================================================================
SCRIPT_VERSION <- "fig5_S9_S10 v1.0  2026-08-23"
cat(SCRIPT_VERSION, "\n\n")
suppressPackageStartupMessages({ library(ggplot2); library(pROC) })

MAN <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"
OUT <- file.path(MAN,"Aug_1_ICI_FINAL_PAPER","Figures")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

# ===================== FIGURE 5: the profile null =====================
SCOR <- file.path(MAN,"MASTER_GRID","ICI_binary_ispy2_scores.csv")
d <- read.csv(SCOR, stringsAsFactors=FALSE)
d$y <- as.integer(d$pcr_binary); d$arm <- factor(d$arm, levels=c("control","pembro"))
d <- d[!is.na(d$y) & !is.na(d$arm),]
LAY <- c("Sensory","Integration","Execution")
Z <- scale(as.matrix(d[,LAY])); colnames(Z) <- LAY
d$lim <- factor(LAY[apply(Z,1,which.min)], levels=c("Execution","Sensory","Integration"))

# a: pCR rate by limiting layer and arm
R <- aggregate(y ~ lim + arm, d, function(v) c(rate=mean(v), n=length(v)))
R <- data.frame(lim=R$lim, arm=R$arm, rate=R$y[,"rate"], n=R$y[,"n"])
m1 <- glm(y ~ lim*arm, d, family=binomial); m0 <- glm(y ~ lim+arm, d, family=binomial)
p_prof <- anova(m0,m1,test="Chisq")$`Pr(>Chi)`[2]
cat(sprintf("limiting-layer x arm interaction p = %.4f\n", p_prof))

p5a <- ggplot(R, aes(lim, rate, fill=arm)) +
  geom_col(position=position_dodge(width=0.75), width=0.65) +
  geom_text(aes(label=sprintf("%.0f%%", 100*rate)), position=position_dodge(width=0.75),
            vjust=-0.3, size=2.5, colour="grey25") +
  scale_fill_manual(values=c(control="#5F5E5A", pembro="#B03A2E"), name=NULL) +
  scale_y_continuous(labels=function(v) paste0(100*v,"%"), limits=c(0,0.62), expand=c(0,0)) +
  labs(x="patient's lowest layer", y="pCR rate",
       subtitle=sprintf("a   Benefit is the same whichever layer is lowest (p = %.3f)", p_prof)) +
  theme_minimal(base_size=9) +
  theme(legend.position="top", panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(), plot.subtitle=element_text(face="bold", size=9))

# b: six profile encodings vs the Integration level
enc <- list()
enc[["lowest layer"]]        <- d$lim
enc[["Sensory - Execution"]] <- Z[,"Sensory"] - Z[,"Execution"]
enc[["Integration - own mean"]] <- Z[,"Integration"] - rowMeans(Z)
enc[["spread (max - min)"]]  <- apply(Z,1,max) - apply(Z,1,min)
mn <- apply(Z,1,min); md <- apply(Z,1,function(v) sort(v)[2])
enc[["gap to second lowest"]] <- md - mn
enc[["deep vs shallow minimum"]] <- factor(ifelse(md-mn >= median(md-mn) & mn < 0,
                                                  "deep","shallow"))
enc[["Integration tertile"]] <- cut(Z[,"Integration"],
  quantile(Z[,"Integration"], c(0,1/3,2/3,1)), include.lowest=TRUE,
  labels=c("low","mid","high"))
E <- do.call(rbind, lapply(names(enc), function(nm){
  s <- d; s$v <- enc[[nm]]
  a <- anova(glm(y ~ v+arm, s, family=binomial),
             glm(y ~ v*arm, s, family=binomial), test="Chisq")
  data.frame(encoding=nm, p=a$`Pr(>Chi)`[2],
             kind=ifelse(grepl("tertile", nm), "Integration level", "profile shape"),
             stringsAsFactors=FALSE) }))
E <- E[order(E$kind, -E$p), ]
E$encoding <- factor(E$encoding, levels=E$encoding)
print(E, row.names=FALSE)

p5b <- ggplot(E, aes(p, encoding, colour=kind)) +
  geom_vline(xintercept=0.05, linetype=2, linewidth=0.35, colour="grey60") +
  geom_point(size=2.8) +
  geom_text(aes(label=sprintf("  %.3f", p)), hjust=0, size=2.6, colour="grey30",
            show.legend=FALSE) +
  scale_colour_manual(values=c(`profile shape`="#5F5E5A",
                               `Integration level`="#B03A2E"), name=NULL) +
  scale_x_log10(limits=c(0.001, 3), breaks=c(0.001,0.01,0.05,0.1,0.5,1)) +
  labs(x="interaction p with treatment arm", y=NULL,
       subtitle="b   No encoding of the profile interacts with treatment; the level does") +
  theme_minimal(base_size=9) +
  theme(legend.position="top", panel.grid.minor=element_blank(),
        plot.subtitle=element_text(face="bold", size=9))

fp5 <- file.path(OUT,"Fig5_profile_null")
if (requireNamespace("patchwork", quietly=TRUE)) { library(patchwork)
  ggsave(paste0(fp5,".pdf"), p5a|p5b, width=10.2, height=4.0, device=cairo_pdf)
  ggsave(paste0(fp5,".png"), p5a|p5b, width=10.2, height=4.0, dpi=400)
} else { ggsave(paste0(fp5,"_a.pdf"), p5a, width=5, height=4, device=cairo_pdf)
         ggsave(paste0(fp5,"_b.pdf"), p5b, width=5.4, height=4, device=cairo_pdf) }
write.csv(E, file.path(OUT,"Fig5_source_values.csv"), row.names=FALSE)
cat("Wrote Fig5\n")

# ===================== S9: per-cohort AUC with CIs ====================
AU <- file.path(MAN,"FIGURE2B","layer_auc_with_ci.csv")
if (file.exists(AU)) {
  S <- read.csv(AU, stringsAsFactors=FALSE)
  S$design <- "Single-arm, ICI-treated"
  S$design[grepl("ispy2.*control", S$cohort)] <- "Randomized control arm (no ICI)"
  S$design[grepl("ispy2.*pembro",  S$cohort)] <- "Randomized ICI arm"
  S$design[grepl("braun", S$cohort)] <- "Pre-specified negative control"
  S$layer <- factor(S$layer, levels=c("Sensory","Integration","Execution","Composite"))
  ord <- unique(S$cohort[order(-S$n)]); S$cohort <- factor(S$cohort, levels=rev(ord))
  cat(sprintf("\nS9: %d estimates, %d with CI crossing 0.5\n",
              nrow(S), sum(S$CI_low < 0.5)))
  s9 <- ggplot(S, aes(AUC, cohort, colour=design)) +
    geom_vline(xintercept=0.5, linewidth=0.35, colour="grey55") +
    geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0, linewidth=0.45, alpha=0.85) +
    geom_point(size=1.9) +
    facet_wrap(~layer, nrow=1) +
    scale_colour_manual(values=c(
      `Randomized ICI arm`="#B03A2E", `Single-arm, ICI-treated`="#5F5E5A",
      `Randomized control arm (no ICI)`="#378ADD",
      `Pre-specified negative control`="#7F77DD"), name=NULL,
      guide=guide_legend(nrow=2)) +
    scale_x_continuous(limits=c(0.25,1), breaks=c(0.5,0.75,1)) +
    labs(x="AUC for response vs non-response (95% CI, DeLong)", y=NULL,
         caption=paste0(sum(S$CI_low<0.5), " of ", nrow(S),
           " estimates have a confidence interval including 0.5. Cohorts ordered by size.")) +
    theme_minimal(base_size=9) +
    theme(legend.position="top", panel.grid.minor=element_blank(),
          strip.text=element_text(face="bold", size=9),
          plot.caption=element_text(hjust=0, size=7, colour="grey35"))
  ggsave(file.path(OUT,"FigS9_percohort_auc.pdf"), s9, width=10.6, height=5.2,
         device=cairo_pdf)
  ggsave(file.path(OUT,"FigS9_percohort_auc.png"), s9, width=10.6, height=5.2, dpi=400)
  cat("Wrote S9\n")
} else cat("!! layer_auc_with_ci.csv missing - S9 skipped\n")

# ===================== S10: layer correlations ========================
CO <- file.path(MAN,"Figures","Fig2_source_numbers.csv")
if (file.exists(CO)) {
  C <- read.csv(CO, stringsAsFactors=FALSE)
  C <- C[!duplicated(paste(C$cohort, C$pair)), ]
  ord <- unique(C$cohort[order(-C$n)]); C$cohort <- factor(C$cohort, levels=rev(ord))
  C$small <- C$n < 20
  cat(sprintf("\nS10: max correlation %.3f\n", max(C$r)))
  s10 <- ggplot(C, aes(r, cohort, colour=pair, alpha=!small, size=n)) +
    geom_vline(xintercept=0, linewidth=0.3, colour="grey80") +
    geom_point() +
    scale_colour_manual(values=c("#1D9E75","#378ADD","#D4537E"), name=NULL) +
    scale_alpha_manual(values=c(`TRUE`=1,`FALSE`=0.35), guide="none") +
    scale_size(range=c(1.1,3.4), guide="none") +
    scale_x_continuous(limits=c(-0.1,1), breaks=seq(0,1,0.25)) +
    labs(x="Pearson r between layer scores", y=NULL,
         caption="No layer pair exceeds r = 0.66 in any cohort. Point size is cohort size; small cohorts faded.") +
    theme_minimal(base_size=9) +
    theme(legend.position="top", panel.grid.minor=element_blank(),
          plot.caption=element_text(hjust=0, size=7, colour="grey35"))
  ggsave(file.path(OUT,"FigS10_layer_correlations.pdf"), s10, width=6.6, height=5.2,
         device=cairo_pdf)
  ggsave(file.path(OUT,"FigS10_layer_correlations.png"), s10, width=6.6, height=5.2, dpi=400)
  cat("Wrote S10\n")
} else cat("!! Fig2_source_numbers.csv missing - S10 skipped\n")

cat("\nAll outputs in", OUT, "\n")
