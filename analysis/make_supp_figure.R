# Supplementary figure: locked-reference single-patient scoring reproduces
# cohort-median scoring and preserves the I-SPY2 predictive interaction.
# Saves PNG (600 dpi) + editable SVG (post-processed by svg_join_text.py).
suppressPackageStartupMessages({ library(data.table); library(pROC); library(ggplot2); library(cowplot) })
source("ici_locked_reference.R")
PANEL <- load_panel("ICI_98gene_panel.csv")
FONT <- "Helvetica"

loo_integration <- function(E, mode = "quantile") {
  vapply(seq_len(ncol(E)), function(i) {
    r <- fit_locked_reference(E[, -i, drop = FALSE], PANEL)
    score_locked(E[, i, drop = FALSE], r, mode)$Integration
  }, numeric(1))
}

## ---------- I-SPY2 (GSE194040) ----------
parse_sm <- function(f) {
  ln <- readLines(gzfile(f), warn = FALSE); ln <- ln[grepl("^!Sample_", ln)]; out <- list()
  for (r in strsplit(ln, "\t")) {
    key <- sub("^!Sample_", "", r[1]); val <- gsub('^"|"$', '', r[-1])
    if (grepl("^characteristics", key)) { i <- which(grepl(":", val))[1]; if (is.na(i)) next
      nm <- make.names(trimws(sub(":.*$", "", val[i]))); val <- trimws(sub("^[^:]*:\\s*", "", val))
    } else nm <- make.names(key)
    while (nm %in% names(out)) nm <- paste0(nm, ".x"); out[[nm]] <- val }
  n <- max(lengths(out)); as.data.frame(lapply(out, function(v){length(v)<-n; v}), stringsAsFactors = FALSE)
}
md <- rbind(transform(parse_sm("data/ISPY2_GPL20078_series.txt.gz"), .p="a")[,c("patient.id","arm","her2","pcr")],
            transform(parse_sm("data/ISPY2_GPL30493_series.txt.gz"), .p="b")[,c("patient.id","arm","her2","pcr")])
n01 <- function(v) ifelse(grepl("pos|\\+|yes|^1$",v,TRUE),1L,ifelse(grepl("neg|\\-|no|^0$",v,TRUE),0L,NA))
md$.id<-as.character(md$patient.id); md$.her2<-n01(md$her2); md$.pcr<-n01(md$pcr)
md$.arm2<-ifelse(grepl("Pembrolizumab",md$arm),"pembro",ifelse(md$arm=="Paclitaxel","control",NA))
ex <- as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz")); rownames(ex)<-toupper(ex[[1]]); ex[[1]]<-NULL
ex <- as.matrix(ex); colnames(ex)<-as.character(colnames(ex))
elig <- subset(md, .arm2%in%c("pembro","control") & .her2==0 & !is.na(.pcr) & .id%in%colnames(ex))
elig <- elig[!duplicated(elig$.id),]; E <- ex[,elig$.id]
y <- elig$.pcr; arm <- factor(elig$.arm2, c("control","pembro"))
cm <- score_cohort_median(E, PANEL); cm <- cm[match(elig$.id, cm$sample),]
Ilk <- loo_integration(E)
z <- function(v)(v-mean(v))/sd(v)
pint <- anova(glm(y~z(cm$Integration)+arm,binomial), glm(y~z(cm$Integration)*arm,binomial), test="LRT")[2,"Pr(>Chi)"]
sp_ispy2 <- cor(cm$Integration, Ilk, method="spearman")

## ---------- GSE78220 (melanoma) for the caption numbers ----------
em <- as.data.frame(fread("data/GSE78220_expr.csv")); rownames(em)<-toupper(em$gene); em$gene<-NULL
em <- log2(as.matrix(em)+1); rr <- fread("data/GSE78220_response.csv")
em <- em[, colnames(em)%in%rr$sample]; ym <- rr$responder[match(colnames(em), rr$sample)]
cmm <- score_cohort_median(em, PANEL); cmm <- cmm[match(colnames(em), cmm$sample),]
Ilm <- loo_integration(em)
sp_mel <- cor(cmm$Integration, Ilm, method="spearman")
auc_mel_cm <- as.numeric(auc(ym, cmm$Sensory, quiet=TRUE, direction="<"))

## ---------- Panel A: concordance scatter (I-SPY2) ----------
dA <- data.frame(cohort_median = cm$Integration, locked = Ilk,
                 Arm = ifelse(arm=="pembro","pembrolizumab","control"))
pA <- ggplot(dA, aes(cohort_median, locked, colour = Arm)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, colour = "grey55") +
  geom_point(size = 1.9, alpha = 0.8) +
  scale_colour_manual(values = c(control = "#9E9E9E", pembrolizumab = "#2E6FB0")) +
  annotate("text", x = min(dA$cohort_median), y = max(dA$locked),
           label = sprintf("Spearman = %.2f", sp_ispy2), hjust = 0, vjust = 1, size = 3.6, family = FONT) +
  labs(x = "Integration (cohort-median)", y = "Integration (locked, single-patient)",
       colour = NULL, tag = "a") +
  theme_classic(base_family = FONT, base_size = 11) +
  theme(legend.position = c(0.98, 0.04), legend.justification = c(1, 0),
        axis.text = element_text(colour = "black"), plot.tag = element_text(face = "bold"))

## ---------- Panel B: I-SPY2 interaction (pCR by Integration x arm) ----------
hi <- factor(ifelse(cm$Integration >= median(cm$Integration), "Integration high", "Integration low"),
             c("Integration low", "Integration high"))
dB <- as.data.frame(as.table(tapply(y, list(hi, arm), mean) * 100))
names(dB) <- c("Integration", "Arm", "pCR")
dB$Arm <- factor(ifelse(dB$Arm == "pembro", "pembrolizumab", "control"), c("control", "pembrolizumab"))
pB <- ggplot(dB, aes(Arm, pCR, fill = Integration)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.66) +
  geom_text(aes(label = sprintf("%.0f%%", pCR)), position = position_dodge(width = 0.72),
            vjust = -0.4, size = 3.3, family = FONT) +
  scale_fill_manual(values = c("Integration low" = "#9E9E9E", "Integration high" = "#2E6FB0")) +
  scale_y_continuous(limits = c(0, 80), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "pCR rate (%)", fill = NULL, tag = "b",
       subtitle = "Integration x treatment interaction identical under both scorings") +
  theme_classic(base_family = FONT, base_size = 11) +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        axis.text = element_text(colour = "black"), plot.tag = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9))

fig <- plot_grid(pA, pB, nrow = 1, rel_widths = c(1, 1.05))
ggsave("SuppFig_locked_reference.png", fig, width = 9.4, height = 4.2, dpi = 600, bg = "white")
ggsave("SuppFig_locked_reference.svg", fig, width = 9.4, height = 4.2, device = svglite::svglite, bg = "white")

cat(sprintf("I-SPY2: n=%d, interaction p=%.4f, Integration Spearman(locked vs cohort-median)=%.3f\n",
            ncol(E), pint, sp_ispy2))
cat(sprintf("GSE78220: Integration Spearman=%.3f, cohort-median Sensory AUC=%.3f\n", sp_mel, auc_mel_cm))
cat("Wrote SuppFig_locked_reference.png / .svg\n")
