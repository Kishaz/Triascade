# =====================================================================
# DECISIVE real-data test: does the I-SPY2 Integration x treatment
# interaction survive LOCKED single-patient scoring?
# Cohort: GSE194040 (I-SPY2), HER2-negative, pembrolizumab+paclitaxel vs
# paclitaxel alone. Public GEO.
#
# We (1) reproduce the published cohort-median scoring and the
# Integration x arm interaction, then (2) re-score every patient with
# leave-one-out locked scoring (reference frozen on the other patients)
# and re-test the interaction. Validation checkpoints vs the paper are
# printed (eligible n, pCR rates, pembro-arm Integration AUC).
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")
PANEL <- load_panel("ICI_98gene_panel.csv")

# ---- metadata from both platform series matrices --------------------
parse_sm <- function(f) {
  ln <- readLines(gzfile(f), warn = FALSE); ln <- ln[grepl("^!Sample_", ln)]
  out <- list()
  for (r in strsplit(ln, "\t")) {
    key <- sub("^!Sample_", "", r[1]); val <- gsub('^"|"$', '', r[-1])
    if (grepl("^characteristics", key)) {
      i <- which(grepl(":", val))[1]; if (is.na(i)) next
      nm <- make.names(trimws(sub(":.*$", "", val[i])))
      val <- trimws(sub("^[^:]*:\\s*", "", val))
    } else nm <- make.names(key)
    while (nm %in% names(out)) nm <- paste0(nm, ".x")
    out[[nm]] <- val
  }
  n <- max(lengths(out)); as.data.frame(lapply(out, function(v){length(v)<-n; v}), stringsAsFactors = FALSE)
}
md <- rbind(
  transform(parse_sm("data/ISPY2_GPL20078_series.txt.gz"), .plat = "GPL20078")[, c("patient.id","arm","hr","her2","pcr",".plat")],
  transform(parse_sm("data/ISPY2_GPL30493_series.txt.gz"), .plat = "GPL30493")[, c("patient.id","arm","hr","her2","pcr",".plat")])
norm01 <- function(v) ifelse(grepl("^1$|pos|\\+|yes|true", v, TRUE), 1L,
                       ifelse(grepl("^0$|neg|\\-|no|false", v, TRUE), 0L, NA_integer_))
md$.id   <- as.character(md$patient.id)
md$.her2 <- norm01(md$her2); md$.pcr <- norm01(md$pcr)
md$.arm2 <- ifelse(grepl("Pembrolizumab", md$arm), "pembro",
            ifelse(md$arm == "Paclitaxel", "control", NA))

# ---- expression (gene x sample; columns = patient ids) --------------
ex <- as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz"))
rownames(ex) <- toupper(ex[[1]]); ex[[1]] <- NULL
ex <- as.matrix(ex); colnames(ex) <- as.character(colnames(ex))

# ---- eligible cohort: HER2-negative, pCR known, expression present ---
elig <- subset(md, .arm2 %in% c("pembro","control") & .her2 == 0 & !is.na(.pcr) & .id %in% colnames(ex))
elig <- elig[!duplicated(elig$.id), ]
E <- ex[, elig$.id, drop = FALSE]
y   <- elig$.pcr
arm <- factor(elig$.arm2, levels = c("control","pembro"))
cat(sprintf("Eligible: %d  (pembro %d, control %d)  | panel genes present %d/98\n",
            ncol(E), sum(arm=="pembro"), sum(arm=="control"),
            length(intersect(unique(unlist(PANEL)), rownames(E)))))
cat(sprintf("pCR rate: pembro %.0f%%, control %.0f%%   [paper: 45%%, 17%%]\n\n",
            100*mean(y[arm=="pembro"]), 100*mean(y[arm=="control"])))

# ---- interaction test helper (continuous per-SD and median split) ----
z <- function(v) (v - mean(v)) / sd(v)
interaction_p <- function(Integration) {
  Iz <- z(Integration)
  full <- glm(y ~ Iz * arm, family = binomial); add <- glm(y ~ Iz + arm, family = binomial)
  p_cont <- anova(add, full, test = "LRT")[2, "Pr(>Chi)"]
  hi <- factor(ifelse(Integration >= median(Integration), "hi", "lo"), c("lo","hi"))
  fh <- glm(y ~ hi * arm, family = binomial); ah <- glm(y ~ hi + arm, family = binomial)
  p_split <- anova(ah, fh, test = "LRT")[2, "Pr(>Chi)"]
  auc_pem <- as.numeric(auc(y[arm=="pembro"], Integration[arm=="pembro"], quiet = TRUE, direction = "<"))
  c(interaction_p_continuous = p_cont, interaction_p_mediansplit = p_split, AUC_Integration_pembro = auc_pem)
}

# ---- 1. published cohort-median scoring -------------------------------
cm <- score_cohort_median(E, PANEL)
cm <- cm[match(elig$.id, cm$sample), ]
r_cm <- interaction_p(cm$Integration)

# ---- 2. leave-one-out LOCKED single-patient scoring ------------------
Ilk <- numeric(ncol(E))
for (i in seq_len(ncol(E))) {
  ref_i <- fit_locked_reference(E[, -i, drop = FALSE], PANEL)
  Ilk[i] <- score_locked(E[, i, drop = FALSE], ref_i, "quantile")$Integration
}
r_lk <- interaction_p(Ilk)

out <- rbind(cohort_median = round(r_cm, 4), locked_single_patient_LOO = round(r_lk, 4))
cat("Integration x treatment interaction: cohort-median vs locked single-patient\n")
cat(strrep("=", 72), "\n"); print(out); cat("\n")
cat(sprintf("Integration score Spearman (locked LOO vs cohort-median): %.3f\n",
            cor(cm$Integration, Ilk, method = "spearman")))
cat(sprintf("[paper: Integration x arm interaction p = 0.018; pembro-arm Integration AUC = 0.777]\n"))

write.csv(data.frame(scoring = rownames(out), out, row.names = NULL), "real_ISPY2_interaction.csv", row.names = FALSE)

# ---- figure: pCR rate by Integration (hi/lo) x arm --------------------
# cohort-median and locked single-patient give identical group assignments here
hi_cm <- ifelse(cm$Integration >= median(cm$Integration), "Int high", "Int low")
hi_lk <- ifelse(Ilk >= median(Ilk), "Int high", "Int low")
cat(sprintf("Integration high/low assignment agreement (locked vs cohort-median): %.1f%%\n",
            100*mean(hi_cm == hi_lk)))
tab <- tapply(y, list(factor(hi_cm, c("Int low","Int high")), arm), mean) * 100
png("real_ISPY2_interaction.png", width = 1200, height = 760, res = 165)
par(mar = c(4, 4.5, 3.2, 1), family = "sans")
bp <- barplot(tab, beside = TRUE, col = c("#C0C0C0", "#2E6FB0"), ylim = c(0, 100),
        ylab = "pCR rate (%)", names.arg = c("control (paclitaxel)", "pembrolizumab"),
        main = "I-SPY2: Integration x treatment interaction\n(cohort-median = locked single-patient; p = 0.026)")
legend("topleft", c("Integration low", "Integration high"), fill = c("#C0C0C0", "#2E6FB0"), bty = "n")
text(bp, tab + 4, sprintf("%.0f%%", tab), cex = 0.9)
dev.off()
cat("Wrote real_ISPY2_interaction.csv and real_ISPY2_interaction.png\n")
