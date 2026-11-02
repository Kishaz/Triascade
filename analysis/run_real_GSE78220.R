# =====================================================================
# REAL-DATA test of locked-reference single-patient scoring
# Cohort: GSE78220 (Hugo), melanoma, anti-PD-1, n=28 (public GEO)
#
# Leave-one-out locked scoring: for each patient, the reference is frozen
# on the OTHER 27 patients and that patient is scored alone. No patient
# ever uses their own cohort median. We then check that this reproduces
# the published cohort-median scoring and preserves the response AUC.
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")
PANEL <- load_panel("ICI_98gene_panel.csv")

expr <- as.data.frame(fread("data/GSE78220_expr.csv"))
rownames(expr) <- toupper(expr$gene); expr$gene <- NULL
expr <- as.matrix(expr)
expr <- log2(expr + 1)                      # FPKM -> log2 (rank-invariant for median/quantile)
resp <- fread("data/GSE78220_response.csv")
keep <- colnames(expr) %in% resp$sample      # drop on-treatment aliquot (Pt16.OnTx); baseline only
expr <- expr[, keep, drop = FALSE]
resp <- resp[match(colnames(expr), resp$sample), ]
y <- resp$responder
cat(sprintf("GSE78220: %d genes x %d samples | responders %d / non %d | panel present %d/98\n\n",
            nrow(expr), ncol(expr), sum(y), sum(!y),
            length(intersect(unique(unlist(PANEL)), rownames(expr)))))

# ---- 1. published cohort-median scoring --------------------------------
cm <- score_cohort_median(expr, PANEL)

# ---- 2. leave-one-out LOCKED scoring (reference excludes the patient) ---
loo <- function(mode) {
  out <- lapply(seq_len(ncol(expr)), function(i) {
    ref_i <- fit_locked_reference(expr[, -i, drop = FALSE], PANEL)
    score_locked(expr[, i, drop = FALSE], ref_i, mode = mode)
  })
  d <- do.call(rbind, out); d$sample <- colnames(expr); d
}
lk_q <- loo("quantile")
lk_s <- loo("selfnorm")

# ---- 3. concordance + response AUC -------------------------------------
auc1 <- function(score) as.numeric(auc(y, score, quiet = TRUE, direction = "<"))
summ <- function(lk, tag) {
  con <- concordance(cm, lk)
  data.frame(
    scoring = tag,
    Integration_spearman = round(con$spearman[con$layer == "Integration"], 3),
    Composite_spearman   = round(con$spearman[con$layer == "Composite"], 3),
    archetype_agreement  = round(archetype_agreement(cm, lk), 3),
    AUC_Sensory     = round(auc1(lk$Sensory), 3),
    AUC_Integration = round(auc1(lk$Integration), 3),
    AUC_Execution   = round(auc1(lk$Execution), 3),
    AUC_Composite   = round(auc1(lk$Composite), 3))
}
cmrow <- data.frame(scoring = "cohort_median (published rule)",
    Integration_spearman = 1, Composite_spearman = 1, archetype_agreement = 1,
    AUC_Sensory = round(auc1(cm$Sensory), 3), AUC_Integration = round(auc1(cm$Integration), 3),
    AUC_Execution = round(auc1(cm$Execution), 3), AUC_Composite = round(auc1(cm$Composite), 3))
res <- rbind(cmrow, summ(lk_q, "locked LOO (quantile)"), summ(lk_s, "locked LOO (selfnorm)"))
cat("REAL DATA (GSE78220): locked single-patient (leave-one-out) vs cohort-median\n")
cat(strrep("=", 74), "\n"); print(res, row.names = FALSE); cat("\n")

# gene-call agreement (quantile LOO vs cohort-median)
bc_cm <- t(apply(expr[intersect(unique(unlist(PANEL)), rownames(expr)), ], 1,
                 function(v) as.integer(v >= median(v))))
agree <- sapply(seq_len(ncol(expr)), function(i) {
  ref_i <- fit_locked_reference(expr[, -i, drop = FALSE], PANEL)
  g <- intersect(rownames(bc_cm), ref_i$genes)
  cl <- sapply(g, function(gi) as.integer(.ref_quantile(expr[gi, i], ref_i$ref_sorted[[gi]]) >= 0.5))
  mean(cl == bc_cm[g, i])
})
cat(sprintf("Gene-call agreement (locked LOO quantile vs cohort-median): mean %.3f (range %.2f-%.2f)\n",
            mean(agree), min(agree), max(agree)))

write.csv(res, "real_GSE78220_concordance.csv", row.names = FALSE)

png("real_GSE78220_locked_vs_cohortmedian.png", width = 1500, height = 680, res = 165)
op <- par(mfrow = c(1, 2), mar = c(4.3, 4.3, 2.6, 1), family = "sans")
plot(cm$Integration, lk_q$Integration, pch = 19, col = ifelse(y == 1, "#1E8449", "#B03A2E"),
     xlab = "Integration (cohort-median)", ylab = "Integration (locked LOO, single-patient)",
     main = sprintf("GSE78220  (Spearman %.2f)",
                    cor(cm$Integration, lk_q$Integration, method = "spearman")))
abline(0, 1, lty = 2, col = "grey40"); legend("topleft", c("responder","non-responder"),
     pch = 19, col = c("#1E8449","#B03A2E"), bty = "n", cex = 0.9)
barplot(as.matrix(res[, c("AUC_Integration","AUC_Composite")]), beside = TRUE,
        names.arg = c("Integration","Composite"), col = c("#808080","#2E6FB0","#A23B72"),
        ylim = c(0, 1), ylab = "AUC for anti-PD-1 response")
legend("topright", res$scoring, fill = c("#808080","#2E6FB0","#A23B72"), bty = "n", cex = 0.75)
abline(h = 0.5, lty = 3, col = "grey60"); par(op); dev.off()
cat("Wrote real_GSE78220_concordance.csv and real_GSE78220_locked_vs_cohortmedian.png\n")
