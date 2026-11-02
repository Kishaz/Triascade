# =====================================================================
# (1) Proportion-normalized layer scores  (coverage-robust: each layer
#     is scored as a FRACTION of its measured genes, so a cohort missing
#     a gene stays on the same scale) and
# (2) Leave-one-gene-out stability of the Integration x treatment result
#     in I-SPY2 (is the signal driven by any single gene?)
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")
PANEL <- load_panel("ICI_98gene_panel.csv")

# proportion-normalized median scoring (fraction of measured genes coded 1)
score_prop <- function(expr, PANEL) {
  binf <- function(gs) { g <- intersect(gs, rownames(expr)); if (!length(g)) return(rep(NA, ncol(expr)))
    b <- t(apply(expr[g, , drop = FALSE], 1, function(v) as.integer(v >= median(v, na.rm = TRUE))))
    colMeans(b, na.rm = TRUE) }                       # proportion high in [0,1]
  S <- binf(PANEL$SENSORY)
  I <- binf(PANEL$INTPOS) - binf(PANEL$INTNEG)         # [-1,1]
  E <- binf(PANEL$EXECPOS) - binf(PANEL$EXECNEG)
  z <- function(v) (v - mean(v)) / sd(v)
  data.frame(sample = colnames(expr), Sensory = S, Integration = I, Execution = E,
             Composite = rowMeans(cbind(z(S), z(I), z(E))), row.names = NULL)
}
auc1 <- function(y, s) as.numeric(auc(y, s, quiet = TRUE, direction = "<"))

## ---- load I-SPY2 ----
parse_sm <- function(f) { ln <- readLines(gzfile(f), warn = FALSE); ln <- ln[grepl("^!Sample_", ln)]; out <- list()
  for (r in strsplit(ln, "\t")) { k <- sub("^!Sample_", "", r[1]); v <- gsub('"', '', r[-1])
    if (grepl("^characteristics", k)) { i <- which(grepl(":", v))[1]; if (is.na(i)) next
      nm <- make.names(trimws(sub(":.*$", "", v[i]))); v <- trimws(sub("^[^:]*: *", "", v)) } else nm <- make.names(k)
    while (nm %in% names(out)) nm <- paste0(nm, ".x"); out[[nm]] <- v }
  n <- max(lengths(out)); as.data.frame(lapply(out, function(x){length(x)<-n; x}), stringsAsFactors = FALSE) }
md <- rbind(parse_sm("data/ISPY2_GPL20078_series.txt.gz")[,c("patient.id","arm","her2","pcr")],
            parse_sm("data/ISPY2_GPL30493_series.txt.gz")[,c("patient.id","arm","her2","pcr")])
n01 <- function(v) ifelse(grepl("pos|[+]|yes|^1$", v, TRUE), 1L, ifelse(grepl("neg|[-]|no|^0$", v, TRUE), 0L, NA))
md$.id <- as.character(md$patient.id); md$.her2 <- n01(md$her2); md$.pcr <- n01(md$pcr)
md$.arm2 <- ifelse(grepl("Pembrolizumab", md$arm), "pembro", ifelse(md$arm == "Paclitaxel", "control", NA))
ex <- as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz")); rownames(ex) <- toupper(ex[[1]]); ex[[1]] <- NULL
ex <- as.matrix(ex); colnames(ex) <- as.character(colnames(ex))
elig <- subset(md, .arm2 %in% c("pembro","control") & .her2 == 0 & !is.na(.pcr) & .id %in% colnames(ex))
elig <- elig[!duplicated(elig$.id), ]; E <- ex[, elig$.id]
y <- elig$.pcr; arm <- factor(elig$.arm2, c("control","pembro"))
z <- function(v) scale(v)[,1]
int_p <- function(sc) anova(glm(y~z(sc)+arm, binomial), glm(y~z(sc)*arm, binomial), test="LRT")[2,"Pr(>Chi)"]

## ---- (1) proportion-normalized vs raw-count: same within-cohort result ----
cm  <- score_cohort_median(E, PANEL); cm <- cm[match(elig$.id, cm$sample), ]
pr  <- score_prop(E, PANEL);          pr <- pr[match(elig$.id, pr$sample), ]
cat("I-SPY2 Integration x treatment interaction:\n")
cat(sprintf("  raw-count score        p = %.3f\n", int_p(cm$Integration)))
cat(sprintf("  proportion-normalized  p = %.3f   (Spearman vs raw = %.3f)\n\n",
            int_p(pr$Integration), cor(cm$Integration, pr$Integration, method="spearman")))

## ---- (2) leave-one-gene-out stability of the Integration result ----
intgenes <- c(PANEL$INTPOS, PANEL$INTNEG)
full_p   <- int_p(pr$Integration); full_auc <- auc1(y[arm=="pembro"], pr$Integration[arm=="pembro"])
loo <- lapply(intgenes, function(drop) {
  P2 <- PANEL; P2$INTPOS <- setdiff(PANEL$INTPOS, drop); P2$INTNEG <- setdiff(PANEL$INTNEG, drop)
  s <- score_prop(E, P2); s <- s[match(elig$.id, s$sample), ]
  data.frame(dropped = drop, interaction_p = int_p(s$Integration),
             AUC_pembro = auc1(y[arm=="pembro"], s$Integration[arm=="pembro"]))
})
loo <- do.call(rbind, loo)
cat(sprintf("Leave-one-gene-out (dropping each of %d Integration genes):\n", length(intgenes)))
cat(sprintf("  interaction p : full %.3f | LOO range %.3f - %.3f | all < 0.05 ? %s\n",
            full_p, min(loo$interaction_p), max(loo$interaction_p), all(loo$interaction_p < 0.05)))
cat(sprintf("  pembro AUC    : full %.3f | LOO range %.3f - %.3f\n", full_auc, min(loo$AUC_pembro), max(loo$AUC_pembro)))
cat("  most influential (largest p when dropped):\n")
print(head(loo[order(-loo$interaction_p), ], 5), row.names = FALSE)
write.csv(loo, "stability_leavegeneout.csv", row.names = FALSE)

## ---- figure ----
png("stability_leavegeneout.png", width = 1500, height = 640, res = 165)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2.6, 1), family = "sans")
plot(sort(loo$interaction_p), pch = 19, col = "#2E6FB0", ylim = c(0, max(0.06, max(loo$interaction_p))),
     xlab = "Integration genes dropped (one at a time)", ylab = "Integration x treatment p",
     main = "Leave-one-gene-out: interaction p")
abline(h = 0.05, lty = 2, col = "#B03A2E"); abline(h = full_p, lty = 3, col = "grey40")
legend("topleft", c("p = 0.05", "full-panel p"), lty = c(2,3), col = c("#B03A2E","grey40"), bty="n", cex=0.85)
hist(loo$AUC_pembro, breaks = 12, col = "#A23B72", border = "white",
     xlab = "pembrolizumab-arm Integration AUC", main = "Leave-one-gene-out: pembro AUC")
abline(v = full_auc, lty = 3, col = "grey30"); par(op); dev.off()
cat("Wrote stability_leavegeneout.csv and stability_leavegeneout.png\n")
