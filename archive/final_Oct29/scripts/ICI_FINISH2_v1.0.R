# =====================================================================
SCRIPT_VERSION <- "ICI_FINISH2  v1.0  2026-10-16"
# 1. Supp Table S4: your ispy2_categorical_check.R on the 248 cohort, minus
#    the 2 controls with no CXCL9 (Integration NA) - same rule as everywhere else.
# 2. Supp Fig S2: your per-cohort AUC table (FIGURE2B/layer_auc_with_ci.csv)
#    redrawn with the I-SPY2 rows taken from the n=248 run (Figure3AB_source.csv).
# =====================================================================
suppressPackageStartupMessages(library(data.table))
MAN <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"; GRID <- file.path(MAN, "MASTER_GRID")
OCT <- file.path(MAN, "Oct 16"); SD <- file.path(OCT, "source_data")
cat("\n", SCRIPT_VERSION, "\n", sep = "")

# ---------------- 1. S4 ------------------------------------------------
d <- fread(file.path(GRID, "ICI_binary_ispy2_scores_248.csv"))
d246 <- d[!is.na(Integration)]
f246 <- file.path(GRID, "ICI_binary_ispy2_scores_246complete.csv"); fwrite(d246, f246)
cat("S4 input:", nrow(d246), "patients (2 with no CXCL9 removed)\n")
src <- readLines(file.path(MAN, "ispy2_categorical_check.R"), warn = FALSE)
src <- gsub("ICI_binary_ispy2_scores\\.csv", "ICI_binary_ispy2_scores_246complete.csv", src)
src <- gsub("CATEGORICAL", "CATEGORICAL_248", src)
tryCatch(eval(parse(text = src), envir = new.env()), error = function(e) cat("ERROR:", conditionMessage(e), "\n"))
for (o in list.files(file.path(MAN, "ISPY2/CATEGORICAL_248"), pattern = "\\.csv$", full.names = TRUE)) {
  cat("\n--", basename(o), "\n"); print(fread(o), digits = 3)
  file.copy(o, file.path(OCT, "supplementary_tables", paste0("Supp_Table_S4_n248__", basename(o))), overwrite = TRUE) }
old <- list.files(file.path(OCT, "supplementary_tables"), pattern = "n189_RERUN", full.names = TRUE); file.remove(old)

# ---------------- 2. S2 ------------------------------------------------
T2 <- fread(file.path(MAN, "FIGURE2B/layer_auc_with_ci.csv"), data.table = FALSE)
AB <- fread(file.path(SD, "Figure3AB_source.csv"), data.table = FALSE)
for (a in c("control", "drug")) {
  arm <- if (a == "drug") "pembro" else "control"
  for (i in seq_len(nrow(AB))) {
    L <- AB$layer[i]; k <- which(grepl("ispy2", T2$cohort) & T2$arm == arm & T2$layer == L)
    T2$n[k] <- AB[[paste0("n_", a)]][i]; T2$n_pos[k] <- AB[[paste0("ev_", a)]][i]
    T2$n_neg[k] <- T2$n[k] - T2$n_pos[k]
    T2$AUC[k] <- round(AB[[paste0("AUC_", a)]][i], 3); T2$CI_low[k] <- round(AB[[paste0("AUClo_", a)]][i], 3)
    T2$CI_high[k] <- round(AB[[paste0("AUChi_", a)]][i], 3)
  }
}
T2$crosses_0.5 <- T2$CI_low <= 0.5
fwrite(T2, file.path(SD, "SuppFigS2_source_n248.csv"))
print(T2[grepl("ispy2", T2$cohort), c("cohort", "layer", "n", "n_pos", "AUC", "CI_low", "CI_high")])
LAY <- c("Sensory", "Integration", "Execution", "Composite"); col <- c("#8E44AD", "#1F6FB4", "#2E8B57", "black")
coh <- unique(T2$cohort)
drawS2 <- function() {
  par(mar = c(4.5, 13, 1, 1)); y0 <- rev(seq_along(coh)) * 5
  plot(NA, xlim = c(0.25, 1), ylim = c(2, max(y0) + 3), yaxt = "n", ylab = "", xlab = "AUC for response (DeLong 95% CI)")
  abline(v = 0.5, lty = 2, col = "grey60")
  for (j in seq_along(coh)) for (l in seq_along(LAY)) {
    r <- T2[T2$cohort == coh[j] & T2$layer == LAY[l], ]; if (!nrow(r)) next
    yy <- y0[j] + 1.5 - l
    segments(r$CI_low, yy, r$CI_high, yy, col = col[l], lwd = 1.5); points(r$AUC, yy, pch = 19, col = col[l], cex = 0.8)
  }
  nn <- sapply(coh, function(c) T2$n[T2$cohort == c][1])
  axis(2, at = y0, labels = sprintf("%s (n=%d)", coh, nn), las = 1, cex.axis = 0.7, tick = FALSE)
  legend("bottomleft", LAY, col = col, pch = 19, bty = "n", cex = 0.75, horiz = TRUE)
}
h <- 1.5 + 0.42 * length(coh)
pdf(file.path(OCT, "supplementary_figures/SuppFigS2_percohort_AUC.pdf"), 7.5, h); drawS2(); invisible(dev.off())
png(file.path(OCT, "supplementary_figures/SuppFigS2_percohort_AUC.png"), 7.5, h, units = "in", res = 300); drawS2(); invisible(dev.off())
cat("\nSupp S2 redrawn with I-SPY2 on n=248. DONE - paste the output to Claude.\n")
