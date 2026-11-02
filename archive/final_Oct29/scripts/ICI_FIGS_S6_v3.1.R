# =====================================================================
SCRIPT_VERSION <- "ICI_FIGS_S6  v3.1  2026-10-27"
# v3.1: pooled rows removed - cohorts are different analyses and are not combined.
# Supp Fig S6 drawn from the numbers ALREADY REPORTED in the manuscript.
# No analysis is run. Every value below is quoted from the v12 text
# ("Integration stratifies response within already-inflamed tumors"),
# except the Braun and Riaz 95% CIs, which come from today's run of the
# 07 code and reproduced the text's ORs exactly (0.85, 1.17).
# Rose and Gide-Pembro CIs are not in the text, so they plot as points only.
# =====================================================================
OCT <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript/Oct 16"
cat("\n", SCRIPT_VERSION, "\n", sep = "")

S0 <- data.frame(
  label = c("IMvigor210, TIS-high", "IMvigor210, CD8A-high", "Rose BLCA, CD8A-high",
            "Riaz melanoma, CD8A-high", "Gide-Pembro, CD8A-high",
            "Pooled validation (Rose, Riaz, Gide)", "Pooled incl. IMvigor210",
            "Braun ccRCC, CD8A-high (negative control)"),
  n     = c(152, 149, 44, 25, 16, NA, NA, 91),
  OR    = c(1.63, 1.53, 1.62, 1.17, 7.55, 1.62, 1.56, 0.85),
  lo    = c(1.11, 1.04, NA, 0.44, NA, 0.93, 1.13, 0.55),
  hi    = c(2.44, 2.28, NA, 3.02, NA, 2.83, 2.15, 1.32),
  p     = c(0.015, 0.033, NA, NA, NA, 0.086, 0.006, 0.48),
  type  = c("cohort", "cohort", "cohort", "cohort", "cohort", "pooled", "pooled", "control"),
  stringsAsFactors = FALSE)
S <- S0[S0$type != "pooled", ]

draw <- function() {
  k <- nrow(S); y <- rev(seq_len(k)); y[S$type == "control"] <- 0.3
  par(mar = c(4.5, 17, 1.5, 6))
  plot(NA, log = "x", xlim = c(0.3, 12), ylim = c(0, k + 0.5), yaxt = "n", ylab = "",
       xlab = "Integration OR for objective response per SD (95% CI)")
  abline(v = 1, lty = 2, col = "grey60")
  col <- c(cohort = "#1F6FB4", pooled = "black", control = "grey50")[S$type]
  ok <- !is.na(S$lo)
  segments(S$lo[ok], y[ok], pmin(S$hi[ok], 12), y[ok], col = col[ok], lwd = 2)
  points(pmin(S$OR, 12), y, pch = ifelse(S$type == "pooled", 18, 19), cex = ifelse(S$type == "pooled", 1.8, 1.2), col = col)
  lab <- ifelse(is.na(S$n), S$label, sprintf("%s (n=%d)", S$label, S$n))
  axis(2, at = y, labels = lab, las = 1, cex.axis = 0.75, tick = FALSE)
  axis(4, at = y, labels = ifelse(is.na(S$p), "", paste0("p=", S$p)), las = 1, cex.axis = 0.7, tick = FALSE)
  mtext("Gide-Pembro unstable (near-complete separation); Rose and Gide CIs not reported", side = 1, line = 3.4, cex = 0.6, adj = 1)
}
pdf(file.path(OCT, "supplementary_figures/SuppFigS6_inflamed_strata.pdf"), 8, 4.6); draw(); invisible(dev.off())
png(file.path(OCT, "supplementary_figures/SuppFigS6_inflamed_strata.png"), 8, 4.6, units = "in", res = 300); draw(); invisible(dev.off())
write.csv(S, file.path(OCT, "source_data/SuppFigS6_source.csv"), row.names = FALSE)
cat("Supp S6 written from the manuscript's reported values.\n")
