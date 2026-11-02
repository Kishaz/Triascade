# =====================================================================
# Locked-reference single-patient scoring  -  SYNTHETIC DEMONSTRATION
# =====================================================================
# This runs on SIMULATED data only, to show the mechanics and quantify
# how well locked-reference single-patient scoring reproduces the
# published cohort-median scoring. It produces no real results.
# On the real cohorts, point `expr_ref` / `expr_new` at the pipeline's
# expression matrices (genes x samples) and read concordance the same way.
# =====================================================================
set.seed(2026)
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")
PANEL <- load_panel("ICI_98gene_panel.csv")
genes <- unique(unlist(PANEL))

# ---- generative model: immune activity + a stromal/vascular barrier ----
# Per-gene baselines are SHARED across cohorts (same biology / same platform);
# the "shift" cohort additionally applies a per-gene affine transform to mimic
# a different expression platform or normalization.
MU   <- setNames(rnorm(length(genes), 6, 1.5), genes)   # shared gene baselines
LOAD <- 1.1
sim_cohort <- function(n, platform = c("ref", "shift")) {
  platform <- match.arg(platform)
  imm     <- rnorm(n)                      # overall immune activity
  barrier <- 0.3 * imm + rnorm(n, 0, 0.95) # stromal/TGF-beta barrier
  E <- sapply(genes, function(g) {
    lat <- if (g %in% c(PANEL$INTNEG, PANEL$EXECNEG)) barrier else imm
    MU[g] + LOAD * lat + rnorm(n, 0, 0.8)
  })
  E <- t(E); rownames(E) <- genes; colnames(E) <- paste0("S", seq_len(n))
  if (platform == "shift") {               # mimic a different platform: per-gene affine
    a <- runif(length(genes), 0.5, 2.0); b <- rnorm(length(genes), 0, 3)
    E <- E * a + b
  }
  attr(E, "imm") <- imm; attr(E, "barrier") <- barrier
  E
}

expr_ref   <- sim_cohort(300, "ref")
expr_same  <- sim_cohort(250, "ref")       # new cohort, SAME platform as reference
expr_shift <- sim_cohort(250, "shift")     # new cohort, DIFFERENT platform
# response: higher immune access (imm - barrier) -> higher response probability
mk_y <- function(E) { acc <- attr(E, "imm") - 0.8 * attr(E, "barrier")
                      rbinom(ncol(E), 1, plogis(1.2 * scale(acc)[, 1] - 0.3)) }
y_same  <- mk_y(expr_same); y_shift <- mk_y(expr_shift)

ref <- fit_locked_reference(expr_ref, PANEL, label = "sim_reference_n300")
cat(sprintf("Locked reference frozen on %d genes (%d samples).\n\n",
            length(ref$genes), 300L))

# ---- helper: binary gene calls (gene x sample) for agreement ----------
`%||%` <- function(a, b) if (is.null(a)) b else a
bincalls <- function(expr, ref = NULL, mode = NULL) {
  g <- intersect(ref$genes %||% rownames(expr), rownames(expr))
  if (is.null(ref)) {   # cohort-median
    b <- t(apply(expr[g, , drop = FALSE], 1, function(v) as.integer(v >= median(v, na.rm = TRUE))))
  } else {
    sc <- expr
    if (mode == "selfnorm") sc[g, ] <- apply(sc[g, , drop = FALSE], 2,
      function(c) { s <- sd(c); if (s == 0) c * 0 else (c - mean(c)) / s })
    b <- t(sapply(g, function(gi) {
      v <- as.numeric(sc[gi, ])
      if (mode == "direct")   as.integer(v >= ref$ref_median[[gi]])
      else if (mode == "selfnorm") as.integer(v >= ref$std_median[[gi]])
      else as.integer(vapply(v, .ref_quantile, numeric(1), sorted_ref = ref$ref_sorted[[gi]]) >= 0.5)
    }))
  }
  rownames(b) <- g; colnames(b) <- colnames(expr); b
}

evaluate <- function(expr, y, platform_label) {
  cm <- score_cohort_median(expr, PANEL)
  bc_cm <- bincalls(expr)
  rows <- list()
  for (mode in c("quantile", "direct", "selfnorm")) {
    lk <- score_locked(expr, ref, mode = mode)
    con <- concordance(cm, lk)
    arche <- archetype_agreement(cm, lk)
    bc_lk <- bincalls(expr, ref, mode)
    gg <- intersect(rownames(bc_cm), rownames(bc_lk))
    gene_agree <- mean(bc_cm[gg, ] == bc_lk[gg, ], na.rm = TRUE)
    auc_cm <- as.numeric(auc(y, cm$Integration, quiet = TRUE, direction = "<"))
    auc_lk <- as.numeric(auc(y, lk$Integration, quiet = TRUE, direction = "<"))
    rows[[mode]] <- data.frame(
      platform = platform_label, locked_mode = mode,
      gene_call_agreement = round(gene_agree, 3),
      Integration_spearman = round(con$spearman[con$layer == "Integration"], 3),
      Composite_spearman   = round(con$spearman[con$layer == "Composite"], 3),
      archetype_agreement  = round(arche, 3),
      AUC_Integration_cohortmedian = round(auc_cm, 3),
      AUC_Integration_locked       = round(auc_lk, 3))
  }
  do.call(rbind, rows)
}

res <- rbind(evaluate(expr_same, y_same, "same_platform"),
             evaluate(expr_shift, y_shift, "shifted_platform"))
rownames(res) <- NULL
cat("CONCORDANCE: locked single-patient scoring vs cohort-median scoring\n")
cat(strrep("=", 70), "\n"); print(res); cat("\n")

# ---- single-patient property: one patient scored ALONE == in a matrix --
one <- expr_same[, 1, drop = FALSE]
alone  <- score_locked(one, ref, "selfnorm")
inbatch <- score_locked(expr_same, ref, "selfnorm")[1, ]
same_one <- isTRUE(all.equal(alone$Integration, inbatch$Integration)) &&
            isTRUE(all.equal(alone$Sensory, inbatch$Sensory)) &&
            isTRUE(all.equal(alone$Execution, inbatch$Execution))
cat(sprintf("Single-patient check: patient scored alone == in cohort?  %s\n\n",
            if (same_one) "YES (scores identical)" else "NO"))

write.csv(res, "demo_concordance_table.csv", row.names = FALSE)

# ---- figure: locked (selfnorm) vs cohort-median Integration, same platform
cm <- score_cohort_median(expr_same, PANEL)
lk <- score_locked(expr_same, ref, "selfnorm")
png("demo_locked_vs_cohortmedian.png", width = 1500, height = 700, res = 170)
op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.4, 1), family = "sans")
plot(cm$Integration, lk$Integration, pch = 19, col = "#2E6FB055",
     xlab = "Integration (cohort-median)", ylab = "Integration (locked, single-patient)",
     main = sprintf("Same platform  (Spearman %.2f)",
                    cor(cm$Integration, lk$Integration, method = "spearman")))
abline(0, 1, col = "grey40", lty = 2)
sub <- res[res$platform == "shifted_platform", ]
barplot(sub$gene_call_agreement, names.arg = sub$locked_mode, ylim = c(0, 1),
        col = "#A23B72", ylab = "Gene-call agreement", main = "Shifted platform: agreement by mode")
abline(h = 1, col = "grey70", lty = 3)
par(op); dev.off()
cat("Wrote demo_concordance_table.csv and demo_locked_vs_cohortmedian.png\n")
