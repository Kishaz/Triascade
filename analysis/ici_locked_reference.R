# =====================================================================
# ICI locked-reference single-patient scoring  (prototype module)
# =====================================================================
# Problem: the published ICI layer scores binarize each gene at the
# COHORT median (gene = 1 if expr >= median across the cohort), so a
# single patient cannot be scored without a cohort around them.
#
# This module freezes a per-gene reference distribution ONCE (from a
# locked reference set) and then scores any new sample, one patient at a
# time, against those frozen references. It also reproduces the original
# cohort-median scoring so the two can be compared head to head.
#
# Scoring rule (unchanged from the manuscript):
#   Sensory     = sum of its 21 genes coded 1
#   Integration = (21 positive coded 1) - (20 negative coded 1)
#   Execution   = (19 positive coded 1) - (19 negative coded 1)
#   Composite   = mean of the three layer scores, each z-scored
#   IFNG and TGFB1 sit in two layers and count in both (98 distinct genes)
#
# Expression matrices are genes (rows) x samples (columns), as in the
# pipeline's score_ici(). A gene absent from a matrix is dropped from
# that layer's counts, exactly as in the cohort-median code.
# =====================================================================

suppressPackageStartupMessages({ library(data.table) })

# ---- panel ----------------------------------------------------------
load_panel <- function(csv = "ICI_98gene_panel.csv") {
  p <- fread(csv, data.table = FALSE)
  comp <- function(layer, dir) p$gene[p$layer == layer & p$direction == dir]
  PANEL <- list(
    SENSORY = p$gene[p$layer == "Sensory"],
    INTPOS  = comp("Integration_positive", "promotes"),
    INTNEG  = comp("Integration_negative", "suppresses"),
    EXECPOS = comp("Execution_positive",   "promotes"),
    EXECNEG = comp("Execution_negative",   "suppresses"))
  stopifnot(identical(unname(lengths(PANEL)), c(21L, 21L, 20L, 19L, 19L)),
            length(unique(unlist(PANEL))) == 98L)
  PANEL
}

# ---- layer scores from a binary gene x sample matrix ----------------
.layer_scores <- function(binlist) {
  # binlist: named list of binary matrices for SENSORY/INTPOS/.../EXECNEG
  csum <- function(b) if (is.null(b) || !nrow(b)) rep(0, 0) else colSums(b, na.rm = TRUE)
  S  <- csum(binlist$SENSORY)
  I  <- csum(binlist$INTPOS) - csum(binlist$INTNEG)
  E  <- csum(binlist$EXECPOS) - csum(binlist$EXECNEG)
  data.frame(sample = names(S), Sensory = as.numeric(S),
             Integration = as.numeric(I), Execution = as.numeric(E),
             row.names = NULL, stringsAsFactors = FALSE)
}

# =====================================================================
# 1. ORIGINAL cohort-median scoring (reproduces score_ici)
# =====================================================================
score_cohort_median <- function(expr, PANEL) {
  binset <- function(gs) {
    g <- intersect(gs, rownames(expr))
    if (!length(g)) return(NULL)
    t(apply(expr[g, , drop = FALSE], 1,
            function(v) as.integer(v >= median(v, na.rm = TRUE))))
  }
  bl <- lapply(PANEL, binset)
  for (k in names(bl)) if (!is.null(bl[[k]])) colnames(bl[[k]]) <- colnames(expr)
  sc <- .layer_scores(bl)
  z <- function(x) { s <- sd(x, na.rm = TRUE); if (is.na(s) || s == 0) x * 0 else (x - mean(x, na.rm = TRUE)) / s }
  sc$Composite <- rowMeans(cbind(z(sc$Sensory), z(sc$Integration), z(sc$Execution)))
  sc$scoring <- "cohort_median"
  sc
}

# =====================================================================
# 2. LOCKED reference: freeze per-gene thresholds + layer z-score stats
# =====================================================================
fit_locked_reference <- function(expr_ref, PANEL, label = "reference") {
  genes <- intersect(unique(unlist(PANEL)), rownames(expr_ref))
  # store the sorted reference values (for quantile mapping) and median
  ref_sorted <- lapply(genes, function(g) sort(as.numeric(expr_ref[g, ]), na.last = NA))
  names(ref_sorted) <- genes
  ref_median <- vapply(ref_sorted, function(v) median(v), numeric(1))
  # per-sample-standardized per-gene median (for the affine-robust selfnorm mode):
  # z-score each reference sample ACROSS the panel genes, then take each gene's median
  M <- as.matrix(expr_ref[genes, , drop = FALSE])
  Z <- apply(M, 2, function(col) { s <- sd(col, na.rm = TRUE)
             if (is.na(s) || s == 0) col * 0 else (col - mean(col, na.rm = TRUE)) / s })
  std_median <- apply(Z, 1, median, na.rm = TRUE); names(std_median) <- genes
  # freeze the layer-score mean/sd on the reference (for Composite z-scoring)
  ref_scores <- score_cohort_median(expr_ref, PANEL)   # ref is its own cohort
  zstats <- lapply(c("Sensory", "Integration", "Execution"),
                   function(L) c(mean = mean(ref_scores[[L]]), sd = sd(ref_scores[[L]])))
  names(zstats) <- c("Sensory", "Integration", "Execution")
  structure(list(label = label, genes = genes, ref_sorted = ref_sorted,
                 ref_median = ref_median, std_median = std_median,
                 zstats = zstats, PANEL = PANEL),
            class = "ici_locked_ref")
}

# quantile of x within a sorted reference vector (fraction of ref <= x)
.ref_quantile <- function(x, sorted_ref) {
  if (is.na(x)) return(NA_real_)
  (findInterval(x, sorted_ref) ) / length(sorted_ref)
}

# Score any sample(s) against a locked reference. Works on one column
# (single patient) or a matrix. mode:
#   "quantile" (default) : code 1 if the value's quantile in the frozen
#                          reference is >= 0.5  (rank-aligned, platform-robust)
#   "direct"             : code 1 if value >= frozen reference median
#                          (same-platform / same-normalization only)
# mode:
#   "quantile" (default): code 1 if the value's quantile in the frozen reference
#                         is >= 0.5  (rank-aligned; same expression platform)
#   "direct"            : code 1 if value >= frozen reference median
#                         (same platform AND same normalization only)
#   "selfnorm"          : z-score the sample ACROSS the panel genes, then code 1
#                         if that standardized value >= the frozen standardized
#                         reference median. Removes per-sample affine (platform)
#                         shifts, so it is the most portable single-patient mode.
score_locked <- function(expr_new, ref, mode = c("quantile", "direct", "selfnorm")) {
  mode <- match.arg(mode)
  if (is.null(dim(expr_new))) expr_new <- matrix(expr_new, ncol = 1,
      dimnames = list(names(expr_new), "patient"))
  PANEL <- ref$PANEL
  # for selfnorm: standardize each new sample across the panel genes present
  if (mode == "selfnorm") {
    gp <- intersect(ref$genes, rownames(expr_new))
    sub <- as.matrix(expr_new[gp, , drop = FALSE])
    expr_new[gp, ] <- apply(sub, 2, function(col) { s <- sd(col, na.rm = TRUE)
                        if (is.na(s) || s == 0) col * 0 else (col - mean(col, na.rm = TRUE)) / s })
  }
  call_gene <- function(g, vals) {
    if (!(g %in% ref$genes)) return(rep(NA_integer_, length(vals)))
    if (mode == "direct")   return(as.integer(vals >= ref$ref_median[[g]]))
    if (mode == "selfnorm") return(as.integer(vals >= ref$std_median[[g]]))
    as.integer(vapply(vals, .ref_quantile, numeric(1), sorted_ref = ref$ref_sorted[[g]]) >= 0.5)
  }
  binset <- function(gs) {
    g <- intersect(gs, rownames(expr_new))
    g <- intersect(g, ref$genes)
    if (!length(g)) return(NULL)
    ns <- ncol(expr_new)
    res <- vapply(g, function(gi) call_gene(gi, as.numeric(expr_new[gi, ])), numeric(ns))
    b <- if (ns == 1) matrix(res, nrow = length(g), ncol = 1) else t(res)   # genes x samples
    rownames(b) <- g; colnames(b) <- colnames(expr_new); b
  }
  bl <- lapply(PANEL, binset)
  sc <- .layer_scores(bl)
  zf <- function(x, L) { m <- ref$zstats[[L]]["mean"]; s <- ref$zstats[[L]]["sd"]
                         if (is.na(s) || s == 0) x * 0 else (x - m) / s }
  sc$Composite <- rowMeans(cbind(zf(sc$Sensory, "Sensory"),
                                 zf(sc$Integration, "Integration"),
                                 zf(sc$Execution, "Execution")))
  sc$scoring <- paste0("locked_", mode)
  sc
}

# =====================================================================
# 3. Concordance between two score tables on the SAME samples
# =====================================================================
concordance <- function(sc_a, sc_b) {
  m <- merge(sc_a, sc_b, by = "sample", suffixes = c("_a", "_b"))
  layers <- c("Sensory", "Integration", "Execution", "Composite")
  out <- lapply(layers, function(L) {
    a <- m[[paste0(L, "_a")]]; b <- m[[paste0(L, "_b")]]
    hi_a <- a >= median(a); hi_b <- b >= median(b)   # median split per scoring
    data.frame(layer = L,
               pearson = cor(a, b, use = "complete.obs"),
               spearman = cor(a, b, method = "spearman", use = "complete.obs"),
               hi_lo_agreement = mean(hi_a == hi_b))
  })
  do.call(rbind, out)
}

# 3-bit archetype (S/I/E high-low) agreement between two scorings
archetype_agreement <- function(sc_a, sc_b) {
  m <- merge(sc_a, sc_b, by = "sample", suffixes = c("_a", "_b"))
  code <- function(df, s) paste0(
    ifelse(df[[paste0("Sensory", s)]]     >= median(df[[paste0("Sensory", s)]]),     "S+", "S-"),
    ifelse(df[[paste0("Integration", s)]] >= median(df[[paste0("Integration", s)]]), "I+", "I-"),
    ifelse(df[[paste0("Execution", s)]]   >= median(df[[paste0("Execution", s)]]),   "E+", "E-"))
  mean(code(m, "_a") == code(m, "_b"))
}
