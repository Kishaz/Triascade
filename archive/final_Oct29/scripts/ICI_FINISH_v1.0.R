# =====================================================================
SCRIPT_VERSION <- "ICI_FINISH  v1.0  2026-10-16"
# Everything still open for the paper, in one run. Nothing new is analysed;
# each part either reads a file you already made, runs one of your own
# scripts with only the input file swapped to the n=248 cohort, or plots
# numbers already stated in the manuscript.
#   1. print Fig 3C/3D/3E source numbers (tertile pCR, TNBC vs HR+, lowest layer)
#   2. Supp Table S4: YOUR ispy2_categorical_check.R on the 248 file
#   3. Supp Table S3: DeLong CIs for the two NSCLC studies
#   4. archetype check: Coordinated vs Multi-Layer Ataxia response, per cohort
#   5. Supp Fig S2: per-cohort AUC table with I-SPY2 rows from the 248 run
#   6. Supp Fig S1: regrouped AUCs, drawn from the values in the text
# Paste the whole console output back to Claude.
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
MAN  <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"
GRID <- file.path(MAN, "MASTER_GRID")
OCT  <- file.path(MAN, "Oct 16")
SD   <- file.path(OCT, "source_data")
cat("\n", SCRIPT_VERSION, "\n", sep = "")
hdr <- function(t) cat("\n", strrep("=", 62), "\n", t, "\n", strrep("=", 62), "\n", sep = "")

# ---------------- 1 ---------------------------------------------------
hdr("1. Figure 3C / 3D / 3E source numbers (n=248)")
for (f in c("Figure3C_source.csv", "Figure3D_source.csv", "Figure3E_source.csv")) {
  p <- file.path(SD, f); if (!file.exists(p)) p <- file.path(MAN, "Figures_v12", f)
  cat("\n--", f, "\n"); if (file.exists(p)) print(fread(p), digits = 3) else cat("   not found\n")
}

# ---------------- 2 ---------------------------------------------------
hdr("2. Supp Table S4: your ispy2_categorical_check.R on the 248 file")
cc <- list.files(MAN, pattern = "^ispy2_categorical_check.*\\.R$", recursive = TRUE, full.names = TRUE)
cc <- cc[!grepl("/Oct 16/", cc)]
if (length(cc)) {
  f <- cc[order(file.mtime(cc), decreasing = TRUE)][1]; cat("using", f, "\n")
  src <- readLines(f, warn = FALSE)
  n_hit <- sum(grepl("ICI_binary_ispy2_scores\\.csv", src))
  src <- gsub("ICI_binary_ispy2_scores\\.csv", "ICI_binary_ispy2_scores_248.csv", src)
  src <- gsub("CATEGORICAL", "CATEGORICAL_248", src)
  cat("score-file references swapped:", n_hit, "\n")
  if (n_hit == 0) cat("   script does not name the score file directly - S4 NOT rerun\n") else
    tryCatch(eval(parse(text = src), envir = new.env()), error = function(e) cat("   ERROR:", conditionMessage(e), "\n"))
  out <- list.files(file.path(MAN, "ISPY2/CATEGORICAL_248"), pattern = "\\.csv$", full.names = TRUE)
  for (o in out) { cat("\n--", basename(o), "\n"); print(fread(o), digits = 3)
    file.copy(o, file.path(OCT, "supplementary_tables", paste0("Supp_Table_S4_n248__", basename(o))), overwrite = TRUE) }
} else cat("ispy2_categorical_check.R not found\n")

# ---------------- 3 ---------------------------------------------------
hdr("3. Supp Table S3: NSCLC split with DeLong 95% CIs")
ns <- fread(file.path(GRID, "ICI_binary_nsclc_scores.csv"), data.table = FALSE)
cat("columns:", paste(names(ns), collapse = ", "), "\n")
rc  <- intersect(c("response_bin", "responder", "response"), names(ns))[1]
cmp <- intersect(c("composite", "ICI_composite", "Composite"), names(ns))[1]
L3 <- c("Sensory", "Integration", "Execution", cmp)
S3 <- do.call(rbind, lapply(split(ns, ns$cohort), function(s) {
  s <- s[!is.na(s[[rc]]), ]
  do.call(rbind, lapply(L3, function(L) {
    r <- roc(s[[rc]], s[[L]], levels = c(0, 1), direction = "<", quiet = TRUE); ci <- as.numeric(ci.auc(r, method = "delong"))
    data.frame(study = s$cohort[1], n = nrow(s), responders = sum(s[[rc]] == 1), layer = L,
               AUC = round(ci[2], 3), lo = round(ci[1], 3), hi = round(ci[3], 3)) })) }))
print(S3)
fwrite(S3, file.path(OCT, "supplementary_tables/Supp_Table_S3_NSCLC_split.csv"))

# ---------------- 4 ---------------------------------------------------
hdr("4. Archetypes: Coordinated vs Multi-Layer Ataxia response rate per cohort")
fs <- list.files(GRID, pattern = "^ICI_binary_.*_scores\\.csv$", full.names = TRUE)
fs <- fs[!grepl("_248|primecut|TCGA|PAAD", fs, ignore.case = TRUE)]
A <- do.call(rbind, lapply(fs, function(f) {
  d <- fread(f, data.table = FALSE)
  rc <- intersect(c("response_bin", "pcr_binary", "responder", "response"), names(d))[1]
  if (is.na(rc) || !all(c("Sensory", "Integration", "Execution") %in% names(d))) return(NULL)
  hi <- sapply(c("Sensory", "Integration", "Execution"), function(L) d[[L]] > median(d[[L]], na.rm = TRUE))
  d$arch <- ifelse(rowSums(hi) == 3, "Coordinated", ifelse(rowSums(!hi) == 3, "Ataxia", "other"))
  grp <- if ("arm" %in% names(d)) d$arm else rep("all", nrow(d))
  do.call(rbind, lapply(split(seq_len(nrow(d)), grp), function(ix) {
    s <- d[ix, ]; s <- s[!is.na(s[[rc]]), ]
    co <- s[which(s$arch == "Coordinated"), rc]; at <- s[which(s$arch == "Ataxia"), rc]
    data.frame(cohort = sub("ICI_binary_(.*)_scores\\.csv", "\\1", basename(f)), arm = grp[ix][1],
               n = nrow(s), n_coord = length(co), resp_coord = round(100 * mean(co), 1),
               n_ataxia = length(at), resp_ataxia = round(100 * mean(at), 1)) }))
}))
A$coord_higher <- A$resp_coord > A$resp_ataxia
print(A)
cat("\nCoordinated > Ataxia in", sum(A$coord_higher, na.rm = TRUE), "of", sum(!is.na(A$coord_higher)), "cohort-by-arm analyses\n")
fwrite(A, file.path(SD, "archetype_by_cohort_check.csv"))

# ---------------- 5 ---------------------------------------------------
hdr("5. Supp Fig S2: per-cohort AUC table, I-SPY2 rows replaced by the 248 run")
t2 <- list.files(MAN, pattern = "^layer_auc_with_ci\\.csv$", recursive = TRUE, full.names = TRUE)
t2 <- t2[!grepl("/Oct 16/", t2)]
ab <- file.path(SD, "Figure3AB_source.csv"); if (!file.exists(ab)) ab <- file.path(MAN, "Figures_v12/Figure3AB_source.csv")
if (length(t2) && file.exists(ab)) {
  T2 <- fread(t2[1], data.table = FALSE); cat("table:", t2[1], "\ncolumns:", paste(names(T2), collapse = ", "), "\n")
  print(head(T2[grepl("ispy", apply(T2, 1, paste, collapse = " "), ignore.case = TRUE), ], 10))
  fwrite(fread(ab), file.path(SD, "SuppFigS2_ISPY2_rows_n248.csv"))
  cat("I-SPY2 n=248 AUC rows saved to source_data/SuppFigS2_ISPY2_rows_n248.csv\n")
} else cat("layer_auc_with_ci.csv or Figure3AB_source.csv not found\n")

# ---------------- 6 ---------------------------------------------------
hdr("6. Supp Fig S1: regrouped AUCs (values as stated in the text)")
S1 <- data.frame(
  grouping = c(rep("Indication", 6), rep("Agent", 4), rep("Target", 3)),
  group = c("Urothelial", "Melanoma", "Esophageal adeno.", "NSCLC", "Renal", "Breast (pembrolizumab arm)",
            "Atezolizumab", "Nivolumab", "Pembrolizumab", "PD-(L)1 overall",
            "PD-1", "PD-1 or PD-L1", "PD-L1"),
  top_layer = c("Integration", "Sensory", "Sensory", "Execution", "Execution", "Integration",
                "Integration", "Execution", "Sensory", "Sensory",
                "Sensory", "Sensory", "Integration"),
  AUC = c(0.63, 0.71, 0.80, 0.76, 0.56, 0.78, 0.62, 0.57, 0.78, 0.67, 0.64, 0.67, 0.62))
col <- c(Sensory = "#8E44AD", Integration = "#1F6FB4", Execution = "#2E8B57")
drawS1 <- function() {
  par(mar = c(4.5, 13, 1, 1)); y <- rev(seq_len(nrow(S1))) + c(rep(2, 6), rep(1, 4), rep(0, 3))
  plot(S1$AUC, y, xlim = c(0.45, 0.9), pch = 19, cex = 1.6, col = col[S1$top_layer], yaxt = "n",
       ylab = "", xlab = "AUC of the best-performing layer")
  abline(v = 0.5, lty = 2, col = "grey60"); segments(0.5, y, S1$AUC, y, col = col[S1$top_layer])
  axis(2, at = y, labels = paste0(S1$group, "  [", S1$grouping, "]"), las = 1, cex.axis = 0.75, tick = FALSE)
  legend("bottomright", names(col), col = col, pch = 19, bty = "n", cex = 0.8)
}
dir.create(file.path(OCT, "supplementary_figures"), showWarnings = FALSE)
pdf(file.path(OCT, "supplementary_figures/SuppFigS1_regrouped_AUC.pdf"), 7, 5); drawS1(); invisible(dev.off())
png(file.path(OCT, "supplementary_figures/SuppFigS1_regrouped_AUC.png"), 7, 5, units = "in", res = 300); drawS1(); invisible(dev.off())
fwrite(S1, file.path(SD, "SuppFigS1_source.csv"))
cat("Supp S1 written. Breast uses the pembrolizumab arm only (n=248 run).\n")
cat("\nDONE - paste everything above to Claude.\n")
