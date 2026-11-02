# =====================================================================
SCRIPT_VERSION <- "ICI_SPLIT_RERUN  v1.0  2026-10-28"
# Reruns three of YOUR scripts so the paper treats the two NSCLC studies
# separately and uses I-SPY2 at n = 248 (246 with Integration):
#   1. build_binary_KM.R             -> Supp Fig S6A-D (Braun) and S6E (NSCLC GSE135222)
#   2. build_archetype_concordance.R -> archetype numbers, Figure 2B by analysis, Supp Fig S5
#   3. build_waterfall_ridgeline.R   -> Cohen's d, Supp Fig S8A
# Your scripts are read and run with only these changes: the NSCLC entry split
# into GSE126044 and GSE135222 (per-study scores file), the I-SPY2 file set to
# the n = 248 file, output folders moved to Oct 16, and two wrong labels fixed
# (the NSCLC curves are progression-free survival, not overall survival; the
# old "n = 971" subtitle). Scripts 2 and 3 are first run unchanged except the
# output folder, to check they reproduce the numbers now in the paper.
# Your MASTER_GRID files are not changed. Old S5 and S8A files are moved to
# supplementary_figures/_removed. Full output of each run: source_data/rerun_split_v1
# =====================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("ICI_BASE", "/Users/dbedi/Desktop/R data FINAL")
MAN  <- file.path(BASE, "ICI_Manuscript"); GRID <- file.path(MAN, "MASTER_GRID")
OCT  <- file.path(MAN, "Oct 16"); SF <- file.path(OCT, "supplementary_figures"); FG <- file.path(OCT, "figures")
RUN  <- file.path(OCT, "source_data", "rerun_split_v1"); dir.create(RUN, showWarnings = FALSE, recursive = TRUE)
OLD  <- file.path(SF, "_removed"); dir.create(OLD, showWarnings = FALSE)
cat("\n", SCRIPT_VERSION, "\n", sep = "")

ISPY   <- "ICI_binary_ispy2_scores_246complete.csv"     # I-SPY2 n = 248, 246 with Integration
NSPLIT <- "ICI_binary_nsclc_scores_by_study.csv"         # NSCLC scored per study (ICI_NSCLC_SPLIT v1.2)
for (f in c(ISPY, NSPLIT)) if (!file.exists(file.path(GRID, f))) stop("not found in MASTER_GRID: ", f)
stopifnot(setequal(unique(fread(file.path(GRID, NSPLIT), data.table = FALSE)$cohort), c("GSE126044", "GSE135222")))

# run one of your scripts with listed text changes; each change must match exactly n times
run_yours <- function(file, subs, tag) {
  txt <- paste(readLines(file.path(MAN, file), warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (BASE != "/Users/dbedi/Desktop/R data FINAL")                       # test runs only
    subs <- c(list(c('"/Users/dbedi/Desktop/R data FINAL"', deparse(BASE), 1)), subs)
  for (s in subs) {
    k <- lengths(regmatches(txt, gregexpr(s[1], txt, fixed = TRUE)))
    if (k != as.integer(s[3])) { cat(sprintf("  STOP: %s - expected %s match(es) of [%s], found %d\n", file, s[3], s[1], k)); return(NULL) }
    txt <- gsub(s[1], s[2], txt, fixed = TRUE)
  }
  log <- file.path(RUN, paste0(tag, ".log")); zz <- file(log, open = "wt")
  sink(zz); sink(zz, type = "message")
  e <- new.env()
  res <- tryCatch({ eval(parse(text = txt, keep.source = FALSE, encoding = "UTF-8"), envir = e); NULL },
                  error = function(err) err)
  sink(type = "message"); sink(); close(zz)
  if (!is.null(res)) {
    cat("  ERROR in", file, "(", tag, "):", conditionMessage(res), "\n  last lines of its output:\n")
    cat(paste("   ", tail(readLines(log, warn = FALSE), 15)), sep = "\n"); return(NULL) }
  cat("  ran", file, "-", tag, "\n")
  e
}
retire <- function(stem) {                    # move an old figure aside (never deleted)
  for (f in list.files(SF, pattern = paste0("^", stem, "\\.(pdf|png)$"), full.names = TRUE)) {
    to <- file.path(OLD, sub("\\.(pdf|png)$", "__NSCLC_combined.\\1", basename(f)))
    if (file.rename(f, to)) cat("  moved old", basename(f), "to _removed\n") }
}
put <- function(from_dir, from_stem, to_dir, to_stem) {
  for (x in c("pdf", "png")) { src <- file.path(from_dir, paste0(from_stem, ".", x))
    if (file.exists(src)) file.copy(src, file.path(to_dir, paste0(to_stem, ".", x)), overwrite = TRUE) }
  cat("  written:", to_stem, "(.pdf/.png) in", basename(to_dir), "\n")
}
outdir <- function(old, tag) c(old, sprintf('OUT_DIR  <- file.path(MAN_DIR, "Oct 16", "source_data", "rerun_split_v1", "%s")', tag), 1)

# ======================= 1. survival (Supp Fig S6A-E) =======================
cat("\n==== 1. survival curves: Braun (S6A-D) and NSCLC GSE135222 (S6E) ====\n")
km <- run_yours("build_binary_KM.R", list(
  outdir('OUT_DIR  <- file.path(GRID_DIR, "Figures_v2", "KM_panels")', "KM"),
  c('write.csv(hr_df, file.path(GRID_DIR, "HR_summary_binary.csv")', 'write.csv(hr_df, file.path(OUT_DIR, "HR_summary_binary.csv")', 1),
  c('"NSCLC", "PFS Time (months)"', '"NSCLC GSE135222", "PFS Time (months)"', 4),
  c('"NSCLC — Binary KM by ICI layers (PFS)"', '"NSCLC GSE135222 — Binary KM by ICI layers (PFS)"', 1),
  c('"Braun ccRCC", "NSCLC", "EAC")', '"Braun ccRCC", "NSCLC GSE135222", "EAC")', 1),
  c('ylab = "Overall Survival")', 'ylab = if (grepl("PFS", km_result$time_label)) "Progression-free survival" else "Overall survival")', 1)),
  "KM")
if (!is.null(km)) {
  H <- km$hr_df[km$hr_df$cohort %in% c("Braun ccRCC", "NSCLC GSE135222"),
                c("cohort", "layer", "hr", "ci_lo", "ci_hi", "p", "n_high", "n_low")]
  H[, 3:6] <- lapply(H[, 3:6], signif, 3); print(H, row.names = FALSE)
  cat("  paper: Braun HR 0.77-0.85, all p > 0.3 | NSCLC Integration HR 0.49 (0.21-1.18), p = 0.11\n")
  put(file.path(RUN, "KM"), "KM_Braun_binary", SF, "SuppFigS6AD_Braun_KM")
  put(file.path(RUN, "KM"), "KM_NSCLC_binary", SF, "SuppFigS6E_NSCLC_GSE135222_KM")
}

# ======================= 2. archetypes (Figure 2B, Supp Fig S5) =============
cat("\n==== 2. archetypes ====\n")
A_OUT <- 'OUT_DIR  <- file.path(GRID_DIR, "Figures_v2", "archetypes")'
A_SUBS <- list(
  c('"ICI_binary_ispy2_scores.csv"', paste0('"', ISPY, '"'), 2),
  c('subtitle = "All 12 cohorts pooled (n = 971)",', 'subtitle = sprintf("All patients (n = %d)", nrow(archetype_df)),', 1),
  c('file.path(GRID_DIR, "archetype_', 'file.path(OUT_DIR, "archetype_', 3))
A_SPLIT <- list(
  c('"NSCLC"           = list(file = "ICI_binary_nsclc_scores.csv", tumor = "NSCLC",',
    paste0('"NSCLC GSE126044" = list(file = "', NSPLIT, '", tumor = "NSCLC", cont_by_study = TRUE,\n',
           '                           filter = function(d) d[d$cohort == "GSE126044", ],\n',
           '                           cont_csv = file.path(MAN_DIR, "bulk_scores", "NSCLC_ICI_scores.csv")),\n',
           '  "NSCLC GSE135222" = list(file = "', NSPLIT, '", tumor = "NSCLC", cont_by_study = TRUE,\n',
           '                           filter = function(d) d[d$cohort == "GSE135222", ],'), 1),
  c('"Riaz SKCM", "NSCLC", "EAC", "Friedrich")', '"Riaz SKCM", "NSCLC GSE126044", "NSCLC GSE135222", "EAC", "Friedrich")', 1),
  c('if (is.na(cont_id)) next',
    'if (is.na(cont_id)) next\n  if (isTRUE(spec$cont_by_study)) cont <- cont[cont[[cont_id]] %in% archetype_df$id[archetype_df$cohort == cn], , drop = FALSE]', 1))

LEV <- c("Coordinated", "Multi-Layer Ataxia", "Sensory-Blind", "Spatially-Excluded", "Execution-Failed", "Mixed")
arch_numbers <- function(e) {
  A <- e$archetype_df; a <- factor(as.character(A$archetype_binary), levels = LEV)
  N <- data.frame(archetype = LEV, n = as.vector(table(a)), pct = round(100 * as.vector(prop.table(table(a))), 1),
                  resp_pct = round(100 * as.vector(tapply(A$resp, a, mean)), 1))
  if (!is.null(e$conc_df)) { C <- e$conc_df; ca <- factor(C$archetype_continuous, levels = LEV)
    N$agree_pct <- round(100 * as.vector(tapply(as.character(C$archetype_binary) == C$archetype_continuous, ca, mean))) }
  N }
per_analysis <- function(e) {
  A <- e$archetype_df; A$a <- as.character(A$archetype_binary)
  do.call(rbind, lapply(split(A, factor(A$cohort, levels = unique(A$cohort))), function(d) {
    co <- d$resp[d$a == "Coordinated"]; at <- d$resp[d$a == "Multi-Layer Ataxia"]
    data.frame(analysis = d$cohort[1], n = nrow(d),
               coord_n = length(co), coord_resp = if (length(co)) round(100 * mean(co)) else NA,
               ataxia_n = length(at), ataxia_resp = if (length(at)) round(100 * mean(at)) else NA,
               coord_pct = round(100 * mean(d$a == "Coordinated")), ataxia_pct = round(100 * mean(d$a == "Multi-Layer Ataxia"))) })) }
conc_line <- function(e) if (!is.null(e$conc_df)) {
  C <- e$conc_df; m <- as.character(C$archetype_binary) == C$archetype_continuous
  cat(sprintf("  continuous vs binary agreement: %.1f%% of %d patients |", 100 * mean(m), nrow(C)))
  cat(paste(sprintf(" %s %.0f%% (n=%d)", names(tapply(m, C$cohort, mean)), 100 * tapply(m, C$cohort, mean),
                    as.vector(table(C$cohort)[names(tapply(m, C$cohort, mean))])), collapse = ";"), "\n") }

a12 <- run_yours("build_archetype_concordance.R", c(A_SUBS, list(outdir(A_OUT, "archetypes_check"))), "archetypes_check")
if (!is.null(a12)) {
  N <- arch_numbers(a12)
  N$paper_pct  <- c(30.1, 21.0, 9.0, 7.0, 6.8, 26.1); N$paper_resp <- c(38.8, 13.4, 32.3, 20.8, 37.1, 24.3)
  if (!is.null(N$agree_pct)) N$paper_agree <- c(71, 70, 41, 50, 21, 41)
  cat("\n  CHECK - NSCLC as one cohort, as in the paper now (patients:", nrow(a12$archetype_df), "; paper 1,028):\n")
  print(N, row.names = FALSE); conc_line(a12); cat("  paper: agreement 49.6% of 224 patients\n")
  P <- per_analysis(a12)
  cat("  Coordinated response higher than Ataxia in", sum(P$coord_resp > P$ataxia_resp, na.rm = TRUE), "of", nrow(P), "analyses\n")
}
a13 <- run_yours("build_archetype_concordance.R", c(A_SUBS, list(outdir(A_OUT, "archetypes_13")), A_SPLIT), "archetypes_13")
if (!is.null(a13)) {
  cat("\n  NEW - NSCLC as two studies (patients:", nrow(a13$archetype_df), "):\n")
  print(arch_numbers(a13), row.names = FALSE); conc_line(a13)
  P <- per_analysis(a13)
  cat("\n  per analysis (response % in Coordinated vs Multi-Layer Ataxia; prevalence %):\n"); print(P, row.names = FALSE)
  hi <- P$coord_resp > P$ataxia_resp
  cat("  Coordinated response higher than Ataxia in", sum(hi, na.rm = TRUE), "of", nrow(P), "analyses;",
      "equal in", sum(P$coord_resp == P$ataxia_resp, na.rm = TRUE), "; lower in", sum(P$coord_resp < P$ataxia_resp, na.rm = TRUE),
      "; a group empty in", sum(is.na(hi)), "\n")
  fwrite(P, file.path(RUN, "archetypes_by_analysis_13.csv"))
  put(file.path(RUN, "archetypes_13"), "archetype_x_cohort_response_binary", FG, "Figure2B_archetype_by_analysis_13")
  retire("SuppFigS5_archetype_concordance")
  put(file.path(RUN, "archetypes_13"), "archetype_continuous_vs_binary_confusion", SF, "SuppFigS5_archetype_concordance")
}

# ======================= 3. Cohen's d (Supp Fig S8A) ========================
cat("\n==== 3. Cohen's d ====\n")
W_OUT  <- 'OUT_DIR  <- file.path(GRID_DIR, "Figures_v2", "waterfall_ridgeline")'
W_SUBS <- list(c('write.csv(delta_df, file.path(GRID_DIR, "waterfall_delta_summary.csv")',
                 'write.csv(delta_df, file.path(OUT_DIR, "waterfall_delta_summary.csv")', 1))
W_SPLIT <- list(
  c('"ICI_binary_ispy2_scores.csv"', paste0('"', ISPY, '"'), 2),
  c('"NSCLC"           = list(file = "ICI_binary_nsclc_scores.csv", tumor = "NSCLC"),',
    paste0('"NSCLC GSE126044" = list(file = "', NSPLIT, '", tumor = "NSCLC", filter = function(d) d[d$cohort == "GSE126044", ]),\n',
           '  "NSCLC GSE135222" = list(file = "', NSPLIT, '", tumor = "NSCLC", filter = function(d) d[d$cohort == "GSE135222", ]),'), 1))
dtab <- function(e) { D <- e$delta_df[e$delta_df$layer == "Composite", c("cohort", "n", "n_R", "cohens_d", "p")]
  D$cohort <- as.character(D$cohort); D$cohens_d <- round(D$cohens_d, 2); D$p <- signif(D$p, 2); D[order(-D$cohens_d), ] }
w12 <- run_yours("build_waterfall_ridgeline.R", c(W_SUBS, list(outdir(W_OUT, "cohens_d_check"))), "cohens_d_check")
if (!is.null(w12)) {
  D <- dtab(w12)
  pd <- c("Gide-Pembro" = 1.31, "I-SPY2 pembro" = 1.00, "NSCLC" = 0.97, "EAC" = 0.69, "Rose BLCA" = 0.53, "IMvigor210" = 0.36, "Braun ccRCC" = 0.14)
  pp <- c("Gide-Pembro" = 7.5e-5, "I-SPY2 pembro" = 7.1e-6, "NSCLC" = 0.003, "EAC" = 0.04, "Rose BLCA" = 0.05, "IMvigor210" = 0.012, "Braun ccRCC" = 0.29)
  D$paper_d <- pd[D$cohort]; D$paper_p <- pp[D$cohort]
  cat("\n  CHECK - script as it is (Composite):\n"); print(D, row.names = FALSE)
}
w13 <- run_yours("build_waterfall_ridgeline.R", c(W_SUBS, list(outdir(W_OUT, "cohens_d_13")), W_SPLIT), "cohens_d_13")
if (!is.null(w13)) {
  cat("\n  NEW - I-SPY2 n = 248, NSCLC as two studies (Composite):\n"); print(dtab(w13), row.names = FALSE)
  fwrite(w13$delta_df, file.path(RUN, "cohens_d_all_layers_13.csv"))
  retire("SuppFigS8A_cohens_d")
  put(file.path(RUN, "cohens_d_13"), "waterfall_cohens_d_binary", SF, "SuppFigS8A_cohens_d")
}
cat("\nDONE - paste everything above to Claude.\n")
