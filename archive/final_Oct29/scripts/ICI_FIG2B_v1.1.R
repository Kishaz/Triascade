# =====================================================================
SCRIPT_VERSION <- "ICI_FIG2B  v1.1  2026-10-29"
# Redraws Figure 2B (response rate by archetype in each analysis) with YOUR
# build_archetype_concordance.R, run exactly as for the current Figure 2B
# (I-SPY2 n = 248, 246 with Integration; the two NSCLC studies separate).
# Layout changes only; the numbers in the cells do not change:
#   - rows in the same order and with the same names and n as Figure 2A (v2.5);
#     the names are your script's own (Rose BLCA, Gide-Pembro, Hugo SKCM ...)
#   - archetype names written across the top, so they no longer run into the title
#   - a cell with no patients is grey with "-" (white is the colour for ~30% response)
#   - plain title
# Old Figure 2B moved to figures/_removed. MASTER_GRID files are not changed.
# =====================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("ICI_BASE", "/Users/dbedi/Desktop/R data FINAL")
MAN  <- file.path(BASE, "ICI_Manuscript"); GRID <- file.path(MAN, "MASTER_GRID")
OCT  <- file.path(MAN, "Oct 16"); FG <- file.path(OCT, "figures")
W    <- file.path(OCT, "source_data", "rerun_split_v1", "fig2b_v1"); dir.create(W, showWarnings = FALSE, recursive = TRUE)
cat("\n", SCRIPT_VERSION, "\n", sep = "")
ISPY   <- "ICI_binary_ispy2_scores_246complete.csv"
NSPLIT <- "ICI_binary_nsclc_scores_by_study.csv"
SRC    <- file.path(MAN, "build_archetype_concordance.R")
for (f in c(file.path(GRID, c(ISPY, NSPLIT)), SRC)) if (!file.exists(f)) stop("not found: ", f)

lab_code <- paste(
  'LAB2 <- c("IMvigor210", "Rose BLCA", "I-SPY2 control", "I-SPY2 pembro", "Gide-Pembro", "Gide-IpiPembro",',
  '          "Hugo SKCM", "Riaz SKCM", "NSCLC GSE126044", "NSCLC GSE135222", "EAC", "Braun ccRCC", "Friedrich")',
  'names(LAB2) <- LAB2',
  'stopifnot(setequal(names(LAB2), unique(archetype_df$cohort)))',
  'LAB2[] <- paste0(LAB2, " (n=", as.vector(table(archetype_df$cohort)[names(LAB2)]), ")")',
  'arch_cohort_resp <- tidyr::complete(arch_cohort_resp, cohort = names(LAB2), archetype_binary, fill = list(n = 0L))',
  sep = "\n")
subs <- list(
  # the run behind the current Figure 2B (ICI_SPLIT_RERUN v1.0)
  c('"ICI_binary_ispy2_scores.csv"', paste0('"', ISPY, '"'), 2),
  c('subtitle = "All 12 cohorts pooled (n = 971)",', 'subtitle = sprintf("All patients (n = %d)", nrow(archetype_df)),', 1),
  c('file.path(GRID_DIR, "archetype_', 'file.path(OUT_DIR, "archetype_', 3),
  c('OUT_DIR  <- file.path(GRID_DIR, "Figures_v2", "archetypes")', paste0("OUT_DIR  <- ", deparse(W)), 1),
  c('"NSCLC"           = list(file = "ICI_binary_nsclc_scores.csv", tumor = "NSCLC",',
    paste0('"NSCLC GSE126044" = list(file = "', NSPLIT, '", tumor = "NSCLC", cont_by_study = TRUE,\n',
           '                           filter = function(d) d[d$cohort == "GSE126044", ],\n',
           '                           cont_csv = file.path(MAN_DIR, "bulk_scores", "NSCLC_ICI_scores.csv")),\n',
           '  "NSCLC GSE135222" = list(file = "', NSPLIT, '", tumor = "NSCLC", cont_by_study = TRUE,\n',
           '                           filter = function(d) d[d$cohort == "GSE135222", ],'), 1),
  c('"Riaz SKCM", "NSCLC", "EAC", "Friedrich")', '"Riaz SKCM", "NSCLC GSE126044", "NSCLC GSE135222", "EAC", "Friedrich")', 1),
  c('if (is.na(cont_id)) next',
    'if (is.na(cont_id)) next\n  if (isTRUE(spec$cont_by_study)) cont <- cont[cont[[cont_id]] %in% archetype_df$id[archetype_df$cohort == cn], , drop = FALSE]', 1),
  # layout of the heatmap only
  c('arch_cohort_resp$cohort <- factor(arch_cohort_resp$cohort,',
    paste0(lab_code, '\narch_cohort_resp$cohort <- factor(LAB2[as.character(arch_cohort_resp$cohort)],'), 1),
  c('levels = rev(cohort_order_lvl))', 'levels = rev(unname(LAB2)))', 1),
  c('100 * response_rate, n)),', '100 * response_rate, n) %>% replace(n == 0, "-")),', 1),
  c('midpoint = 0.3, limits = c(0, 1),', 'midpoint = 0.3, limits = c(0, 1), na.value = "grey90",', 1),
  c('scale_x_discrete(position = "top") +',
    paste0('scale_x_discrete(position = "top", labels = c("Coordinated" = "Coordinated", "Spatially-Excluded" = "Spatially-\\nExcluded", ',
           '"Sensory-Blind" = "Sensory-\\nBlind", "Execution-Failed" = "Execution-\\nFailed", "Multi-Layer Ataxia" = "Multi-Layer\\nAtaxia", "Mixed" = "Mixed")) +'), 1),
  c('axis.text.x.top = element_text(angle = 35, hjust = 0, face = "bold"),',
    'axis.text.x.top = element_text(angle = 0, hjust = 0.5, face = "bold", lineheight = 0.9),', 1),
  c('ggsave(file.path(OUT_DIR, "archetype_x_cohort_response_binary.pdf"), p3,',
    paste0('p3 <- p3 + labs(title = "Response rate by archetype in each analysis",\n',
           '                subtitle = "Each cell: response rate (number of patients). Grey: no patients with that archetype.")\n',
           'ggsave(file.path(OUT_DIR, "archetype_x_cohort_response_binary.pdf"), p3,'), 1))
if (BASE != "/Users/dbedi/Desktop/R data FINAL")                          # test runs only
  subs <- c(list(c('"/Users/dbedi/Desktop/R data FINAL"', deparse(BASE), 1)), subs)

txt <- paste(readLines(SRC, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
for (s in subs) {
  k <- lengths(regmatches(txt, gregexpr(s[1], txt, fixed = TRUE)))
  if (k != as.integer(s[3])) stop(sprintf("expected %s match(es) of [%s], found %d - paste this to Claude", s[3], s[1], k))
  txt <- gsub(s[1], s[2], txt, fixed = TRUE) }
log <- file.path(W, "fig2b.log"); zz <- file(log, open = "wt"); sink(zz); sink(zz, type = "message")
e <- new.env()
res <- tryCatch({ eval(parse(text = txt, keep.source = FALSE, encoding = "UTF-8"), envir = e); NULL }, error = function(err) err)
sink(type = "message"); sink(); close(zz)
if (!is.null(res)) { cat("ERROR:", conditionMessage(res), "\n"); cat(tail(readLines(log, warn = FALSE), 15), sep = "\n"); stop("paste this to Claude") }

# ---- the numbers (unchanged) and the new figure into Oct 16 ---------------------
A <- e$archetype_df; cat(sprintf("patients: %d in %d analyses\n", nrow(A), length(unique(A$cohort))))
cell <- function(r) if (length(r)) sprintf("%.0f%% (n=%d)", 100 * mean(r), length(r)) else "-"
P <- do.call(rbind, lapply(names(e$LAB2), function(cn) { d <- A[A$cohort == cn, ]; a <- as.character(d$archetype_binary)
  data.frame(analysis = e$LAB2[[cn]], Coordinated = cell(d$resp[a == "Coordinated"]), Multi_Layer_Ataxia = cell(d$resp[a == "Multi-Layer Ataxia"])) }))
cat("response rate, Coordinated vs Multi-Layer Ataxia (unchanged from the current Figure 2B):\n")
print(P, row.names = FALSE)
dir.create(file.path(FG, "_removed"), showWarnings = FALSE)
for (x in c("pdf", "png")) {
  old <- file.path(FG, paste0("Figure2B_archetype_by_analysis_13.", x))
  dest <- file.path(FG, "_removed", paste0("Figure2B_archetype_by_analysis_13__old_layout.", x)); i <- 1
  while (file.exists(dest)) { i <- i + 1; dest <- file.path(FG, "_removed", sprintf("Figure2B_archetype_by_analysis_13__old_layout_%d.%s", i, x)) }
  if (file.exists(old)) file.rename(old, dest)
  file.copy(file.path(W, paste0("archetype_x_cohort_response_binary.", x)), old) }
cat("\nwritten: Oct 16/figures/Figure2B_archetype_by_analysis_13 (.pdf/.png); old one in figures/_removed\nDONE - paste everything above to Claude.\n")
