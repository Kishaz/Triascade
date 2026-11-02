# =====================================================================
SCRIPT_VERSION <- "ICI_S1_CORRELATIONS  v1.0  2026-10-28"
# Supplementary Figure S1 (layer correlations) redrawn for the paper's 13
# cohort-by-arm analyses with YOUR two scripts:
#   fig2_layer_separability.R      computes Pearson r per analysis (its own rules:
#                                  fewer than 10 patients skipped, under 20 faded)
#   Oct4_ICI/scripts/fig5_S9_S10.R (S10 part) draws the figure
# Changes: only the paper's cohorts (TCGA-PAAD, PRIMECUT and Huang out); I-SPY2
# at n = 248 by arm, Gide by arm, NSCLC by study; patients with a response (the
# same 1,028 as Figure 2A); readable names; and the caption calculated from the
# data instead of the typed-in "r = 0.66". The old S1 is moved to
# supplementary_figures/_removed. Your MASTER_GRID files are not changed.
# =====================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("ICI_BASE", "/Users/dbedi/Desktop/R data FINAL")
MAN  <- file.path(BASE, "ICI_Manuscript"); GRID <- file.path(MAN, "MASTER_GRID")
OCT  <- file.path(MAN, "Oct 16"); SF <- file.path(OCT, "supplementary_figures")
W    <- file.path(OCT, "source_data", "rerun_split_v1", "S1_correlations"); IN <- file.path(W, "scores")
dir.create(IN, showWarnings = FALSE, recursive = TRUE); unlink(list.files(IN, full.names = TRUE))
dir.create(file.path(SF, "_removed"), showWarnings = FALSE)
cat("\n", SCRIPT_VERSION, "\n", sep = "")
F2 <- file.path(MAN, "fig2_layer_separability.R"); F5 <- file.path(MAN, "Oct4_ICI", "scripts", "fig5_S9_S10.R")
for (f in c(F2, F5)) if (!file.exists(f)) stop("not found: ", f)

# ---- the 13 analyses, patients with a response --------------------------------
rd   <- function(f) fread(file.path(GRID, f), data.table = FALSE)
resp <- function(d) d[[intersect(c("response_bin", "pcr_binary"), names(d))[1]]]
LABS <- integer(0)
put  <- function(d, lab) {
  d <- d[!is.na(resp(d)) & complete.cases(d[, c("Sensory", "Integration", "Execution")]), ]
  fwrite(d, file.path(IN, paste0("ICI_binary_", lab, "_scores.csv"))); LABS[lab] <<- nrow(d) }
is <- rd("ICI_binary_ispy2_scores_246complete.csv"); gi <- rd("ICI_binary_gide_scores.csv")
ns <- rd("ICI_binary_nsclc_scores_by_study.csv")
put(rd("ICI_binary_imvigor_scores.csv"), "IMvigor210");      put(rd("ICI_binary_rose_scores.csv"), "Rose")
put(is[is$arm == "pembro", ], "I-SPY2 pembrolizumab");       put(is[is$arm == "control", ], "I-SPY2 paclitaxel only")
put(gi[gi$arm == "Pembro", ], "Gide pembrolizumab");          put(gi[gi$arm == "Ipi+Pembro", ], "Gide ipi+pembro")
put(rd("ICI_binary_hugo_scores.csv"), "Hugo");               put(rd("ICI_binary_riaz_scores.csv"), "Riaz")
put(ns[ns$cohort == "GSE126044", ], "NSCLC GSE126044");       put(ns[ns$cohort == "GSE135222", ], "NSCLC GSE135222")
put(rd("ICI_binary_eac_scores.csv"), "EAC");                 put(rd("ICI_binary_friedrich_scores.csv"), "Friedrich")
put(rd("ICI_binary_braun_scores.csv"), "Braun")
cat("analyses:", length(LABS), "| patients:", sum(LABS), "\n")

# ---- run one of your scripts' text with exact replacements ----------------------
run_txt <- function(txt, subs, tag, env = new.env()) {
  if (BASE != "/Users/dbedi/Desktop/R data FINAL" && grepl('"/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"', txt, fixed = TRUE))
    subs <- c(list(c('"/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"', deparse(MAN), 1)), subs)   # test runs only
  for (s in subs) {
    k <- lengths(regmatches(txt, gregexpr(s[1], txt, fixed = TRUE)))
    if (k != as.integer(s[3])) stop(sprintf("%s: expected %s match(es) of [%s], found %d - paste this to Claude", tag, s[3], s[1], k))
    txt <- gsub(s[1], s[2], txt, fixed = TRUE) }
  log <- file.path(W, paste0(tag, ".log")); zz <- file(log, open = "wt"); sink(zz); sink(zz, type = "message")
  res <- tryCatch({ eval(parse(text = txt, keep.source = FALSE, encoding = "UTF-8"), envir = env); NULL }, error = function(e) e)
  sink(type = "message"); sink(); close(zz)
  if (!is.null(res)) { cat("ERROR in", tag, ":", conditionMessage(res), "\n"); cat(tail(readLines(log, warn = FALSE), 15), sep = "\n"); stop("paste this to Claude") }
  env }
rl <- function(f) paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")

# 1. your correlation script on the 13 analyses
invisible(run_txt(rl(F2), list(
  c('GRID   <- file.path(MAN, "MASTER_GRID")', paste0("GRID   <- ", deparse(IN)), 1),
  c('OUTDIR <- file.path(MAN, "Figures")', paste0("OUTDIR <- ", deparse(W)), 1)), "fig2_correlations"))
C <- read.csv(file.path(W, "Fig2_source_numbers.csv"), stringsAsFactors = FALSE, check.names = FALSE)
C <- C[!duplicated(paste(C$cohort, C$pair)), c("cohort", "n", "pair", "r")]
wide <- reshape(C, idvar = c("cohort", "n"), timevar = "pair", direction = "wide")
names(wide) <- sub("^r\\.", "", names(wide)); wide <- wide[order(-wide$n), ]
wide[, -(1:2)] <- round(wide[, -(1:2)], 2)
cat("\nPearson r between layer scores:\n"); print(wide, row.names = FALSE)
mx <- C[which.max(C$r), ]; gone <- setdiff(names(LABS), unique(C$cohort))
cat(sprintf("\nhighest: r = %.3f (%s, %s, n = %d)\n", mx$r, mx$pair, mx$cohort, mx$n))
if (length(gone)) cat("not shown (your script skips fewer than 10 patients):", paste0(gone, " (n = ", LABS[gone], ")", collapse = ", "), "\n")

# 2. your S10 drawing code, caption calculated
CAP <- sprintf("Highest correlation r = %.2f (%s, %s).\nPoint size is the number of patients; analyses under 20 patients faded.%s",
               mx$r, mx$pair, mx$cohort,
               if (length(gone)) paste0("\n", paste0(gone, " (n = ", LABS[gone], ")", collapse = ", "), " not shown: fewer than 10 patients.") else "")
f5 <- rl(F5); a <- regexpr("# ===================== S10", f5, fixed = TRUE); b <- regexpr('cat("\\nAll outputs in"', f5, fixed = TRUE)
if (a < 0 || b < 0) stop("S10 block not found in fig5_S9_S10.R - paste this to Claude")
e <- new.env(); e$OUT <- W; e$CAP <- CAP
invisible(run_txt(substr(f5, a, b - 1), list(
  c('CO <- file.path(MAN,"Figures","Fig2_source_numbers.csv")', paste0("CO <- ", deparse(file.path(W, "Fig2_source_numbers.csv"))), 1),
  c('caption="No layer pair exceeds r = 0.66 in any cohort. Point size is cohort size; small cohorts faded.")', "caption=CAP)", 1),
  c('scale_x_continuous(limits=c(-0.1,1), breaks=seq(0,1,0.25))', 'scale_x_continuous(limits=c(min(-0.1, floor(10*min(C$r))/10), 1), breaks=seq(-1,1,0.25))', 1)),
  "S10_figure", e))

# 3. into Oct 16
for (x in c("pdf", "png")) {
  old <- file.path(SF, paste0("SuppFigS1_layer_correlations.", x))
  if (file.exists(old)) file.rename(old, file.path(SF, "_removed", paste0("SuppFigS1_layer_correlations__old_cohorts.", x)))
  file.copy(file.path(W, paste0("FigS10_layer_correlations.", x)), old, overwrite = TRUE) }
cat("\nwritten: SuppFigS1_layer_correlations (.pdf/.png); old one moved to _removed\nDONE - paste everything above to Claude.\n")
