# =====================================================================
SCRIPT_VERSION <- "ICI_PURITY_SPLIT  v1.0  2026-10-28"
# Reruns YOUR MASTER_GRID/05_tumor_purity_sensitivity.R with the two Gide arms
# (and the two NSCLC studies) as separate cohorts, as everywhere else in the
# paper. Your script runs unchanged except its data folder, which is a copy of
# your MASTER_GRID score files plus the split files, made in
# Oct 16/source_data/rerun_split_v1/purity. The combined Gide file stays in the
# copy as a check: it must reproduce the paper (n = 57; Sensory OR 8.74,
# Integration 2.44, Composite 3.16). Your MASTER_GRID files are not changed.
# =====================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("ICI_BASE", "/Users/dbedi/Desktop/R data FINAL")
MAN  <- file.path(BASE, "ICI_Manuscript"); GRID <- file.path(MAN, "MASTER_GRID")
PD   <- file.path(MAN, "Oct 16", "source_data", "rerun_split_v1", "purity"); dir.create(PD, showWarnings = FALSE, recursive = TRUE)
cat("\n", SCRIPT_VERSION, "\n", sep = "")
SCRIPT <- file.path(GRID, "05_tumor_purity_sensitivity.R")
for (f in c(SCRIPT, file.path(GRID, "per_cohort_tumor_purity.csv"))) if (!file.exists(f)) stop("not found: ", f)

# copy of the score files your script reads, plus the split ones
src <- list.files(GRID, pattern = "^ICI_binary_.*_scores\\.csv$", full.names = TRUE)
invisible(file.copy(c(src, file.path(GRID, "per_cohort_tumor_purity.csv")), PD, overwrite = TRUE))
g <- fread(file.path(GRID, "ICI_binary_gide_scores.csv"), data.table = FALSE)
cat("Gide arms:", paste(names(table(g$arm)), table(g$arm), sep = " n=", collapse = ", "), "\n")
fwrite(g[g$arm == "Pembro", ], file.path(PD, "ICI_binary_gidePembro_scores.csv"))
fwrite(g[g$arm == "Ipi+Pembro", ], file.path(PD, "ICI_binary_gideIpiPembro_scores.csv"))
ns <- fread(file.path(GRID, "ICI_binary_nsclc_scores_by_study.csv"), data.table = FALSE)
for (s in c("GSE126044", "GSE135222")) fwrite(ns[ns$cohort == s, ], file.path(PD, paste0("ICI_binary_nsclc", s, "_scores.csv")))

txt <- paste(readLines(SCRIPT, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
old <- 'DATA_DIR <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript/MASTER_GRID/"'
if (lengths(regmatches(txt, gregexpr(old, txt, fixed = TRUE))) != 1) stop("DATA_DIR line not found as expected - paste this to Claude")
txt <- gsub(old, paste0("DATA_DIR <- ", deparse(paste0(PD, "/"))), txt, fixed = TRUE)

log <- file.path(PD, "purity_run.log"); zz <- file(log, open = "wt"); sink(zz); sink(zz, type = "message")
res <- tryCatch({ eval(parse(text = txt, keep.source = FALSE, encoding = "UTF-8"), envir = new.env()); NULL },
                error = function(e) e)
sink(type = "message"); sink(); close(zz)
L <- readLines(log, warn = FALSE)
if (!is.null(res)) { cat("ERROR:", conditionMessage(res), "\n"); cat(tail(L, 20), sep = "\n"); stop("paste this to Claude") }

# print the Gide and NSCLC blocks of your script's output
hdr <- grep("^=== ", L)
for (h in hdr) {
  nm <- sub("^=== (\\S+) .*", "\\1", L[h])
  if (nm %in% c("gide", "gidePembro", "gideIpiPembro", "nsclc", "nsclcGSE126044", "nsclcGSE135222")) {
    end <- c(hdr[hdr > h], grep("^Saved |^Skipped", L))[1] - 1; if (is.na(end)) end <- length(L)
    blk <- L[h:end]; cat("\n"); cat(blk[nzchar(blk) & !grepl("^\\[", blk)], sep = "\n") }
}
cat("\npaper (both Gide arms together, n = 57, 31 responders): Sensory OR 8.74 (p < 0.001), Integration 2.44 (p = 0.013), Composite 3.16 (p = 0.005)\n")
cat("\nfull output:", log, "\nDONE - paste everything above to Claude.\n")
