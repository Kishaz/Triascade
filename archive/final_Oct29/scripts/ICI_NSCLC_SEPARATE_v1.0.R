# =====================================================================
SCRIPT_VERSION <- "ICI_NSCLC_SEPARATE  v1.0  2026-10-28"
# Reports the two NSCLC studies as separate cohorts - no combined n=43 row.
#  1. per-cohort AUC table with 13 analyses:
#     Oct 16/source_data/percohort_AUC_13.csv
#     (the NSCLC n=43 rows replaced by GSE126044 and GSE135222 from
#      Supplementary Table S3; every other row unchanged)
#  2. prints how many of the 52 estimates have a 95% CI that includes 0.5
#     (and the same count for the old 12-analysis table, which should be the
#     paper's 27 of 48), and the dominant layer of each analysis
#  3. redraws the per-cohort AUC supplementary figure with 13 rows
#     (file named S5 if you already ran the renumbering script, else S2)
#  4. lists your scripts that still treat NSCLC as one cohort (archetypes /
#     Figure 2B, Cohen's d / Supp Fig S9, NSCLC survival) - upload those
# =====================================================================
suppressPackageStartupMessages(library(data.table))
OCT <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript/Oct 16"
OCT <- Sys.getenv("ICI_OCT16", OCT)
SD  <- file.path(OCT, "source_data"); SF <- file.path(OCT, "supplementary_figures")
cat("\n", SCRIPT_VERSION, "\n", sep = "")

T2 <- fread(file.path(SD, "SuppFigS2_source_FINAL.csv"), data.table = FALSE)
T2 <- T2[, c("cohort", "layer", "n", "n_pos", "AUC", "CI_low", "CI_high")]
T2$cohort <- trimws(T2$cohort)
if (sum(tolower(T2$cohort) == "nsclc") != 4) stop("expected 4 combined NSCLC rows in SuppFigS2_source_FINAL.csv - paste this to Claude")
cat("\ncombined NSCLC rows removed:\n"); print(T2[tolower(T2$cohort) == "nsclc", ], row.names = FALSE)

S3 <- fread(file.path(OCT, "supplementary_tables", "Supp_Table_S3_NSCLC_split.csv"), data.table = FALSE)
stopifnot(nrow(S3) == 8, setequal(S3$study, c("GSE126044", "GSE135222")))
add <- data.frame(cohort = paste("nsclc", S3$study), layer = S3$layer, n = S3$n, n_pos = S3$responders,
                  AUC = S3$AUC, CI_low = S3$CI_low, CI_high = S3$CI_high)
T13 <- rbind(T2[tolower(T2$cohort) != "nsclc", ], add)

ORD <- c("imvigor", "rose", "ispy2 [pembro]", "ispy2 [control]", "gide [pembro]", "gide [ipi+pembro]",
         "hugo", "eac", "riaz", "nsclc gse126044", "nsclc gse135222", "friedrich", "braun")
LAY <- c("Sensory", "Integration", "Execution", "Composite")
T13 <- T13[order(match(tolower(T13$cohort), ORD), match(T13$layer, LAY)), ]
stopifnot(!any(is.na(match(tolower(T13$cohort), ORD))))
fwrite(T13, file.path(SD, "percohort_AUC_13.csv"))

cat("\nanalyses:", length(unique(T13$cohort)), "| patients:", sum(tapply(T13$n, T13$cohort, max)), "\n")
inc <- function(d) sum(d$CI_low <= 0.5 & d$CI_high >= 0.5)
cat("estimates whose 95% CI includes 0.5 - old 12-analysis table:", inc(T2), "of", nrow(T2),
    "| 13 analyses:", inc(T13), "of", nrow(T13), "\n")

dom <- do.call(rbind, lapply(split(T13, factor(tolower(T13$cohort), levels = ORD)), function(d) {
  a <- setNames(d$AUC, d$layer)[c("Sensory", "Integration", "Execution")]
  data.frame(analysis = d$cohort[1], Sensory = a[1], Integration = a[2], Execution = a[3],
             dominant = if (max(a) >= 0.55) names(a)[which.max(a)] else "None") }))
cat("\ndominant layer (highest of the three, if AUC >= 0.55):\n"); print(dom, row.names = FALSE)
print(table(factor(dom$dominant, levels = c("Sensory", "Integration", "Execution", "None"))))
cat("\nNSCLC rows now:\n"); print(T13[grepl("^nsclc", T13$cohort), ], row.names = FALSE)

# ---------------- per-cohort AUC figure, 13 rows -----------------------------
LAB <- c("imvigor" = "IMvigor210", "rose" = "Rose BLCA", "ispy2 [pembro]" = "I-SPY2 pembrolizumab",
         "ispy2 [control]" = "I-SPY2 paclitaxel only", "gide [pembro]" = "Gide pembrolizumab",
         "gide [ipi+pembro]" = "Gide ipi+pembro", "hugo" = "Hugo", "eac" = "Esophageal (EAC)",
         "riaz" = "Riaz nivolumab", "nsclc gse126044" = "NSCLC GSE126044", "nsclc gse135222" = "NSCLC GSE135222",
         "friedrich" = "Friedrich myeloma", "braun" = "Braun ccRCC")
col <- c("#8E44AD", "#1F6FB4", "#2E8B57", "black")
coh <- unique(T13$cohort)
draw <- function() {
  par(mar = c(4.5, 13, 1, 1)); y0 <- rev(seq_along(coh)) * 5
  plot(NA, xlim = c(0.2, 1), ylim = c(0, max(y0) + 3), yaxt = "n", ylab = "", xlab = "AUC for response (DeLong 95% CI)")
  abline(v = 0.5, lty = 2, col = "grey60")
  for (j in seq_along(coh)) for (l in seq_along(LAY)) {
    r <- T13[T13$cohort == coh[j] & T13$layer == LAY[l], ]; if (!nrow(r)) next
    yy <- y0[j] + 1.5 - l
    segments(r$CI_low, yy, r$CI_high, yy, col = col[l], lwd = 1.5); points(r$AUC, yy, pch = 19, col = col[l], cex = 0.8) }
  nn <- sapply(coh, function(c) max(T13$n[T13$cohort == c]))
  axis(2, at = y0, labels = sprintf("%s (n=%d)", LAB[tolower(coh)], nn), las = 1, cex.axis = 0.7, tick = FALSE)
  legend("bottomleft", LAY, col = col, pch = 19, bty = "n", cex = 0.75, horiz = TRUE) }
renumbered <- file.exists(file.path(SF, "SuppFigS2_PRISMA_flow.png")) || file.exists(file.path(SF, "SuppFigS5_percohort_AUC.png"))
stem <- file.path(SF, if (renumbered) "SuppFigS5_percohort_AUC" else "SuppFigS2_percohort_AUC")
h <- 1.5 + 0.42 * length(coh)
pdf(paste0(stem, ".pdf"), 7.5, h); draw(); invisible(dev.off())
png(paste0(stem, ".png"), 7.5, h, units = "in", res = 300); draw(); invisible(dev.off())
cat("\nper-cohort AUC figure redrawn:", basename(stem), "(.pdf/.png)",
    if (renumbered) "" else "- old numbering; the renumbering script will make it S5", "\n")
# ---------------- 4. your scripts that still use NSCLC as one cohort -----------
BASE <- dirname(dirname(OCT))                               # .../R data FINAL
rs <- list.files(BASE, pattern = "\\.R$", recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
rs <- rs[!grepl("/Oct 16/|/datasets/", rs)]
txt <- lapply(rs, function(f) tryCatch(readLines(f, warn = FALSE), error = function(e) character(0)))
has <- function(p) vapply(txt, function(s) any(grepl(p, s, ignore.case = TRUE)), logical(1))
mine <- vapply(txt, function(s) any(grepl("^SCRIPT_VERSION <-", s)), logical(1))   # scripts Claude wrote
find <- data.frame(use = ifelse(has("Ataxia"), "archetypes / Fig 2B",
                         ifelse(has("cohen"), "Cohen's d / Supp S9",
                         ifelse(has("GSE135222|nsclc") & has("Surv\\("), "NSCLC survival / Supp S7E", NA))),
                   script = sub(paste0(BASE, "/"), "", rs, fixed = TRUE), modified = format(file.mtime(rs), "%Y-%m-%d %H:%M"))
find <- find[!is.na(find$use) & !mine, ]; find <- find[order(find$use, find$modified, decreasing = c(FALSE, TRUE), method = "radix"), ]
cat("\n==== your scripts that still treat NSCLC as one cohort - please upload the newest of each ====\n")
if (nrow(find)) print(find, row.names = FALSE) else cat("none found under", BASE, "\n")

cat("\nDONE - paste everything above to Claude.\n")
