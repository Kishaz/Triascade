# =====================================================================
SCRIPT_VERSION <- "ICI_FIG3_FIG4_LRT  v1.0  2026-10-28"
# Redraws Figure 3 (I-SPY2, n = 248) and Figure 4 (NeoTRIP) with
# likelihood-ratio p-values, so the figures and the text use one test.
# Runs the script that drew them, ICI_FIGS_v12_v1.1.R, with these exact changes:
#   1. panels A and B: likelihood-ratio p instead of Wald p (both figures)
#   2. panel C, tertile x arm: likelihood-ratio test on the three tertiles,
#      the same test as Supplementary Table S4
#   3. Figure 3E, lowest layer: layers z-scored over the 246 patients with all
#      three layers, as in your ispy2_categorical_check.R (Table S4)
#   4. Figure 3D: wider left margin, so the "HR+/HER2-" label is not cut off
#   5. output folder: Oct 16/source_data/rerun_split_v1/fig3_fig4_LRT
# Everything else (data file, ORs, AUCs, gating) is unchanged.
# New figures go to Oct 16/figures (old ones to figures/_removed); source CSVs
# to Oct 16/source_data (old ones to source_data/_removed). Prints each number
# the text uses next to its current value. MASTER_GRID files are not changed.
# =====================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("ICI_BASE", "/Users/dbedi/Desktop/R data FINAL")
MAN  <- file.path(BASE, "ICI_Manuscript"); OCT <- file.path(MAN, "Oct 16")
FIG  <- file.path(OCT, "figures"); SD <- file.path(OCT, "source_data")
W    <- file.path(SD, "rerun_split_v1", "fig3_fig4_LRT"); dir.create(W, showWarnings = FALSE, recursive = TRUE)
cat("\n", SCRIPT_VERSION, "\n", sep = "")
SRC <- file.path(MAN, "ICI_FIGS_v12_v1.1.R")
if (!file.exists(SRC)) {
  s <- list.files(MAN, pattern = "^ICI_FIGS_v12_v1\\.1\\.R$", recursive = TRUE, full.names = TRUE)
  s <- s[!grepl("/Oct 16/", s, fixed = TRUE)]
  if (!length(s)) stop("ICI_FIGS_v12_v1.1.R not found - paste this to Claude")
  SRC <- s[1] }
cat("figure script:", SRC, "\n")

txt <- paste(readLines(SRC, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
k4 <- gregexpr("# FIGURE 4", txt, fixed = TRUE)[[1]]
if (length(k4) != 1 || k4[1] < 0) stop("Figure 4 section not found as expected - paste this to Claude")
part <- list(fig3 = substr(txt, 1, k4 - 1), fig4 = substr(txt, k4, nchar(txt)))

subs3 <- list(
  c('OUT       <- file.path(BASE, "Figures_v12")', paste0("OUT       <- ", deparse(W)), 1),
  c('fmtp(S$p_wald)', 'fmtp(S$p_lrt)', 2),
  c('p <- coef(summary(glm(y ~ Tn * A, binomial, d)))["Tn:Adrug", 4]',
    'p <- anova(glm(y ~ T + A, binomial, d), glm(y ~ T * A, binomial, d), test = "LRT")[2, "Pr(>Chi)"]', 1),
  c('z <- scale(d[, c("Sensory", "Integration", "Execution")])',
    paste0('cc3 <- complete.cases(d[, c("Sensory", "Integration", "Execution")]); z <- matrix(NA_real_, nrow(d), 3); ',
           'z[cc3, ] <- scale(as.matrix(d[cc3, c("Sensory", "Integration", "Execution")]))'), 1),
  c('  panel_or(rec, "Paclitaxel + pembrolizumab", "D", "Integration OR per SD (95% CI)")',
    '  par(mar = c(4.2, 7.5, 2.6, 1)); panel_or(rec, "Paclitaxel + pembrolizumab", "D", "Integration OR per SD (95% CI)"); par(mar = c(4.2, 6, 2.6, 1))', 1))
if (BASE != "/Users/dbedi/Desktop/R data FINAL")                          # test runs only
  subs3 <- c(list(c('BASE      <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"', paste0("BASE      <- ", deparse(MAN)), 1)), subs3)
for (s in subs3) {
  k <- lengths(regmatches(part$fig3, gregexpr(s[1], part$fig3, fixed = TRUE)))
  if (k != as.integer(s[3])) stop(sprintf("expected %s match(es) of [%s], found %d - paste this to Claude", s[3], s[1], k))
  part$fig3 <- gsub(s[1], s[2], part$fig3, fixed = TRUE) }

e <- new.env()
run <- function(code, tag) {
  log <- file.path(W, paste0(tag, ".log")); zz <- file(log, open = "wt"); sink(zz); sink(zz, type = "message")
  res <- tryCatch({ eval(parse(text = code, keep.source = FALSE, encoding = "UTF-8"), envir = e); NULL }, error = function(err) err)
  sink(type = "message"); sink(); close(zz)
  if (!is.null(res)) { cat("ERROR in", tag, ":", conditionMessage(res), "\n"); cat(tail(readLines(log, warn = FALSE), 12), sep = "\n") }
  is.null(res) }
if (!run(part$fig3, "figure3")) stop("Figure 3 not redrawn - paste this to Claude")
ok4 <- run(part$fig4, "figure4")
if (!ok4) cat("Figure 4 NOT redrawn (error above); Figure 3 is done.\n")

# ---- into Oct 16 (old versions kept in _removed) --------------------------------
swap <- function(new, dir) {
  dir.create(file.path(dir, "_removed"), showWarnings = FALSE)
  old <- file.path(dir, basename(new))
  if (file.exists(old)) file.rename(old, file.path(dir, "_removed", sub("(\\.[A-Za-z]+)$", "__wald\\1", basename(new))))
  file.copy(new, old) }
fig_new <- list.files(W, pattern = paste0("^(Figure3_ISPY2_n248", if (ok4) "|Figure4_NeoTRIP", ")\\.(pdf|png)$"), full.names = TRUE)
src_new <- list.files(W, pattern = paste0("^Figure3[A-F]+_source\\.csv$", if (ok4) "|^Figure4[A-D]+_source\\.csv$"), full.names = TRUE)
invisible(lapply(fig_new, swap, dir = FIG)); invisible(lapply(src_new, swap, dir = SD))

# ---- numbers for the text ---------------------------------------------------------
f3 <- function(x) sprintf("%.3f", x); f2 <- function(x) sprintf("%.2f", x); pc <- function(x, n) sprintf("%.0f%%", 100 * x / n)
d <- e$d; S3 <- e$S3
cat(sprintf("\n==== FIGURE 3, I-SPY2: %d patients (control %d, pembro %d) ====\n", nrow(d), sum(d$arm == "control"), sum(d$arm == "drug")))
cat("receptor:", paste(names(table(d$receptor)), table(d$receptor), collapse = ", "),
    sprintf("| pCR pembro %s, control %s   (text: 114 TNBC, 134 HR+; 45%%, 17%%)\n",
            pc(sum(d$y[d$arm == "drug"]), sum(d$arm == "drug")), pc(sum(d$y[d$arm == "control"]), sum(d$arm == "control"))))
now <- c(Sensory = "3.10 vs 1.54; p = 0.068", Integration = "3.16 vs 1.31; p = 0.025",
         Execution = "2.24 vs 1.94; p = 0.71", Composite = "4.00 vs 1.61; p = 0.032")
cat("\nper SD, layer x arm (3B):  n | OR pembro (95% CI) | OR control (95% CI) | LRT p | text now\n")
for (i in seq_len(nrow(S3))) with(S3[i, ], cat(sprintf("  %-11s %3d | %s (%s-%s) | %s (%s-%s) | %s | %s\n", layer, n_control + n_drug,
  f2(OR_drug), f2(ORlo_drug), f2(ORhi_drug), f2(OR_control), f2(ORlo_control), f2(ORhi_control), f3(p_lrt), now[layer])))
cat("AUC pembro / control (3A):", paste0(S3$layer, " ", f3(S3$AUC_drug), "/", f3(S3$AUC_control), collapse = ", "),
    "\n   (text: 0.754/0.616, 0.777/0.574, 0.687/0.673, 0.804/0.635)\n")
tt <- e$T3$tab
cat(sprintf("\ntertiles (3C): pembro %s, control %s | LRT p = %s   (text: 11/46/68%%, 14/16/20%%; p = 0.018; Table S4 0.034)\n",
            paste(pc(tt$x[tt$arm == "drug"], tt$n[tt$arm == "drug"]), collapse = "/"),
            paste(pc(tt$x[tt$arm == "control"], tt$n[tt$arm == "control"]), collapse = "/"), f3(e$T3$p)))
rc <- e$rec
cat("receptor (3D):\n")
for (i in seq_len(nrow(rc))) with(rc[i, ], cat(sprintf("  %-18s interaction OR %s (95%% CI %s-%s) | LRT p = %s\n", layer, f2(int_OR), f2(int_lo), f2(int_hi), f3(p_lrt))))
cat("   (text: TNBC n = 112, OR 4.45 (1.32-14.98), p = 0.016; HR+ n = 134, OR 1.22, p = 0.71)\n")
lt <- e$low_tab
cat(sprintf("\nlowest layer (3E), Sensory/Integration/Execution: control %s, pembro %s | n = %d | LRT p = %s\n",
            paste(pc(lt$x[lt$arm == "control"], lt$n[lt$arm == "control"]), collapse = "/"),
            paste(pc(lt$x[lt$arm == "drug"], lt$n[lt$arm == "drug"]), collapse = "/"), sum(lt$n), f3(e$p_low)))
cat("   (text: control 16/23/11%, pembro 41/33/59%, p = 0.062; Table S4 0.065)\n")
g <- e$gate_tab
cat("\ngating, pembrolizumab arm (3F):\n")
for (i in seq_len(nrow(g))) with(g[i, ], cat(sprintf("  %-20s %2d/%2d = %.0f%% (%.0f-%.0f%%)\n", gate, pcr, n, pct, lo, hi)))
cat("   (text: all 45% (34-57); Sensory-high 63% of 35 (46-77); + Integration-high 75% of 28 (57-87);\n          Sensory-low 26% of 34 (15-43); + Execution-high 70% of 23 (49-84))\n")

P <- data.frame(figure = "3", panel = c(rep("A,B", nrow(S3)), "C", rep("D", nrow(rc)), "E"),
                test = c(paste(S3$layer, "x arm, per SD"), "Integration tertile x arm", paste("Integration x arm,", rc$layer), "lowest layer x arm"),
                LRT_p = signif(c(S3$p_lrt, e$T3$p, rc$p_lrt, e$p_low), 3))
if (ok4) {
  S4 <- e$S4
  cat(sprintf("\n==== FIGURE 4, NeoTRIP: %d patients (chemotherapy %d, atezolizumab %d) ====\n", nrow(e$ph), sum(e$ph$arm == "control"), sum(e$ph$arm == "drug")))
  cat("layer       | OR atezo | OR chemo | interaction OR (95% CI) | LRT p | AUC atezo/chemo\n")
  for (i in seq_len(nrow(S4))) with(S4[i, ], cat(sprintf("  %-11s | %s | %s | %s (%s-%s) | %s | %s/%s\n", layer,
    f2(OR_drug), f2(OR_control), f2(int_OR), f2(int_lo), f2(int_hi), f3(p_lrt), f3(AUC_drug), f3(AUC_control))))
  cat(sprintf("tertiles (4C): LRT p = %s\n", f3(e$T4$p)))
  cat("   (text: Integration 2.03 vs 1.15, 1.76 (1.02-3.05), p = 0.040; AUC 0.689 vs 0.547;\n          Sensory 0.94, p = 0.82; Execution 1.13, p = 0.64; Composite 1.28, p = 0.35)\n")
  P <- rbind(P, data.frame(figure = "4", panel = c(rep("A,B,D", nrow(S4)), "C"),
                           test = c(paste(S4$layer, "x arm, per SD"), "Integration tertile x arm"), LRT_p = signif(c(S4$p_lrt, e$T4$p), 3)))
}
fwrite(P, file.path(SD, "Figure3_Figure4_pvalues_LRT.csv"))
cat("\nwritten:", paste(basename(fig_new), collapse = ", "), "-> Oct 16/figures (old ones in figures/_removed)\n")
cat("p-values: Oct 16/source_data/Figure3_Figure4_pvalues_LRT.csv\nDONE - paste everything above to Claude.\n")
