# =====================================================================
SCRIPT_VERSION <- "ICI_FIGS_IMVIGOR  v1.1  2026-10-16"
# Builds the two IMvigor210 panels the v12 text cites but that do not exist:
#   Figure 5F    mean z-scored layer by IHC immune phenotype, Integration split by response
#   Supp Fig S6  Integration vs CD8A / TIS / Teff / Composite within TIS-high and CD8A-high strata
# Then searches the ICI_Manuscript tree for the three old supplementary
# figures that were not found (S1 regrouped AUC, S4A-E Braun/NSCLC KM,
# S5B dominant-layer Cohen's d) and copies candidates for you to pick.
# Column names are detected, not assumed; the console says what it used.
# Expected (from the text): phenotypes desert 69 / excluded 113 / inflamed 62;
#   Sensory mean z -0.53 / -0.12 / +0.93; TIS-high n=152 Integration OR 1.63;
#   CD8A-high n=149 Integration OR 1.53.
# =====================================================================
suppressPackageStartupMessages(library(data.table))
BASE  <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"
OCT   <- file.path(BASE, "Oct 16")
SC    <- file.path(BASE, "MASTER_GRID/ICI_binary_imvigor_scores.csv")
CLIN  <- c(file.path(BASE, "Oct4_ICI/ICI-framework/data/imvigor_clinical_full.csv"),
           file.path(BASE, "Aug_1_ICI_FINAL_PAPER/ICI-framework/data/imvigor_clinical_full.csv"),
           file.path(BASE, "Aug_1_ICI_FINAL_PAPER/imvigor_clinical_full.csv"))
BENCH <- c(file.path(BASE, "benchmark_figures/IMvigor210_benchmark_scores.csv"),
           file.path(BASE, "IMvigor210_ICI_results/benchmark_figures/IMvigor210_benchmark_scores.csv"))
for (d in c("figures", "supplementary_figures", "supplementary_figures/_candidates", "source_data"))
  dir.create(file.path(OCT, d), showWarnings = FALSE, recursive = TRUE)
cat("\n", SCRIPT_VERSION, "\n", strrep("=", 62), "\n", sep = "")

first_existing <- function(v) { v <- v[file.exists(v)]; if (!length(v)) stop("none of these exist:\n", paste(v, collapse = "\n")); v[1] }
idcol <- function(df, ids) {           # the column whose values best match the score-file sample ids
  ov <- sapply(df, function(x) sum(as.character(x) %in% ids)); names(which.max(ov)) }
pick <- function(df, pat, what) {
  hit <- grep(pat, names(df), ignore.case = TRUE, value = TRUE)
  cat(sprintf("  %-11s -> %s\n", what, if (length(hit)) hit[1] else "NOT FOUND")); if (length(hit)) hit[1] else NA }

# ---------------- load and join -------------------------------------
sc <- fread(SC, data.table = FALSE); sc$sample <- as.character(sc$sample)
comp <- intersect(names(sc), c("composite", "ICI_composite", "Composite"))[1]; sc$Composite <- sc[[comp]]
cl_f <- first_existing(CLIN); bn_f <- first_existing(BENCH)
cl <- fread(cl_f, data.table = FALSE); bn <- fread(bn_f, data.table = FALSE)
cat("scores:", nrow(sc), "| clinical:", basename(cl_f), nrow(cl), "| benchmarks:", basename(bn_f), nrow(bn), "\n")
cid <- idcol(cl, sc$sample); bid <- idcol(bn, sc$sample)
cat("join keys: clinical$", cid, " (", sum(cl[[cid]] %in% sc$sample), " match) | benchmarks$", bid,
    " (", sum(bn[[bid]] %in% sc$sample), " match)\n", sep = "")
cat("columns used:\n")
# v1.1: phenotype found by VALUES (a column containing desert/excluded/inflamed), searched in
# the clinical file, then IMvigor_benchmark_covariates.csv, then the benchmark file
COV <- c(file.path(BASE, "Oct4_ICI/ICI-framework/data/IMvigor_benchmark_covariates.csv"),
         file.path(BASE, "Aug_1_ICI_FINAL_PAPER/ICI-framework/data/IMvigor_benchmark_covariates.csv"))
find_pheno <- function(df) {
  h <- names(df)[sapply(df, function(x) sum(tolower(as.character(x)) %in% c("desert", "excluded", "inflamed")) > 50)]
  if (length(h)) h[1] else NA }
ph_src <- NULL
for (cand in list(list(cl, cl_f), if (any(file.exists(COV))) { f <- COV[file.exists(COV)][1]; list(fread(f, data.table = FALSE), f) }, list(bn, bn_f))) {
  if (is.null(cand)) next
  col <- find_pheno(cand[[1]])
  if (!is.na(col)) { ph_src <- cand; c_ph <- col; break }
}
if (is.null(ph_src)) stop("No column with desert/excluded/inflamed values in clinical, covariates or benchmark files.")
cat(sprintf("  %-11s -> %s  (in %s)\n", "phenotype", c_ph, basename(ph_src[[2]])))
pdf_ <- ph_src[[1]]; pid <- idcol(pdf_, sc$sample)
b_tis <- pick(bn, "^TIS|GEP|tis_|_tis", "TIS/GEP")
b_tef <- pick(bn, "teff|t_eff", "Teff")
b_cd8 <- pick(bn, "^CD8A$|cd8a", "CD8A")
b_ifn <- pick(bn, "ifn", "IFN-gamma")

d <- merge(sc, unique(data.frame(sample = as.character(pdf_[[pid]]), phenotype = tolower(as.character(pdf_[[c_ph]])))), by = "sample", all.x = TRUE)
bsub <- data.frame(sample = as.character(bn[[bid]]))
bmap <- c(TIS = b_tis, Teff = b_tef, CD8A = b_cd8, IFNG = b_ifn); bmap <- bmap[!is.na(bmap)]
for (k in names(bmap)) bsub[[k]] <- suppressWarnings(as.numeric(bn[[bmap[[k]]]]))
d <- merge(d, bsub, by = "sample", all.x = TRUE)
d$y <- as.integer(d$response_bin)
cat("merged:", nrow(d), "patients | responders:", sum(d$y, na.rm = TRUE), "\n")

COL3 <- c(desert = "#C9C9C9", excluded = "#8FB3D9", inflamed = "#1F6FB4")
tag <- function(t) mtext(t, side = 3, adj = -0.1, font = 2, line = 0.6, cex = 1.1)

# =====================================================================
# FIGURE 5F — layers by IHC phenotype
# =====================================================================
cat("\n== FIGURE 5F ==\n")
LAY <- c("Sensory", "Integration", "Execution", "Composite")
for (L in LAY) d[[paste0("z_", L)]] <- as.numeric(scale(d[[L]]))
dp <- d[d$phenotype %in% c("desert", "excluded", "inflamed"), ]
dp$phenotype <- factor(dp$phenotype, levels = c("desert", "excluded", "inflamed"))
print(table(dp$phenotype))
M  <- sapply(LAY, function(L) tapply(dp[[paste0("z_", L)]], dp$phenotype, mean, na.rm = TRUE))
SE <- sapply(LAY, function(L) tapply(dp[[paste0("z_", L)]], dp$phenotype, function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v)))))
cat("mean z by phenotype:\n"); print(round(M, 2))
kp <- sapply(LAY, function(L) kruskal.test(dp[[paste0("z_", L)]] ~ dp$phenotype)$p.value)
dp$grp <- factor(paste(dp$phenotype, ifelse(dp$y == 1, "R", "NR")),
                 levels = c("desert NR", "desert R", "excluded NR", "excluded R", "inflamed NR", "inflamed R"))
gI  <- tapply(dp$z_Integration, dp$grp, mean, na.rm = TRUE)
gIs <- tapply(dp$z_Integration, dp$grp, function(v) sd(v, na.rm = TRUE) / sqrt(sum(!is.na(v))))
gn  <- table(dp$grp)
cat("Integration mean z by phenotype x response:\n"); print(round(gI, 2)); print(gn)
p_infl <- wilcox.test(z_Integration ~ y, data = dp[dp$phenotype == "inflamed", ])$p.value

draw5F <- function() {
  par(mfrow = c(1, 2), mar = c(5, 4.5, 2.6, 1))
  bp <- barplot(M, beside = TRUE, col = COL3, border = NA, ylim = range(c(M - SE, M + SE, 0)) * 1.25,
                ylab = "Mean layer score (z, ±SE)", las = 1)
  arrows(bp, M - SE, bp, M + SE, angle = 90, code = 3, length = 0.03, col = "grey30"); abline(h = 0)
  mtext(paste0("p=", formatC(kp, format = "g", digits = 2)), side = 1, at = colMeans(bp), line = 2.4, cex = 0.7)
  legend("topleft", c("Desert", "Excluded", "Inflamed"), fill = COL3, border = NA, bty = "n", cex = 0.8)
  tag("F (left)")
  cols <- rep(COL3, each = 2); dens <- rep(c(25, NA), 3)
  bp2 <- barplot(gI, col = cols, density = dens, border = cols, ylim = range(c(gI - gIs, gI + gIs, 0)) * 1.3,
                 ylab = "Integration (z, ±SE)", las = 2, names.arg = sub(" ", "\n", names(gI)), cex.names = 0.75)
  arrows(bp2, gI - gIs, bp2, gI + gIs, angle = 90, code = 3, length = 0.03, col = "grey30"); abline(h = 0)
  text(bp2, 0, paste0("n=", gn), pos = ifelse(gI >= 0, 1, 3), cex = 0.65)
  mtext(paste0("inflamed R vs NR p = ", formatC(p_infl, format = "g", digits = 2)), side = 3, line = -1.2, adj = 0.98, cex = 0.7)
  tag("F (right)")
}
pdf(file.path(OCT, "figures/Figure5F_IMvigor_phenotype.pdf"), 10, 4.8); draw5F(); invisible(dev.off())
png(file.path(OCT, "figures/Figure5F_IMvigor_phenotype.png"), 10, 4.8, units = "in", res = 300); draw5F(); invisible(dev.off())
fwrite(data.frame(phenotype = rownames(M), M, n = as.vector(table(dp$phenotype))), file.path(OCT, "source_data/Figure5F_left_source.csv"))
fwrite(data.frame(group = names(gI), mean_z = gI, se = gIs, n = as.vector(gn)), file.path(OCT, "source_data/Figure5F_right_source.csv"))

# =====================================================================
# SUPP FIG S6 — Integration within inflamed strata
# =====================================================================
cat("\n== SUPP FIG S6 ==\n")
SIGS <- intersect(c("Integration", "CD8A", "TIS", "Teff", "Composite"), names(d))
strata <- list()
if ("TIS" %in% names(d))  strata[["TIS-high"]]  <- d[!is.na(d$TIS)  & d$TIS  > median(d$TIS,  na.rm = TRUE), ]
if ("CD8A" %in% names(d)) strata[["CD8A-high"]] <- d[!is.na(d$CD8A) & d$CD8A > median(d$CD8A, na.rm = TRUE), ]
if (!length(strata)) stop("No TIS or CD8A column found in the benchmark file - send Claude names(bn).")
S6 <- do.call(rbind, lapply(names(strata), function(st) {
  s <- strata[[st]]; s <- s[!is.na(s$y), ]
  do.call(rbind, lapply(SIGS, function(g) {
    x <- as.numeric(scale(s[[g]])); m <- glm(s$y ~ x, binomial); ci <- confint.default(m)["x", ]
    data.frame(stratum = st, n = nrow(s), responders = sum(s$y), signature = g,
               OR = exp(coef(m)["x"]), lo = exp(ci[1]), hi = exp(ci[2]), p = coef(summary(m))["x", 4], row.names = NULL) }))
}))
print(S6, digits = 3)
drawS6 <- function() {
  par(mfrow = c(1, length(strata)), mar = c(4.5, 7, 2.6, 1))
  for (st in names(strata)) {
    s <- S6[S6$stratum == st, ]; y <- rev(seq_len(nrow(s)))
    col <- ifelse(s$signature == "Integration", "#1F6FB4", "grey45")
    plot(NA, log = "x", xlim = c(0.5, 4), ylim = c(0.5, nrow(s) + 0.5), yaxt = "n", ylab = "",
         xlab = "OR for response per SD (95% CI)", main = sprintf("%s (n=%d)", st, s$n[1]), font.main = 1)
    abline(v = 1, lty = 2, col = "grey60")
    segments(s$lo, y, s$hi, y, lwd = 2, col = col); points(s$OR, y, pch = 19, col = col)
    axis(2, at = y, labels = s$signature, las = 1)
    text(3.9, y, paste0("p=", formatC(s$p, format = "f", digits = 3)), adj = 1, cex = 0.75)
  }
}
pdf(file.path(OCT, "supplementary_figures/SuppFigS6_inflamed_strata.pdf"), 10, 4.2); drawS6(); invisible(dev.off())
png(file.path(OCT, "supplementary_figures/SuppFigS6_inflamed_strata.png"), 10, 4.2, units = "in", res = 300); drawS6(); invisible(dev.off())
fwrite(S6, file.path(OCT, "source_data/SuppFigS6_source.csv"))

# =====================================================================
# Search for the three old supplementary figures not yet located
# =====================================================================
cat("\n== SEARCHING FOR OLD S1 / S4A-E / S5B ==\n")
all_f <- list.files(BASE, recursive = TRUE, full.names = TRUE, pattern = "\\.(pdf|png|tiff?)$", ignore.case = TRUE)
all_f <- all_f[!grepl("/Oct 16/|ARCHIVE_pre_Aug|/_archive/", all_f)]
want <- list(SuppFigS1_regrouped  = "indication|by_agent|by_target|regroup|FigureS2|Supp.*S2",
             SuppFigS4AE_Braun_NSCLC_KM = "braun.*(km|surv)|nsclc.*(km|surv)|FigureS6|Supp.*S6",
             SuppFigS5B_dominant_d = "dominant|FigureS7B|S7B|cohens.*dom|best_layer")
for (nm in names(want)) {
  hit <- all_f[grepl(want[[nm]], basename(all_f), ignore.case = TRUE)]
  cat(sprintf("%-28s %d candidate(s)\n", nm, length(hit)))
  for (h in hit) { cat("    ", sub(BASE, "", h, fixed = TRUE), "\n")
    file.copy(h, file.path(OCT, "supplementary_figures/_candidates", paste0(nm, "__", basename(h))), overwrite = TRUE) }
}
cat("\nDone. New panels in Oct 16/figures and Oct 16/supplementary_figures; candidates in supplementary_figures/_candidates\n")
