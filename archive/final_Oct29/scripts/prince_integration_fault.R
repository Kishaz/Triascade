# =====================================================================
# PRINCE — where inside Integration is the fault?
# =====================================================================
# Integration = INTPOS - INTNEG, so a low score has two opposite causes
# with two opposite treatments:
#
#   INTPOS low   nothing is recruiting immune cells. No chemokine
#                gradient, no adhesion signal. Fix = induce recruitment
#                (STING agonist, interferon, epigenetic therapy)
#
#   INTNEG high  the calls are being made but a barrier blocks entry.
#                TGFB, collagen, FAP, VEGF. Fix = dismantle the barrier
#                (TGFB inhibitor, anti-angiogenic, stromal targeting)
#
# Both give the same Integration number. Only the halves distinguish
# them, and only the distinction implies a drug.
#
# Runs in the A1 (anti-PD-1) arm where Integration tracked survival at
# HR 0.43, then per gene to name the specific culprits.
#
# n = 17 with 14 deaths. Per-gene work across 41 genes is exploratory —
# the script prints the Bonferroni threshold and expects roughly two
# genes below 0.05 by chance alone.
# =====================================================================

suppressPackageStartupMessages({ library(dplyr); library(survival) })

MAN <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"
OUT <- file.path(MAN, "PRINCE")
ROOT <- file.path(MAN, "data", "prince-trial-data")

INTPOS <- c("CXCL9","CXCL10","CXCL11","CXCL13","CCL5","CCL4","CCL19","CCL21",
            "ICAM1","VCAM1","MADCAM1","CXCR3","CCR5","CCR7",
            "CD80","CD86","CD40","BATF3","IRF8","CR2","TNFRSF17")
INTNEG <- c("VEGFA","VEGFB","VEGFC","FLT1","KDR","TGFB1","TGFB2","TGFB3",
            "COL1A1","COL1A2","COL3A1","COL4A1","FAP","ACTA2","FN1",
            "MMP2","MMP9","MMP14","CXCL12","EDNRB")

# functional sub-blocks of the barrier, since they imply different drugs
SUB <- list(
  RECRUIT_chemokine = c("CXCL9","CXCL10","CXCL11","CXCL13","CCL5","CCL4","CCL19","CCL21"),
  RECRUIT_adhesion  = c("ICAM1","VCAM1","MADCAM1","CXCR3","CCR5","CCR7"),
  RECRUIT_APC       = c("CD80","CD86","CD40","BATF3","IRF8","CR2","TNFRSF17"),
  BARRIER_TGFB      = c("TGFB1","TGFB2","TGFB3"),
  BARRIER_matrix    = c("COL1A1","COL1A2","COL3A1","COL4A1","FN1","FAP","ACTA2"),
  BARRIER_MMP       = c("MMP2","MMP9","MMP14"),
  BARRIER_vascular  = c("VEGFA","VEGFB","VEGFC","FLT1","KDR","EDNRB"),
  BARRIER_CXCL12    = "CXCL12")

# ---- rebuild the filtered matrix -------------------------------------
files <- list.files(file.path(ROOT, "RNAseq", "Data"), pattern = "\\.tsv$", full.names = TRUE)
raw <- lapply(files, function(f) read.delim(f, check.names = FALSE, stringsAsFactors = FALSE))
sids <- sapply(raw, function(x) as.character(x$sample[1]))
vecs <- lapply(raw, function(x) {
  v <- setNames(as.numeric(x$TPM), toupper(x$`Gene Symbol`)); v[!duplicated(names(v))] })
g <- Reduce(intersect, lapply(vecs, names))
EALL <- sapply(vecs, function(v) v[g]); rownames(EALL) <- g; colnames(EALL) <- sids
rm(raw, vecs); invisible(gc())

md <- read.csv(file.path(ROOT, "RNAseq", "NatureMed_GX_ph2_metadata.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
cl <- read.csv(file.path(ROOT, "Clinical", "PICI0002_ph2_clinical.csv"),
               check.names = FALSE, stringsAsFactors = FALSE)
sel <- md[!is.na(md$gdc.anatomic.site.name) & md$gdc.anatomic.site.name == "Liver", ]
cols <- intersect(sel$sample.id, colnames(EALL))
E <- EALL[, cols, drop = FALSE]
if (max(E, na.rm = TRUE) > 50) E <- log2(E + 1)

sum_binary <- function(gs) {
  gs <- intersect(gs, rownames(E)); if (!length(gs)) return(rep(NA_real_, ncol(E)))
  Reduce(`+`, lapply(gs, function(x) {
    v <- as.numeric(E[x, ]); as.integer(v >= median(v, na.rm = TRUE)) }))
}

d <- data.frame(sample.id = colnames(E), INTPOS_sum = sum_binary(INTPOS),
                INTNEG_sum = sum_binary(INTNEG), stringsAsFactors = FALSE)
d$Integration <- d$INTPOS_sum - d$INTNEG_sum
for (nm in names(SUB)) d[[nm]] <- sum_binary(SUB[[nm]])

d <- d %>%
  left_join(sel[, c("sample.id", "Deidentified.ID")], by = "sample.id") %>%
  left_join(cl, by = "Deidentified.ID")
d$arm <- d$`Actual Arm`
d$os_time <- suppressWarnings(as.numeric(d$clinical.observation.os))
ev <- d$clinical.observation.os.event
d$os_event <- if (is.numeric(ev)) as.integer(ev > 0) else
  as.integer(grepl("true|dead|death|^1$|yes", tolower(as.character(ev))))

A <- d[d$arm == "A1" & !is.na(d$os_time) & !is.na(d$os_event), ]
cat(sprintf("A1 (anti-PD-1): n=%d, %d deaths\n\n", nrow(A), sum(A$os_event)))

cox1 <- function(dat, v) {
  x <- as.numeric(scale(dat[[v]]))
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(NULL)
  s <- tryCatch(summary(coxph(Surv(os_time, os_event) ~ x, data = dat)),
                error = function(e) NULL)
  if (is.null(s)) return(NULL)
  c(HR = s$conf.int[1], lo = s$conf.int[3], hi = s$conf.int[4], p = s$coef[1, 5])
}

# ---- 1. which half carries it? ---------------------------------------
cat("=== 1. Which half of Integration drives the signal? ===\n")
cat("    HR below 1 = higher value, longer survival\n\n")
for (v in c("Integration", "INTPOS_sum", "INTNEG_sum")) {
  r <- cox1(A, v); if (is.null(r)) next
  cat(sprintf("  %-14s HR %.2f (%.2f-%.2f)  p=%.4f%s\n", v, r[1], r[2], r[3], r[4],
              ifelse(v == "INTNEG_sum" && r[1] > 1, "   <- barrier: more is worse", "")))
}
cat("\n  INTPOS carries it  -> recruitment failure; induce chemokines\n")
cat("  INTNEG carries it  -> barrier; TGFB / stromal / anti-angiogenic\n")
cat("  neither alone      -> the contrast matters, not either half\n")

# ---- 2. which functional block? --------------------------------------
cat("\n=== 2. Which sub-block? ===\n\n")
res <- list()
for (nm in names(SUB)) {
  r <- cox1(A, nm); if (is.null(r)) next
  cat(sprintf("  %-18s (%2d genes) HR %.2f (%.2f-%.2f)  p=%.4f\n",
              nm, length(intersect(SUB[[nm]], rownames(E))), r[1], r[2], r[3], r[4]))
  res[[nm]] <- data.frame(block = nm, HR = r[1], lo = r[2], hi = r[3], p = r[4])
}

# ---- 3. per gene ------------------------------------------------------
allg <- intersect(c(INTPOS, INTNEG), rownames(E))
cat(sprintf("\n=== 3. Per gene (%d genes, Bonferroni %.4f) ===\n", length(allg), 0.05/length(allg)))
cat("    expect ~2 below 0.05 by chance; treat as a ranking\n\n")
G <- do.call(rbind, lapply(allg, function(gg) {
  A$.g <- as.numeric(E[gg, A$sample.id])
  r <- cox1(A, ".g"); if (is.null(r)) return(NULL)
  data.frame(gene = gg, side = ifelse(gg %in% INTPOS, "INTPOS", "INTNEG"),
             HR = r[1], lo = r[2], hi = r[3], p = r[4], stringsAsFactors = FALSE)
}))
G <- G[order(G$p), ]
print(utils::head(G, 12), row.names = FALSE)

write.csv(G, file.path(OUT, "PRINCE_A1_pergene_survival.csv"), row.names = FALSE)
if (length(res)) write.csv(do.call(rbind, res),
    file.path(OUT, "PRINCE_A1_subblocks.csv"), row.names = FALSE)

cat("\nCross-check: the trial's own analysis found a higher TGF-beta\n")
cat("signature predicted shorter survival after nivo/chemo. If\n")
cat("BARRIER_TGFB comes out with HR above 1 here, the panel recovers\n")
cat("their result independently.\n")
cat("\nWrote", file.path(OUT, "PRINCE_A1_pergene_survival.csv"), "\n")
