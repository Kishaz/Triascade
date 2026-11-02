# =====================================================================
SCRIPT_VERSION <- "ICI_NSCLC_SPLIT  v1.2  2026-10-27"
# v1.2: the old combined NSCLC scores were not correct (one median across both studies).
#       The comparison with them is now for information only and never blocks the
#       correct per-study scoring. Defaults: log2(CPM+1) for GSE126044 (as in your NSCLC
#       scripts) and the >= median rule used for every other cohort.
# Scores the two NSCLC studies separately (GSE126044, GSE135222), each at its
# own cohort median - the same binary method as every other cohort.
# v1.1: reads the expression files from "R data FINAL/datasets", prepared the
#       way your NSCLC scripts prepare them (GSE135222: Ensembl -> symbol with
#       org.Hs.eg.db, log2(TPM+1); GSE126044: log2(CPM+1)).
#  1. for information: compares a both-studies-together scoring with your old
#     MASTER_GRID/ICI_binary_nsclc_scores.csv (the incorrect combined version)
#  2. scores each study at its own median (the correct analysis), writes
#     MASTER_GRID/ICI_binary_nsclc_scores_by_study.csv (old file untouched)
#     and the per-study AUCs (DeLong 95% CI) as Supplementary Table S3
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
BASE <- "/Users/dbedi/Desktop/R data FINAL"
MAN  <- file.path(BASE, "ICI_Manuscript"); GRID <- file.path(MAN, "MASTER_GRID"); OCT <- file.path(MAN, "Oct 16")
DS   <- file.path(BASE, "datasets")
F135 <- file.path(DS, "GSE135222", "GSE135222", "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz")
F126 <- file.path(DS, "GSE126044", "GSE126044", "GSE126044_counts.txt.gz")
cat("\n", SCRIPT_VERSION, "\n", sep = "")
for (f in c(F135, F126)) if (!file.exists(f)) stop("not found: ", f)

PANEL <- list(
  Sensory = c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","B2M","TAP1","TAP2","TAPBP","PSMB8","PSMB9",
              "PSMB10","IRF1","STAT1","STAT2","JAK1","JAK2","IFNG","IFNGR1","IFNGR2"),
  IntPos  = c("CXCL9","CXCL10","CXCL11","CXCL13","CCL4","CCL5","CCL19","CCL21","ICAM1","VCAM1","MADCAM1",
              "CXCR3","CCR5","CCR7","CD80","CD86","CD40","BATF3","IRF8","CR2","TNFRSF17"),
  IntNeg  = c("VEGFA","VEGFB","VEGFC","FLT1","KDR","TGFB1","TGFB2","TGFB3","COL1A1","COL1A2","COL3A1",
              "COL4A1","FAP","ACTA2","FN1","MMP2","MMP9","MMP14","CXCL12","EDNRB"),
  ExecPos = c("GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","GNLY","NKG7","FASLG","IFNG","TNF","LTA","CD8A",
              "CD8B","TBX21","EOMES","RUNX3","IL2","IL12A"),
  ExecNeg = c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","BTLA","VSIR","CD274","PDCD1LG2","CD160","CD244",
              "TOX","NR4A1","NR4A2","NR4A3","ENTPD1","IL10","TGFB1","IDO1"))
ALLG <- unique(unlist(PANEL)); stopifnot(length(ALLG) == 98)

old <- fread(file.path(GRID, "ICI_binary_nsclc_scores.csv"), data.table = FALSE)
old$sample_id <- as.character(old$sample_id)
cat("your NSCLC file:", nrow(old), "patients;", paste(names(table(old$cohort)), table(old$cohort), sep = " n=", collapse = ", "), "\n")
for (s in unique(old$cohort)) cat("  ", s, "IDs e.g.:", paste(head(old$sample_id[old$cohort == s], 3), collapse = ", "), "\n")

fix_hla <- function(m) {                     # HLA.A -> HLA-A where needed
  for (g in ALLG) { alt <- gsub("-", ".", g)
    if (!(g %in% rownames(m)) && alt %in% rownames(m)) rownames(m)[rownames(m) == alt] <- g }
  m[rownames(m) %in% ALLG, , drop = FALSE]
}

# ---------------- load GSE135222 (as in your NSCLC scripts) --------------------
x <- read.delim(gzfile(F135), check.names = FALSE)
ens <- sub("\\..*", "", as.character(x[[1]]))
m135 <- as.matrix(x[, -1])
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("org.Hs.eg.db is needed (it was used by your NSCLC scripts)")
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db, keys = unique(ens),
                                              keytype = "ENSEMBL", columns = "SYMBOL"))
map <- map[!duplicated(map$ENSEMBL) & !is.na(map$SYMBOL), ]
idx <- match(ens, map$ENSEMBL); has <- !is.na(idx)
m135 <- m135[has, , drop = FALSE]; rownames(m135) <- make.unique(map$SYMBOL[idx[has]])
m135 <- fix_hla(log2(m135 + 1))

# ---------------- load GSE126044 -----------------------------------------------
y <- read.delim(gzfile(F126), check.names = FALSE)
m126raw <- as.matrix(y[, -1]); rownames(m126raw) <- make.unique(as.character(y[[1]]))
m126cpm <- log2(sweep(m126raw, 2, colSums(m126raw, na.rm = TRUE), "/") * 1e6 + 1)
m126 <- list(cpm = fix_hla(m126cpm), raw = fix_hla(m126raw))
cat("panel genes found: GSE135222", nrow(m135), "| GSE126044", nrow(m126$cpm), "\n")
cat("expression columns e.g.: GSE135222", paste(head(colnames(m135), 3), collapse = ", "),
    "| GSE126044", paste(head(colnames(m126raw), 3), collapse = ", "), "\n")

# ---------------- match expression columns to your sample IDs ------------------
gsm_table <- function(gse, to_col) {          # expression column <-> GSM via local series matrix
  sm <- list.files(file.path(DS, gse), pattern = "series_matrix.*\\.txt(\\.gz)?$", recursive = TRUE, full.names = TRUE)[1]
  if (is.na(sm) || !requireNamespace("GEOquery", quietly = TRUE)) return(NULL)
  pd <- Biobase::pData(suppressMessages(GEOquery::getGEO(filename = sm, getGPL = FALSE)))
  data.frame(col = to_col(pd$title), gsm = pd$geo_accession, title = pd$title, stringsAsFactors = FALSE)
}
align <- function(m, gse, to_col) {
  want <- old$sample_id[old$cohort == gse]
  if (sum(colnames(m) %in% want) >= length(want) - 1) return(m[, colnames(m) %in% want, drop = FALSE])
  tb <- gsm_table(gse, to_col); if (is.null(tb)) return(m[, 0, drop = FALSE])
  for (key in c("gsm", "title")) {
    nm <- tb[[key]][match(colnames(m), tb$col)]
    if (sum(nm %in% want) >= length(want) - 1) { colnames(m) <- nm; return(m[, colnames(m) %in% want, drop = FALSE]) }
  }
  m[, 0, drop = FALSE]
}
m135 <- align(m135, "GSE135222", function(t) gsub("NSCLC ", "NSCLC", t))
m126 <- lapply(m126, align, gse = "GSE126044", to_col = function(t) gsub("RNA-seq_", "", t))
cat("samples matched to your file: GSE135222", ncol(m135), "of", sum(old$cohort == "GSE135222"),
    "| GSE126044", ncol(m126$cpm), "of", sum(old$cohort == "GSE126044"), "\n")
if (ncol(m135) == 0 || ncol(m126$cpm) == 0) stop("sample IDs could not be matched - paste this output to Claude")

# ---------------- scoring (cohort-median binarization) -------------------------
score <- function(m, tie) {
  b <- t(apply(m, 1, function(v) { md <- median(v, na.rm = TRUE); as.integer(if (tie == ">=") v >= md else v > md) }))
  colnames(b) <- colnames(m)
  s <- function(gs) colSums(b[intersect(gs, rownames(b)), , drop = FALSE], na.rm = TRUE)
  S <- s(PANEL$Sensory); I <- s(PANEL$IntPos) - s(PANEL$IntNeg); E <- s(PANEL$ExecPos) - s(PANEL$ExecNeg)
  z <- function(v) (v - mean(v)) / sd(v)
  data.frame(sample_id = colnames(m), Sensory = S, Integration = I, Execution = E,
             Sensory_z = z(S), Integration_z = z(I), Execution_z = z(E),
             composite = (z(S) + z(I) + z(E)) / 3, stringsAsFactors = FALSE)
}
cmp <- function(new) {
  mm <- merge(old, new, by = "sample_id", suffixes = c("_old", "_new"))
  c(n = nrow(mm), S = max(abs(mm$Sensory_old - mm$Sensory_new)), I = max(abs(mm$Integration_old - mm$Integration_new)),
    E = max(abs(mm$Execution_old - mm$Execution_new)), C = max(abs(mm$composite_old - mm$composite_new)))
}
combine <- function(a, b, genes) {
  g <- if (genes == "intersect") intersect(rownames(a), rownames(b)) else union(rownames(a), rownames(b))
  pad <- function(m) { o <- m[match(g, rownames(m)), , drop = FALSE]; rownames(o) <- g; o }
  cbind(pad(a), pad(b))
}

# ---------------- 1. reproduce your combined scores ----------------------------
cat("\n==== for information: old combined scores vs scoring both studies together ====\n")
hit <- NULL; tried <- list()
for (norm in c("cpm", "raw")) for (genes in c("union", "intersect")) for (tie in c(">=", ">")) {
  d <- cmp(score(combine(m135, m126[[norm]], genes), tie))
  tried[[length(tried) + 1]] <- data.frame(GSE126044 = norm, genes = genes, tie = tie, t(round(d, 4)))
  if (d["n"] == nrow(old) && all(d[c("S", "I", "E")] == 0) && d["C"] < 1e-6 && is.null(hit)) hit <- list(norm = norm, genes = genes, tie = tie)
}
print(do.call(rbind, tried), row.names = FALSE)
if (is.null(hit)) {
  hit <- list(norm = "cpm", genes = "union", tie = ">=")
  cat("\nThe old combined file is not reproduced exactly by these preparations (it was made elsewhere).\n",
      "Per-study scoring below uses log2(CPM+1) for GSE126044 and the >= median rule.\n", sep = "")
} else cat("\nOld combined file reproduced exactly with: GSE126044", hit$norm, "| median rule", hit$tie,
           "- the same preparation is used for the per-study scoring below.\n")

# ---------------- 2. score each study at its own median ------------------------
by_study <- rbind(transform(score(m135, hit$tie), cohort = "GSE135222"),
                  transform(score(m126[[hit$norm]], hit$tie), cohort = "GSE126044"))
keep <- intersect(c("sample_id", "responder", "response_bin"), names(old))
by_study <- merge(by_study, old[, keep], by = "sample_id")
fwrite(by_study, file.path(GRID, "ICI_binary_nsclc_scores_by_study.csv"))
cat("\nwritten: MASTER_GRID/ICI_binary_nsclc_scores_by_study.csv (", nrow(by_study), "patients )\n")

auc_tab <- do.call(rbind, lapply(split(by_study, by_study$cohort), function(d) {
  d <- d[!is.na(d$response_bin), ]
  do.call(rbind, lapply(c("Sensory", "Integration", "Execution", "composite"), function(L) {
    r <- roc(d$response_bin, d[[L]], levels = c(0, 1), direction = "<", quiet = TRUE)
    ci <- as.numeric(ci.auc(r, method = "delong"))
    data.frame(study = d$cohort[1], n = nrow(d), responders = sum(d$response_bin == 1),
               layer = if (L == "composite") "Composite" else L,
               AUC = round(ci[2], 3), CI_low = round(max(0, ci[1]), 3), CI_high = round(min(1, ci[3]), 3)) })) }))
cat("\n==== per-study AUC, each study scored at its own median ====\n"); print(auc_tab, row.names = FALSE)
fwrite(auc_tab, file.path(OCT, "supplementary_tables", "Supp_Table_S3_NSCLC_split.csv"))
cat("\nSupplementary Table S3 updated (Oct 16/supplementary_tables/Supp_Table_S3_NSCLC_split.csv)\n")
cat("\nDONE - paste everything above to Claude.\n")
