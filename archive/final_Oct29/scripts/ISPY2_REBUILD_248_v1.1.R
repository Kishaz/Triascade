# =====================================================================
SCRIPT_VERSION <- "ISPY2_REBUILD_248  v1.1  2026-10-04"
# Rebuild the I-SPY2 ICI score file on the full eligible cohort (n=248)
# from GEO metadata, replacing the stored n=189 (which was the
# intersection with the old SRKI files; 59 paclitaxel controls on
# GPL30493 were never parsed).
#
# Panel, expression parsing, metadata parsing and scoring are copied
# unchanged from ici_ispy2_hrpos_v1_2.R. Nothing new is estimated.
#
# Output: MASTER_GRID/ICI_binary_ispy2_scores_248.csv
#   same columns as the stored file (pid_str, Sensory, Integration,
#   Execution, composite column under its stored name, arm, pcr_binary)
#   plus hr, her2, receptor, platform.
# The old file is NOT overwritten.
# =====================================================================

suppressPackageStartupMessages({ library(data.table) })

# ----------------------------- CONFIG --------------------------------
BASE   <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"
EXPR   <- file.path(BASE, "ISPY2/GSE194040",
                    "GSE194040_ISPY2ResID_AgilentGeneExp_990_FrshFrzn_meanCol_geneLevel_n988.txt.gz")
SERIES <- file.path(BASE, "ISPY2", c("GSE194040-GPL20078_series_matrix.txt.gz",
                                     "GSE194040-GPL30493_series_matrix.txt.gz"))
STORED <- file.path(BASE, "MASTER_GRID/ICI_binary_ispy2_scores.csv")
OUTFILE <- file.path(BASE, "MASTER_GRID/ICI_binary_ispy2_scores_248.csv")
AUDIT   <- file.path(BASE, "ISPY2/REBUILD_248")
COL_ID <- "patient.id"; COL_ARM <- "arm"; COL_HR <- "hr"; COL_HER2 <- "her2"; COL_PCR <- "pcr"
# ---------------------------------------------------------------------

dir.create(AUDIT, showWarnings = FALSE, recursive = TRUE)
cat("\n", SCRIPT_VERSION, "\n", strrep("=", 62), "\n", sep = "")

# ============ STAGE 0: frozen panel (unchanged from v1.2) =============
panel_src <- list.files(BASE, pattern = "^build_scRNA_celltype\\.R$",
                        recursive = TRUE, full.names = TRUE)[1]
pe <- new.env(); lines <- readLines(panel_src, warn = FALSE)
keep <- grep("^\\s*(SENSORY|INTPOS|INTNEG|EXECPOS|EXECNEG)\\s*<-", lines)
blk <- c()
for (k in keep) {
  j <- k
  repeat {
    blk <- c(blk, lines[j]); s <- paste(lines[k:j], collapse = "")
    if (lengths(regmatches(s, gregexpr("\\(", s))) ==
        lengths(regmatches(s, gregexpr("\\)", s)))) break
    j <- j + 1; if (j > k + 25) stop("Runaway parsing panel at line ", k)
  }
}
eval(parse(text = blk), envir = pe)
PANEL <- list(SENSORY = get("SENSORY", pe), INTPOS = get("INTPOS", pe),
              INTNEG = get("INTNEG", pe), EXECPOS = get("EXECPOS", pe),
              EXECNEG = get("EXECNEG", pe))
stopifnot(identical(unname(lengths(PANEL)), c(21L,21L,20L,19L,19L)),
          length(unique(unlist(PANEL))) == 98L)
cat("Panel verified: 98 unique genes from", basename(panel_src), "\n")

# ============ STAGE 1: expression =====================================
ex <- suppressWarnings(fread(EXPR, data.table = FALSE))
rownames(ex) <- ex[[1]]; ex <- as.matrix(ex[, -1, drop = FALSE])
offp <- unlist(lapply(PANEL, function(g) setdiff(g, rownames(ex))))
cat("expression:", nrow(ex), "x", ncol(ex),
    "| off-platform:", paste(offp, collapse = ", "), "\n")   # expect VSIR

# ============ STAGE 2: metadata, BOTH platforms =======================
parse_sm <- function(f) {
  ln <- readLines(gzfile(f), warn = FALSE); ln <- ln[grepl("^!Sample_", ln)]
  out <- list()
  for (r in strsplit(ln, "\t")) {
    key <- sub("^!Sample_", "", r[1]); val <- gsub('^"|"$', '', r[-1])
    if (grepl("^characteristics", key)) {
      i <- which(grepl(":", val))[1]; if (is.na(i)) next
      nm <- make.names(trimws(sub(":.*$", "", val[i])))
      val <- trimws(sub("^[^:]*:\\s*", "", val))
    } else nm <- make.names(key)
    while (nm %in% names(out)) nm <- paste0(nm, ".x")
    out[[nm]] <- val
  }
  n <- max(lengths(out))
  as.data.frame(lapply(out, function(v) { length(v) <- n; v }), stringsAsFactors = FALSE)
}
stopifnot(all(file.exists(SERIES)))                 # both platforms required
mds <- lapply(SERIES, parse_sm)
for (i in seq_along(mds)) mds[[i]]$.platform <- sub("^GSE194040-(GPL\\d+).*$", "\\1", basename(SERIES[i]))
cmn <- Reduce(intersect, lapply(mds, names))
md <- do.call(rbind, lapply(mds, function(d) d[, cmn, drop = FALSE]))
md$.id <- as.character(md[[COL_ID]])
cat("metadata:", nrow(md), "samples |", paste(names(table(md$.platform)), table(md$.platform), collapse = "  "), "\n")

norm01 <- function(v) {
  v <- as.character(v)
  ifelse(grepl("^1$|pos|\\+|yes|true", v, ignore.case = TRUE), 1L,
  ifelse(grepl("^0$|neg|\\-|no|false", v, ignore.case = TRUE), 0L, NA_integer_))
}
md$.hr <- norm01(md[[COL_HR]]); md$.her2 <- norm01(md[[COL_HER2]])
md$.pcr <- norm01(md[[COL_PCR]]); md$.arm <- as.character(md[[COL_ARM]])
md$.arm2 <- ifelse(grepl("Pembrolizumab", md$.arm), "pembro",
            ifelse(md$.arm == "Paclitaxel", "control", NA))
md$.receptor <- ifelse(md$.her2 == 1, "HER2pos",
                ifelse(md$.hr == 1, "HRpos", "TNBC"))

# ============ STAGE 3: the eligible cohort ============================
# pembro (Paclitaxel + Pembrolizumab) or paclitaxel-alone control,
# HER2-negative, pCR known, expression present.
elig <- subset(md, .arm2 %in% c("pembro","control") & .her2 == 0 &
                   !is.na(.pcr) & .id %in% colnames(ex))
stopifnot(!any(duplicated(elig$.id)))
cat("\neligible cohort:", nrow(elig), "\n")
print(table(arm = elig$.arm2, receptor = elig$.receptor))
print(table(arm = elig$.arm2, platform = elig$.platform))
cat("pCR by arm:\n"); print(table(arm = elig$.arm2, pcr = elig$.pcr))
if (nrow(elig) != 248) warning("Expected 248 eligible patients, got ", nrow(elig),
                               " - check metadata before using the output.")

# ============ STAGE 4: scoring at the 248 median (unchanged logic) ====
score_ici <- function(m) {
  bin <- function(gs) {
    g <- intersect(gs, rownames(m))
    t(apply(m[g, , drop = FALSE], 1,
            function(v) as.integer(v >= median(v, na.rm = TRUE))))
  }
  S <- colSums(bin(PANEL$SENSORY))
  IP <- colSums(bin(PANEL$INTPOS)); IN <- colSums(bin(PANEL$INTNEG))
  EP <- colSums(bin(PANEL$EXECPOS)); EN <- colSums(bin(PANEL$EXECNEG))
  d <- data.frame(sample = colnames(m), Sensory = S, Integration = IP - IN,
                  Execution = EP - EN, row.names = NULL)
  d$ICI_composite <- rowMeans(scale(d[, c("Sensory","Integration","Execution")]))
  d
}
sc <- score_ici(ex[, elig$.id, drop = FALSE])
out <- merge(sc, elig[, c(".id",".arm2",".pcr",".hr",".her2",".receptor",".platform")],
             by.x = "sample", by.y = ".id")

# mirror the stored file's column names so downstream scripts need no edits
st <- fread(STORED, data.table = FALSE)
comp_name <- intersect(names(st), c("ICI_composite","composite","Composite"))[1]
if (is.na(comp_name)) comp_name <- "ICI_composite"
new <- data.frame(pid_str = out$sample,
                  Sensory = out$Sensory, Integration = out$Integration,
                  Execution = out$Execution, row.names = NULL)
new[[comp_name]] <- out$ICI_composite
new$arm <- out$.arm2; new$pcr_binary <- out$.pcr
new$hr <- out$.hr; new$her2 <- out$.her2
new$receptor <- out$.receptor; new$platform <- out$.platform
fwrite(new, OUTFILE)
cat("\nwritten:", OUTFILE, " (", nrow(new), "rows )\n")

# v1.1: report NA scores (a panel gene NA in some expression columns)
na_rows <- new[!complete.cases(new[, c("Sensory","Integration","Execution",comp_name)]), ]
cat("patients with an NA layer score:", nrow(na_rows), "\n")
if (nrow(na_rows)) {
  print(na_rows[, c("pid_str","arm","receptor","platform","Sensory","Integration","Execution")])
  g <- unlist(PANEL); nag <- g[colSums(is.na(ex[intersect(g, rownames(ex)), na_rows$pid_str, drop = FALSE])) > 0]
  cat("panel genes with NA expression in those patients:",
      paste(unique(names(which(rowSums(is.na(ex[intersect(g, rownames(ex)), na_rows$pid_str, drop = FALSE])) > 0))), collapse = ", "), "\n")
  fwrite(na_rows, file.path(AUDIT, "ispy2_248_NA_scores.csv"))
}

# ============ STAGE 5: audit against the stored 189 ===================
cat("\n", strrep("=", 62), "\nAUDIT vs stored n=189\n", strrep("=", 62), "\n", sep = "")
st$.id <- as.character(st$pid_str)
cat("stored:", nrow(st), " | stored ids inside the 248:", sum(st$.id %in% new$pid_str), "\n")
ov <- merge(st, new, by.x = ".id", by.y = "pid_str", suffixes = c("_189","_248"))
cat("arm agreement:", all(ov$arm_189 == ov$arm_248),
    "| pCR agreement:", all(ov$pcr_binary_189 == ov$pcr_binary_248), "\n")
cat("correlation of layer scores on the 189 overlap (threshold moves a little, expect >0.95):\n")
for (L in c("Sensory","Integration","Execution"))
  cat(sprintf("  %-12s r = %.3f\n", L, cor(ov[[paste0(L,"_189")]], ov[[paste0(L,"_248")]])))

# headline check (v1.1: NA scores dropped): Integration x arm on the 248 (should be ~p 0.025 perSD / 0.026 binarized)
fit <- function(dat, L, mode) {
  x <- dat[[L]]
  dat$X <- if (mode == "perSD") as.numeric(scale(x)) else as.integer(x >= median(x, na.rm = TRUE))
  dat <- dat[!is.na(dat$X), ]                      # v1.1: drop NA scores
  dat$A <- factor(dat$arm, levels = c("control","pembro"))
  m <- glm(pcr_binary ~ X * A, binomial, dat)
  orarm <- function(a) unname(exp(coef(glm(pcr_binary ~ X, binomial, dat[dat$arm == a, ]))[2]))
  data.frame(layer = L, mode = mode, OR_pembro = orarm("pembro"), OR_control = orarm("control"),
             p_interaction = coef(summary(m))["X:Apembro", 4])
}
chk <- do.call(rbind, lapply(c("Sensory","Integration","Execution",comp_name), function(L)
  rbind(fit(new, L, "perSD"), fit(new, L, "binarized"))))
cat("\nlayer x arm interaction on n =", nrow(new), "\n"); print(chk, digits = 3)
fwrite(chk, file.path(AUDIT, "ispy2_248_interaction_check.csv"))
fwrite(as.data.frame(table(arm = new$arm, receptor = new$receptor, pcr = new$pcr_binary)),
       file.path(AUDIT, "ispy2_248_composition.csv"))
cat("audit files in", AUDIT, "\n")

cat("\n", strrep("=", 62), "\nNEXT\n", strrep("=", 62), "\n",
    "Point these at ICI_binary_ispy2_scores_248.csv and rerun unchanged:\n",
    "  rebuild_layer_auc_table.R   (Figure 2B I-SPY2 rows)\n",
    "  fig2_ispy2_predictive.R\n",
    "  ispy2_categorical_check.R\n",
    "  layer_dominance_test.R      (pembro arm only in the Breast row)\n",
    "  the gating / min-vs-mean block and the purity script\n",
    "Keep the old 189 file for the supplementary row.\n", sep = "")
