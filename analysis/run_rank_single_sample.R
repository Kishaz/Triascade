# =====================================================================
# Rank-based SINGLE-SAMPLE scoring (singscore-style) vs cohort-median
# For each tumour, rank all measured genes WITHIN that sample; a layer
# score is the normalized mean rank of its up genes plus (0.5 - mean
# rank) of its down genes. No cohort is used: each patient is scored
# independently, continuously, and platform-robustly (ranks, not values).
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")
PANEL <- load_panel("ICI_98gene_panel.csv")

# ---- rank-based single-sample layer scores --------------------------
ss_layers <- function(expr) {                       # expr: genes x samples
  N <- nrow(expr)
  P <- apply(expr, 2, function(x) (rank(x, ties.method = "average") - 0.5) / N)  # [0,1] within sample
  up   <- function(g) { g <- intersect(g, rownames(expr)); colMeans(P[g, , drop = FALSE]) - 0.5 }
  down <- function(g) { g <- intersect(g, rownames(expr)); 0.5 - colMeans(P[g, , drop = FALSE]) }
  S <- up(PANEL$SENSORY)
  I <- up(PANEL$INTPOS) + down(PANEL$INTNEG)
  E <- up(PANEL$EXECPOS) + down(PANEL$EXECNEG)
  z <- function(v) (v - mean(v)) / sd(v)
  data.frame(sample = colnames(expr), Sensory = S, Integration = I, Execution = E,
             Composite = rowMeans(cbind(z(S), z(I), z(E))), row.names = NULL)
}

auc1 <- function(y, s) as.numeric(auc(y, s, quiet = TRUE, direction = "<"))

## ============ GSE78220 (melanoma, anti-PD-1) ============
em <- as.data.frame(fread("data/GSE78220_expr.csv")); rownames(em) <- toupper(em$gene); em$gene <- NULL
em <- log2(as.matrix(em) + 1); rr <- fread("data/GSE78220_response.csv")
em <- em[, colnames(em) %in% rr$sample]; ym <- rr$responder[match(colnames(em), rr$sample)]
cm_m <- score_cohort_median(em, PANEL); cm_m <- cm_m[match(colnames(em), cm_m$sample), ]
ss_m <- ss_layers(em);                 ss_m <- ss_m[match(colnames(em), ss_m$sample), ]
mel <- data.frame(layer = c("Sensory","Integration","Execution","Composite"),
  AUC_median = sapply(c("Sensory","Integration","Execution","Composite"), function(L) round(auc1(ym, cm_m[[L]]), 3)),
  AUC_ranksingle = sapply(c("Sensory","Integration","Execution","Composite"), function(L) round(auc1(ym, ss_m[[L]]), 3)))
cat("GSE78220 (n=", ncol(em), ") response AUC: cohort-median vs rank-single-sample\n", sep=""); print(mel, row.names = FALSE); cat("\n")

## ============ I-SPY2 ============
parse_sm <- function(f) { ln <- readLines(gzfile(f), warn = FALSE); ln <- ln[grepl("^!Sample_", ln)]; out <- list()
  for (r in strsplit(ln, "\t")) { k <- sub("^!Sample_", "", r[1]); v <- gsub('"', '', r[-1])
    if (grepl("^characteristics", k)) { i <- which(grepl(":", v))[1]; if (is.na(i)) next
      nm <- make.names(trimws(sub(":.*$", "", v[i]))); v <- trimws(sub("^[^:]*: *", "", v)) } else nm <- make.names(k)
    while (nm %in% names(out)) nm <- paste0(nm, ".x"); out[[nm]] <- v }
  n <- max(lengths(out)); as.data.frame(lapply(out, function(x){length(x)<-n; x}), stringsAsFactors = FALSE) }
md <- rbind(parse_sm("data/ISPY2_GPL20078_series.txt.gz")[,c("patient.id","arm","her2","pcr")],
            parse_sm("data/ISPY2_GPL30493_series.txt.gz")[,c("patient.id","arm","her2","pcr")])
n01 <- function(v) ifelse(grepl("pos|[+]|yes|^1$", v, TRUE), 1L, ifelse(grepl("neg|[-]|no|^0$", v, TRUE), 0L, NA))
md$.id <- as.character(md$patient.id); md$.her2 <- n01(md$her2); md$.pcr <- n01(md$pcr)
md$.arm2 <- ifelse(grepl("Pembrolizumab", md$arm), "pembro", ifelse(md$arm == "Paclitaxel", "control", NA))
ex <- as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz")); rownames(ex) <- toupper(ex[[1]]); ex[[1]] <- NULL
ex <- as.matrix(ex); colnames(ex) <- as.character(colnames(ex))
elig <- subset(md, .arm2 %in% c("pembro","control") & .her2 == 0 & !is.na(.pcr) & .id %in% colnames(ex))
elig <- elig[!duplicated(elig$.id), ]; E <- ex[, elig$.id]
y <- elig$.pcr; arm <- factor(elig$.arm2, c("control","pembro"))
cm <- score_cohort_median(E, PANEL); cm <- cm[match(elig$.id, cm$sample), ]
ssp <- ss_layers(E);                 ssp <- ssp[match(elig$.id, ssp$sample), ]
z <- function(v) scale(v)[,1]
row_for <- function(score) {
  iz <- z(score)
  p <- anova(glm(y~iz+arm, binomial), glm(y~iz*arm, binomial), test="LRT")[2,"Pr(>Chi)"]
  c(interaction_p = round(p,3), AUC_pembro = round(auc1(y[arm=="pembro"], score[arm=="pembro"]),3),
    AUC_control = round(auc1(y[arm=="control"], score[arm=="control"]),3))
}
tab <- function(sc) do.call(rbind, lapply(c("Sensory","Integration","Execution","Composite"), function(L) data.frame(layer=L, t(row_for(sc[[L]])))))
cat("I-SPY2 (n=", ncol(E), ")  cohort-median scoring:\n", sep=""); print(tab(cm), row.names = FALSE)
cat("\nI-SPY2  RANK single-sample scoring:\n"); print(tab(ssp), row.names = FALSE); cat("\n")

# single-sample invariance: one patient scored alone == within the matrix
one <- ss_layers(E[, 1, drop = FALSE])
cat(sprintf("Single-sample invariance (patient scored alone == in cohort)? %s\n",
            isTRUE(all.equal(one$Integration, ssp$Integration[ssp$sample == colnames(E)[1]]))))
write.csv(tab(ssp), "ranksingle_ISPY2.csv", row.names = FALSE)
