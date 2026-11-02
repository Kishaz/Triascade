# =====================================================================
# NR TAXONOMY on the LEGACY 98-gene panel (proof of method, I-SPY2).
# Legacy 98 only (study decision 2026-11-02). Modules = ICI_98_module_assignment.csv
# (M1-M6); layer structure matches Fig 2B:
#   Sense       = M1 (IFN/antigen)                           [z of prop-high]
#   Integration = M2 activating  -  mean(M3 stromal, M4 angio) suppressive
#   Execution   = M5 effector     -  M6 checkpoint/exhaustion         <- legacy - arm
# Classify every tumour by its LIMITING layer from BASELINE expression
# (outcome-blind), then ask:
#   (1) among non-responders (no pCR), which layer dominates?
#   (2) does each class get the ICB benefit (pembro - control pCR)?
#   (3) within Integration-suppressed NR, which suppressive MODULE
#       (stromal/ECM vs angiogenic/CAF-TGFb) is highest -> matched agent.
# Scoring is cohort-median (same regime as Fig 2C). Numbers never typed.
# =====================================================================
suppressPackageStartupMessages({ library(data.table) })

asg <- read.csv("ICI_98_module_assignment.csv", stringsAsFactors = FALSE); asg$gene <- toupper(asg$gene)
MOD <- list(
  Sense       = asg$gene[asg$module == "M1"],
  Traffic     = asg$gene[asg$module == "M2"],
  Stromal_ECM = asg$gene[asg$module == "M3"],
  Angio_CAF   = asg$gene[asg$module == "M4"],
  Execute     = asg$gene[asg$module == "M5"],
  Checkpoint  = asg$gene[asg$module == "M6"])
cat("module sizes (legacy 98):\n"); print(lengths(MOD))

## ---------- load I-SPY2 ----------
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

## ---------- module scores: proportion of module genes high (>= cohort median), z-scored ----------
prop_high <- function(genes) { g <- intersect(genes, rownames(E)); if (!length(g)) return(rep(NA, ncol(E)))
  b <- t(apply(E[g,,drop=FALSE], 1, function(v) as.integer(v >= median(v, na.rm=TRUE)))); colMeans(b, na.rm=TRUE) }
Z <- function(v) (v - mean(v, na.rm=TRUE)) / sd(v, na.rm=TRUE)
M <- sapply(MOD, prop_high); Mz <- apply(M, 2, Z)

## ---------- layer scores (outcome-blind), matching Fig 2B legacy structure ----------
Sense <- Mz[,"Sense"]
Integ <- Z(M[,"Traffic"] - (M[,"Stromal_ECM"] + M[,"Angio_CAF"])/2)   # activating - mean suppressive (2 modules)
Exec  <- Z(M[,"Execute"] - M[,"Checkpoint"])                           # effector - exhaustion (legacy - arm)

## limiting-layer class = WEAKEST LINK (outcome-blind): lowest standardized score,
## called only if meaningfully deficient (z < FLOOR); else cascade intact.
FLOOR <- -0.5
tiers <- cbind(Sense = Sense, Integration = Integ, Execution = Exec)
minz  <- apply(tiers, 1, min); which_min <- colnames(tiers)[apply(tiers, 1, which.min)]
klass <- ifelse(minz >= FLOOR, "Cascade-intact",
         ifelse(which_min == "Sense",       "Sense-low (cold)",
         ifelse(which_min == "Integration", "Integration-suppressed", "Execution-low")))
klass <- factor(klass, c("Sense-low (cold)","Integration-suppressed","Execution-low","Cascade-intact"))

cat(sprintf("\neligible n=%d; responders (pCR) %d, non-responders %d\n", length(y), sum(y==1), sum(y==0)))
cat("\n(1) baseline layer-class distribution, ALL patients:\n"); print(table(klass))
cat("\n   among NON-RESPONDERS only:\n"); print(table(klass[y==0]))
cat("\n   among RESPONDERS only:\n");     print(table(klass[y==1]))

cat("\n(2) pCR rate by layer class x arm (does each class get the ICB benefit?):\n")
tb <- round(tapply(y, list(klass, arm), function(z) 100*mean(z)), 0)
nn <- tapply(y, list(klass, arm), length)
show <- data.frame(class=rownames(tb), control_pCR=tb[,"control"], pembro_pCR=tb[,"pembro"],
                   benefit_pts = tb[,"pembro"]-tb[,"control"],
                   n_ctrl=nn[,"control"], n_pembro=nn[,"pembro"])
print(show, row.names = FALSE)

cat("\n(3) within Integration-suppressed NON-RESPONDERS, mean suppressive-module z (which agent?):\n")
iss <- which(klass=="Integration-suppressed" & y==0)
supp <- colMeans(Mz[iss, c("Stromal_ECM","Angio_CAF"), drop=FALSE])
print(round(sort(supp, decreasing=TRUE), 2))
cat(sprintf("   (n = %d Integration-suppressed non-responders)\n", length(iss)))

## ---------- save the NR composition (for a Fig 2C explanatory element) ----------
## Fig 2C uses the ICB(pembro)-arm denominator: the framework's job is to explain
## ICB failure, so pooled-arm NR counts are not the honest "why ICB failed" figure.
nr_comp <- as.data.frame(table(klass[y==0])); names(nr_comp) <- c("class","n_NR")
nr_comp$pct_NR <- round(100*nr_comp$n_NR/sum(nr_comp$n_NR), 1)
pnr <- arm=="pembro" & y==0
nr_comp_pembro <- as.data.frame(table(klass[pnr])); names(nr_comp_pembro) <- c("class","n_NR")
nr_comp_pembro$pct_NR <- round(100*nr_comp_pembro$n_NR/sum(nr_comp_pembro$n_NR), 1)
cat("\n(4) ICB(pembro)-arm NR composition (the honest 'why ICB failed' denominator):\n")
print(nr_comp_pembro, row.names=FALSE)
saveRDS(list(nr_comp=nr_comp, nr_comp_pembro=nr_comp_pembro, by_class=show,
             n_R=sum(y==1), n_NR=sum(y==0), n_NR_pembro=sum(pnr),
             supp_module=round(sort(supp, decreasing=TRUE),2)),
        "nr_taxonomy_legacy98_ISPY2.rds")
cat("\nsaved nr_taxonomy_legacy98_ISPY2.rds\n")
