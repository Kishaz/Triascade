# =====================================================================
# NR TAXONOMY (proof of method, I-SPY2)
# Classify every tumour by its LIMITING cascade tier from BASELINE
# expression (outcome-blind), then ask:
#   (1) among non-responders (no pCR), which block dominates?
#   (2) does the Integration-suppressed block get the ICB benefit?
#   (3) within Integration-suppressed NR, which suppressive MODULE
#       (Treg/checkpoint, stromal/ECM, angiogenic/CAF-TGFb) is highest
#       -> names the matched targeted agent.
# Tiers/modules from the DERIVED panel (ICI_PANEL_v1_locked.csv).
# =====================================================================
suppressPackageStartupMessages({ library(data.table) })

DER <- "../ICI-study-redesign/provenance-review-tight/derivation/ICI_PANEL_v1_locked.csv"
dp  <- fread(DER, data.table = FALSE); dp$gene <- toupper(dp$gene)
MOD <- list(
  Sense        = dp$gene[grepl("^Sensory", dp$layer)],
  Traffic      = dp$gene[dp$layer == "Integration (chemokine trafficking)"],
  Treg_ckpt    = dp$gene[dp$layer == "Integration (Treg / checkpoint regulation)"],
  Stromal_ECM  = dp$gene[dp$layer == "Resistance (stromal / ECM-coagulation)"],
  Angio_CAF    = dp$gene[dp$layer == "Resistance (angiogenic / CAF-TGFb)"],
  Execute      = dp$gene[grepl("^Execution", dp$layer)])
cat("module sizes:\n"); print(lengths(MOD))

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

## ---------- tier scores (outcome-blind), adequacy = tier z >= 0 (cohort median) ----------
Sense <- Mz[,"Sense"]
Integ <- Z(M[,"Traffic"] - (M[,"Treg_ckpt"] + M[,"Stromal_ECM"] + M[,"Angio_CAF"])/3)  # net: activating - mean suppressive
Exec  <- Mz[,"Execute"]

## limiting-tier class = WEAKEST LINK (outcome-blind): the tier with the
## lowest standardized score, called only if meaningfully deficient (z < FLOOR);
## otherwise the cascade is intact. One limiting tier per tumour.
FLOOR <- -0.5
tiers <- cbind(Sense = Sense, Integration = Integ, Execution = Exec)
minz  <- apply(tiers, 1, min); which_min <- colnames(tiers)[apply(tiers, 1, which.min)]
klass <- ifelse(minz >= FLOOR, "Cascade-intact",
         ifelse(which_min == "Sense",       "Sense-low (cold)",
         ifelse(which_min == "Integration", "Integration-suppressed", "Execution-low")))
klass <- factor(klass, c("Sense-low (cold)","Integration-suppressed","Execution-low","Cascade-intact"))

cat(sprintf("\neligible n=%d; responders (pCR) %d, non-responders %d\n", length(y), sum(y==1), sum(y==0)))
cat("\n(1) baseline block class distribution, ALL patients:\n"); print(table(klass))
cat("\n   among NON-RESPONDERS only:\n"); print(table(klass[y==0]))
cat("\n   among RESPONDERS only:\n");     print(table(klass[y==1]))

cat("\n(2) pCR rate by block class x arm (does Integration-suppressed get the ICB benefit?):\n")
tb <- round(tapply(y, list(klass, arm), function(z) 100*mean(z)), 0)
nn <- tapply(y, list(klass, arm), length)
show <- data.frame(class=rownames(tb), control_pCR=tb[,"control"], pembro_pCR=tb[,"pembro"],
                   benefit_pts = tb[,"pembro"]-tb[,"control"],
                   n_ctrl=nn[,"control"], n_pembro=nn[,"pembro"])
print(show, row.names = FALSE)

cat("\n(3) within Integration-suppressed NON-RESPONDERS, mean suppressive-module z (which agent?):\n")
iss <- which(klass=="Integration-suppressed" & y==0)
supp <- colMeans(Mz[iss, c("Treg_ckpt","Stromal_ECM","Angio_CAF"), drop=FALSE])
print(round(sort(supp, decreasing=TRUE), 2))
cat(sprintf("   (n = %d Integration-suppressed non-responders)\n", length(iss)))
