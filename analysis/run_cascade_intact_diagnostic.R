# =====================================================================
# DIAGNOSTIC: why is 35% of the cohort "cascade-intact"? (legacy 98, I-SPY2)
# Tests three explanations for the large cascade-intact fraction:
#   (A) arm split - are cascade-intact NRs mostly CONTROL-arm (never got ICB)?
#       The framework's job is to explain ICB failure, so the honest denominator
#       is ICB(pembro)-arm NRs.
#   (B) threshold geometry - with 3 z-scored layers + FLOOR, how much cascade-intact
#       is expected by construction? Sweep FLOOR.
#   (C) borderline - are cascade-intact tumours just above the floor (near-misses)?
# =====================================================================
suppressPackageStartupMessages({ library(data.table) })
asg <- read.csv("ICI_98_module_assignment.csv", stringsAsFactors = FALSE); asg$gene <- toupper(asg$gene)
MOD <- list(Sense=asg$gene[asg$module=="M1"], Traffic=asg$gene[asg$module=="M2"],
            Stromal_ECM=asg$gene[asg$module=="M3"], Angio_CAF=asg$gene[asg$module=="M4"],
            Execute=asg$gene[asg$module=="M5"], Checkpoint=asg$gene[asg$module=="M6"])
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
prop_high <- function(genes) { g <- intersect(genes, rownames(E)); if (!length(g)) return(rep(NA, ncol(E)))
  b <- t(apply(E[g,,drop=FALSE], 1, function(v) as.integer(v >= median(v, na.rm=TRUE)))); colMeans(b, na.rm=TRUE) }
Z <- function(v) (v - mean(v, na.rm=TRUE)) / sd(v, na.rm=TRUE)
M <- sapply(MOD, prop_high)
Sense <- Z(M[,"Sense"]); Integ <- Z(M[,"Traffic"] - (M[,"Stromal_ECM"]+M[,"Angio_CAF"])/2); Exec <- Z(M[,"Execute"]-M[,"Checkpoint"])
tiers <- cbind(Sense=Sense, Integration=Integ, Execution=Exec)
minz  <- apply(tiers, 1, min); which_min <- colnames(tiers)[apply(tiers, 1, which.min)]
FLOOR <- -0.5
intact <- minz >= FLOOR

cat("=== (A) cascade-intact split by ARM and outcome ===\n")
cat(sprintf("total cascade-intact: %d of %d (%.0f%%)\n", sum(intact), length(intact), 100*mean(intact)))
tab <- table(arm = arm, pCR = ifelse(y==1,"R","NR"))[, c("R","NR")]
tabi<- table(arm = arm[intact], pCR = ifelse(y[intact]==1,"R","NR"))[, c("R","NR")]
cat("\nAll eligible by arm x outcome:\n"); print(addmargins(tab))
cat("\nCascade-intact only, by arm x outcome:\n"); print(addmargins(tabi))
cat(sprintf("\nICB(pembro)-arm cascade-intact pCR rate = %.0f%%  (immune-competent -> responds to ICB)\n",
            100*mean(y[intact & arm=="pembro"])))

cat("\n=== framework's real job: explain ICB(pembro)-arm non-responders ===\n")
pnr <- arm=="pembro" & y==0
cat(sprintf("pembro-arm NRs: n=%d\n", sum(pnr)))
cat(sprintf("  of these, cascade-intact (unexplained): %d (%.0f%%)\n", sum(pnr & intact), 100*mean(intact[pnr])))
cat(sprintf("  with an identifiable limiting layer:     %d (%.0f%%)\n", sum(pnr & !intact), 100*mean(!intact[pnr])))
cat("  pembro-NR limiting-layer composition:\n")
kl <- ifelse(intact,"Cascade-intact", ifelse(which_min=="Sense","Sense-low",ifelse(which_min=="Integration","Integration-suppressed","Execution-low")))
print(table(factor(kl[pnr], c("Sense-low","Integration-suppressed","Execution-low","Cascade-intact"))))

cat("\n=== (B) threshold geometry: cascade-intact fraction vs FLOOR ===\n")
for (f in c(0, -0.25, -0.5, -0.75, -1)) {
  obs <- mean(minz >= f)
  exp_indep <- (1-pnorm(f))^3   # 3 independent standard-normal layers
  cat(sprintf("  FLOOR=%+.2f : observed intact %.0f%% | expected if 3 indep layers %.0f%%\n", f, 100*obs, 100*exp_indep)) }

cat("\n=== (C) borderline: how deep below the floor are the limiting tumours? ===\n")
cat(sprintf("cascade-intact min-layer z: median %.2f (all sit just above the %.2f floor)\n", median(minz[intact]), FLOOR))
cat(sprintf("limiting tumours min-layer z: median %.2f (clearly deficient)\n", median(minz[!intact])))
near <- intact & minz < (FLOOR+0.25)
cat(sprintf("cascade-intact within 0.25 SD of the floor (near-miss): %d of %d (%.0f%%)\n", sum(near), sum(intact), 100*mean(near[intact])))
