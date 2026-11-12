# =============================================================================
# Reference and scoring sensitivity analyses (reviewers P6 / M14; the numeric
# results CY's edit moved out of the Methods now have a reproducible home).
# In IMvigor210 (full transcriptomes available): (i) random gene dropout of 10%
# and 20% of the 98-gene panel -> Pearson correlation of the single-patient
# Integration score vs full coverage, response AUC, and archetype concordance;
# (ii) reference portability -> score IMvigor210 against the frozen NeoTRIP
# reference vs a reference refit from IMvigor210 itself.
# Single-patient selfnorm scoring (score_locked, mode="selfnorm").
# =============================================================================
suppressPackageStartupMessages({ library(pROC); library(data.table) })
source("ici_locked_reference.R")
assets <- readRDS("../app/references/ICI_v1.rds"); REF <- assets$ref; PANEL <- assets$panel
o <- readRDS("data/imvigor210/imvigor210_clean.rds"); E <- as.matrix(o$E); rownames(E) <- toupper(rownames(o$E))
y01 <- ifelse(o$pd$binaryResponse=="CR/PR",1L, ifelse(o$pd$binaryResponse=="SD/PD",0L,NA))
keep <- !is.na(y01); E <- E[, keep]; y <- y01[keep]
PG <- intersect(REF$genes, rownames(E))                     # panel genes present
cat(sprintf("IMvigor210: %d samples with response (R %d / NR %d); %d/%d panel genes present\n",
            ncol(E), sum(y), sum(y==0), length(PG), length(REF$genes)))

layers_of <- function(Emat, ref=REF) score_locked(Emat, ref, "selfnorm")
auc1 <- function(s) as.numeric(auc(y, s, direction="<", quiet=TRUE, levels=c(0,1)))
arch <- function(sc){ Z <- apply(sc[,c("Sensory","Integration","Execution")], 2, function(v) (v-mean(v))/sd(v))
  w <- c("Sensory","Integration","Execution")[max.col(-Z, ties.method="first")]
  ifelse(apply(Z,1,min) >= -0.5, "intact",
    ifelse(w=="Sensory","cold", ifelse(w=="Integration","blocked","execution-impaired"))) }

full <- layers_of(E); int_full <- full$Integration; arch_full <- arch(full)
cat(sprintf("Full coverage: Integration AUC = %.3f\n\n", auc1(int_full)))

# ---- (i) random gene dropout (mean over 25 seeds) ---------------------------
dropout <- function(frac, nseed=25){
  st <- t(sapply(seq_len(nseed), function(s){ set.seed(s)
    drop <- sample(PG, round(frac*length(PG)))
    Ed <- E[setdiff(rownames(E), drop), , drop=FALSE]
    sc <- layers_of(Ed)
    c(r = cor(int_full, sc$Integration, use="complete.obs"),
      auc = auc1(sc$Integration),
      conc = mean(arch(sc) == arch_full)) }))
  colMeans(st) }
d10 <- dropout(0.10); d20 <- dropout(0.20)

# ---- (ii) reference portability: NeoTRIP vs IMvigor210-derived reference ----
ALT <- fit_locked_reference(E[PG, , drop=FALSE], PANEL, "IMvigor210-derived")
alt <- score_locked(E, ALT, "selfnorm")
port <- c(r = cor(int_full, alt$Integration, use="complete.obs"),
          auc_neotrip = auc1(int_full), auc_alt = auc1(alt$Integration),
          conc = mean(arch(alt) == arch_full))

TAB <- data.frame(
  analysis = c("10% gene dropout","20% gene dropout","Reference portability (IMvigor210-derived)"),
  Integration_r = round(c(d10["r"], d20["r"], port["r"]), 2),
  response_AUC  = round(c(d10["auc"], d20["auc"], port["auc_alt"]), 2),
  archetype_concordance_pct = round(100*c(d10["conc"], d20["conc"], port["conc"])),
  note = c("vs full-coverage AUC 0.61","vs full-coverage AUC 0.61",
           sprintf("NeoTRIP AUC %.2f vs IMvigor210-ref AUC %.2f", port["auc_neotrip"], port["auc_alt"])))
out <- "../results/tables/ICI_SuppTable_scoring_sensitivity.csv"
write.csv(TAB, out, row.names=FALSE)
cat("wrote", out, "\n\n"); print(TAB[,1:4], row.names=FALSE)
cat(sprintf("\nFull-coverage Integration AUC = %.2f; portability NeoTRIP %.2f vs IMvigor-ref %.2f\n",
            auc1(int_full), port["auc_neotrip"], port["auc_alt"]))
