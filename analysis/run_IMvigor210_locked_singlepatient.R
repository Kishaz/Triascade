# =====================================================================
# IMvigor210 single-patient locked scoring (legacy 98-gene panel).
# Does the ICI framework discriminate ICB responders in a PURE-ICB,
# non-chemo cohort at the SINGLE-PATIENT level?
#   - atezolizumab MONOTHERAPY, urothelial carcinoma (no chemo backbone
#     -> removes the I-SPY2 pCR chemo/proliferation confound)
#   - n=298 evaluable (CR/PR vs SD/PD) -> finally powers the single-patient claim
#   - different cancer -> generalizability
# Three scorings compared, per layer (S/I/E), by AUC for CR/PR vs SD/PD:
#   (a) cohort-median        (reference / upper bound; needs the cohort)
#   (b) single-patient self  (selfnorm, reference frozen on IMvigor210 itself)
#   (c) single-patient ext   (selfnorm, reference frozen on EXTERNAL NeoTRIP
#                             breast baseline -> cross-cancer deployability)
# Also: NR limiting-layer taxonomy (why non-responders fail).  Numbers from data.
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")
PANEL <- load_panel("ICI_98gene_panel.csv")

obj <- readRDS("data/imvigor210/imvigor210_clean.rds"); E <- obj$E; pd <- obj$pd
rownames(E) <- toupper(rownames(E)); E <- E[!duplicated(rownames(E)), ]
resp <- pd[["binaryResponse"]]
keep <- !is.na(resp) & resp %in% c("CR/PR","SD/PD")
E <- as.matrix(E[, keep]); y <- as.integer(resp[keep]=="CR/PR")
cat(sprintf("IMvigor210 (atezo mono, urothelial): evaluable n=%d  (R=%d / NR=%d, %.0f%% response)\n",
            length(y), sum(y), sum(y==0), 100*mean(y)))
cat(sprintf("panel genes present in IMvigor210: %d / %d\n",
            length(intersect(unique(unlist(PANEL)), rownames(E))), length(unique(unlist(PANEL)))))

ord <- function(sc) sc[match(colnames(E), sc$sample), ]
auc1 <- function(score) tryCatch(as.numeric(auc(y, score, quiet=TRUE, direction="<")), error=function(e) NA)

# (a) cohort-median
cm <- ord(score_cohort_median(E, PANEL))
# (b) single-patient, reference frozen on IMvigor210 itself
ref_self <- fit_locked_reference(E, PANEL, "IMvigor210_self")
sp_self  <- ord(score_locked(E, ref_self, mode="selfnorm"))
# (c) single-patient, reference frozen on EXTERNAL NeoTRIP breast baseline
nt <- as.data.frame(fread("data/NeoTRIP_expr.txt.gz")); rownames(nt) <- toupper(nt$HGNC_approved_symbol)
ntb <- as.matrix(nt[, grep("_Baseline$", colnames(nt))])
ref_ext <- fit_locked_reference(ntb, PANEL, "NeoTRIP_baseline")
sp_ext   <- ord(score_locked(E, ref_ext, mode="selfnorm"))

cat("\n== AUC for ICB response (CR/PR vs SD/PD), per layer ==\n")
cat(sprintf("%-12s  cohort-median   single-pt(self)   single-pt(NeoTRIP ext ref)\n", "layer"))
L3 <- c("Sensory","Integration","Execution")
aucs <- sapply(L3, function(L) c(cm=auc1(cm[[L]]), sp_self=auc1(sp_self[[L]]), sp_ext=auc1(sp_ext[[L]])))
for (L in L3) cat(sprintf("%-12s     %.3f           %.3f            %.3f\n",
                          L, aucs["cm",L], aucs["sp_self",L], aucs["sp_ext",L]))

# ---- NR limiting-layer taxonomy (cohort-median z, weakest link) ----
Z <- function(v) (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE)
tiers <- cbind(Sensory=Z(cm$Sensory), Integration=Z(cm$Integration), Execution=Z(cm$Execution))
minz <- apply(tiers,1,min); wmin <- colnames(tiers)[apply(tiers,1,which.min)]
FLOOR <- -0.5
klass <- ifelse(minz>=FLOOR,"Cascade-intact",
         ifelse(wmin=="Sensory","Sense-low (cold)",
         ifelse(wmin=="Integration","Integration-suppressed","Execution-low")))
klass <- factor(klass, c("Sense-low (cold)","Integration-suppressed","Execution-low","Cascade-intact"))
cat(sprintf("\n== NR limiting-layer composition (n=%d non-responders) ==\n", sum(y==0)))
nrc <- table(klass[y==0]); print(nrc)
cat("\n== response rate (%) by limiting layer ==\n")
print(round(tapply(y, klass, function(z) 100*mean(z)),0))

saveRDS(list(cohort="IMvigor210 (atezo mono, urothelial)", n=length(y), n_R=sum(y), n_NR=sum(y==0),
             pct_resp=100*mean(y), auc=aucs, nr_comp=as.data.frame(nrc), klass=klass, y=y),
        "imvigor210_locked_singlepatient.rds")
cat("\nsaved imvigor210_locked_singlepatient.rds\ndone\n")
