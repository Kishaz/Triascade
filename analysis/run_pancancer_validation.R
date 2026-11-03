# =====================================================================
# Pan-cancer ICB validation of the legacy-98 ICI framework.
# ONE evaluator applied identically to every cohort (airtight by construction):
#   - gene IDs upper-cased + de-duplicated; panel coverage reported
#   - BASELINE / pre-treatment samples only (enforced per-cohort loader)
#   - response = CR/PR (or durable benefit) vs the rest, binary y
#   - scoring: (a) cohort-median  (b) single-patient selfnorm, reference frozen
#     on the cohort itself  (c) single-patient selfnorm, EXTERNAL NeoTRIP-breast
#     reference (cross-cancer deployability)
#   - AUC with PRE-SPECIFIED direction (higher layer score = responder; never flipped)
#   - Integration AUC 95% CI (DeLong) for cohort-median and single-patient
#   - NR limiting-layer taxonomy
# Numbers are computed here, never typed.
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")
PANEL <- load_panel("ICI_98gene_panel.csv"); PANEL_GENES <- unique(unlist(PANEL))

# external reference: NeoTRIP breast baseline, fit once
nt <- as.data.frame(fread("data/NeoTRIP_expr.txt.gz")); rownames(nt) <- toupper(nt$HGNC_approved_symbol)
REF_EXT <- fit_locked_reference(as.matrix(nt[, grep("_Baseline$", colnames(nt))]), PANEL, "NeoTRIP_baseline")

# ---------------------------------------------------------------------
# Competitor RNA biomarkers, scored on the SAME E / SAME patients as the ICI
# panel (apples-to-apples). Each = mean of within-cohort z-scored signature genes
# (platform-agnostic, matches the selfnorm spirit); direction PRE-SPECIFIED
# "immune-hot = responder" (higher = responder, never flipped), as for the panel.
#   GEP  = T-cell-inflamed Gene Expression Profile (Ayers 2017, 18-gene / TIS)
#   IFNG = IFN-gamma signature (Ayers 2017, 6-gene)
#   CYT  = cytolytic activity (Rooney 2015: GZMA, PRF1)
#   CD274= PD-L1 mRNA (single-gene RNA proxy for PD-L1 IHC)
# TMB and PD-L1 IHC are NOT RNA-derivable and are absent from these cohorts ->
# cannot be scored on the same patients (reported as a limitation, never faked).
BENCH_SETS <- list(
  GEP  = c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6",
           "HLA.DQA1","HLA-DQA1","HLA.DRB1","HLA-DRB1","HLA.E","HLA-E",
           "IDO1","LAG3","NKG7","PDCD1LG2","PSMB10","STAT1","TIGIT"),
  IFNG = c("IDO1","CXCL10","CXCL9","HLA.DRA","HLA-DRA","STAT1","IFNG"),
  CYT  = c("GZMA","PRF1"),
  CD274= c("CD274"))
score_bench <- function(E){                        # E: genes x samples (rownames upper-cased)
  zmean <- function(genes){ g<-intersect(toupper(genes), rownames(E))
    if(length(g)==0) return(list(s=rep(NA_real_, ncol(E)), cov=0L, k=0L))
    M <- E[g,,drop=FALSE]
    Z <- t(apply(M, 1, function(v){ s<-sd(v); if(is.na(s)||s==0) rep(0,length(v)) else (v-mean(v))/s }))
    list(s=colMeans(Z, na.rm=TRUE), cov=length(g), k=length(unique(toupper(genes)))) }
  lapply(BENCH_SETS, zmean) }

icb_eval <- function(E, y, cohort, cancer, drug) {
  E <- as.matrix(E); storage.mode(E) <- "numeric"
  rownames(E) <- toupper(rownames(E)); E <- E[!duplicated(rownames(E)), , drop=FALSE]
  ok <- !is.na(y); E <- E[, ok, drop=FALSE]; y <- as.integer(y[ok])
  cov <- length(intersect(PANEL_GENES, rownames(E)))
  ord <- function(sc) sc[match(colnames(E), sc$sample), ]
  cm  <- ord(score_cohort_median(E, PANEL))
  sp  <- ord(score_locked(E, fit_locked_reference(E, PANEL, "self"), mode="selfnorm"))
  spx <- ord(score_locked(E, REF_EXT, mode="selfnorm"))
  auc1 <- function(s) tryCatch(as.numeric(auc(y, s, quiet=TRUE, direction="<")), error=function(e) NA_real_)
  cii  <- function(s){ r<-roc(y,s,quiet=TRUE,direction="<"); c(as.numeric(r$auc), ci.auc(r)[1], ci.auc(r)[3]) }
  layers <- c("Sensory","Integration","Execution")
  A <- sapply(layers, function(L) c(cohort=auc1(cm[[L]]), single_self=auc1(sp[[L]]), single_ext=auc1(spx[[L]])))
  ## ---- competitor RNA biomarkers on the SAME patients (same E, same y) ----
  bset <- score_bench(E)                                   # list per biomarker: s, cov, k
  bciF <- function(s) if(all(is.na(s)) || length(unique(s[!is.na(s)]))<2) c(NA,NA,NA) else cii(s)
  BEN  <- t(sapply(names(bset), function(b){ c(bciF(bset[[b]]$s), cov=bset[[b]]$cov, k=bset[[b]]$k) }))
  colnames(BEN) <- c("auc","lo","hi","cov","k")
  Z <- function(v) (v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE)
  tiers <- cbind(Sensory=Z(cm$Sensory),Integration=Z(cm$Integration),Execution=Z(cm$Execution))
  wmin <- colnames(tiers)[apply(tiers,1,which.min)]; minz<-apply(tiers,1,min)
  kl <- ifelse(minz>=-0.5,"intact",ifelse(wmin=="Sensory","cold",ifelse(wmin=="Integration","blocked","no-kill")))
  cat(sprintf("\n### %-16s %-11s %-12s n=%d (R %d/NR %d, %2.0f%%)  panel %d/98\n",
              cohort, cancer, drug, length(y), sum(y), sum(y==0), 100*mean(y), cov))
  cat(sprintf("  %-12s cohort  single(self)  single(ext)\n",""))
  for (L in layers) cat(sprintf("  %-12s %.3f     %.3f        %.3f\n", L, A["cohort",L],A["single_self",L],A["single_ext",L]))
  ci<-cii(cm$Integration); cs<-cii(sp$Integration)
  cat(sprintf("  Integration CI: cohort %.3f[%.2f-%.2f]  single %.3f[%.2f-%.2f]\n", ci[1],ci[2],ci[3],cs[1],cs[2],cs[3]))
  cat("  NR limiting layer: "); print(table(factor(kl[y==0],c("cold","blocked","no-kill","intact"))))
  cat("  bench AUC (same patients): ");
  cat(paste(sprintf("%s %.2f", rownames(BEN), BEN[,"auc"]), collapse="  "), "\n")
  list(cohort=cohort,cancer=cancer,drug=drug,n=length(y),nR=sum(y),nNR=sum(y==0),cov=cov,
       auc=A,int_ci=ci,int_ci_single=cs,bench=BEN,
       nr_comp=table(factor(kl[y==0],c("cold","blocked","no-kill","intact"))),
       scores=data.frame(sample=as.character(sp$sample), y=y,
                         sens_single =as.numeric(sp$Sensory),     sens_cohort =as.numeric(cm$Sensory),
                         integ_single=as.numeric(sp$Integration), integ_cohort=as.numeric(cm$Integration),
                         exec_single =as.numeric(sp$Execution),   exec_cohort =as.numeric(cm$Execution),
                         GEP=bset$GEP$s, IFNG=bset$IFNG$s, CYT=bset$CYT$s, CD274=bset$CD274$s,
                         limiting    =kl,                         # per-patient limiting-layer / archetype class
                         row.names=NULL))
}

RES <- list()
## ---- IMvigor210 (urothelial, atezolizumab monotherapy) ----
obj<-readRDS("data/imvigor210/imvigor210_clean.rds")
rsp<-obj$pd[["binaryResponse"]]; kp<-rsp%in%c("CR/PR","SD/PD")
RES$IMvigor210 <- icb_eval(obj$E[,kp], as.integer(rsp[kp]=="CR/PR"), "IMvigor210","urothelial","atezo-mono")

## ---- GSE78220 (melanoma, Hugo 2016, anti-PD1; baseline biopsies) ----
ex<-read.csv("data/GSE78220_expr.csv",row.names=1,check.names=FALSE)
rp<-read.csv("data/GSE78220_response.csv",stringsAsFactors=FALSE); rp<-rp[match(colnames(ex),rp$sample),]
RES$GSE78220 <- icb_eval(ex, rp$responder, "GSE78220","melanoma","anti-PD1")

saveRDS(RES, "pancancer_validation.rds")
cat("\nsaved pancancer_validation.rds\n")
