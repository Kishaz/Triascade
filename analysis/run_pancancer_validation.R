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
  list(cohort=cohort,cancer=cancer,drug=drug,n=length(y),nR=sum(y),nNR=sum(y==0),cov=cov,
       auc=A,int_ci=ci,int_ci_single=cs,nr_comp=table(factor(kl[y==0],c("cold","blocked","no-kill","intact"))),
       scores=data.frame(sample=as.character(sp$sample), y=y,
                         sens_single =as.numeric(sp$Sensory),     sens_cohort =as.numeric(cm$Sensory),
                         integ_single=as.numeric(sp$Integration), integ_cohort=as.numeric(cm$Integration),
                         exec_single =as.numeric(sp$Execution),   exec_cohort =as.numeric(cm$Execution),
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
