# =====================================================================
# Benchmark: ICI Integration (single-patient) vs established RNA biomarkers,
# scored on the SAME patients / SAME expression per cohort (icb_eval $bench + $scores).
#   competitors: GEP (T-cell-inflamed, Ayers), IFNG (Ayers 6-gene),
#                CYT (Rooney), CD274 (PD-L1 mRNA)
# Two airtight comparisons across the VALIDATED (anti-PD-(L)1) bucket:
#   (1) pooled AUC per biomarker — inverse-variance DerSimonian-Laird random
#       effects on the logit-AUC scale (identical math to run_pooled_auc.R);
#   (2) per-cohort PAIRED DeLong test (ICI vs each competitor, same patients),
#       summarised as wins/ties/losses + pooled paired dAUC.
# Numbers computed here, never typed. Writes benchmark tables for the figure.
# =====================================================================
suppressPackageStartupMessages({ library(pROC) })
R <- readRDS("pancancer_validation.rds")
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")

BIOM <- c(ICI="ICI Integration", GEP="T-cell-inflamed GEP", IFNG="IFN-γ signature",
          CYT="Cytolytic activity", CD274="PD-L1 mRNA (CD274)")

## ---- per-cohort AUC + 95% CI for every biomarker (from stored fields) -------
# ICI uses int_ci_single (DeLong); competitors use the $bench matrix (DeLong, same).
ac <- function(k, b){ if(b=="ICI") R[[k]]$int_ci_single else R[[k]]$bench[b, c("auc","lo","hi")] }
PC <- do.call(rbind, lapply(WORKS, function(k) do.call(rbind, lapply(names(BIOM), function(b){
  v <- ac(k,b); data.frame(cohort=k, cancer=R[[k]]$cancer, n=R[[k]]$n, biom=b,
                           auc=v[1], lo=v[2], hi=v[3], row.names=NULL) }))))
write.csv(PC, "../results/tables/ICI_benchmark_per_cohort_auc.csv", row.names=FALSE)

## ---- pooled AUC per biomarker: DerSimonian-Laird random effects (logit) -----
lg<-function(a) log(a/(1-a)); ilg<-function(x) 1/(1+exp(-x)); se_logit<-function(a,se) se/(a*(1-a))
pool <- function(b){
  d <- PC[PC$biom==b, ]; a<-d$auc; se_a<-(d$hi-d$lo)/(2*1.96)
  d <- data.frame(y=lg(a), sey=se_logit(a, se_a)); d<-d[is.finite(d$y)&is.finite(d$sey)&d$sey>0,]
  k<-nrow(d); wf<-1/d$sey^2; mf<-sum(wf*d$y)/sum(wf)
  Q<-sum(wf*(d$y-mf)^2); C<-sum(wf)-sum(wf^2)/sum(wf); tau2<-max(0,(Q-(k-1))/C)
  wr<-1/(d$sey^2+tau2); mr<-sum(wr*d$y)/sum(wr); ser<-sqrt(1/sum(wr))
  I2<-ifelse(Q>k-1, 100*(Q-(k-1))/Q, 0)
  c(auc=ilg(mr), lo=ilg(mr-1.96*ser), hi=ilg(mr+1.96*ser), I2=I2, k=k) }
POOL <- as.data.frame(t(sapply(names(BIOM), pool))); POOL$biom<-rownames(POOL); POOL$label<-BIOM[POOL$biom]
POOL <- POOL[order(-POOL$auc), ]
write.csv(POOL[,c("biom","label","auc","lo","hi","I2","k")],
          "../results/tables/ICI_benchmark_pooled_auc.csv", row.names=FALSE)

## ---- per-cohort PAIRED DeLong: ICI vs each competitor (same patients) -------
paircol <- c(GEP="GEP", IFNG="IFNG", CYT="CYT", CD274="CD274")
PAIR <- do.call(rbind, lapply(WORKS, function(k){ s<-R[[k]]$scores
  r1 <- roc(s$y, s$integ_single, quiet=TRUE, direction="<")
  do.call(rbind, lapply(names(paircol), function(b){
    r2 <- roc(s$y, s[[paircol[b]]], quiet=TRUE, direction="<")
    p  <- tryCatch(roc.test(r1, r2, method="delong", paired=TRUE)$p.value, error=function(e) NA_real_)
    data.frame(cohort=k, comp=b, dAUC=as.numeric(r1$auc)-as.numeric(r2$auc), p=p, row.names=NULL) })) }))
write.csv(PAIR, "../results/tables/ICI_benchmark_paired_delong.csv", row.names=FALSE)

## ---- summary ----------------------------------------------------------------
cat("=== Pooled AUC across", length(WORKS), "validated cohorts (random effects) ===\n")
for(i in seq_len(nrow(POOL))) cat(sprintf("  %-22s %.3f [%.3f-%.3f]  I2=%2.0f%%\n",
  POOL$label[i], POOL$auc[i], POOL$lo[i], POOL$hi[i], POOL$I2[i]))
cat("\n=== ICI vs competitor: per-cohort wins (dAUC>0) / ties / losses ===\n")
for(b in names(paircol)){ d<-PAIR[PAIR$comp==b,]
  cat(sprintf("  ICI vs %-6s  win %2d  loss %2d  (mean dAUC %+.3f, %d/%d cohorts favour ICI)\n",
    b, sum(d$dAUC>0), sum(d$dAUC<0), mean(d$dAUC), sum(d$dAUC>0), nrow(d))) }
cat("\ndone — benchmark tables written\n")
