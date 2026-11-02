# =====================================================================
# Pooled Integration AUC across the pan-cancer ICB validation.
# Reads pancancer_validation.rds (15 cohorts from run_pancancer_extra.R).
# Pools the SINGLE-PATIENT selfnorm Integration AUC (the paper's thesis metric;
# cohort-median also reported) by inverse-variance on the logit-AUC scale,
# DerSimonian-Laird random effects + fixed effect, with I^2 and 95% CI.
# Per-cohort SE backed out of the stored DeLong 95% CI.
# Buckets are declared explicitly (not parsed) so the pooling is auditable.
# =====================================================================
R <- readRDS("pancancer_validation.rds")

# ---- bucket membership (explicit & prespecified) --------------------
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose")   # pure PD-(L)1, immunogenic
CTLA4 <- c("VanAllen2015")                                    # CTLA4 monotherapy (reported apart)
NULLB <- c("GSE165252","PRINCE","Braun")                      # cold / chemo+ICB / neg-control
EXCL  <- c("Friedrich")                                       # bispecific, n=9, not ICB -> excluded

# ---- helpers --------------------------------------------------------
se_from_ci <- function(ci) (ci[3]-ci[2])/(2*qnorm(0.975))     # SE on AUC scale
lg  <- function(a) log(a/(1-a))                               # logit
ilg <- function(x) 1/(1+exp(-x))                              # inv-logit
se_logit <- function(a, se_a) se_a/(a*(1-a))                  # delta method

pool <- function(cohorts, metric=c("single","cohort")){
  metric <- match.arg(metric)
  fld <- if(metric=="single") "int_ci_single" else "int_ci"
  rows <- lapply(cohorts, function(k){
    x <- R[[k]]; ci <- x[[fld]]
    a <- ci[1]; sa <- se_from_ci(ci)
    data.frame(cohort=x$cohort, cancer=x$cancer, n=x$n, nR=x$nR, nNR=x$nNR,
               auc=a, lo=ci[2], hi=ci[3],
               y=lg(a), sey=se_logit(a, sa), stringsAsFactors=FALSE)
  })
  d <- do.call(rbind, rows)
  w  <- 1/d$sey^2
  # fixed effect
  yf <- sum(w*d$y)/sum(w); sef <- sqrt(1/sum(w))
  # DerSimonian-Laird random effects
  Q  <- sum(w*(d$y-yf)^2); k <- nrow(d)
  C  <- sum(w) - sum(w^2)/sum(w)
  tau2 <- max(0, (Q-(k-1))/C)
  wr <- 1/(d$sey^2 + tau2)
  yr <- sum(wr*d$y)/sum(wr); ser <- sqrt(1/sum(wr))
  I2 <- max(0, 100*(Q-(k-1))/Q)
  list(d=d, k=k, n=sum(d$n), nR=sum(d$nR), nNR=sum(d$nNR),
       fixed = c(auc=ilg(yf), lo=ilg(yf-1.96*sef), hi=ilg(yf+1.96*sef)),
       random= c(auc=ilg(yr), lo=ilg(yr-1.96*ser), hi=ilg(yr+1.96*ser)),
       tau2=tau2, I2=I2, Q=Q)
}

report <- function(tag, cohorts, metric){
  p <- pool(cohorts, metric)
  cat(sprintf("\n=== %s  [%s-patient Integration AUC]  k=%d, n=%d (R %d / NR %d) ===\n",
              tag, metric, p$k, p$n, p$nR, p$nNR))
  for(i in seq_len(nrow(p$d))) with(p$d[i,],
    cat(sprintf("  %-12s %-11s n=%-4d  %.3f [%.2f-%.2f]\n", cohort, cancer, n, auc, lo, hi)))
  cat(sprintf("  ---- fixed  : %.3f [%.3f-%.3f]\n", p$fixed["auc"], p$fixed["lo"], p$fixed["hi"]))
  cat(sprintf("  ---- random : %.3f [%.3f-%.3f]  (I2=%.0f%%, tau2=%.3f)\n",
              p$random["auc"], p$random["lo"], p$random["hi"], p$I2, p$tau2))
  invisible(p)
}

cat("############ POOLED INTEGRATION AUC ############\n")
pw_s <- report("WORKS (pure PD-(L)1)", WORKS, "single")
pw_c <- report("WORKS (pure PD-(L)1)", WORKS, "cohort")

# contrast buckets (not pooled into the works estimate)
cat("\n=== CTLA4 monotherapy (separate) ===\n")
for(k in CTLA4) with(R[[k]], cat(sprintf("  %-12s %-11s n=%-4d  single %.3f  cohort %.3f\n",
     cohort,cancer,n,int_ci_single[1],int_ci[1])))
cat("\n=== NULL bucket (cold / chemo+ICB / neg-control) ===\n")
for(k in NULLB) with(R[[k]], cat(sprintf("  %-12s %-11s n=%-4d  single %.3f  cohort %.3f\n",
     cohort,cancer,n,int_ci_single[1],int_ci[1])))
cat("\n=== excluded: ", paste(EXCL,collapse=", "), " (n<10, not checkpoint) ===\n", sep="")

# ---- write a tidy table ---------------------------------------------
bucket_of <- function(k) ifelse(k %in% WORKS,"works",
                         ifelse(k %in% CTLA4,"ctla4",
                         ifelse(k %in% NULLB,"null","excluded")))
TAB <- do.call(rbind, lapply(names(R), function(k){ x<-R[[k]]
  data.frame(cohort=x$cohort, cancer=x$cancer, drug=x$drug, bucket=bucket_of(k),
             n=x$n, nR=x$nR, nNR=x$nNR,
             integ_cohort=x$int_ci[1],        cohort_lo=x$int_ci[2],        cohort_hi=x$int_ci[3],
             integ_single=x$int_ci_single[1], single_lo=x$int_ci_single[2], single_hi=x$int_ci_single[3],
             row.names=NULL) }))
pooled_tab <- data.frame(
  bucket="works_pooled", metric=c("single","cohort"),
  k=pw_s$k, n=pw_s$n, nR=pw_s$nR, nNR=pw_s$nNR,
  auc_fixed =c(pw_s$fixed["auc"],  pw_c$fixed["auc"]),
  fixed_lo  =c(pw_s$fixed["lo"],   pw_c$fixed["lo"]),
  fixed_hi  =c(pw_s$fixed["hi"],   pw_c$fixed["hi"]),
  auc_random=c(pw_s$random["auc"], pw_c$random["auc"]),
  random_lo =c(pw_s$random["lo"],  pw_c$random["lo"]),
  random_hi =c(pw_s$random["hi"],  pw_c$random["hi"]),
  I2=c(pw_s$I2, pw_c$I2), row.names=NULL)
dir.create("../results/tables", showWarnings=FALSE, recursive=TRUE)
write.csv(TAB,        "../results/tables/pancancer_cohort_table.csv", row.names=FALSE)
write.csv(pooled_tab, "../results/tables/pooled_integration_auc.csv", row.names=FALSE)
cat("\nwrote ../results/tables/{pancancer_cohort_table,pooled_integration_auc}.csv\n")
