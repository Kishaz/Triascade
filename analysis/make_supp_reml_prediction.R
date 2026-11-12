# =============================================================================
# REML sensitivity analysis + random-effects prediction interval for the pooled
# single-patient Integration response AUC (reviewers P5 / M33). Complements the
# primary DerSimonian-Laird (method-of-moments) estimate by (i) re-estimating the
# between-cohort variance tau^2 by restricted maximum likelihood and (ii)
# reporting the 95% prediction interval (the range for a new cohort's AUC).
# Per-cohort single-patient Integration AUC + DeLong SE, pooled on the logit scale.
# =============================================================================
suppressPackageStartupMessages(library(pROC))
R <- readRDS("pancancer_validation.rds")
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
lg <- function(a) log(a/(1-a)); ilg <- function(x) 1/(1+exp(-x))

# per-cohort logit-AUC and its SE (from the DeLong CI)
d <- do.call(rbind, lapply(WORKS, function(k){
  s <- R[[k]]$scores; r <- roc(s$y, s$integ_single, direction="<", quiet=TRUE, levels=c(0,1))
  ci <- as.numeric(ci.auc(r, method="delong"))               # lo, auc, hi
  se_auc <- (ci[3]-ci[1])/(2*qnorm(0.975))
  data.frame(cohort=k, auc=ci[2], y=lg(ci[2]), se=se_auc/(ci[2]*(1-ci[2]))) }))
y <- d$y; s2 <- d$se^2; k <- nrow(d)

# ---- DerSimonian-Laird (primary) -------------------------------------------
w <- 1/s2; yf <- sum(w*y)/sum(w); Q <- sum(w*(y-yf)^2)
C <- sum(w)-sum(w^2)/sum(w); tau2_dl <- max(0,(Q-(k-1))/C)
pool <- function(tau2){ wr<-1/(s2+tau2); yr<-sum(wr*y)/sum(wr); se<-sqrt(1/sum(wr)); list(yr=yr,se=se,wr=wr) }
DL <- pool(tau2_dl); I2 <- max(0,100*(Q-(k-1))/Q)

# ---- REML tau^2 (iterative) -------------------------------------------------
tau2 <- tau2_dl
for (it in 1:200){ wr<-1/(s2+tau2); yr<-sum(wr*y)/sum(wr)
  tau2_new <- sum(wr^2*((y-yr)^2 - s2))/sum(wr^2) + 1/sum(wr)
  tau2_new <- max(0, tau2_new); if (abs(tau2_new-tau2)<1e-8) break; tau2 <- tau2_new }
RE <- pool(tau2)

# ---- 95% prediction interval (Higgins-Thompson-Spiegelhalter, REML) ---------
tcrit <- qt(0.975, df=k-2)
pi_lo <- RE$yr - tcrit*sqrt(tau2 + RE$se^2); pi_hi <- RE$yr + tcrit*sqrt(tau2 + RE$se^2)

TAB <- data.frame(
  model = c("DerSimonian-Laird (primary)","REML (sensitivity)"),
  tau2  = round(c(tau2_dl, tau2), 4),
  pooled_AUC = round(ilg(c(DL$yr, RE$yr)), 3),
  CI_low  = round(ilg(c(DL$yr-1.96*DL$se, RE$yr-1.96*RE$se)), 3),
  CI_high = round(ilg(c(DL$yr+1.96*DL$se, RE$yr+1.96*RE$se)), 3),
  I2_pct = round(c(I2, I2), 0))
out <- "../results/tables/ICI_SuppTable_reml_prediction.csv"
write.csv(TAB, out, row.names=FALSE)
cat("wrote", out, "\n\n"); print(TAB, row.names=FALSE)
cat(sprintf("\nREML pooled Integration AUC = %.3f (95%% CI %.3f-%.3f), I2 = %.0f%%\n",
            ilg(RE$yr), ilg(RE$yr-1.96*RE$se), ilg(RE$yr+1.96*RE$se), I2))
cat(sprintf("95%% prediction interval (REML): %.2f to %.2f (AUC scale)\n", ilg(pi_lo), ilg(pi_hi)))
