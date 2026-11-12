# =============================================================================
# Supplementary Table 4 (addendum): pooled single-patient response AUC with 95%
# CI for ALL THREE ICI layers (Sensory / Integration / Execution), so the three
# layers can be compared quantitatively rather than described as weaker
# (reviewer R11 / comment C47). Per-cohort single-patient AUC + DeLong 95% CI
# for each layer, pooled across the 11 anti-PD-(L)1 cohorts by the SAME logit-
# scale DerSimonian-Laird random-effects model used for Integration in Fig 3b.
# Reads pancancer_validation.rds. Editable table; values never typed.
# =============================================================================
suppressPackageStartupMessages({ library(pROC) })
R <- readRDS("pancancer_validation.rds")

WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
DISP <- c(IMvigor210="IMvigor210", GSE78220="Hugo 2016", GSE91061="Riaz 2017",
          GSE126044="Cho 2020", GSE135222="Jung 2019", Gide2019="Gide 2019",
          Liu2019="Liu 2019", Puch="Puch", Kim2018gas="Kim 2018", Rose="Rose",
          ISPY2="I-SPY2 (breast)")
LAYERS <- c(Sensory="sens_single", Integration="integ_single", Execution="exec_single")

# per-cohort single-patient AUC + DeLong 95% CI for one layer (higher = responder)
aucci <- function(y, score){
  ok <- is.finite(score) & is.finite(y)
  r  <- roc(y[ok], score[ok], direction="<", quiet=TRUE, levels=c(0,1))
  ci <- as.numeric(ci.auc(r, method="delong"))          # c(lo, auc, hi)
  c(auc=ci[2], lo=ci[1], hi=ci[3])
}

# logit-scale DerSimonian-Laird random-effects pool (matches make_fig3b_auc.R)
se_a <- function(lo,hi)(hi-lo)/(2*qnorm(0.975)); lg<-function(a)log(a/(1-a)); ilg<-function(x)1/(1+exp(-x))
pool <- function(auc,lo,hi){
  y<-lg(auc); sy<-se_a(lo,hi)/(auc*(1-auc)); w<-1/sy^2
  yf<-sum(w*y)/sum(w); Q<-sum(w*(y-yf)^2); k<-length(auc); C<-sum(w)-sum(w^2)/sum(w)
  tau2<-max(0,(Q-(k-1))/C); wr<-1/(sy^2+tau2); yr<-sum(wr*y)/sum(wr); ser<-sqrt(1/sum(wr))
  list(auc=ilg(yr), lo=ilg(yr-1.96*ser), hi=ilg(yr+1.96*ser), I2=max(0,100*(Q-(k-1))/Q))
}

# ---- per-cohort table (all three layers, AUC + CI) --------------------------
fmtci <- function(v) sprintf("%.2f [%.2f-%.2f]", v["auc"], v["lo"], v["hi"])
rows <- lapply(WORKS, function(k){
  x <- R[[k]]; s <- x$scores; y <- s$y
  L <- lapply(LAYERS, function(col) aucci(y, s[[col]]))
  data.frame(cohort=DISP[[k]], cancer=x$cancer, n=x$n, nR=x$nR, nNR=x$nNR,
             Sensory_AUC=fmtci(L$Sensory), Integration_AUC=fmtci(L$Integration),
             Execution_AUC=fmtci(L$Execution), stringsAsFactors=FALSE,
             .s=L$Sensory["auc"], .sl=L$Sensory["lo"], .sh=L$Sensory["hi"],
             .i=L$Integration["auc"], .il=L$Integration["lo"], .ih=L$Integration["hi"],
             .e=L$Execution["auc"], .el=L$Execution["lo"], .eh=L$Execution["hi"]) })
TAB <- do.call(rbind, rows)

# ---- pooled row per layer ---------------------------------------------------
pS <- pool(TAB$.s, TAB$.sl, TAB$.sh)
pI <- pool(TAB$.i, TAB$.il, TAB$.ih)
pE <- pool(TAB$.e, TAB$.el, TAB$.eh)
pooled <- data.frame(cohort="Pooled (11 anti-PD-(L)1, random-effects)", cancer="-",
  n=sum(TAB$n), nR=sum(TAB$nR), nNR=sum(TAB$nNR),
  Sensory_AUC=fmtci(c(auc=pS$auc,lo=pS$lo,hi=pS$hi)),
  Integration_AUC=fmtci(c(auc=pI$auc,lo=pI$lo,hi=pI$hi)),
  Execution_AUC=fmtci(c(auc=pE$auc,lo=pE$lo,hi=pE$hi)), stringsAsFactors=FALSE)
I2row <- data.frame(cohort="I-squared (%)", cancer="-", n=NA, nR=NA, nNR=NA,
  Sensory_AUC=sprintf("%.0f", pS$I2), Integration_AUC=sprintf("%.0f", pI$I2),
  Execution_AUC=sprintf("%.0f", pE$I2), stringsAsFactors=FALSE)

OUT <- rbind(TAB[order(-TAB$.i), c("cohort","cancer","n","nR","nNR","Sensory_AUC","Integration_AUC","Execution_AUC")],
             pooled, I2row)
out <- "../results/tables/ICI_SuppTable4_pooled_layer_auc.csv"
write.csv(OUT, out, row.names=FALSE)

cat("wrote", out, "\n\n")
cat(sprintf("POOLED single-patient response AUC (11 anti-PD-(L)1 cohorts, %d patients):\n", sum(TAB$n)))
cat(sprintf("  Sensory      %.3f [%.3f-%.3f]  I2=%.0f%%\n", pS$auc, pS$lo, pS$hi, pS$I2))
cat(sprintf("  Integration  %.3f [%.3f-%.3f]  I2=%.0f%%   <- headline (expect ~0.72 [0.66-0.78], I2~42)\n", pI$auc, pI$lo, pI$hi, pI$I2))
cat(sprintf("  Execution    %.3f [%.3f-%.3f]  I2=%.0f%%\n", pE$auc, pE$lo, pE$hi, pE$I2))
