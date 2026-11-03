# =====================================================================
# Added (incremental) value: does the ICI Integration score predict ICB response
# BEYOND each established RNA biomarker, and vice-versa?  Cohort-stratified
# logistic regression over the 11 validated cohorts (cohort as a factor absorbs
# between-cohort response-rate / platform differences); every predictor z-scored
# WITHIN its cohort so ORs are per-SD and comparable. Nested likelihood-ratio test
# for the added term. Symmetric both ways so the reader sees ICI adds signal the
# inflammation markers miss, while they add little on top of ICI.
# Reads per-sample scores from pancancer_validation.rds ($scores). Numbers computed
# here, never typed. Writes the added-value table.
# =====================================================================
R <- readRDS("pancancer_validation.rds")
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
COMP <- c(GEP="T-cell-inflamed GEP", IFNG="IFN-γ signature",
          CYT="Cytolytic activity", CD274="PD-L1 mRNA (CD274)")

## ---- pooled per-patient frame: predictors z-scored WITHIN cohort ------------
zc <- function(v) { s<-sd(v, na.rm=TRUE); if(is.na(s)||s==0) v*0 else (v-mean(v,na.rm=TRUE))/s }
D <- do.call(rbind, lapply(WORKS, function(k){ s<-R[[k]]$scores
  data.frame(cohort=k, y=s$y, ICI=zc(s$integ_single),
             GEP=zc(s$GEP), IFNG=zc(s$IFNG), CYT=zc(s$CYT), CD274=zc(s$CD274)) }))
D$cohort <- factor(D$cohort)
cat(sprintf("pooled n=%d across %d cohorts (R %d / NR %d)\n", nrow(D), nlevels(D$cohort), sum(D$y), sum(D$y==0)))

## ---- nested LRT helper: does `add` improve y ~ cohort + base ? ---------------
# returns added-term OR per SD [95% CI], LRT chi2 (df=1) and p.
lrt_add <- function(base, add){
  f0 <- glm(reformulate(c("cohort", base), "y"), binomial, D)
  f1 <- glm(reformulate(c("cohort", base, add), "y"), binomial, D)
  b <- coef(summary(f1))[add, ]; or<-exp(b["Estimate"]); se<-b["Std. Error"]
  chi <- deviance(f0) - deviance(f1); p <- pchisq(chi, df=1, lower.tail=FALSE)
  c(OR=unname(or), lo=unname(exp(b["Estimate"]-1.96*se)), hi=unname(exp(b["Estimate"]+1.96*se)),
    chi2=unname(chi), p=unname(p)) }

## ---- (1) ICI added on top of each competitor; (2) competitor on top of ICI ---
AV <- do.call(rbind, lapply(names(COMP), function(c){
  a <- lrt_add(base=c,     add="ICI")      # ICI beyond competitor
  b <- lrt_add(base="ICI", add=c)          # competitor beyond ICI
  rbind(
    data.frame(base=COMP[c], added="ICI Integration", added_key="ICI",
               OR=a["OR"], lo=a["lo"], hi=a["hi"], chi2=a["chi2"], p=a["p"], row.names=NULL),
    data.frame(base="ICI Integration", added=COMP[c], added_key=c,
               OR=b["OR"], lo=b["lo"], hi=b["hi"], chi2=b["chi2"], p=b["p"], row.names=NULL)) }))
AV$padj <- p.adjust(AV$p, "BH")
write.csv(AV, "../results/tables/ICI_added_value.csv", row.names=FALSE)

## ---- report ----------------------------------------------------------------
star <- function(p) ifelse(p<1e-3,"***",ifelse(p<1e-2,"**",ifelse(p<0.05,"*","ns")))
cat("\n=== ICI Integration added ON TOP OF each biomarker (cohort-stratified LRT) ===\n")
for(c in names(COMP)){ d<-AV[AV$base==COMP[c] & AV$added_key=="ICI",]
  cat(sprintf("  + ICI beyond %-20s OR/SD %.2f [%.2f-%.2f]  chi2=%5.1f  p=%.1e %s\n",
              COMP[c], d$OR, d$lo, d$hi, d$chi2, d$p, star(d$p))) }
cat("\n=== each biomarker added ON TOP OF ICI Integration (reverse) ===\n")
for(c in names(COMP)){ d<-AV[AV$base=="ICI Integration" & AV$added_key==c,]
  cat(sprintf("  + %-20s beyond ICI   OR/SD %.2f [%.2f-%.2f]  chi2=%5.1f  p=%.2f %s\n",
              COMP[c], d$OR, d$lo, d$hi, d$chi2, d$p, star(d$p))) }
cat("\ndone — added-value table written\n")
