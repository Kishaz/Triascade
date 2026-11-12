# =============================================================================
# Trial-specific and pooled Integration-by-treatment interaction ODDS RATIOS
# with 95% CI (reviewer R23 / comment C70): a significant pooled interaction
# does not imply both trials have similar treatment-by-Integration effects, so
# report the per-trial interaction OR and CI, not only the P value.
# Reuses results/tables/ICI_predictive_checkpoint_patients.csv (per-patient
# within-trial z of each layer score, arm, pCR, trial) written by
# run_predictive_prognostic.R. Same model -> pooled Integration reproduces the
# manuscript OR 1.81 [1.04-3.16].
# =============================================================================
D <- read.csv("../results/tables/ICI_predictive_checkpoint_patients.csv", stringsAsFactors = FALSE)
D$arm <- factor(D$arm, c("ctrl","icb"))

# interaction OR (exp of the layer x armicb coefficient) + Wald 95% CI + P
fit_or <- function(df, L, pooled = FALSE){
  form <- if (pooled) reformulate(c(sprintf("%s*arm", L), "trial"), "y")
          else        reformulate(sprintf("%s*arm", L), "y")
  f  <- glm(form, binomial, df)
  co <- sprintf("%s:armicb", L)
  est <- coef(f)[[co]]; se <- sqrt(diag(vcov(f)))[[co]]       # Wald 95% CI
  c(OR = exp(est), lo = exp(est - 1.96*se), hi = exp(est + 1.96*se),
    p = 2*pnorm(-abs(est/se)))
}

# established likelihood-ratio interaction P values (headline; what the manuscript cites)
PT <- read.csv("../results/tables/ICI_predictive_interaction.csv", stringsAsFactors=FALSE)
PT <- PT[PT$scoring=="locked single-patient", ]
PC <- read.csv("../results/tables/ICI_predictive_interaction_continuous.csv", stringsAsFactors=FALSE)
lrt_p <- function(L, tr){
  if (tr=="Pooled (trial-stratified)") return(PC$interaction_p[PC$layer==L])
  key <- if (tr=="I-SPY2") "I-SPY2 (pembro vs pac)" else "GeparNuevo (durva vs placebo)"
  PT$interaction_p[PT$trial==key & PT$layer==L] }

LAYERS <- c("Sensory","Integration","Execution")
rows <- list()
for (L in LAYERS){
  for (tr in c("I-SPY2","GeparNuevo","Pooled (trial-stratified)")){
    sub <- if (tr=="Pooled (trial-stratified)") D else D[D$trial==tr, ]
    v <- fit_or(sub, L, pooled = (tr=="Pooled (trial-stratified)"))
    rows[[length(rows)+1]] <- data.frame(layer=L, model=tr, n=nrow(sub),
      OR=v["OR"], CI_low=v["lo"], CI_high=v["hi"], interaction_P_LRT=lrt_p(L, tr), row.names=NULL)
  }
}
TAB <- do.call(rbind, rows)
num <- c("OR","CI_low","CI_high"); TAB[num] <- round(TAB[num], 2); TAB$interaction_P_LRT <- signif(TAB$interaction_P_LRT, 2)
out <- "../results/tables/ICI_SuppTable_interaction_OR.csv"
write.csv(TAB, out, row.names = FALSE)
cat("wrote", out, "\n\n")
cat("Integration-by-treatment interaction OR [95% CI], LRT P:\n")
print(TAB[TAB$layer=="Integration", c("model","n","OR","CI_low","CI_high","interaction_P_LRT")], row.names = FALSE)
cat("\n(pooled Integration expected ~ OR 1.81 [1.04-3.16], P=0.035)\n")
