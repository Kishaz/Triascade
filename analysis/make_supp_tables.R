# =============================================================================
# Supplementary tables for the pan-cancer ICB validation.
# Reads pancancer_validation.rds (run_pancancer_extra.R). Writes, per cohort,
# the single-patient AND cohort-median AUC for ALL THREE layers (Sensory /
# Integration / Execution), plus the external-reference single-patient AUC and
# the Integration DeLong 95% CI. This is the source table behind Fig 3b and the
# Sensory/Execution supplementary separation figures.
# =============================================================================
R <- readRDS("pancancer_validation.rds")

WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
CTLA4 <- c("VanAllen2015"); NULLB <- c("GSE165252","PRINCE","Braun")
bkt <- function(k) if(k %in% WORKS) "validated" else if(k %in% CTLA4) "ctla4" else if(k %in% NULLB) "null" else "excluded"

TAB <- do.call(rbind, lapply(names(R), function(k){ x <- R[[k]]; A <- x$auc
  data.frame(cohort=x$cohort, cancer=x$cancer, drug=x$drug, bucket=bkt(k),
             n=x$n, nR=x$nR, nNR=x$nNR, panel_cov=x$cov,
             # single-patient selfnorm (the paper's thesis metric)
             Sensory_single=A["single_self","Sensory"], Integration_single=A["single_self","Integration"], Execution_single=A["single_self","Execution"],
             # cohort-median
             Sensory_cohort=A["cohort","Sensory"], Integration_cohort=A["cohort","Integration"], Execution_cohort=A["cohort","Execution"],
             # single-patient, external NeoTRIP reference
             Sensory_extref=A["single_ext","Sensory"], Integration_extref=A["single_ext","Integration"], Execution_extref=A["single_ext","Execution"],
             # Integration DeLong 95% CI (single-patient)
             Integration_single_lo=x$int_ci_single[2], Integration_single_hi=x$int_ci_single[3],
             row.names=NULL) }))
TAB[ ,sapply(TAB, is.numeric)] <- round(TAB[ ,sapply(TAB, is.numeric)], 3)
TAB <- TAB[order(factor(TAB$bucket, c("validated","ctla4","null","excluded")), -TAB$Integration_single), ]

dir.create("../results/tables", showWarnings=FALSE, recursive=TRUE)
out <- "../results/tables/ICI_SuppTable_layer_auc_by_cohort.csv"
write.csv(TAB, out, row.names=FALSE)
cat("wrote", out, "-", nrow(TAB), "cohorts x 3 layers\n")
print(TAB[, c("cohort","cancer","bucket","n","Sensory_single","Integration_single","Execution_single")], row.names=FALSE)
