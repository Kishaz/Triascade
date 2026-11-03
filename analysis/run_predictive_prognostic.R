# =============================================================================
# PREDICTIVE vs PROGNOSTIC  (Rung 2 of the biomarker ladder)
# Is the ICI signal SPECIFIC to immune-checkpoint blockade, or merely a good
# prognostic marker that tracks outcome regardless of therapy?
#
# Two complementary strands, both with PRE-SPECIFIED direction (higher layer
# score = responder; never flipped) and the project's LOCKED 98-gene panel:
#
#  (A) Randomized biomarker x treatment-arm INTERACTION tests  [gold standard]
#       - I-SPY2    : pembrolizumab+paclitaxel   vs  paclitaxel alone   (HER2-)
#       - GeparNuevo: durvalumab + chemo         vs  placebo + chemo    (TNBC)
#      A significant score x arm interaction (nested LRT) = predictive.
#      AUC is reported separately in the ICB arm and the control arm.
#      Scoring: cohort-median AND locked single-patient (leave-one-out),
#      exactly as elsewhere in the project (ici_locked_reference.R).
#
#  (B) Cross-context SPECIFICITY gradient across all 16 cohorts: Integration
#      single-patient AUC grouped by treatment context. A checkpoint-specific
#      marker should work for anti-PD-(L)1, partially for anti-CTLA4, and FAIL
#      for non-checkpoint immunotherapy (bispecific) and in the ccRCC negative
#      control (where angiogenesis, not inflammation, drives response).
#
# Writes: results/tables/ICI_predictive_interaction.csv
#         results/tables/ICI_predictive_crosscontext.csv
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")
PANEL <- load_panel("ICI_98gene_panel.csv")
OUT   <- "../results/tables"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)

z     <- function(v) scale(v)[,1]
auc1  <- function(yy, s) as.numeric(auc(yy, s, quiet=TRUE, direction="<"))
int_p <- function(sc, y, arm)
  anova(glm(y ~ z(sc) + arm, binomial),
        glm(y ~ z(sc) * arm, binomial), test="LRT")[2, "Pr(>Chi)"]

# leave-one-out locked single-patient scoring of an expression matrix ---------
score_loo <- function(E){
  lk <- do.call(rbind, lapply(seq_len(ncol(E)), function(i){
    r <- fit_locked_reference(E[,-i,drop=FALSE], PANEL)
    score_locked(E[,i,drop=FALSE], r, "quantile") }))
  lk[match(colnames(E), lk$sample), ]
}

# one interaction block, all four layers, both scorings ----------------------
interaction_block <- function(trial, E, y, arm){
  cm <- score_cohort_median(E, PANEL); cm <- cm[match(colnames(E), cm$sample), ]
  lk <- score_loo(E)
  mk <- function(df, scoring){
    do.call(rbind, Map(function(sc, tag) data.frame(
        trial=trial, scoring=scoring, layer=tag,
        n_icb=sum(arm=="icb"), n_ctrl=sum(arm=="ctrl"),
        interaction_p = int_p(sc, y, arm),
        AUC_icb  = auc1(y[arm=="icb"],  sc[arm=="icb"]),
        AUC_ctrl = auc1(y[arm=="ctrl"], sc[arm=="ctrl"]), row.names=NULL),
      list(df$Sensory, df$Integration, df$Execution, df$Composite),
      c("Sensory","Integration","Execution","Composite"))) }
  rbind(mk(cm, "cohort-median"), mk(lk, "locked single-patient"))
}

# ---- (A1) I-SPY2 : pembro+paclitaxel vs paclitaxel alone (HER2-negative) ----
parse_sm <- function(f){
  ln <- readLines(gzfile(f), warn=FALSE); ln <- ln[grepl("^!Sample_", ln)]; out <- list()
  for (r in strsplit(ln, "\t")){ k <- sub("^!Sample_","",r[1]); v <- gsub('"','',r[-1])
    if (grepl("^characteristics", k)){ i <- which(grepl(":", v))[1]; if (is.na(i)) next
      nm <- make.names(trimws(sub(":.*$","",v[i]))); v <- trimws(sub("^[^:]*: *","",v)) } else nm <- make.names(k)
    while (nm %in% names(out)) nm <- paste0(nm,".x"); out[[nm]] <- v }
  n <- max(lengths(out)); as.data.frame(lapply(out, function(x){ length(x)<-n; x }), stringsAsFactors=FALSE) }
md <- rbind(parse_sm("data/ISPY2_GPL20078_series.txt.gz")[,c("patient.id","arm","her2","hr","pcr")],
            parse_sm("data/ISPY2_GPL30493_series.txt.gz")[,c("patient.id","arm","her2","hr","pcr")])
n01 <- function(v) ifelse(grepl("pos|[+]|yes|^1$",v,TRUE),1L, ifelse(grepl("neg|[-]|no|^0$",v,TRUE),0L,NA))
md$.id  <- as.character(md$patient.id); md$.her2 <- n01(md$her2); md$.pcr <- n01(md$pcr)
md$.arm <- ifelse(grepl("Pembrolizumab", md$arm), "icb", ifelse(md$arm=="Paclitaxel","ctrl",NA))
ex <- as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz")); rownames(ex) <- toupper(ex[[1]]); ex[[1]] <- NULL
ex <- as.matrix(ex); colnames(ex) <- as.character(colnames(ex))
el <- subset(md, .arm %in% c("icb","ctrl") & .her2==0 & !is.na(.pcr) & .id %in% colnames(ex))
el <- el[!duplicated(el$.id), ]
A_ispy <- interaction_block("I-SPY2 (pembro vs pac)", ex[, el$.id],
                            el$.pcr, factor(el$.arm, c("ctrl","icb")))

# ---- (A2) GeparNuevo : durvalumab+chemo vs placebo+chemo (TNBC) -------------
G  <- "data/geparnuevo"
gx <- as.data.frame(fread(file.path(G,"G9173839_expr.txt.gz"), header=TRUE))
gg <- toupper(gx[[1]]); gx[[1]] <- NULL
GE <- rowsum(as.matrix(gx), gg); colnames(GE) <- as.character(colnames(GE))
bm <- as.data.frame(fread(file.path(G,"G9173839_biomarkers.csv.gz"))); bm$ResearchID <- as.character(bm$ResearchID)
bm <- bm[bm$ResearchID %in% colnames(GE) & bm$pCR.status %in% c(0,1), ]
GE <- GE[, bm$ResearchID]
gy <- as.integer(bm$pCR.status)
garm <- factor(ifelse(bm$Arm=="control","ctrl","icb"), c("ctrl","icb"))
A_gep <- interaction_block("GeparNuevo (durva vs placebo)", GE, gy, garm)

A <- rbind(A_ispy, A_gep)
A[,c("interaction_p","AUC_icb","AUC_ctrl")] <- round(A[,c("interaction_p","AUC_icb","AUC_ctrl")], 3)
write.csv(A, file.path(OUT,"ICI_predictive_interaction.csv"), row.names=FALSE)

# per-arm pCR rate by Integration median-split (intuitive view of the interaction)
wilson <- function(k,n){ if(n==0) return(c(NA,NA)); z<-1.96; p<-k/n; d<-1+z^2/n
  c(100*((p+z^2/(2*n))/d - z*sqrt(p*(1-p)/n+z^2/(4*n^2))/d),
    100*((p+z^2/(2*n))/d + z*sqrt(p*(1-p)/n+z^2/(4*n^2))/d)) }
arm_rates <- function(trial, E, y, arm){
  sc <- score_loo(E)$Integration
  split <- ifelse(sc > median(sc), "High", "Low")          # within-trial median split
  do.call(rbind, lapply(c("ctrl","icb"), function(a) do.call(rbind, lapply(c("Low","High"), function(g){
    sel <- arm==a & split==g; k <- sum(y[sel]); n <- sum(sel); ci <- wilson(k,n)
    data.frame(trial=trial, arm=ifelse(a=="icb","ICB arm","Control arm"),
               integ=g, k=k, n=n, pct=100*k/n, lo=ci[1], hi=ci[2], row.names=NULL) }))))
}
RT <- rbind(arm_rates("I-SPY2 (pembro vs pac)", ex[, el$.id], el$.pcr, factor(el$.arm, c("ctrl","icb"))),
            arm_rates("GeparNuevo (durva vs placebo)", GE, gy, garm))
RT[,c("pct","lo","hi")] <- round(RT[,c("pct","lo","hi")], 1)
write.csv(RT, file.path(OUT,"ICI_predictive_interaction_rates.csv"), row.names=FALSE)

cat("=============== (A) RANDOMIZED INTERACTION TESTS ===============\n")
for (tr in unique(A$trial)){ cat("\n##", tr, "\n")
  print(A[A$trial==tr & A$scoring=="locked single-patient",
          c("layer","interaction_p","AUC_icb","AUC_ctrl","n_icb","n_ctrl")], row.names=FALSE) }

# ---- (A3) POOLED continuous layer x arm across the TWO checkpoint trials -----
# Max-power version: within-trial z of each cohort-median layer score, then one
# trial-stratified logistic with a score x arm interaction. Confirms the discrete
# per-trial result and pools I-SPY2 + GeparNuevo (n=348).
patient_frame <- function(trial, E, y, arm){
  cm <- score_cohort_median(E, PANEL); cm <- cm[match(colnames(E), cm$sample), ]
  zt <- function(v) scale(v)[,1]
  data.frame(trial=trial, y=y, arm=factor(arm, c("ctrl","icb")),
             Sensory=zt(cm$Sensory), Integration=zt(cm$Integration), Execution=zt(cm$Execution)) }
CKPT <- rbind(patient_frame("I-SPY2", ex[, el$.id], el$.pcr, el$.arm),
              patient_frame("GeparNuevo", GE, gy, as.character(garm)))
write.csv(CKPT, file.path(OUT,"ICI_predictive_checkpoint_patients.csv"), row.names=FALSE)
cont_int <- function(L){
  b <- glm(reformulate(c(L,"arm","trial"),"y"), binomial, CKPT)
  f <- glm(reformulate(c(sprintf("%s*arm",L),"trial"),"y"), binomial, CKPT)
  data.frame(layer=L, n=nrow(CKPT),
             interaction_p = anova(b,f,test="LRT")[2,"Pr(>Chi)"],
             OR_per_SD     = exp(coef(f)[sprintf("%s:armicb",L)]), row.names=NULL) }
CI_CONT <- do.call(rbind, lapply(c("Sensory","Integration","Execution"), cont_int))
CI_CONT[,c("interaction_p","OR_per_SD")] <- round(CI_CONT[,c("interaction_p","OR_per_SD")], 3)
write.csv(CI_CONT, file.path(OUT,"ICI_predictive_interaction_continuous.csv"), row.names=FALSE)
cat("\n## POOLED continuous (I-SPY2 + GeparNuevo, n=348), trial-stratified score x arm:\n")
print(CI_CONT, row.names=FALSE)

# ---- (B) cross-context specificity gradient (Integration single-patient) ----
R <- readRDS("pancancer_validation.rds")
cc <- do.call(rbind, lapply(names(R), function(k){ x <- R[[k]]
  ci <- if (!is.null(x$int_ci_single)) x$int_ci_single else x$int_ci
  data.frame(cohort=k, cancer=x$cancer, drug=x$drug, n=x$n,
             AUC=as.numeric(ci[1]), lo=as.numeric(ci[2]), hi=as.numeric(ci[3])) }))
cc$context <- "anti-PD-(L)1"
cc$context[grepl("CTLA4|ipi", cc$drug, TRUE) & !grepl("PD", cc$drug, TRUE)] <- "anti-CTLA4"
cc$context[grepl("CD40", cc$drug, TRUE)]            <- "chemo+ICB combo"
cc$context[grepl("neg-ctrl",  cc$drug)]             <- "negative control (ccRCC)"
cc$context[grepl("spec-ctrl", cc$drug)]             <- "specificity control (bispecific)"
cc[,c("AUC","lo","hi")] <- round(cc[,c("AUC","lo","hi")], 3)
cc <- cc[order(factor(cc$context, levels=c("anti-PD-(L)1","anti-CTLA4","chemo+ICB combo",
       "negative control (ccRCC)","specificity control (bispecific)")), -cc$AUC), ]
write.csv(cc, file.path(OUT,"ICI_predictive_crosscontext.csv"), row.names=FALSE)

cat("\n=============== (B) CROSS-CONTEXT SPECIFICITY ===============\n")
print(cc[,c("cohort","cancer","context","n","AUC","lo","hi")], row.names=FALSE)
cat("\nmean Integration AUC by context:\n")
print(round(tapply(cc$AUC, cc$context, mean), 3))
cat("\nwrote ICI_predictive_interaction.csv + ICI_predictive_crosscontext.csv\n")
