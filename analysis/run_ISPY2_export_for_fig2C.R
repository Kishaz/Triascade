# =============================================================================
# Export real I-SPY2 performance data for Figure 2 Panel C (ICI in action).
# Reuses the locked-panel (tier-map A) scoring from run_derived_panel_ISPY2.R:
#   - ROC points + AUC for each layer in the ICB (pembro) arm
#   - Integration-high/low x arm pCR benefit table (the predictive headline)
# Saves fig2C_ISPY2_perf.rds (numbers are computed here, never typed into the figure).
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")

PANEL <- load_panel("ICI_98gene_panel.csv")

parse_sm <- function(f){ ln<-readLines(gzfile(f),warn=FALSE); ln<-ln[grepl("^!Sample_",ln)]; out<-list()
  for (r in strsplit(ln,"\t")){ k<-sub("^!Sample_","",r[1]); v<-gsub('"','',r[-1])
    if (grepl("^characteristics",k)){ i<-which(grepl(":",v))[1]; if(is.na(i)) next
      nm<-make.names(trimws(sub(":.*$","",v[i]))); v<-trimws(sub("^[^:]*: *","",v)) } else nm<-make.names(k)
    while (nm %in% names(out)) nm<-paste0(nm,".x"); out[[nm]]<-v }
  n<-max(lengths(out)); as.data.frame(lapply(out,function(x){length(x)<-n;x}),stringsAsFactors=FALSE) }
md <- rbind(parse_sm("data/ISPY2_GPL20078_series.txt.gz")[,c("patient.id","arm","her2","pcr")],
            parse_sm("data/ISPY2_GPL30493_series.txt.gz")[,c("patient.id","arm","her2","pcr")])
n01 <- function(v) ifelse(grepl("pos|[+]|yes|^1$",v,TRUE),1L,ifelse(grepl("neg|[-]|no|^0$",v,TRUE),0L,NA))
md$.id<-as.character(md$patient.id); md$.her2<-n01(md$her2); md$.pcr<-n01(md$pcr)
md$.arm2<-ifelse(grepl("Pembrolizumab",md$arm),"pembro",ifelse(md$arm=="Paclitaxel","control",NA))
ex<-as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz")); rownames(ex)<-toupper(ex[[1]]); ex[[1]]<-NULL
ex<-as.matrix(ex); colnames(ex)<-as.character(colnames(ex))
elig<-subset(md,.arm2 %in% c("pembro","control") & .her2==0 & !is.na(.pcr) & .id %in% colnames(ex))
elig<-elig[!duplicated(elig$.id),]; E<-ex[,elig$.id]
y<-elig$.pcr; arm<-factor(elig$.arm2,c("control","pembro"))

score_cm <- function(expr, PANEL){
  binset<-function(gs){g<-intersect(gs,rownames(expr)); if(!length(g)) return(NULL)
    b<-t(apply(expr[g,,drop=FALSE],1,function(v) as.integer(v>=median(v,na.rm=TRUE)))); colnames(b)<-colnames(expr); b}
  cs<-function(b) if(is.null(b)) rep(0,ncol(expr)) else colSums(b,na.rm=TRUE)
  S<-cs(binset(PANEL$SENSORY)); I<-cs(binset(PANEL$INTPOS))-cs(binset(PANEL$INTNEG))
  Ex<-cs(binset(PANEL$EXECPOS))-cs(binset(PANEL$EXECNEG))
  data.frame(sample=colnames(expr),Sensory=as.numeric(S),Integration=as.numeric(I),Execution=as.numeric(Ex)) }
cm <- score_cm(E, PANEL); cm <- cm[match(elig$.id, cm$sample), ]

## ---- ROC in the ICB (pembro) arm, per layer ----
pem <- arm=="pembro"; yp <- y[pem]
roc_pts <- list(); aucs <- c()
for (L in c("Sensory","Integration","Execution")) {
  r <- roc(yp, cm[[L]][pem], quiet=TRUE, direction="<")
  roc_pts[[L]] <- data.frame(layer=L, fpr=1-r$specificities, tpr=r$sensitivities)
  aucs[L] <- as.numeric(r$auc)
}
roc_df <- do.call(rbind, lapply(roc_pts, function(d) d[order(d$fpr, d$tpr), ]))

## ---- Integration-high/low x arm benefit (predictive headline) ----
hi <- cm$Integration >= median(cm$Integration)
tb <- tapply(y, list(ifelse(hi,"high","low"), arm), mean)*100
benefit <- data.frame(
  grp = c("Integration-low","Integration-high"),
  control = c(tb["low","control"], tb["high","control"]),
  pembro  = c(tb["low","pembro"],  tb["high","pembro"]))
benefit$benefit <- benefit$pembro - benefit$control
int_p <- anova(glm(y~scale(cm$Integration)+arm, binomial),
               glm(y~scale(cm$Integration)*arm, binomial), test="LRT")[2,"Pr(>Chi)"]

# ---- single-patient DEPLOYABILITY: freeze thresholds on an external breast
#      reference (NeoTRIP baseline), score each I-SPY2 patient ALONE (selfnorm),
#      and measure whether Integration's ranking is preserved in the ICB arm. ----
nt <- as.data.frame(fread("data/NeoTRIP_expr.txt.gz")); rownames(nt) <- toupper(nt$HGNC_approved_symbol)
ntb <- as.matrix(nt[, grep("_Baseline$", colnames(nt))])
refb <- fit_locked_reference(ntb, PANEL, "NeoTRIP_baseline")
sp <- score_locked(E, refb, mode="selfnorm")            # E columns are in elig order
auc_sp_int <- as.numeric(auc(y[pem], sp$Integration[pem], quiet=TRUE, direction="<"))

out <- list(roc=roc_df, auc=aucs, benefit=benefit, int_p=int_p,
            n_pembro=sum(pem), n_control=sum(!pem),
            pcr_pembro=100*mean(y[pem]), pcr_control=100*mean(y[!pem]),
            auc_sp_int=auc_sp_int, sp_ref="NeoTRIP breast reference")
saveRDS(out, "fig2C_ISPY2_perf.rds")
cat(sprintf("n pembro=%d control=%d | pCR pembro=%.0f%% control=%.0f%%\n",
            out$n_pembro, out$n_control, out$pcr_pembro, out$pcr_control))
cat(sprintf("TRIAL-LEVEL AUC pembro: S=%.2f I=%.2f E=%.2f | interaction p=%.3f\n",
            aucs["Sensory"],aucs["Integration"],aucs["Execution"], int_p))
cat(sprintf("SINGLE-PATIENT Integration AUC (selfnorm, NeoTRIP ref) = %.2f\n", auc_sp_int))
cat("saved fig2C_ISPY2_perf.rds\n")
