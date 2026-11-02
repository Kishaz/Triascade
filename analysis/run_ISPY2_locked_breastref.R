# Single-patient locked scoring of I-SPY2 with a TISSUE-MATCHED breast reference
# (NeoTRIP, baseline samples). Fair test of the single-patient claim.
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")
PANEL <- load_panel("ICI_98gene_panel.csv")

# ---- NeoTRIP breast reference (baseline samples) ----
nt <- as.data.frame(fread("data/NeoTRIP_expr.txt.gz"))
rownames(nt) <- toupper(nt$HGNC_approved_symbol)
samp <- grep("_Baseline$", colnames(nt), value=TRUE)
ntm <- as.matrix(nt[, samp]);
ref <- fit_locked_reference(ntm, PANEL, label="NeoTRIP_baseline")
cat(sprintf("breast reference frozen on %d panel genes (NeoTRIP, %d baseline samples)\n", length(ref$genes), length(samp)))

# ---- I-SPY2 ----
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

auc1 <- function(yy,s) as.numeric(auc(yy,s,quiet=TRUE,direction="<"))
int_p <- function(sc) anova(glm(y~scale(sc)+arm,binomial), glm(y~scale(sc)*arm,binomial), test="LRT")[2,"Pr(>Chi)"]
ben <- function(sc){ hi<-sc>=median(sc); tb<-tapply(y,list(ifelse(hi,"hi","lo"),arm),mean)*100
  c(low=tb["lo","pembro"]-tb["lo","control"], high=tb["hi","pembro"]-tb["hi","control"]) }
report <- function(tag, cm){ cm<-cm[match(elig$.id,cm$sample),]; cat(sprintf("\n== %s ==\n", tag))
  for (L in c("Sensory","Integration","Execution")) cat(sprintf("  %-12s AUC pembro %.3f | interaction p %.4f\n",
      L, auc1(y[arm=="pembro"], cm[[L]][arm=="pembro"]), int_p(cm[[L]])))
  b<-ben(cm$Integration); cat(sprintf("  Integration benefit: low %+.0f / high %+.0f pts\n", b["low"], b["high"])) }

report("LOCKED single-patient (selfnorm, ref=NeoTRIP breast)", score_locked(E, ref, mode="selfnorm"))
report("LOCKED single-patient (quantile, ref=NeoTRIP breast)", score_locked(E, ref, mode="quantile"))
cat("\ndone\n")
