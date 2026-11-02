suppressPackageStartupMessages({library(data.table); library(pROC)})
dp<-read.csv("ICI_panel_LOCKED_v1.csv"); dp$gene<-toupper(dp$gene)
EXH<-c("PDCD1","LAG3","TIGIT","HAVCR2","ENTPD1","PDCD1LG2")
EXEC_ALL<-toupper(dp$gene[dp$layer=="Execution"])
EFF<-setdiff(EXEC_ALL,EXH)
cat("Execution effector genes (n=",length(EFF),"); exhaustion genes (n=",length(EXH),")\n",sep="")
prophi<-function(expr,g){g<-intersect(g,rownames(expr)); if(!length(g))return(rep(NA,ncol(expr)))
  colMeans(t(apply(expr[g,,drop=FALSE],1,function(v) as.integer(v>=median(v,na.rm=TRUE)))),na.rm=TRUE)}
z<-function(v)scale(v)[,1]; A<-function(y,s)as.numeric(auc(y,s,quiet=TRUE,direction="<"))

## ---------- IMvigor210 ----------
obj<-readRDS("data/imvigor210/imvigor210_clean.rds"); E<-obj$E; pd<-obj$pd
eff<-prophi(E,EFF); exh<-prophi(E,EXH)
resp<-ifelse(pd$binaryResponse=="CR/PR",1L,ifelse(pd$binaryResponse=="SD/PD",0L,NA))
ip<-pd[["Immune phenotype"]]
cat("\n=== IMvigor210 (atezo, n responders usable) ===\n")
ok<-!is.na(resp)
cat(sprintf("effector score  -> response AUC = %.3f\n",A(resp[ok],eff[ok])))
cat(sprintf("exhaustion score-> response AUC = %.3f   (>0.5 means exhaustion predicts RESPONSE)\n",A(resp[ok],exh[ok])))
cat(sprintf("effector score  mean: inflamed %.2f / excluded %.2f / desert %.2f\n",
  mean(eff[ip=="inflamed"],na.rm=T),mean(eff[ip=="excluded"],na.rm=T),mean(eff[ip=="desert"],na.rm=T)))
cat(sprintf("exhaustion score mean: inflamed %.2f / excluded %.2f / desert %.2f\n",
  mean(exh[ip=="inflamed"],na.rm=T),mean(exh[ip=="excluded"],na.rm=T),mean(exh[ip=="desert"],na.rm=T)))
# Among effector-HIGH tumors, does exhaustion separate R vs NR? (the "armed but restrained" question)
hi<-eff>=median(eff,na.rm=T) & ok
cat(sprintf("\nAmong effector-HIGH tumors (n=%d): response rate exhaustion-hi %.0f%% vs exhaustion-lo %.0f%%\n",
  sum(hi), 100*mean(resp[hi & exh>=median(exh[hi])],na.rm=T), 100*mean(resp[hi & exh<median(exh[hi])],na.rm=T)))
lo<-eff<median(eff,na.rm=T) & ok
cat(sprintf("Among effector-LOW tumors  (n=%d): response rate exhaustion-hi %.0f%% vs exhaustion-lo %.0f%%\n",
  sum(lo), 100*mean(resp[lo & exh>=median(exh[lo],na.rm=T)],na.rm=T), 100*mean(resp[lo & exh<median(exh[lo],na.rm=T)],na.rm=T)))

## ---------- I-SPY2: does SUBTRACTING exhaustion help or hurt the predictive score? ----------
parse_sm<-function(f){ln<-readLines(gzfile(f),warn=FALSE);ln<-ln[grepl("^!Sample_",ln)];o<-list()
 for(r in strsplit(ln,"\t")){k<-sub("^!Sample_","",r[1]);v<-gsub('"','',r[-1])
  if(grepl("^characteristics",k)){i<-which(grepl(":",v))[1];if(is.na(i))next;nm<-make.names(trimws(sub(":.*$","",v[i])));v<-trimws(sub("^[^:]*: *","",v))}else nm<-make.names(k)
  while(nm%in%names(o))nm<-paste0(nm,".x");o[[nm]]<-v};n<-max(lengths(o));as.data.frame(lapply(o,function(x){length(x)<-n;x}),stringsAsFactors=FALSE)}
md<-rbind(parse_sm("data/ISPY2_GPL20078_series.txt.gz")[,c("patient.id","arm","her2","pcr")],parse_sm("data/ISPY2_GPL30493_series.txt.gz")[,c("patient.id","arm","her2","pcr")])
n01<-function(v)ifelse(grepl("pos|[+]|yes|^1$",v,TRUE),1L,ifelse(grepl("neg|[-]|no|^0$",v,TRUE),0L,NA))
md$.id<-as.character(md$patient.id);md$.her2<-n01(md$her2);md$.pcr<-n01(md$pcr)
md$.arm<-ifelse(grepl("Pembrolizumab",md$arm),"treat",ifelse(md$arm=="Paclitaxel","control",NA))
ex<-as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz"));rownames(ex)<-toupper(ex[[1]]);ex[[1]]<-NULL;ex<-as.matrix(ex)
el<-subset(md,.arm%in%c("treat","control")&.her2==0&!is.na(.pcr)&.id%in%colnames(ex));el<-el[!duplicated(el$.id),];Ei<-ex[,el$.id]
y<-el$.pcr; arm<-factor(el$.arm,c("control","treat"))
effi<-prophi(Ei,EFF); exhi<-prophi(Ei,EXH)
test_int<-function(score,lab){Iz<-z(score);p<-anova(glm(y~Iz+arm,binomial),glm(y~Iz*arm,binomial),test="LRT")[2,"Pr(>Chi)"]
  cat(sprintf("  Execution = %-28s interaction p = %.3f\n",lab,p))}
cat("\n=== I-SPY2: Execution x treatment interaction under different exhaustion handling ===\n")
test_int(effi,"effector only (drop exh)")
test_int(effi+exhi,"effector + exhaustion (current, all +)")
test_int(effi-exhi,"effector - exhaustion (subtracted)")
