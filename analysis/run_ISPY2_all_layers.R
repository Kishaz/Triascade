suppressPackageStartupMessages({library(data.table); library(pROC)}); source("ici_locked_reference.R")
PANEL<-load_panel("ICI_98gene_panel.csv")
parse_sm<-function(f){ln<-readLines(gzfile(f),warn=FALSE);ln<-ln[grepl("^!Sample_",ln)];out<-list()
 for(r in strsplit(ln,"\t")){k<-sub("^!Sample_","",r[1]);v<-gsub('"','',r[-1])
  if(grepl("^characteristics",k)){i<-which(grepl(":",v))[1];if(is.na(i))next;nm<-make.names(trimws(sub(":.*$","",v[i])));v<-trimws(sub("^[^:]*: *","",v))}else nm<-make.names(k)
  while(nm%in%names(out))nm<-paste0(nm,".x");out[[nm]]<-v};n<-max(lengths(out));as.data.frame(lapply(out,function(x){length(x)<-n;x}),stringsAsFactors=FALSE)}
md<-rbind(parse_sm("data/ISPY2_GPL20078_series.txt.gz")[,c("patient.id","arm","her2","hr","pcr")],parse_sm("data/ISPY2_GPL30493_series.txt.gz")[,c("patient.id","arm","her2","hr","pcr")])
n01<-function(v)ifelse(grepl("pos|[+]|yes|^1$",v,TRUE),1L,ifelse(grepl("neg|[-]|no|^0$",v,TRUE),0L,NA));md$.id<-as.character(md$patient.id);md$.her2<-n01(md$her2);md$.hr<-n01(md$hr);md$.pcr<-n01(md$pcr);md$.arm2<-ifelse(grepl("Pembrolizumab",md$arm),"pembro",ifelse(md$arm=="Paclitaxel","control",NA))
ex<-as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz"));rownames(ex)<-toupper(ex[[1]]);ex[[1]]<-NULL;ex<-as.matrix(ex);colnames(ex)<-as.character(colnames(ex))
elig<-subset(md,.arm2 %in% c("pembro","control") & .her2==0 & !is.na(.pcr) & .id %in% colnames(ex));elig<-elig[!duplicated(elig$.id),];E<-ex[,elig$.id]
y<-elig$.pcr;arm<-factor(elig$.arm2,c("control","pembro"))
cm<-score_cohort_median(E,PANEL);cm<-cm[match(elig$.id,cm$sample),]
# LOO locked for all layers
lk<-do.call(rbind,lapply(seq_len(ncol(E)),function(i){r<-fit_locked_reference(E[,-i,drop=FALSE],PANEL);score_locked(E[,i,drop=FALSE],r,"quantile")}))
lk<-lk[match(elig$.id,lk$sample),]
z<-function(v)scale(v)[,1]
row_for<-function(score,tag){
 iz<-z(score)
 p<-anova(glm(y~iz+arm,binomial),glm(y~iz*arm,binomial),test="LRT")[2,"Pr(>Chi)"]
 a_pem<-as.numeric(auc(y[arm=="pembro"],score[arm=="pembro"],quiet=TRUE,direction="<"))
 a_con<-as.numeric(auc(y[arm=="control"],score[arm=="control"],quiet=TRUE,direction="<"))
 data.frame(layer=tag,interaction_p=round(p,3),AUC_pembro=round(a_pem,3),AUC_control=round(a_con,3))
}
cat("=== COHORT-MEDIAN scoring ===\n")
cmres<-do.call(rbind,Map(row_for,list(cm$Sensory,cm$Integration,cm$Execution,cm$Composite),c("Sensory","Integration","Execution","Composite")))
print(cmres,row.names=FALSE)
cat("\n=== LOCKED single-patient (leave-one-out) scoring ===\n")
lkres<-do.call(rbind,Map(row_for,list(lk$Sensory,lk$Integration,lk$Execution,lk$Composite),c("Sensory","Integration","Execution","Composite")))
print(lkres,row.names=FALSE)
