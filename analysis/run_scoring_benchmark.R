# Scoring-method benchmark: is cohort-median binarization the best way to score
# each layer? Compare on (1) the predictive signal (I-SPY2 Integration x treatment)
# and (2) construct validity (IMvigor210 layer scores separating immune phenotypes).
suppressPackageStartupMessages({ library(data.table); library(pROC) })
DER <- "../ICI-study-redesign/provenance-review-tight/derivation/ICI_PANEL_v1_locked.csv"
dp  <- fread(DER, data.table=FALSE); dp$gene<-toupper(dp$gene)
SENS <- dp$gene[grepl("^Sensory",dp$layer)]
INTP <- dp$gene[dp$layer=="Integration (chemokine trafficking)"]
INTN <- dp$gene[dp$layer=="Integration (Treg / checkpoint regulation)"|grepl("^Resistance",dp$layer)]
EXEC <- dp$gene[grepl("^Execution",dp$layer)]

## ---- per-gene scoring primitives (expr: genes x samples) ----
gmedian <- function(expr,g){g<-intersect(g,rownames(expr)); if(!length(g)) return(rep(0,ncol(expr)))
  m<-expr[g,,drop=FALSE]; m }
# metric for a gene set, by method
score_set <- function(expr, g, method){
  g<-intersect(g,rownames(expr)); if(!length(g)) return(rep(0,ncol(expr)))
  m<-expr[g,,drop=FALSE]
  if(method=="bin")   return(colSums(t(apply(m,1,function(v) as.integer(v>=median(v,na.rm=TRUE)))),na.rm=TRUE))
  if(method=="prop")  return(colMeans(t(apply(m,1,function(v) as.integer(v>=median(v,na.rm=TRUE)))),na.rm=TRUE))
  if(method=="zmean") return(colMeans(t(apply(m,1,function(v){s<-sd(v,na.rm=TRUE); if(is.na(s)||s==0) v*0 else (v-mean(v,na.rm=TRUE))/s})),na.rm=TRUE))
  if(method=="medcen")return(colMeans(t(apply(m,1,function(v) v-median(v,na.rm=TRUE))),na.rm=TRUE))
  if(method=="singscore"){ N<-nrow(expr); P<-apply(expr,2,function(x)(rank(x,ties.method="average")-0.5)/N); return(colMeans(P[g,,drop=FALSE])-0.5) }
}
integ <- function(expr, method){ score_set(expr,INTP,method) - score_set(expr,INTN,method) }
exece <- function(expr, method){ score_set(expr,EXEC,method) }
sens  <- function(expr, method){ score_set(expr,SENS,method) }
z<-function(v)scale(v)[,1]; auc1<-function(y,s)as.numeric(auc(y,s,quiet=TRUE,direction="<"))
methods<-c("bin","prop","zmean","medcen","singscore")

## ================= (1) I-SPY2 Integration x treatment =================
parse_sm<-function(f){ln<-readLines(gzfile(f),warn=FALSE);ln<-ln[grepl("^!Sample_",ln)];o<-list()
 for(r in strsplit(ln,"\t")){k<-sub("^!Sample_","",r[1]);v<-gsub('"','',r[-1])
  if(grepl("^characteristics",k)){i<-which(grepl(":",v))[1];if(is.na(i))next;nm<-make.names(trimws(sub(":.*$","",v[i])));v<-trimws(sub("^[^:]*: *","",v))}else nm<-make.names(k)
  while(nm%in%names(o))nm<-paste0(nm,".x");o[[nm]]<-v};n<-max(lengths(o));as.data.frame(lapply(o,function(x){length(x)<-n;x}),stringsAsFactors=FALSE)}
md<-rbind(parse_sm("data/ISPY2_GPL20078_series.txt.gz")[,c("patient.id","arm","her2","pcr")],parse_sm("data/ISPY2_GPL30493_series.txt.gz")[,c("patient.id","arm","her2","pcr")])
n01<-function(v)ifelse(grepl("pos|[+]|yes|^1$",v,TRUE),1L,ifelse(grepl("neg|[-]|no|^0$",v,TRUE),0L,NA))
md$.id<-as.character(md$patient.id);md$.her2<-n01(md$her2);md$.pcr<-n01(md$pcr)
md$.arm<-ifelse(grepl("Pembrolizumab",md$arm),"treat",ifelse(md$arm=="Paclitaxel","control",NA))
ex<-as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz"));rownames(ex)<-toupper(ex[[1]]);ex[[1]]<-NULL;ex<-as.matrix(ex);colnames(ex)<-as.character(colnames(ex))
el<-subset(md,.arm%in%c("treat","control")&.her2==0&!is.na(.pcr)&.id%in%colnames(ex));el<-el[!duplicated(el$.id),];E<-ex[,el$.id]
y<-el$.pcr; arm<-factor(el$.arm,c("control","treat"))
cat("=== (1) I-SPY2 Integration x treatment (n=248) by scoring method ===\n")
for(mth in methods){ Iz<-z(integ(E,mth))
  m1<-glm(y~Iz+arm,binomial);m2<-glm(y~Iz*arm,binomial);co<-summary(m2)$coefficients;r<-grep(":",rownames(co))
  p<-anova(m1,m2,test="LRT")[2,"Pr(>Chi)"]
  cat(sprintf("  %-10s interaction p=%.3f  OR=%.2f  | AUC pembro %.3f\n",mth,p,exp(co[r,1]),auc1(y[arm=="treat"],integ(E,mth)[arm=="treat"]))) }

## ================= (2) IMvigor210 construct validity =================
obj<-readRDS("data/imvigor210/imvigor210_clean.rds");Ei<-obj$E;pd<-obj$pd
ip<-factor(pd[["Immune phenotype"]],c("desert","excluded","inflamed")); resp<-ifelse(pd$binaryResponse=="CR/PR",1L,ifelse(pd$binaryResponse=="SD/PD",0L,NA))
infdes<-ip%in%c("inflamed","desert"); yid<-as.integer(ip[infdes]=="inflamed")
cat("\n=== (2) IMvigor210: Execution-score separation of inflamed vs desert, by method ===\n")
for(mth in methods){ exsc<-exece(Ei,mth)
  cat(sprintf("  %-10s Execution inflamed-vs-desert AUC=%.3f | Sense AUC=%.3f\n",mth, auc1(yid,exsc[infdes]), auc1(yid, sens(Ei,mth)[infdes]))) }
cat("\n(bin = current cohort-median binary count; prop = fraction high; zmean = mean z; medcen = mean(expr-median); singscore = rank single-sample)\n")
