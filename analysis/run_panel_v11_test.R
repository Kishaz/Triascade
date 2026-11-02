# =====================================================================
# Test a gap-filled panel v1.1 (108 + 14 recurrence>=5 canonical genes)
# against the locked 108, on the predictive signal (I-SPY2 pembro + durva,
# pooled interaction). Adds are all literature-recurrent (>=5 signatures),
# so principled, not hand-curation.
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
L <- fread("ICI_panel_LOCKED_v1.csv", data.table=FALSE); L$gene<-toupper(L$gene)
P108 <- list(SENSORY=L$gene[L$layer=="Sense"],
  INTPOS=L$gene[L$layer=="Integration"&L$direction=="positive"],
  INTNEG=L$gene[L$layer=="Integration"&L$direction=="negative"],
  EXECPOS=L$gene[L$layer=="Execution"], EXECNEG=character(0))
# gap-fill additions (all recurrence>=5, from the 147-gene denovo panel)
add <- list(
  SENSORY = c("HLA-DRA"),
  INTPOS  = c("ICOS","ICOSLG","TNFSF9","TNFSF4","TNFRSF18"),          # costimulation (activating)
  INTNEG  = c("SPP1"),                                                # TAM/myeloid suppression
  EXECPOS = c("KLRK1","KLRD1","CD3G","CD3E","ZAP70","FAS","TCF7"))    # effector / T-NK identity
P122 <- Map(function(a,b) unique(c(a,b)), P108, c(add, EXECNEG=list(character(0))))
cat("v1.1 sizes:", paste(names(P122),lengths(P122),sep="="), " total",length(unique(unlist(P122))),"\n\n")

score_int <- function(E,P){ bs<-function(gs){g<-intersect(gs,rownames(E));if(!length(g))return(rep(0,ncol(E)))
   colSums(t(apply(E[g,,drop=FALSE],1,function(v) as.integer(v>=median(v,na.rm=TRUE)))))}
  bs(P$INTPOS)-bs(P$INTNEG)}
int_est <- function(Iz,arm,y){ m<-glm(y~Iz*arm,binomial);co<-summary(m)$coefficients;r<-grep(":",rownames(co))
  p<-anova(glm(y~Iz+arm,binomial),m,test="LRT")[2,"Pr(>Chi)"]; c(b=co[r,1],se=co[r,2],p=p)}
z<-function(v)scale(v)[,1]

## pembro
parse_sm<-function(f){ln<-readLines(gzfile(f),warn=FALSE);ln<-ln[grepl("^!Sample_",ln)];o<-list()
 for(r in strsplit(ln,"\t")){k<-sub("^!Sample_","",r[1]);v<-gsub('"','',r[-1])
  if(grepl("^characteristics",k)){i<-which(grepl(":",v))[1];if(is.na(i))next;nm<-make.names(trimws(sub(":.*$","",v[i])));v<-trimws(sub("^[^:]*: *","",v))}else nm<-make.names(k)
  while(nm%in%names(o))nm<-paste0(nm,".x");o[[nm]]<-v};n<-max(lengths(o));as.data.frame(lapply(o,function(x){length(x)<-n;x}),stringsAsFactors=FALSE)}
md<-rbind(parse_sm("data/ISPY2_GPL20078_series.txt.gz")[,c("patient.id","arm","her2","pcr")],parse_sm("data/ISPY2_GPL30493_series.txt.gz")[,c("patient.id","arm","her2","pcr")])
n01<-function(v)ifelse(grepl("pos|[+]|yes|^1$",v,TRUE),1L,ifelse(grepl("neg|[-]|no|^0$",v,TRUE),0L,NA))
md$.id<-as.character(md$patient.id);md$.her2<-n01(md$her2);md$.pcr<-n01(md$pcr)
md$.arm<-ifelse(grepl("Pembrolizumab",md$arm),"treat",ifelse(md$arm=="Paclitaxel","control",NA))
ex<-as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz"));rownames(ex)<-toupper(ex[[1]]);ex[[1]]<-NULL;ex<-as.matrix(ex);colnames(ex)<-as.character(colnames(ex))
el<-subset(md,.arm%in%c("treat","control")&.her2==0&!is.na(.pcr)&.id%in%colnames(ex));el<-el[!duplicated(el$.id),]
pembro<-list(E=ex[,el$.id],y=el$.pcr,arm=factor(el$.arm,c("control","treat")))
## durva
G<-"data/geparnuevo";exd<-as.data.frame(fread(file.path(G,"G9173839_expr.txt.gz"),header=TRUE));gg<-toupper(exd[[1]]);exd[[1]]<-NULL
Ed<-rowsum(as.matrix(exd),gg);colnames(Ed)<-as.character(colnames(Ed))
bm<-as.data.frame(fread(file.path(G,"G9173839_biomarkers.csv.gz")));bm$ResearchID<-as.character(bm$ResearchID)
bm<-bm[bm$ResearchID%in%colnames(Ed)&bm$pCR.status%in%c(0,1),]
durva<-list(E=Ed[,bm$ResearchID],y=as.integer(bm$pCR.status),arm=factor(ifelse(bm$Arm=="control","control","treat"),c("control","treat")))

pool <- function(P){ rows<-lapply(c("pembro","durva"),function(co){d<-get(co,envir=.GlobalEnv)
    e<-int_est(z(score_int(d$E,P)),d$arm,d$y); data.frame(cohort=co,logOR=e[1],se=e[2],p=e[3])})
  M<-do.call(rbind,rows); w<-1/M$se^2;b<-sum(w*M$logOR)/sum(w);se<-sqrt(1/sum(w));list(M=M,OR=exp(b),ci=exp(c(b-1.96*se,b+1.96*se)),p=2*pnorm(-abs(b/se)))}

for(tag in c("108 (locked v1)","122 (v1.1 gap-filled)")){
  P <- if(grepl("108",tag)) P108 else P122; r<-pool(P)
  cat(sprintf("== %s ==\n", tag))
  for(i in 1:nrow(r$M)) cat(sprintf("   %-7s OR=%.2f p=%.3f\n", r$M$cohort[i], exp(r$M$logOR[i]), r$M$p[i]))
  cat(sprintf("   POOLED OR=%.2f (95%% CI %.2f-%.2f) p=%.4f\n\n", r$OR, r$ci[1], r$ci[2], r$p))
}
