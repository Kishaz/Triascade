# =====================================================================
# (A) BrighTNess (n=482, chemo-only TNBC): predictive-vs-prognostic
#     control. Does Integration show a treatment INTERACTION with
#     carboplatin? (Expected: NO - the ICB benefit split is ICB-specific.)
# (B) Meta-estimate of the Integration x treatment interaction across the
#     randomized ICB arms (I-SPY2 pembro, I-SPY2 durvalumab+olaparib),
#     inverse-variance pooled, derived + legacy panels.
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")

## ---- panels ----
DER <- "../ICI-study-redesign/provenance-review-tight/derivation/ICI_PANEL_v1_locked.csv"
dp  <- fread(DER, data.table=FALSE); dp$gene<-toupper(dp$gene)
PANEL_DER <- list(SENSORY=dp$gene[grepl("^Sensory",dp$layer)],
  INTPOS=dp$gene[dp$layer=="Integration (chemokine trafficking)"],
  INTNEG=dp$gene[dp$layer=="Integration (Treg / checkpoint regulation)"|grepl("^Resistance",dp$layer)],
  EXECPOS=dp$gene[grepl("^Execution",dp$layer)], EXECNEG=character(0))
PANEL_LEG <- load_panel("ICI_98gene_panel.csv")

score_cm <- function(expr,P){ bs<-function(gs){g<-intersect(gs,rownames(expr));if(!length(g))return(NULL)
   b<-t(apply(expr[g,,drop=FALSE],1,function(v) as.integer(v>=median(v,na.rm=TRUE))));colnames(b)<-colnames(expr);b}
  cs<-function(b) if(is.null(b)) rep(0,ncol(expr)) else colSums(b,na.rm=TRUE)
  Integration<-cs(bs(P$INTPOS))-cs(bs(P$INTNEG)); data.frame(sample=colnames(expr),Integration=as.numeric(Integration),row.names=NULL)}
z<-function(v)scale(v)[,1]
# interaction logOR + SE for Integration x arm  (arm: factor control,treat)
int_est <- function(Iz, arm, y){
  m<-glm(y~Iz*arm, binomial); s<-summary(m)$coefficients
  r<-grep(":", rownames(s)); c(beta=s[r,1], se=s[r,2],
    p=anova(glm(y~Iz+arm,binomial),m,test="LRT")[2,"Pr(>Chi)"]) }

## =========== I-SPY2 pembro (n=248) ===========
parse_sm <- function(f){ln<-readLines(gzfile(f),warn=FALSE);ln<-ln[grepl("^!Sample_",ln)];out<-list()
  for(r in strsplit(ln,"\t")){k<-sub("^!Sample_","",r[1]);v<-gsub('"','',r[-1])
   if(grepl("^characteristics",k)){i<-which(grepl(":",v))[1];if(is.na(i))next
     nm<-make.names(trimws(sub(":.*$","",v[i])));v<-trimws(sub("^[^:]*: *","",v))}else nm<-make.names(k)
   while(nm%in%names(out))nm<-paste0(nm,".x");out[[nm]]<-v};n<-max(lengths(out))
  as.data.frame(lapply(out,function(x){length(x)<-n;x}),stringsAsFactors=FALSE)}
md<-rbind(parse_sm("data/ISPY2_GPL20078_series.txt.gz")[,c("patient.id","arm","her2","pcr")],
          parse_sm("data/ISPY2_GPL30493_series.txt.gz")[,c("patient.id","arm","her2","pcr")])
n01<-function(v)ifelse(grepl("pos|[+]|yes|^1$",v,TRUE),1L,ifelse(grepl("neg|[-]|no|^0$",v,TRUE),0L,NA))
md$.id<-as.character(md$patient.id);md$.her2<-n01(md$her2);md$.pcr<-n01(md$pcr)
md$.arm<-ifelse(grepl("Pembrolizumab",md$arm),"treat",ifelse(md$arm=="Paclitaxel","control",NA))
ex<-as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz"));rownames(ex)<-toupper(ex[[1]]);ex[[1]]<-NULL
ex<-as.matrix(ex);colnames(ex)<-as.character(colnames(ex))
el<-subset(md,.arm%in%c("treat","control")&.her2==0&!is.na(.pcr)&.id%in%colnames(ex));el<-el[!duplicated(el$.id),]
pembro<-list(E=ex[,el$.id], y=el$.pcr, arm=factor(el$.arm,c("control","treat")))

## =========== I-SPY2 durvalumab+olaparib (n=100) ===========
G<-"data/geparnuevo"; exd<-as.data.frame(fread(file.path(G,"G9173839_expr.txt.gz"),header=TRUE))
gg<-toupper(exd[[1]]);exd[[1]]<-NULL; Ed<-rowsum(as.matrix(exd),gg);colnames(Ed)<-as.character(colnames(Ed))
bm<-as.data.frame(fread(file.path(G,"G9173839_biomarkers.csv.gz")));bm$ResearchID<-as.character(bm$ResearchID)
bm<-bm[bm$ResearchID%in%colnames(Ed)&bm$pCR.status%in%c(0,1),]
durva<-list(E=Ed[,bm$ResearchID], y=as.integer(bm$pCR.status), arm=factor(ifelse(bm$Arm=="control","control","treat"),c("control","treat")))

## ---- build the meta table ----
cat("=== (B) Integration x treatment interaction, per randomized ICB arm ===\n")
meta_rows <- list()
for(tag in c("DERIVED","LEGACY")){ P<-if(tag=="DERIVED")PANEL_DER else PANEL_LEG
  for(co in c("pembro","durva")){ d<-get(co)
    sc<-score_cm(d$E,P); sc<-sc[match(colnames(d$E),sc$sample),]
    e<-int_est(z(sc$Integration), d$arm, d$y)
    meta_rows[[paste(tag,co)]]<-data.frame(panel=tag,cohort=co,n=length(d$y),
      logOR=e["beta"],se=e["se"],interaction_p=e["p"])
    cat(sprintf("  %-8s %-7s n=%3d  interaction logOR=%+.2f (SE %.2f)  OR=%.2f  p=%.3f\n",
        tag,co,length(d$y),e["beta"],e["se"],exp(e["beta"]),e["p"])) }}
MT<-do.call(rbind,meta_rows); rownames(MT)<-NULL

## ---- inverse-variance fixed-effect pool (per panel) ----
cat("\n=== pooled interaction (inverse-variance fixed effect) ===\n")
for(tag in c("DERIVED","LEGACY")){ s<-MT[MT$panel==tag,]
  w<-1/s$se^2; b<-sum(w*s$logOR)/sum(w); se<-sqrt(1/sum(w)); p<-2*pnorm(-abs(b/se))
  Q<-sum(w*(s$logOR-b)^2); I2<-max(0,(Q-(nrow(s)-1))/Q)*100
  cat(sprintf("  %-8s pooled OR=%.2f (95%% CI %.2f-%.2f)  p=%.4f  | I2=%.0f%% (k=%d, both I-SPY2)\n",
      tag,exp(b),exp(b-1.96*se),exp(b+1.96*se),p,I2,nrow(s))) }

## =========== (A) BrighTNess control ===========
cat("\n=== (A) BrighTNess (chemo-only, n=482): Integration x CARBOPLATIN interaction ===\n")
e<-new.env(); load("data/geparnuevo/valid_datasets_HTG.RData",envir=e)
Eb<-e$brtn.htg.Rseq; rownames(Eb)<-toupper(rownames(Eb))
cli<-e$brtn.htg.cli.dat
idcol<-names(which.max(sapply(cli,function(c) mean(as.character(c)%in%colnames(Eb)))))
cli<-cli[match(colnames(Eb),as.character(cli[[idcol]])),]
yb<-ifelse(grepl("pCR",cli$pathologic_complete_response,ignore.case=TRUE)|cli$pathologic_complete_response=="pCR",1L,0L)
carbo<-factor(ifelse(grepl("PBO FOR CARBOPLATIN",cli$planned_arm),"control","treat"),c("control","treat"))
ok<-!is.na(yb)&!is.na(carbo)
cat(sprintf("  matched n=%d (carbo+ %d, carbo- %d); pCR carbo+ %.0f%% / carbo- %.0f%%\n",
    sum(ok),sum(carbo[ok]=="treat"),sum(carbo[ok]=="control"),
    100*mean(yb[ok&carbo=="treat"]),100*mean(yb[ok&carbo=="control"])))
for(tag in c("DERIVED","LEGACY")){ P<-if(tag=="DERIVED")PANEL_DER else PANEL_LEG
  sc<-score_cm(Eb[,ok],P); Iz<-z(sc$Integration); a<-carbo[ok]; yy<-yb[ok]
  e2<-int_est(Iz,a,yy)
  hi<-sc$Integration>=median(sc$Integration); tb<-tapply(yy,list(ifelse(hi,"hi","lo"),a),function(z)100*mean(z))
  cat(sprintf("  %-8s Integration x carbo interaction p=%.3f (OR=%.2f) | benefit lo %+.0f / hi %+.0f pts\n",
      tag,e2["p"],exp(e2["beta"]),tb["lo","treat"]-tb["lo","control"],tb["hi","treat"]-tb["hi","control"]))
  cat(sprintf("           Integration->pCR main assoc AUC (all chemo): %.3f\n",
      as.numeric(auc(yy,sc$Integration,quiet=TRUE,direction="<")))) }
write.csv(MT,"meta_interaction_table.csv",row.names=FALSE)
