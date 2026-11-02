# =====================================================================
# I-SPY2 durvalumab+olaparib arm (GSE173839, n=105) - independent-ish
# randomized anti-PD-L1-containing cohort. Test Integration x treatment
# (derived + legacy panels) and benefit-by-limiting-layer.
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
G <- "data/geparnuevo"

ex <- as.data.frame(fread(file.path(G,"G9173839_expr.txt.gz"), header=TRUE))
g  <- toupper(ex[[1]]); ex[[1]] <- NULL
E  <- rowsum(as.matrix(ex), g); colnames(E) <- as.character(colnames(E))
bm <- as.data.frame(fread(file.path(G,"G9173839_biomarkers.csv.gz")))
bm$ResearchID <- as.character(bm$ResearchID)
bm <- bm[bm$ResearchID %in% colnames(E) & bm$pCR.status %in% c(0,1), ]
E  <- E[, bm$ResearchID]
y  <- as.integer(bm$pCR.status)
arm <- factor(ifelse(bm$Arm=="control","control","treat"), c("control","treat"))
cat(sprintf("\nn=%d, treat=%d, control=%d; pCR treat %.0f%% / control %.0f%%\n",
            length(y), sum(arm=="treat"), sum(arm=="control"),
            100*mean(y[arm=="treat"]), 100*mean(y[arm=="control"])))

## panels
DER <- "../ICI-study-redesign/provenance-review-tight/derivation/ICI_PANEL_v1_locked.csv"
dp <- fread(DER, data.table=FALSE); dp$gene<-toupper(dp$gene)
PANEL_DER <- list(
  SENSORY = dp$gene[grepl("^Sensory",dp$layer)],
  INTPOS  = dp$gene[dp$layer=="Integration (chemokine trafficking)"],
  INTNEG  = dp$gene[dp$layer=="Integration (Treg / checkpoint regulation)" | grepl("^Resistance",dp$layer)],
  EXECPOS = dp$gene[grepl("^Execution",dp$layer)], EXECNEG=character(0))
source("ici_locked_reference.R"); PANEL_LEG <- load_panel("ICI_98gene_panel.csv")

score_cm <- function(expr,P){ bs<-function(gs){g<-intersect(gs,rownames(expr));if(!length(g))return(NULL)
   b<-t(apply(expr[g,,drop=FALSE],1,function(v) as.integer(v>=median(v,na.rm=TRUE))));colnames(b)<-colnames(expr);b}
  cs<-function(b) if(is.null(b)) rep(0,ncol(expr)) else colSums(b,na.rm=TRUE)
  S<-cs(bs(P$SENSORY));I<-cs(bs(P$INTPOS))-cs(bs(P$INTNEG));Ex<-cs(bs(P$EXECPOS))-cs(bs(P$EXECNEG))
  z<-function(x){s<-sd(x);if(is.na(s)||s==0)x*0 else (x-mean(x))/s}
  data.frame(sample=colnames(expr),Sensory=as.numeric(S),Integration=as.numeric(I),Execution=as.numeric(Ex),
             Composite=rowMeans(cbind(z(S),z(I),z(Ex))),row.names=NULL)}
z<-function(v)scale(v)[,1]; auc1<-function(yy,s)as.numeric(auc(yy,s,quiet=TRUE,direction="<"))
int_p<-function(sc) anova(glm(y~z(sc)+arm,binomial),glm(y~z(sc)*arm,binomial),test="LRT")[2,"Pr(>Chi)"]

for(tag in c("DERIVED","LEGACY")){
  P <- if(tag=="DERIVED") PANEL_DER else PANEL_LEG
  cm <- score_cm(E,P); cm<-cm[match(bm$ResearchID,cm$sample),]
  cat(sprintf("\n==== %s panel ====\n",tag))
  for(L in c("Sensory","Integration","Execution","Composite"))
    cat(sprintf("  %-12s interaction p=%.3f | AUC treat %.3f / control %.3f\n",L,int_p(cm[[L]]),
        auc1(y[arm=="treat"],cm[[L]][arm=="treat"]), auc1(y[arm=="control"],cm[[L]][arm=="control"])))
  hi<-cm$Integration>=median(cm$Integration)
  tb<-tapply(y,list(ifelse(hi,"hi","lo"),arm),function(z)100*mean(z))
  cat(sprintf("  Integration pCR%%: control lo/hi %.0f/%.0f | treat lo/hi %.0f/%.0f | benefit lo %+.0f / hi %+.0f\n",
      tb["lo","control"],tb["hi","control"],tb["lo","treat"],tb["hi","treat"],
      tb["lo","treat"]-tb["lo","control"],tb["hi","treat"]-tb["hi","control"]))
}
