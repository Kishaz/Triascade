# Archetype distribution on IMvigor210 using the LEGACY 98-gene panel.
# Legacy E has a suppressive arm: E = mean(ExecPos high) - mean(ExecNeg high) in [-1,+1].
# Layer HIGH at principled midpoints: S >= 0.5 ; I >= 0 ; E >= 0.
suppressPackageStartupMessages({ library(data.table) })
source("ici_locked_reference.R")
PAN <- load_panel("ICI_98gene_panel.csv")   # SENSORY21 INTPOS21 INTNEG20 EXECPOS19 EXECNEG19

obj <- readRDS("data/imvigor210/imvigor210_clean.rds"); E <- obj$E; pd <- obj$pd
ip <- factor(pd[["Immune phenotype"]], c("desert","excluded","inflamed"))
up <- function(g) rownames(E)[match(intersect(toupper(g), toupper(rownames(E))), toupper(rownames(E)))]
prop_high <- function(g){ g<-up(g); g<-g[!is.na(g)]; if(!length(g)) return(rep(NA_real_,ncol(E)))
  b<-t(apply(E[g,,drop=FALSE],1,function(v) as.integer(v>=median(v,na.rm=TRUE)))); colMeans(b,na.rm=TRUE) }

S  <- prop_high(PAN$SENSORY)
I  <- prop_high(PAN$INTPOS) - prop_high(PAN$INTNEG)
Ex <- prop_high(PAN$EXECPOS) - prop_high(PAN$EXECNEG)
Sh<-as.integer(S>=0.5); Ih<-as.integer(I>=0); Eh<-as.integer(Ex>=0); pat<-paste0(Sh,Ih,Eh)

NAME <- c("111"="Inflamed","110"="Execution Failed","101"="Excluded (1,0,1)",
          "100"="Excluded (1,0,0)","011"="Sensing Failed","010"="I-only","001"="E-only","000"="Desert")
n<-length(pat)
tab<-as.data.frame(table(factor(pat,names(NAME))),responseName="n"); names(tab)[1]<-"pattern"
tab$archetype<-NAME[as.character(tab$pattern)]; tab$pct<-round(100*tab$n/n,1)
tab<-tab[order(-tab$n),c("pattern","archetype","n","pct")]
cat(sprintf("N=%d\n== LEGACY-panel S/I/E pattern frequencies ==\n", n)); print(tab,row.names=FALSE)
cat(sprintf("\nSensing Failed (0,1,1): %d (%.1f%%)\n", sum(pat=="011"), 100*mean(pat=="011")))
cat("\nSpearman among layer scores (S-E coupling? I orthogonal?):\n")
print(round(cor(cbind(S=S,I=I,E=Ex),method="spearman",use="complete.obs"),2))
cat("\nmean layer score by immune phenotype (sanity):\n")
print(round(rbind(S=tapply(S,ip,mean,na.rm=TRUE),I=tapply(I,ip,mean,na.rm=TRUE),E=tapply(Ex,ip,mean,na.rm=TRUE)),3))
cat("\npattern x phenotype:\n"); print(addmargins(table(archetype=NAME[pat],phenotype=ip)))
saveRDS(data.frame(sample=colnames(E),S=S,I=I,E=Ex,pattern=pat,phenotype=as.character(ip)),
        "archetype_scores_legacy_IMvigor210.rds")
cat("\ndone\n")
