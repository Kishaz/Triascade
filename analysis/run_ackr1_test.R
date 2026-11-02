# =====================================================================
# ACKR1/DARC mechanistic test in IMvigor210 (n=348).
# ACKR1 = atypical chemokine receptor on postcapillary-venule endothelium
# (chemokine scavenging / T-cell entry). Hypothesis: it tracks the
# "excluded" / Integration-suppressed class and the angiogenic/CAF program.
# This is a DESCRIPTIVE/mechanistic test; ACKR1 is NOT in the scoring panel.
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
obj <- readRDS("data/imvigor210/imvigor210_clean.rds"); E<-obj$E; pd<-obj$pd
resp <- ifelse(pd$binaryResponse=="CR/PR",1L,ifelse(pd$binaryResponse=="SD/PD",0L,NA))
ip   <- factor(pd[["Immune phenotype"]], c("desert","excluded","inflamed"))

G <- "ACKR1"
cat("ACKR1 present in IMvigor210:", G %in% rownames(E), "\n")
a <- E[G, ]
cat(sprintf("ACKR1 log2TPM: mean %.2f, range %.2f-%.2f\n\n", mean(a), min(a), max(a)))

## module/layer scores (same as taxonomy)
dp<-fread("../ICI-study-redesign/provenance-review-tight/derivation/ICI_PANEL_v1_locked.csv",data.table=FALSE);dp$gene<-toupper(dp$gene)
MOD<-list(Sense=dp$gene[grepl("^Sensory",dp$layer)],Traffic=dp$gene[dp$layer=="Integration (chemokine trafficking)"],
 Treg_ckpt=dp$gene[dp$layer=="Integration (Treg / checkpoint regulation)"],Stromal_ECM=dp$gene[dp$layer=="Resistance (stromal / ECM-coagulation)"],
 Angio_CAF=dp$gene[dp$layer=="Resistance (angiogenic / CAF-TGFb)"],Execute=dp$gene[grepl("^Execution",dp$layer)])
ph<-function(g){g<-intersect(g,rownames(E));b<-t(apply(E[g,,drop=FALSE],1,function(v) as.integer(v>=median(v,na.rm=TRUE))));colMeans(b,na.rm=TRUE)}
Z<-function(v)(v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE);M<-sapply(MOD,ph);Mz<-apply(M,2,Z)
Sense<-Mz[,"Sense"];Exec<-Mz[,"Execute"];Integ<-Z(M[,"Traffic"]-(M[,"Treg_ckpt"]+M[,"Stromal_ECM"]+M[,"Angio_CAF"])/3)
tiers<-cbind(Sense,Integration=Integ,Execution=Exec);wmin<-colnames(tiers)[apply(tiers,1,which.min)]
klass<-ifelse(apply(tiers,1,min)>=-0.5,"Cascade-intact",ifelse(wmin=="Sense","Sense-low",ifelse(wmin=="Integration","Integration-suppressed","Execution-low")))
klass<-factor(klass,c("Sense-low","Integration-suppressed","Execution-low","Cascade-intact"))

cat("(1) ACKR1 by immune phenotype (Kruskal p):\n")
print(round(tapply(a,ip,mean),2)); cat(sprintf("   Kruskal p = %.2e\n", kruskal.test(a,ip)$p.value))
cat(sprintf("   excluded-vs-inflamed AUC = %.3f | excluded-vs-desert AUC = %.3f\n",
   as.numeric(auc(as.integer(ip[ip%in%c("excluded","inflamed")]=="excluded"), a[ip%in%c("excluded","inflamed")],quiet=TRUE,direction="<")),
   as.numeric(auc(as.integer(ip[ip%in%c("excluded","desert")]=="excluded"), a[ip%in%c("excluded","desert")],quiet=TRUE,direction="<"))))

cat("\n(2) ACKR1 by limiting-layer class:\n")
print(round(tapply(a,klass,mean),2)); cat(sprintf("   Kruskal p = %.2e\n", kruskal.test(a,klass)$p.value))

cat("\n(3) ACKR1 correlation with module/layer scores (Spearman):\n")
cc<-sapply(c("Sense","Traffic","Treg_ckpt","Stromal_ECM","Angio_CAF","Execute"),function(m) cor(a,Mz[,m],method="spearman"))
cc<-c(cc, Integration_net=cor(a,Integ,method="spearman"))
print(round(sort(cc),2))

cat("\n(4) ACKR1 vs response:\n")
cat(sprintf("   mean NR %.2f vs R %.2f | response AUC %.3f | Wilcoxon p %.3f\n",
   mean(a[resp==0],na.rm=T),mean(a[resp==1],na.rm=T),
   as.numeric(auc(resp[!is.na(resp)],a[!is.na(resp)],quiet=TRUE,direction="<")),
   wilcox.test(a~resp)$p.value))

cat("\n(5) within Integration-suppressed: ACKR1 vs the suppressive modules it may mark:\n")
iss<-klass=="Integration-suppressed"
cat(sprintf("   ACKR1~Angio_CAF r=%.2f | ACKR1~Stromal r=%.2f (within Integration-suppressed, n=%d)\n",
   cor(a[iss],Mz[iss,"Angio_CAF"],method="spearman"),cor(a[iss],Mz[iss,"Stromal_ECM"],method="spearman"),sum(iss)))
