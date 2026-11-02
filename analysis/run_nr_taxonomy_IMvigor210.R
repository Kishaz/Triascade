# =====================================================================
# NR TAXONOMY on IMvigor210 (n=348 urothelial, atezolizumab, single-arm)
# Mariathasan et al. 2018. Uses the ALIGNED cds object (imvigor210_clean.rds;
# entrez-keyed counts -> TPM, log2). The earlier mimifp CSVs had scrambled
# feature order and are NOT used.
# Classify each tumour by its limiting cascade LAYER (Sense / Integration /
# Execution), outcome-blind, then:
#   (a) do modules carry the biology (module z by immune phenotype)?
#   (1) limiting-layer distribution; (2) response by layer;
#   (3) external validation vs Mariathasan immune phenotype;
#   (4) dominant suppressive module in Integration-suppressed NR.
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })

obj <- readRDS("data/imvigor210/imvigor210_clean.rds"); E <- obj$E; pd <- obj$pd
resp <- ifelse(pd$binaryResponse=="CR/PR",1L, ifelse(pd$binaryResponse=="SD/PD",0L,NA))
ip   <- factor(pd[["Immune phenotype"]], c("desert","excluded","inflamed"))

DER <- "../ICI-study-redesign/provenance-review-tight/derivation/ICI_PANEL_v1_locked.csv"
dp  <- fread(DER, data.table = FALSE); dp$gene <- toupper(dp$gene)
MOD <- list(
  Sense       = dp$gene[grepl("^Sensory", dp$layer)],
  Traffic     = dp$gene[dp$layer == "Integration (chemokine trafficking)"],
  Treg_ckpt   = dp$gene[dp$layer == "Integration (Treg / checkpoint regulation)"],
  Stromal_ECM = dp$gene[dp$layer == "Resistance (stromal / ECM-coagulation)"],
  Angio_CAF   = dp$gene[dp$layer == "Resistance (angiogenic / CAF-TGFb)"],
  Execute     = dp$gene[grepl("^Execution", dp$layer)])

prop_high <- function(g){g<-intersect(g,rownames(E)); if(!length(g)) return(rep(NA,ncol(E)))
  b<-t(apply(E[g,,drop=FALSE],1,function(v) as.integer(v>=median(v,na.rm=TRUE)))); colMeans(b,na.rm=TRUE)}
Z<-function(v)(v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE)
cat("panel genes measured per module:\n"); print(sapply(MOD, function(g) length(intersect(g,rownames(E)))))
M<-sapply(MOD,prop_high); Mz<-apply(M,2,Z)
Sense<-Mz[,"Sense"]; Exec<-Mz[,"Execute"]
Integ<-Z(M[,"Traffic"]-(M[,"Treg_ckpt"]+M[,"Stromal_ECM"]+M[,"Angio_CAF"])/3)

cat("\n(a) mean module z by IMMUNE PHENOTYPE (sanity: effector/Sense rise toward inflamed):\n")
print(round(t(sapply(colnames(Mz), function(m) tapply(Mz[,m], ip, mean))),2))
cat("   layer scores by immune phenotype:\n")
print(round(rbind(Sense=tapply(Sense,ip,mean),Integration=tapply(Integ,ip,mean),Execution=tapply(Exec,ip,mean)),2))

FLOOR <- -0.5
tiers <- cbind(Sense=Sense, Integration=Integ, Execution=Exec)
minz<-apply(tiers,1,min); wmin<-colnames(tiers)[apply(tiers,1,which.min)]
klass<-ifelse(minz>=FLOOR,"Cascade-intact",
       ifelse(wmin=="Sense","Sense-low (cold)",
       ifelse(wmin=="Integration","Integration-suppressed","Execution-low")))
klass<-factor(klass,c("Sense-low (cold)","Integration-suppressed","Execution-low","Cascade-intact"))

cat("\n(1) limiting-layer distribution (all / NR / R):\n")
print(rbind(all=table(klass), NR=table(klass[resp==0]), R=table(klass[resp==1])))
cat("\n(2) response rate (CR/PR %) by limiting layer:\n")
rr<-tapply(resp,klass,function(z)100*mean(z,na.rm=TRUE)); nn<-tapply(resp,klass,function(z)sum(!is.na(z)))
print(data.frame(class=names(rr),response_pct=round(rr,0),n=as.integer(nn),row.names=NULL))
cat("\n(3) EXTERNAL VALIDATION - limiting layer vs immune phenotype (row proportions):\n")
ct<-table(klass,ip); print(round(prop.table(ct,1),2))
cat("\n(4) Integration-suppressed NON-RESPONDERS: dominant suppressive module (which agent?):\n")
iss<-which(klass=="Integration-suppressed" & resp==0)
print(round(sort(colMeans(Mz[iss,c("Treg_ckpt","Stromal_ECM","Angio_CAF"),drop=FALSE]),decreasing=TRUE),2))
cat(sprintf("   (n=%d)\n", length(iss)))
