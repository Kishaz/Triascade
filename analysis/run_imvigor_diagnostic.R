# Diagnostic: do the panel MODULES carry the established IMvigor210 biology,
# independent of any taxonomy rule? If yes, the null is the classifier, not
# the panel. Checks: (a) module z by immune phenotype (desert/excluded/
# inflamed), (b) module z by response, (c) collinearity of layer scores,
# (d) per-module AUC for response and for inflamed-vs-desert.
suppressPackageStartupMessages({ library(data.table); library(pROC) })
D <- "data/imvigor210"
DER <- "../ICI-study-redesign/provenance-review-tight/derivation/ICI_PANEL_v1_locked.csv"
dp  <- fread(DER, data.table = FALSE); dp$gene <- toupper(dp$gene)
MOD <- list(
  Sense=dp$gene[grepl("^Sensory",dp$layer)],
  Traffic=dp$gene[dp$layer=="Integration (chemokine trafficking)"],
  Treg_ckpt=dp$gene[dp$layer=="Integration (Treg / checkpoint regulation)"],
  Stromal_ECM=dp$gene[dp$layer=="Resistance (stromal / ECM-coagulation)"],
  Angio_CAF=dp$gene[dp$layer=="Resistance (angiogenic / CAF-TGFb)"],
  Execute=dp$gene[grepl("^Execution",dp$layer)])

cnt <- as.data.frame(fread(file.path(D,"exmat_censored_IMvigor210.csv"))); cnt[[1]]<-NULL; cnt<-as.matrix(cnt)
fd  <- as.data.frame(fread(file.path(D,"fData_IMvigor210.csv")))
sym <- toupper(ifelse(is.na(fd$symbol)|fd$symbol=="",fd$Symbol,fd$symbol)); len<-as.numeric(fd$length)
keep<- !is.na(sym)&sym!=""&!is.na(len)&len>0; cnt<-cnt[keep,];sym<-sym[keep];len<-len[keep]
tpm <- sweep(cnt/len,2,colSums(cnt/len),"/")*1e6; tpm<-rowsum(tpm,sym); E<-log2(tpm+1)
pd  <- as.data.frame(fread(file.path(D,"pData_IMvigor210.csv"))); rownames(pd)<-pd[[1]]; pd<-pd[colnames(E),]
resp<- ifelse(pd$binaryResponse=="CR/PR",1L,ifelse(pd$binaryResponse=="SD/PD",0L,NA))
ip  <- factor(pd$Immune.phenotype, c("desert","excluded","inflamed"))

prop_high <- function(g){g<-intersect(g,rownames(E));b<-t(apply(E[g,,drop=FALSE],1,function(v) as.integer(v>=median(v,na.rm=TRUE))));colMeans(b,na.rm=TRUE)}
Z<-function(v)(v-mean(v,na.rm=TRUE))/sd(v,na.rm=TRUE)
M<-sapply(MOD,prop_high); Mz<-apply(M,2,Z)
Sense<-Mz[,"Sense"]; Exec<-Mz[,"Execute"]
Integ<-Z(M[,"Traffic"]-(M[,"Treg_ckpt"]+M[,"Stromal_ECM"]+M[,"Angio_CAF"])/3)

cat("(a) mean module z by IMMUNE PHENOTYPE (should rise desert->excluded->inflamed for effector;\n    suppressive/stroma should peak in EXCLUDED):\n")
by_ip <- t(sapply(colnames(Mz), function(m) tapply(Mz[,m], ip, mean)))
print(round(by_ip,2))
cat("\n   layer scores by immune phenotype:\n")
print(round(rbind(Sense=tapply(Sense,ip,mean), Integration=tapply(Integ,ip,mean),
                  Execution=tapply(Exec,ip,mean)),2))

cat("\n(b) mean module z by RESPONSE (R=CR/PR vs NR):\n")
print(round(t(sapply(colnames(Mz), function(m) tapply(Mz[,m], factor(resp,c(0,1),c("NR","R")), mean))),2))

cat("\n(c) collinearity - Pearson r among layer scores and key modules:\n")
cm <- cbind(Sense=Sense, Integration=Integ, Execution=Exec,
            Traffic=Mz[,"Traffic"], Effector=Mz[,"Execute"], Suppr=Z((M[,"Treg_ckpt"]+M[,"Stromal_ECM"]+M[,"Angio_CAF"])/3))
print(round(cor(cm, use="complete.obs"),2))

cat("\n(d) per-module AUC for RESPONSE and for INFLAMED-vs-DESERT:\n")
a_resp <- sapply(colnames(Mz), function(m){s<-Mz[!is.na(resp),m];as.numeric(auc(resp[!is.na(resp)],s,quiet=TRUE,direction="<"))})
infdes <- ip %in% c("inflamed","desert"); yid <- as.integer(ip[infdes]=="inflamed")
a_id <- sapply(colnames(Mz), function(m) as.numeric(auc(yid, Mz[infdes,m], quiet=TRUE, direction="<")))
print(round(data.frame(AUC_response=a_resp, AUC_inflamed_vs_desert=a_id),3))
cat(sprintf("\n  effector(Execute) module: response AUC %.3f, inflamed-vs-desert AUC %.3f\n",
            a_resp["Execute"], a_id["Execute"]))
