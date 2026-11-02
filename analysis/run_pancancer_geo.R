# =====================================================================
# Add the downloaded GEO ICB cohorts to the pan-cancer validation, using the
# SAME evaluator (icb_eval) as the local cohorts. Airtight handling per cohort:
#   GSE91061  Riaz 2017   melanoma  nivolumab  -> Entrez->symbol; PRE samples only
#   GSE126044 Cho 2020    NSCLC     anti-PD1   -> symbols; all pre-treatment
#   GSE135222 Jung 2019   NSCLC     anti-PD(L)1-> Ensembl->symbol; DCB = PFS>=6mo
# Response/timepoint parsed from GEO SOFT (data/_dl/*_soft.txt).
# =====================================================================
source("run_pancancer_validation.R")   # -> icb_eval, REF_EXT, PANEL (prints+saves 2 local cohorts)
RES <- readRDS("pancancer_validation.rds")
suppressPackageStartupMessages({ library(org.Hs.eg.db); library(AnnotationDbi) })

read_mat <- function(f){ dt<-fread(f); rn<-as.character(dt[[1]])
  M<-as.matrix(dt[,-1,with=FALSE]); rownames(M)<-rn; storage.mode(M)<-"numeric"; M }
to_symbol <- function(M, from){ ids<-rownames(M); if(from=="ENSEMBL") ids<-sub("\\..*$","",ids)
  sym<-suppressMessages(mapIds(org.Hs.eg.db, keys=ids, column="SYMBOL", keytype=from, multiVals="first"))
  ok<-!is.na(sym) & !is.na(rowMeans(M,na.rm=TRUE)); M<-M[ok,,drop=FALSE]; sym<-sym[ok]
  o<-order(rowMeans(M,na.rm=TRUE),decreasing=TRUE); M<-M[o,,drop=FALSE]; sym<-sym[o]   # keep highest-expressed per symbol
  d<-!duplicated(sym); M<-M[d,,drop=FALSE]; rownames(M)<-sym[d]; M }
soft_df <- function(f){ ln<-readLines(f); idx<-cumsum(grepl("^\\^SAMPLE",ln))
  o<-lapply(split(ln,idx), function(b){ ti<-sub("^!Sample_title = ","",grep("^!Sample_title ",b,value=TRUE))[1]
    ch<-sub("^!Sample_characteristics_ch1 = ","",grep("^!Sample_characteristics_ch1 ",b,value=TRUE))
    list(title=ti, kv=setNames(trimws(sub("^[^:]*: *","",ch)), make.names(trimws(sub(":.*$","",ch))))) })
  ti<-sapply(o,function(x)x$title); keys<-unique(unlist(lapply(o,function(x)names(x$kv))))
  d<-data.frame(title=ti, stringsAsFactors=FALSE)
  for(k in keys) d[[k]]<-sapply(o, function(x) if(k%in%names(x$kv)) x$kv[[k]] else NA_character_); d }

## ---- GSE91061 (Riaz, melanoma, nivolumab) : PRE only, CR/PR vs SD/PD ----
ri <- to_symbol(read_mat("data/_dl/GSE91061_fpkm.csv.gz"), "ENTREZID")
m <- soft_df("data/_dl/GSE91061_soft.txt")
vis <- m[[grep("visit", names(m), ignore.case=TRUE, value=TRUE)[1]]]; rsp <- m[["response"]]
pre <- m$title[vis=="Pre" & rsp %in% c("PRCR","SD","PD")]   # Riaz codes responders as "PRCR"; UNK excluded
Eri <- ri[, intersect(colnames(ri), pre), drop=FALSE]
yri <- setNames(ifelse(rsp=="PRCR",1L,0L), m$title)[colnames(Eri)]
RES$GSE91061 <- icb_eval(Eri, yri, "GSE91061","melanoma","nivo(Pre)")

## ---- GSE126044 (Cho, NSCLC, anti-PD1) : all pre-treatment ----
cc <- read_mat("data/_dl/GSE126044_counts.txt.gz")       # symbols already
m2 <- soft_df("data/_dl/GSE126044_soft.txt")
yco <- setNames(ifelse(grepl("non", m2[["patient.response"]], ignore.case=TRUE),0L,1L),
                sub("^RNA-seq_","",m2$title))[colnames(cc)]
RES$GSE126044 <- icb_eval(cc, yco, "GSE126044","NSCLC","anti-PD1")

## ---- GSE135222 (Jung, NSCLC) : DCB = PFS>=180d ; progressed<180 = NDB ; censored<180 excluded ----
ju <- to_symbol(read_mat("data/_dl/GSE135222_exp.tsv.gz"), "ENSEMBL")
m3 <- soft_df("data/_dl/GSE135222_soft.txt")
pfs <- as.numeric(m3[[grep("^progression.free", names(m3), value=TRUE)[1]]])
pft <- as.numeric(m3[["pfs.time"]])
dcb <- ifelse(pft>=180, 1L, ifelse(pfs==1 & pft<180, 0L, NA_integer_))
yju <- setNames(dcb, gsub(" ","",m3$title))[colnames(ju)]
RES$GSE135222 <- icb_eval(ju, yju, "GSE135222","NSCLC","aPD(L)1 DCB>=6mo")

saveRDS(RES, "pancancer_validation.rds")
cat("\n=== SUMMARY: Integration AUC across cohorts ===\n")
cat(sprintf("%-11s %-10s %-16s n   R/NR   cohort  single(self)\n","cohort","cancer","drug"))
for(r in RES) cat(sprintf("%-11s %-10s %-16s %-4d %d/%-4d %.3f   %.3f\n",
  r$cohort, r$cancer, r$drug, r$n, r$nR, r$nNR, r$auc["cohort","Integration"], r$auc["single_self","Integration"]))
cat("\ndone\n")
