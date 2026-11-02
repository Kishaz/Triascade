# =====================================================================
# Extra public ICB cohorts added to the pan-cancer validation, SAME evaluator
# (icb_eval). Airtight, baseline-only, response from GEO SOFT.
#   GSE165252  PERFECT (van den Ende 2021)  esophageal adenocarcinoma (EAC)
#              neoadjuvant atezolizumab + chemoradiotherapy  -> BASELINE biopsies only.
#              NOTE: chemo+ICB (not pure ICB) -> label as such, like I-SPY2/NeoTRIP.
#              Expression = Ensembl norm counts; sample codes live in Sample_description;
#              timepoint in Sample_title; response = responder/non-responder.
# =====================================================================
source("run_pancancer_geo.R")                 # RES = all current cohorts (IMvigor,Hugo,Riaz,Cho,Jung)

## ---- GSE165252 EAC (PERFECT), baseline only ----
parse_perfect <- function(f){
  ln<-readLines(f); idx<-cumsum(grepl("^\\^SAMPLE",ln))
  do.call(rbind, lapply(split(ln,idx), function(b){
    ti<-sub("^!Sample_title = ","",grep("^!Sample_title ",b,value=TRUE))[1]
    de<-sub("^!Sample_description = ","",grep("^!Sample_description ",b,value=TRUE))[1]
    rs<-sub(".*response: *","",grep("response:",b,value=TRUE)[1])
    data.frame(code=de, title=ti, resp=rs, stringsAsFactors=FALSE) }))
}
pf <- parse_perfect("data/_dl/GSE165252_full.txt")
pf <- pf[!is.na(pf$code) & pf$code!="NONE", ]
base <- pf[grepl("baseline", pf$title, ignore.case=TRUE) & pf$resp %in% c("responder","non-responder"), ]
Eac <- to_symbol(read_mat("data/_dl/GSE165252_norm.txt.gz"), "ENSEMBL")
Eac <- Eac[, intersect(colnames(Eac), base$code), drop=FALSE]
yac <- setNames(ifelse(base$resp=="responder",1L,0L), base$code)[colnames(Eac)]
RES$GSE165252 <- icb_eval(Eac, yac, "GSE165252","EAC","atezo+CRT(neoadj)")

## ---- xmuyulab/ims_gene_signature processed ICB cohorts --------------------
# Gene SYMBOLS; quantile-normalised log2 (rank/median-robust -> fine for cohort-median
# + single-patient selfnorm scoring). response: 1 = CR/PR (responder); 0 = PD, -1 = SD
# -> non-responder. Direction OS-verified (resp=1 long OS/few deaths; 0 short/high;
# -1 intermediate = SD). PRE-treatment biopsies only. cli/exp files in data/_dl/<stem>_*.
xmu <- function(stem, key, cancer, drug){
  cli <- read.csv(sprintf("data/_dl/%s_cli.csv", stem), row.names=1)
  exp <- read_mat(sprintf("data/_dl/%s_exp.csv", stem))
  pre <- rownames(cli)[cli$treatment=="pre"]
  E   <- exp[, intersect(colnames(exp), pre), drop=FALSE]
  y   <- setNames(ifelse(cli$response==1, 1L, 0L), rownames(cli))[colnames(E)]
  icb_eval(E, y, key, cancer, drug)
}
RES$Gide2019     <- xmu("gide",       "Gide2019",    "melanoma", "aPD1/ipi+pembro(Pre)")  # PRJEB23709
RES$Liu2019      <- xmu("mel_liu19",  "Liu2019",     "melanoma", "aPD1(Pre)")
RES$VanAllen2015 <- xmu("mel_van",    "VanAllen2015","melanoma", "ipi/aCTLA4(Pre)")        # CTLA4
RES$Puch         <- xmu("mel_puch",   "Puch",        "melanoma", "aPD1(Pre)")
RES$Kim2018gas   <- xmu("gas_korean", "Kim2018",     "gastric",  "pembro(Pre)")            # new tumour type

## ---- PRINCE / PICI0002 (pancreatic; gem/nab-paclitaxel + nivolumab +/- CD40), BASELINE ----
# Source: ParkerICI/prince-trial-data (public). Per-sample TPM reports merged across the
# NIVOLUMAB-containing arms only (A1 chemo+nivo, C2 chemo+nivo+CD40); arm B2 (CD40 only, no
# checkpoint) excluded. response = RECIST Best Overall Response CR/PR=1 vs SD/PD=0.
# chemo+ICB (not pure ICB), pancreatic = ICB-resistant -> interpret like EAC, separate bucket.
pr_gx <- read.csv("data/_dl/prince/gx_meta.csv", check.names=FALSE)
pr_cl <- read.csv("data/_dl/prince/clinical.csv", check.names=FALSE)
pr_cl$y <- ifelse(toupper(pr_cl$`Best Overall Response`) %in% c("COMPLETE RESPONSE","PARTIAL RESPONSE"),1L,
            ifelse(toupper(pr_cl$`Best Overall Response`) %in% c("STABLE DISEASE","PROGRESSIVE DISEASE"),0L,NA))
g <- merge(pr_gx[pr_gx$timepoint.id=="Baseline",], pr_cl[,c("Deidentified.ID","Actual Arm","y")], by="Deidentified.ID")
g <- g[g$`Actual Arm` %in% c("A1","C2") & !is.na(g$y), ]
g <- g[!duplicated(g$Deidentified.ID), ]
rd <- lapply(g$Filename, function(f){ d<-read.delim(file.path("data/_dl/prince/Data",f), check.names=FALSE)
             d<-d[!duplicated(d$`Gene Symbol`),]; setNames(d$TPM, d$`Gene Symbol`) })
genes <- Reduce(intersect, lapply(rd, names))
Epr <- sapply(rd, function(v) v[genes]); rownames(Epr) <- genes; colnames(Epr) <- g$Deidentified.ID
ypr <- setNames(g$y, g$Deidentified.ID)
RES$PRINCE <- icb_eval(Epr, ypr, "PRINCE","pancreatic","chemo+nivo±CD40")

saveRDS(RES, "pancancer_validation.rds")
cat("\n=== SUMMARY: Integration AUC across cohorts ===\n")
cat(sprintf("%-11s %-10s %-18s n   R/NR   cohort  single(self)\n","cohort","cancer","drug"))
for(r in RES) cat(sprintf("%-11s %-10s %-18s %-4d %d/%-4d %.3f   %.3f\n",
  r$cohort, r$cancer, r$drug, r$n, r$nR, r$nNR, r$auc["cohort","Integration"], r$auc["single_self","Integration"]))
cat("\ndone\n")
