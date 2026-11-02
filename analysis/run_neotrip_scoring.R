suppressPackageStartupMessages({library(data.table)}); source("ici_locked_reference.R"); PANEL<-load_panel("ICI_98gene_panel.csv")
x<-as.data.frame(fread("data/NeoTRIP_expr.txt.gz"))
g<-toupper(x[["HGNC_approved_symbol"]]); x<-x[,grep("_Baseline$",colnames(x)),drop=FALSE]   # pretreatment only
x<-as.matrix(x); rownames(x)<-g; x<-x[!is.na(g)&g!="",]; x<-log2(x+1)
cat(sprintf("NeoTRIP baseline: %d genes x %d samples | panel present %d/98 (missing: %s)\n",
  nrow(x),ncol(x),length(intersect(unique(unlist(PANEL)),rownames(x))),
  paste(setdiff(unique(unlist(PANEL)),rownames(x)),collapse=",")))
cm<-score_cohort_median(x,PANEL); cm<-cm[match(colnames(x),cm$sample),]
lk<-do.call(rbind,lapply(seq_len(ncol(x)),function(i){r<-fit_locked_reference(x[,-i,drop=FALSE],PANEL);score_locked(x[,i,drop=FALSE],r,"quantile")}))
lk<-lk[match(colnames(x),lk$sample),]
con<-concordance(cm,lk)
cat("\nLocked single-patient (LOO) vs cohort-median on NeoTRIP (FFPE RNA-seq TPM):\n"); print(con,row.names=FALSE)
cat(sprintf("\nArchetype agreement: %.3f\n",archetype_agreement(cm,lk)))
# gene-call agreement
bc<-t(apply(x[intersect(unique(unlist(PANEL)),rownames(x)),],1,function(v)as.integer(v>=median(v))))
ag<-sapply(seq_len(ncol(x)),function(i){r<-fit_locked_reference(x[,-i,drop=FALSE],PANEL);gg<-intersect(rownames(bc),r$genes)
  cl<-sapply(gg,function(gi)as.integer(.ref_quantile(x[gi,i],r$ref_sorted[[gi]])>=0.5)); mean(cl==bc[gg,i])})
cat(sprintf("Gene-call agreement (locked LOO vs cohort-median): mean %.3f (range %.2f-%.2f)\n",mean(ag),min(ag),max(ag)))
one<-score_locked(x[,1,drop=FALSE],fit_locked_reference(x[,-1,drop=FALSE],PANEL),"quantile")
cat(sprintf("Single-patient property holds: %s\n", isTRUE(all.equal(one$Integration,lk$Integration[1]))))
