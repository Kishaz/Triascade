# =====================================================================
# PANEL GAP ANALYSIS. Audit canonical ICB genes (by layer + sub-function)
# against: locked 108, the broader derived panels (133, 147), and the full
# literature recurrence universe (2280 genes). Classify each canonical gene:
#   PRESENT (in locked) | RECOVERABLE (recurrence>=5, not in locked) |
#   SUBTHRESHOLD (recurrence 1-4) | ABSENT (not in substrate at all).
# This is an audit for coverage, NOT hand-curation: it flags where the
# literature-derived panel under-covers established biology.
# =====================================================================
suppressPackageStartupMessages(library(data.table))
B <- "../ICI-study-redesign/provenance-review-tight/derivation"
locked <- toupper(fread("ICI_panel_LOCKED_v1.csv",data.table=FALSE)$gene)
p133 <- toupper(fread(file.path(B,"FINAL_DERIVED_PANEL.csv"),data.table=FALSE)$gene)
p147 <- toupper(fread(file.path(B,"ICI_PANEL_denovo.csv"),data.table=FALSE)$gene)
rec  <- fread(file.path(B,"core_results_v4/gene_recurrence.csv"),data.table=FALSE); rec$gene<-toupper(rec$gene)
recn <- setNames(rec$recurrence, rec$gene)

CANON <- list(
 "Sense: antigen-presentation" = c("B2M","HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","TAP1","TAP2","TAPBP","NLRC5","PSMB8","PSMB9","CALR","CIITA","HLA-DRA","HLA-DRB1","CD74"),
 "Sense: IFN signalling"       = c("STAT1","STAT2","IRF1","IRF9","JAK1","JAK2","IFNGR1","IFNGR2","ISG15","MX1","OAS1","GBP1"),
 "Sense: innate sensing"       = c("CGAS","MB21D1","STING1","TMEM173","DDX58","IFIH1","TLR3","TLR4","TLR9"),
 "Sense: recruit chemokine"    = c("CXCL9","CXCL10","CXCL11"),
 "Integr(+): costimulation"    = c("CD28","CD80","CD86","ICOS","ICOSLG","CD40","CD40LG","CD27","CD70","TNFRSF9","TNFSF9","TNFRSF4","TNFSF4","TNFRSF18"),
 "Integr(-): coinhibition"     = c("PDCD1","CD274","PDCD1LG2","CTLA4","LAG3","HAVCR2","TIGIT","BTLA","VSIR","CD276","VTCN1","IDO1","LGALS9","CEACAM1","ENTPD1","NT5E","CD96"),
 "Integr(-): Treg"             = c("FOXP3","IL2RA","IKZF2","CCR8","IL10"),
 "Integr: trafficking"         = c("CCL5","CCL2","CCL4","CCL19","CCL21","CCR7","CXCL13","SELL","CXCR3","CX3CL1"),
 "Integr(-): myeloid/stroma"   = c("CD68","CD163","CSF1R","CSF1","ARG1","MRC1","NOS2","MARCO","SIGLEC1","VEGFA","TGFB1","TGFB2","TGFB3","TGFBR1","TGFBR2","FAP","ACTA2","COL1A1","COL1A2","COL3A1","CXCL12","PDGFRB","IL6","LGALS1","SPP1"),
 "Execution: cytotoxic"        = c("GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","GNLY","NKG7","KLRD1","KLRK1","FASLG","FAS","IFNG","TNF"),
 "Execution: T/NK identity"    = c("CD8A","CD8B","CD3D","CD3E","CD3G","EOMES","TBX21","KLRG1","ZAP70","CX3CR1"),
 "Execution: exhaustion"       = c("TOX","LAYN","TCF7"))

classify <- function(g){ if(g %in% locked) "PRESENT"
  else if(!is.na(recn[g]) && recn[g]>=5) "RECOVERABLE(rec>=5)"
  else if(!is.na(recn[g])) "SUBTHRESHOLD(rec1-4)"
  else "ABSENT(not in substrate)" }

cat("GENE                 status                    rec  in133 in147\n")
summ <- list()
for(grp in names(CANON)){
  cat(sprintf("\n--- %s ---\n", grp))
  st <- sapply(CANON[[grp]], classify)
  for(g in CANON[[grp]]) cat(sprintf("  %-18s %-24s %3s   %-3s  %-3s\n", g, classify(g),
      ifelse(is.na(recn[g]),"-",recn[g]), ifelse(g%in%p133,"y","."), ifelse(g%in%p147,"y",".")))
  summ[[grp]] <- table(factor(st, c("PRESENT","RECOVERABLE(rec>=5)","SUBTHRESHOLD(rec1-4)","ABSENT(not in substrate)")))
}
cat("\n================ SUMMARY (counts per sub-function) ================\n")
S <- do.call(rbind, summ)
print(S)
cat("\nRECOVERABLE genes (recurrent >=5 in literature but below the Grade A+B gate) - the defensible add candidates:\n")
allc <- unique(unlist(CANON))
r <- allc[sapply(allc, classify)=="RECOVERABLE(rec>=5)"]
r <- r[order(-recn[r])]
for(g in r) cat(sprintf("  %-10s rec=%d  %s\n", g, recn[g], ifelse(g%in%p147,"(in 147-gene denovo)", ifelse(g%in%p133,"(in 133)","(not in any panel)"))))
