# =====================================================================
# ESCAPE PROBE: why do the cascade-intact ICB-arm non-responders fail,
# when the 3 ICI layers see no defect?  (legacy 98, I-SPY2)
# Contrast within the cascade-intact, pembro-treated group:
#   NR (n~12)  vs  R (n~23)
# on resistance programs ORTHOGONAL to the 98-gene panel (panel genes removed),
# since the panel's own axes are high in these tumours by definition.
# Exploratory (small n): report mean-z gap + Wilcoxon p, FDR across programs.
# =====================================================================
suppressPackageStartupMessages({ library(data.table) })
asg <- read.csv("ICI_98_module_assignment.csv", stringsAsFactors = FALSE); PANEL98 <- toupper(asg$gene)
MOD <- list(Sense=asg$gene[asg$module=="M1"], Traffic=asg$gene[asg$module=="M2"],
            Stromal_ECM=asg$gene[asg$module=="M3"], Angio_CAF=asg$gene[asg$module=="M4"],
            Execute=asg$gene[asg$module=="M5"], Checkpoint=asg$gene[asg$module=="M6"])

# resistance programs OUTSIDE the panel (literature ICB-resistance axes)
ESCAPE <- list(
  Proliferation   = c("MKI67","TOP2A","CCNB1","CDK1","BIRC5","PCNA","AURKA","CCNB2","CDC20"),
  WNT_bcatenin    = c("CTNNB1","LEF1","AXIN2","NOTUM","RNF43","DKK1","TCF7L2"),
  Hypoxia         = c("VEGFA","CA9","SLC2A1","LDHA","HIF1A","PGK1","ADM"),
  Immsupp_metab   = c("IDO1","ARG1","NOS2","CD274","PDCD1LG2","ENTPD1","NT5E"),
  TGFb_EMT        = c("TGFB1","TWIST1","ZEB1","SNAI1","VIM","SERPINE1","LOX"))
# drop any panel genes so the probe is orthogonal to the framework
ESCAPE <- lapply(ESCAPE, function(g){ gg<-setdiff(toupper(g), PANEL98); gg })
cat("escape programs (panel genes removed):\n"); print(lengths(ESCAPE))

parse_sm <- function(f) { ln <- readLines(gzfile(f), warn = FALSE); ln <- ln[grepl("^!Sample_", ln)]; out <- list()
  for (r in strsplit(ln, "\t")) { k <- sub("^!Sample_", "", r[1]); v <- gsub('"', '', r[-1])
    if (grepl("^characteristics", k)) { i <- which(grepl(":", v))[1]; if (is.na(i)) next
      nm <- make.names(trimws(sub(":.*$", "", v[i]))); v <- trimws(sub("^[^:]*: *", "", v)) } else nm <- make.names(k)
    while (nm %in% names(out)) nm <- paste0(nm, ".x"); out[[nm]] <- v }
  n <- max(lengths(out)); as.data.frame(lapply(out, function(x){length(x)<-n; x}), stringsAsFactors = FALSE) }
md <- rbind(parse_sm("data/ISPY2_GPL20078_series.txt.gz")[,c("patient.id","arm","her2","pcr")],
            parse_sm("data/ISPY2_GPL30493_series.txt.gz")[,c("patient.id","arm","her2","pcr")])
n01 <- function(v) ifelse(grepl("pos|[+]|yes|^1$", v, TRUE), 1L, ifelse(grepl("neg|[-]|no|^0$", v, TRUE), 0L, NA))
md$.id <- as.character(md$patient.id); md$.her2 <- n01(md$her2); md$.pcr <- n01(md$pcr)
md$.arm2 <- ifelse(grepl("Pembrolizumab", md$arm), "pembro", ifelse(md$arm == "Paclitaxel", "control", NA))
ex <- as.data.frame(fread("data/ISPY2_geneExp_n988.txt.gz")); rownames(ex) <- toupper(ex[[1]]); ex[[1]] <- NULL
ex <- as.matrix(ex); colnames(ex) <- as.character(colnames(ex))
elig <- subset(md, .arm2 %in% c("pembro","control") & .her2 == 0 & !is.na(.pcr) & .id %in% colnames(ex))
elig <- elig[!duplicated(elig$.id), ]; E <- ex[, elig$.id]
y <- elig$.pcr; arm <- factor(elig$.arm2, c("control","pembro"))

prop_high <- function(genes) { g <- intersect(genes, rownames(E)); if (!length(g)) return(rep(NA, ncol(E)))
  b <- t(apply(E[g,,drop=FALSE], 1, function(v) as.integer(v >= median(v, na.rm=TRUE)))); colMeans(b, na.rm=TRUE) }
Z <- function(v) (v - mean(v, na.rm=TRUE)) / sd(v, na.rm=TRUE)
M <- sapply(MOD, prop_high)
Sense <- Z(M[,"Sense"]); Integ <- Z(M[,"Traffic"]-(M[,"Stromal_ECM"]+M[,"Angio_CAF"])/2); Exec <- Z(M[,"Execute"]-M[,"Checkpoint"])
minz <- apply(cbind(Sense,Integ,Exec), 1, min); intact <- minz >= -0.5

# signature score = mean of per-gene z across the program (genes present)
sig_z <- function(genes){ g <- intersect(genes, rownames(E)); if(!length(g)) return(rep(NA,ncol(E)))
  zz <- t(apply(E[g,,drop=FALSE],1,Z)); colMeans(zz, na.rm=TRUE) }
SIG <- sapply(ESCAPE, sig_z)

grp <- intact & arm=="pembro"
nr <- grp & y==0; rr <- grp & y==1
cat(sprintf("\ncascade-intact + pembro:  NR n=%d,  R n=%d\n", sum(nr), sum(rr)))
cat("\nprogram            meanZ_NR  meanZ_R   gap(NR-R)   Wilcoxon_p\n")
res <- data.frame()
for (p in colnames(SIG)) { v <- SIG[,p]
  a <- v[nr]; b <- v[rr]; pw <- tryCatch(wilcox.test(a,b)$p.value, error=function(e) NA)
  res <- rbind(res, data.frame(program=p, meanZ_NR=mean(a,na.rm=T), meanZ_R=mean(b,na.rm=T),
                               gap=mean(a,na.rm=T)-mean(b,na.rm=T), p=pw))
  cat(sprintf("%-16s   %+7.2f  %+7.2f    %+7.2f     %8.3f\n", p, mean(a,na.rm=T), mean(b,na.rm=T), mean(a,na.rm=T)-mean(b,na.rm=T), pw)) }
res$fdr <- p.adjust(res$p, "BH")
cat("\nFDR (BH) across the 5 programs:\n"); print(data.frame(program=res$program, gap=round(res$gap,2), p=round(res$p,3), fdr=round(res$fdr,3)), row.names=FALSE)
cat(sprintf("\nNOTE: exploratory, n=%d vs %d; a positive gap = program HIGHER in the unexplained NRs.\n", sum(nr), sum(rr)))
saveRDS(res, "escape_probe_intact_ICB_NR.rds")
