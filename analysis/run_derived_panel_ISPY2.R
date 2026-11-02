# =====================================================================
# MAKE-OR-BREAK CHECK
# Does the DERIVED panel (6 systematic-review modules nested into 3
# cascade tiers) reproduce the Integration x treatment predictive result
# that the legacy 98-gene panel gave in I-SPY2?
#
# Derived panel = ICI_PANEL_v1_locked.csv (108 genes, 6 modules).
# Tier mapping (the hybrid 6->3 nesting):
#   SENSE    = Sensory (IFN / antigen-presentation)
#   INTEGRATE(+) = Integration (chemokine trafficking)
#   INTEGRATE(-) = Integration (Treg/checkpoint) + Resistance (stromal/ECM)
#                  + Resistance (angiogenic/CAF-TGFb)      <- suppressive arm
#   EXECUTE(+)   = Execution (cytotoxic effector)
#   EXECUTE(-)   = any Execution-module gene coded 'suppresses'
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })
source("ici_locked_reference.R")

DER <- "../ICI-study-redesign/provenance-review-tight/derivation/ICI_PANEL_v1_locked.csv"
dp  <- fread(DER, data.table = FALSE)
dp$gene <- toupper(dp$gene)

sens  <- dp$gene[grepl("^Sensory", dp$layer)]
intpos<- dp$gene[dp$layer == "Integration (chemokine trafficking)"]
intneg<- dp$gene[dp$layer == "Integration (Treg / checkpoint regulation)" | grepl("^Resistance", dp$layer)]
execg <- dp$gene[grepl("^Execution", dp$layer)]
execpos <- execg[dp$direction[match(execg, dp$gene)] != "suppresses"]
execneg <- execg[dp$direction[match(execg, dp$gene)] == "suppresses"]

PANEL_DER <- list(SENSORY = sens, INTPOS = intpos, INTNEG = intneg,
                  EXECPOS = execpos, EXECNEG = execneg)
cat("DERIVED panel tier sizes:\n"); print(lengths(PANEL_DER))
cat(sprintf("distinct genes: %d\n\n", length(unique(unlist(PANEL_DER)))))

# legacy panel for side-by-side
PANEL_LEG <- load_panel("ICI_98gene_panel.csv")

## ---------- load I-SPY2 (same pipeline as run_stability.R) ----------
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
cat(sprintf("eligible n=%d (pembro %d, control %d); pCR pembro %.0f%% / control %.0f%%\n\n",
            ncol(E), sum(arm=="pembro"), sum(arm=="control"),
            100*mean(y[arm=="pembro"]), 100*mean(y[arm=="control"])))

# robust cohort-median scorer (handles an empty layer arm as zeros)
score_cm <- function(expr, PANEL) {
  binset <- function(gs){ g<-intersect(gs,rownames(expr)); if(!length(g)) return(NULL)
    b<-t(apply(expr[g,,drop=FALSE],1,function(v) as.integer(v>=median(v,na.rm=TRUE))))
    colnames(b)<-colnames(expr); b }
  cs <- function(b) if(is.null(b)) rep(0, ncol(expr)) else colSums(b, na.rm=TRUE)
  S<-cs(binset(PANEL$SENSORY)); I<-cs(binset(PANEL$INTPOS))-cs(binset(PANEL$INTNEG))
  Ex<-cs(binset(PANEL$EXECPOS))-cs(binset(PANEL$EXECNEG))
  zz<-function(x){s<-sd(x);if(is.na(s)||s==0) x*0 else (x-mean(x))/s}
  data.frame(sample=colnames(expr),Sensory=as.numeric(S),Integration=as.numeric(I),
             Execution=as.numeric(Ex),Composite=rowMeans(cbind(zz(S),zz(I),zz(Ex))),row.names=NULL)
}
z <- function(v) scale(v)[,1]
auc1 <- function(yy, s) as.numeric(auc(yy, s, quiet = TRUE, direction = "<"))
int_p <- function(sc) anova(glm(y~z(sc)+arm, binomial), glm(y~z(sc)*arm, binomial), test="LRT")[2,"Pr(>Chi)"]

report <- function(tag, PANEL) {
  cm <- score_cm(E, PANEL); cm <- cm[match(elig$.id, cm$sample), ]
  cat(sprintf("==== %s ====\n", tag))
  for (L in c("Sensory","Integration","Execution","Composite")) {
    s <- cm[[L]]
    cat(sprintf("  %-12s interaction p = %.4f | AUC pembro %.3f / control %.3f\n",
                L, int_p(s), auc1(y[arm=="pembro"], s[arm=="pembro"]),
                auc1(y[arm=="control"], s[arm=="control"])))
  }
  # pCR-by-Integration x arm (the headline table)
  hi <- cm$Integration >= median(cm$Integration)
  tb <- tapply(y, list(ifelse(hi,"Int-high","Int-low"), arm), mean)*100
  cat("  Integration pCR%:  control low/high ", sprintf("%.0f/%.0f", tb["Int-low","control"], tb["Int-high","control"]),
      " | pembro low/high ", sprintf("%.0f/%.0f", tb["Int-low","pembro"], tb["Int-high","pembro"]), "\n", sep="")
  cat(sprintf("  drug benefit (pembro-control): low %+.0f pts | high %+.0f pts\n\n",
              tb["Int-low","pembro"]-tb["Int-low","control"], tb["Int-high","pembro"]-tb["Int-high","control"]))
  invisible(cm)
}

report("LEGACY 98-gene panel", PANEL_LEG)
report("DERIVED tier-map A: Integration = trafficking(+) vs [Treg/checkpoint + ALL resistance](-)", PANEL_DER)

# ---- Variant B (more standard, pre-specified on biology): Integration =
#      chemokine-trafficking(+) vs Treg/checkpoint(-) ONLY; resistance held
#      OUT as its own axis (a resistance modifier), not folded into Integration.
intneg_B <- dp$gene[dp$layer == "Integration (Treg / checkpoint regulation)"]
resist   <- dp$gene[grepl("^Resistance", dp$layer)]
PANEL_DER_B <- list(SENSORY = sens, INTPOS = intpos, INTNEG = intneg_B,
                    EXECPOS = execpos, EXECNEG = execneg)
cat("Variant B tier sizes (resistance excluded from Integration):\n"); print(lengths(PANEL_DER_B))
report("DERIVED tier-map B: Integration = trafficking(+) vs Treg/checkpoint(-) only", PANEL_DER_B)

# resistance as its own axis (high = more stromal/angiogenic suppression)
rb <- score_cm(E, list(SENSORY=resist, INTPOS=character(0), INTNEG=character(0),
                       EXECPOS=character(0), EXECNEG=character(0)))
rb <- rb[match(elig$.id, rb$sample),]
cat(sprintf("RESISTANCE axis (%d genes): interaction p = %.4f | AUC pembro %.3f / control %.3f\n\n",
            length(resist), int_p(rb$Sensory),
            auc1(y[arm=="pembro"], rb$Sensory[arm=="pembro"]), auc1(y[arm=="control"], rb$Sensory[arm=="control"])))

# correlation of the two Integration scores on the same patients
cmL <- score_cm(E, PANEL_LEG); cmL <- cmL[match(elig$.id, cmL$sample),]
cmD <- score_cm(E, PANEL_DER); cmD <- cmD[match(elig$.id, cmD$sample),]
cat(sprintf("Integration score Spearman (legacy vs derived) = %.3f\n",
            cor(cmL$Integration, cmD$Integration, method="spearman")))
cat(sprintf("Composite   score Spearman (legacy vs derived) = %.3f\n",
            cor(cmL$Composite,   cmD$Composite,   method="spearman")))
