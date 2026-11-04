# Build shipped assets for the ICI scoring app (run once, from project root):
#   Rscript app/make_app_assets.R
#   -> app/ici_reference.rds   (frozen locked reference + panel)
#   -> app/example_cohort.csv  (a realistic example upload)
suppressPackageStartupMessages({ library(data.table) })
source("analysis/ici_locked_reference.R")
PANEL <- load_panel("analysis/ICI_98gene_panel.csv")
PG <- unique(toupper(unlist(PANEL)))

## 1. freeze the locked reference (NeoTRIP breast baseline = validated external ref) ----
nt <- as.data.frame(fread("analysis/data/NeoTRIP_expr.txt.gz"))
rownames(nt) <- toupper(nt$HGNC_approved_symbol)
refmat <- as.matrix(nt[, grep("_Baseline$", colnames(nt)), drop=FALSE])
REF <- fit_locked_reference(refmat, PANEL, "ICI_locked_v1")
saveRDS(list(ref=REF, panel=PANEL, panel_genes=PG,
             n_ref=ncol(refmat), ref_label="NeoTRIP breast baseline (n=241)",
             built=Sys.Date()), "app/ici_reference.rds")
cat(sprintf("froze reference: %d panel genes, %d reference samples -> app/ici_reference.rds\n",
            length(PG), ncol(refmat)))

## 2. example cohort (IMvigor210 subset: all panel genes + random extras, realistic upload) ----
obj <- readRDS("analysis/data/imvigor210/imvigor210_clean.rds"); E <- obj$E
rownames(E) <- toupper(rownames(E))
set.seed(1)
samp <- sample(colnames(E), 60)
extra <- sample(setdiff(rownames(E), PG), 400)
genes <- intersect(c(PG, extra), rownames(E))
ex <- E[genes, samp, drop=FALSE]
out <- data.frame(gene=rownames(ex), round(ex, 3), check.names=FALSE)
fwrite(out, "app/example_cohort.csv")
cat(sprintf("example cohort: %d genes x %d samples -> app/example_cohort.csv\n", nrow(out), ncol(ex)))
