# Build shipped assets for the ICI scoring app (run once, from project root):
#   Rscript app/make_app_assets.R
#   -> app/ici_reference.rds   (frozen locked reference + panel)
#   -> app/example_cohort.csv  (a realistic example upload)
suppressPackageStartupMessages({ library(data.table) })
source("analysis/ici_locked_reference.R"); source("app/score_cohort.R")
PANEL <- load_panel("analysis/ICI_98gene_panel.csv")
PG <- unique(toupper(unlist(PANEL)))

## 1. freeze the locked reference (NeoTRIP breast baseline = validated external ref) ----
nt <- as.data.frame(fread("analysis/data/NeoTRIP_expr.txt.gz"))
rownames(nt) <- toupper(nt$HGNC_approved_symbol)
refmat <- as.matrix(nt[, grep("_Baseline$", colnames(nt)), drop=FALSE])
REF <- fit_locked_reference(refmat, PANEL, "ICI_locked_v1")

## alias / Ensembl -> canonical-symbol map for the panel genes (org.Hs.eg.db) ----
suppressMessages(library(org.Hs.eg.db))
eg <- AnnotationDbi::mapIds(org.Hs.eg.db, PG, "ENTREZID", "SYMBOL"); eg <- eg[!is.na(eg)]
al <- AnnotationDbi::mapIds(org.Hs.eg.db, eg, "ALIAS",   "ENTREZID", multiVals="list")
en <- AnnotationDbi::mapIds(org.Hs.eg.db, eg, "ENSEMBL", "ENTREZID", multiVals="list")
mp <- do.call(rbind, lapply(names(eg), function(sym){ e<-eg[[sym]]
  syns <- unique(toupper(c(sym, unlist(al[[e]]), unlist(en[[e]])))); syns <- syns[!is.na(syns) & syns!=""]
  data.frame(syn=syns, canon=sym, row.names=NULL) }))
mp <- mp[!(mp$syn %in% names(which(table(mp$syn) > 1))), ]   # drop ambiguous synonyms
alias_map <- setNames(mp$canon, mp$syn)
cat(sprintf("alias/Ensembl map: %d synonyms -> %d panel genes\n", length(alias_map), length(unique(alias_map))))

## reference sub-module score distribution (for the per-patient drill-down) -----
SMref <- submod_scores(refmat, REF)
submod_stats <- data.frame(submodule=colnames(SMref), mean=colMeans(SMref, na.rm=TRUE),
                           sd=apply(SMref, 2, sd, na.rm=TRUE), row.names=NULL)
cat("submodule stats:", paste(submod_stats$submodule, collapse=", "), "\n")

dir.create("app/references", showWarnings=FALSE, recursive=TRUE)
saveRDS(list(ref=REF, panel=PANEL, panel_genes=PG, alias_map=alias_map, submod_stats=submod_stats,
             version="ICI panel v1", panel_version="v1", n_genes=length(PG),
             layers=names(PANEL), n_ref=ncol(refmat),
             ref_label="NeoTRIP breast baseline (n=241)", built=Sys.Date()),
        "app/references/ICI_v1.rds")
cat(sprintf("froze reference 'ICI panel v1': %d panel genes, %d reference samples -> app/references/ICI_v1.rds\n",
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
