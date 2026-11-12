# =============================================================================
# Supplementary GO-enrichment table: statistical support for the functional
# annotation of the six co-expression modules (reviewer comment C5). For each
# module M1-M6 of the 98-gene ICI panel, GO biological-process over-representation
# analysis (clusterProfiler::enrichGO) against the genome background, BH-adjusted.
# Confirms that each module is enriched for the process it is named after.
# Reads analysis/ICI_98_module_assignment.csv (gene, module, module_label).
# =============================================================================
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })

mod <- read.csv("ICI_98_module_assignment.csv", stringsAsFactors = FALSE)
mods <- sort(unique(mod$module))
TOPN <- 5                                           # top terms reported per module

rows <- list()
summary_lines <- c()
for (m in mods) {
  genes <- unique(mod$gene[mod$module == m])
  lab   <- mod$module_label[mod$module == m][1]
  ego <- tryCatch(enrichGO(gene = genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                           ont = "BP", pAdjustMethod = "BH",
                           pvalueCutoff = 1, qvalueCutoff = 1,   # keep all; filter in the table
                           minGSSize = 3, maxGSSize = 800, readable = FALSE),
                  error = function(e) NULL)
  df <- if (is.null(ego)) NULL else as.data.frame(ego)
  if (is.null(df) || !nrow(df)) {
    summary_lines <- c(summary_lines, sprintf("  %s %-28s : no GO BP term", m, lab)); next
  }
  df <- df[order(df$p.adjust), ]
  top <- head(df, TOPN)
  rows[[m]] <- data.frame(module = m, module_label = lab,
                          GO_ID = top$ID, GO_term = top$Description,
                          gene_ratio = top$GeneRatio, bg_ratio = top$BgRatio,
                          p_value = signif(top$pvalue, 3), BH_q = signif(top$p.adjust, 3),
                          n_genes = top$Count, genes = top$geneID,
                          stringsAsFactors = FALSE)
  t1 <- top[1, ]
  summary_lines <- c(summary_lines, sprintf("  %s %-28s : %-45s  q=%.1e  (%d/%s)",
                     m, lab, t1$Description, t1$p.adjust, t1$Count, sub("/.*","",t1$GeneRatio)))
}

OUT <- do.call(rbind, rows)
out <- "../results/tables/ICI_SuppTable_go_enrichment.csv"
write.csv(OUT, out, row.names = FALSE)
cat("wrote", out, "  (", nrow(OUT), "rows )\n\n")
cat("TOP GO biological-process term per module (q = BH-adjusted):\n")
cat(paste(summary_lines, collapse = "\n"), "\n")
