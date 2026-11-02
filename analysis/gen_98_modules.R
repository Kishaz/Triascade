# =============================================================================
# Fit the official 98-gene panel into the 6 systematic-review modules, nested
# into the 3-layer cascade (S/I/E). Writes ICI_98_module_assignment.csv so the
# derivation figure's module/layer counts are read from data, never typed.
#   M1 IFN / antigen-presentation   -> Sense        (+)
#   M2 chemokine / trafficking      -> Integration  (+)
#   M3 stromal / ECM                -> Integration  (-)
#   M4 angiogenic / CAF-TGFb        -> Integration  (-)
#   M5 cytotoxic effector           -> Execution    (+)
#   M6 checkpoint / exhaustion      -> Execution    (-)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
leg <- fread("ICI_98gene_panel.csv"); leg$gene <- toupper(leg$gene)

STROMAL <- c("COL1A1","COL1A2","COL3A1","COL4A1","FAP","ACTA2","FN1","MMP2","MMP9","MMP14","CXCL12")

assign_module <- function(gene, layer){
  if (layer=="Sensory")               return(c("M1","IFN / antigen-presentation","Sense","positive"))
  if (layer=="Integration_positive")  return(c("M2","chemokine / trafficking","Integration","positive"))
  if (layer=="Integration_negative")  return(if (gene %in% STROMAL)
                                          c("M3","stromal / ECM","Integration","negative")
                                        else c("M4","angiogenic / CAF-TGFb","Integration","negative"))
  if (layer=="Execution_positive")    return(c("M5","cytotoxic effector","Execution","positive"))
  if (layer=="Execution_negative")    return(c("M6","checkpoint / exhaustion","Execution","negative"))
}
asg <- t(mapply(assign_module, leg$gene, leg$layer))
out <- data.frame(gene=leg$gene, module=asg[,1], module_label=asg[,2],
                  layer=asg[,3], direction=asg[,4], stringsAsFactors=FALSE)
out <- out[order(out$module, out$gene), ]
fwrite(out, "ICI_98_module_assignment.csv")

cat("=== module counts (rows) ===\n")
print(table(factor(out$module, paste0("M",1:6))))
cat(sprintf("\ndistinct genes: %d | rows: %d (IFNG, TGFB1 span two arms)\n",
            length(unique(out$gene)), nrow(out)))
cat("\n=== layer x direction ===\n"); print(table(out$layer, out$direction))
cat("\nwrote ICI_98_module_assignment.csv\n")
