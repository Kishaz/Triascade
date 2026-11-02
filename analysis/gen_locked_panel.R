# =====================================================================
# LOCK the ICI panel (Phase 1). Canonicalise ICI_PANEL_v1_locked.csv
# (108 genes, Grade A+B) into the frozen 6-module -> 3-layer definition
# with per-gene layer + direction, module GO read, recurrence and grade.
# Output: ICI_panel_LOCKED_v1.csv (source of truth for all scoring).
# =====================================================================
suppressPackageStartupMessages(library(data.table))
src <- "../ICI-study-redesign/provenance-review-tight/derivation/ICI_PANEL_v1_locked.csv"
p <- fread(src, data.table=FALSE); p$gene <- toupper(p$gene)

# module -> (layer, direction) lock map  + enrichment read (g:Profiler GO/Reactome)
MAP <- list(
 "Sensory (IFN / antigen-presentation)"        = c("Sense","positive","response to biotic stimulus / defense response (GO p=6e-59)"),
 "Integration (chemokine trafficking)"          = c("Integration","positive","inflammatory response / response to chemical (GO p=4e-46)"),
 "Integration (Treg / checkpoint regulation)"   = c("Integration","negative","leukocyte / cell activation; regulatory (GO p=1e-28)"),
 "Execution (cytotoxic effector)"               = c("Execution","positive","immune response / lymphocyte-mediated cytotoxicity (GO p=1e-64)"),
 "Resistance (stromal / ECM-coagulation)"       = c("Integration","negative","proteolysis / coagulation / ECM organization (GO p=1e-33)"),
 "Resistance (angiogenic / CAF-TGFb)"           = c("Integration","negative","vasculature development / cell migration (GO p=2e-25)"))

# Module numbers run top-to-bottom in cascade order (Sense -> Integration + ->
# Integration - -> Execution) for consistent flow with Figure 1.
mod_short <- c(
 "Sensory (IFN / antigen-presentation)"="M1_IFN_antigen",
 "Integration (chemokine trafficking)"="M2_chemokine_trafficking",
 "Integration (Treg / checkpoint regulation)"="M3_Treg_checkpoint",
 "Resistance (stromal / ECM-coagulation)"="M4_stromal_ECM",
 "Resistance (angiogenic / CAF-TGFb)"="M5_angiogenic_CAF_TGFb",
 "Execution (cytotoxic effector)"="M6_cytotoxic_effector")

out <- data.frame(
  gene        = p$gene,
  module       = mod_short[p$layer],
  module_label = p$layer,
  layer        = sapply(p$layer, function(L) MAP[[L]][1]),
  direction    = sapply(p$layer, function(L) MAP[[L]][2]),
  recurrence   = p$recurrence,
  grade        = p$grade,
  legacy_status= p$legacy_status,
  module_enrichment = sapply(p$layer, function(L) MAP[[L]][3]),
  stringsAsFactors = FALSE, row.names=NULL)
out <- out[order(factor(out$layer,c("Sense","Integration","Execution")),
                 out$direction!="positive", out$module, -out$recurrence), ]

dest1 <- "ICI_panel_LOCKED_v1.csv"
dest2 <- "../ICI-study-redesign/provenance-review-tight/derivation/ICI_panel_LOCKED_v1.csv"
write.csv(out, dest1, row.names=FALSE); write.csv(out, dest2, row.names=FALSE)

cat("LOCKED PANEL  ICI_panel_LOCKED_v1.csv  (", nrow(out), "genes )\n\n")
cat("by layer x direction:\n"); print(table(out$layer, out$direction))
cat("\n3-layer scoring structure:\n")
cat(sprintf("  Sense       = %d genes (all positive)\n", sum(out$layer=="Sense")))
cat(sprintf("  Integration = %d (%d activating - %d suppressive)\n",
  sum(out$layer=="Integration"), sum(out$layer=="Integration"&out$direction=="positive"),
  sum(out$layer=="Integration"&out$direction=="negative")))
cat(sprintf("  Execution   = %d genes (all positive)\n", sum(out$layer=="Execution")))
cat(sprintf("  distinct genes: %d | grades A/B: %s\n", length(unique(out$gene)),
            paste(names(table(out$grade)),table(out$grade),sep="=",collapse=" ")))
cat("\nmodule sizes + enrichment read:\n")
for(m in unique(out$module)){ s<-out[out$module==m,][1,]
  cat(sprintf("  %-26s n=%2d  %-11s %-8s  %s\n", m, sum(out$module==m), s$layer, s$direction, s$module_enrichment)) }
