#!/usr/bin/env bash
# Assemble the canonical submission figure set (non-destructive): copy the current
# editable sources to Fig#/FigS# names. Working stems, the _c Inkscape edits, and the
# build scripts are left untouched. Re-run after regenerating any source figure.
# Figures needing Inkscape assembly (Fig 3, Fig 4, FigS1/S4/S5 multi-panel) are NOT
# copied here; see MANIFEST.md.
set -euo pipefail
cd "$(dirname "$0")"              # results/figures/final
C=../composite
S=../supplementary

copy() {  # copy <src-stem-no-ext> <dest-stem>
  local src="$1" dst="$2" ok=0
  for ext in svg png; do
    if [[ -f "$src.$ext" ]]; then cp "$src.$ext" "$dst.$ext"; ok=1; fi
  done
  [[ $ok -eq 1 ]] && echo "  $dst  <-  $src" || echo "  MISSING: $src (.svg/.png)"
}

echo "Main figures:"
copy "$C/edits/ICI_Fig1_full"          Fig1          # hand-edited composite
copy "$C/ICI_Fig2"                     Fig2          # 2C trimmed
copy "$C/ICI_Fig3_composite"           Fig3          # raster preview (A banner + B-E); editable assembly in Inkscape
copy "$C/ICI_Fig4_composite"           Fig4          # raster preview; editable assembly in Inkscape
copy "$C/ICI_Fig_predictive"           Fig5
copy "$C/ICI_Fig_cascade_biology"      Fig6
copy "$C/ICI_Fig_archetype_biology"    Fig7

echo "Supplementary figures:"
copy "$C/edits/ICI_SuppFig_coexpression_modules" FigS1   # co-expression network + module heatmap (hand-assembled)
copy "$C/edits/ICI_SuppFig_prisma"               FigS2   # PRISMA study-selection flow (hand-assembled)
copy "$C/ICI_SuppFig_scoring_flow"               FigS3
copy "$C/ICI_Fig3b_auc"                          FigS4
# FigS5: merge separation sensory+execution; FigS6: merge added-value bars+rates
copy "$S/ICI_SuppFig_stability"                  FigS7
copy "$C/ICI_Fig_io_landscape"                   FigS8
copy "$S/ICI_SuppFig_cascadescope_tool"          FigS9
echo "done"
