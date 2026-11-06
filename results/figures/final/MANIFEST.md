# Canonical figure set — manifest

Non-destructive submission set. `assemble_final.sh` copies the current editable sources
to these canonical names; working stems, `_c` Inkscape edits, and build scripts are
untouched. Re-run the script after regenerating any source. Numbering follows
`manuscript/FIGURE_RECONCILIATION_v1.md` (decisions D1–D5, 2026-11-06).

## Main figures

| Canonical | Source asset | Builder | In set? |
|---|---|---|---|
| **Fig1** | `composite/edits/ICI_Fig1_full` (hand-edited) | `make_fig1_full.R` | yes |
| **Fig2** | `composite/ICI_Fig2` (2C panel trimmed) | `make_fig2_combined.R` | yes |
| **Fig3** | `composite/ICI_Fig3_composite` (a flow banner + b atlas + c violins + d dose-response + e alluvial) | `make_fig3_composite.R` (PNG) + `compose_svgs.py` (editable SVG) | yes (SVG + PNG) |
| **Fig4** | `composite/ICI_Fig4_composite` (a ranking + b heatmap + c grids) | `make_fig4_composite.R` (PNG) + `compose_svgs.py` (editable SVG) | yes (SVG + PNG) |

Panel letters are **lowercase** (Nature style: **a**, **b**, **c**…), matching the manuscript
legends. `Fig3.svg` / `Fig4.svg` are composed from the panel SVGs by `compose_svgs.py` and keep
editable text (only the atlas anatomy icons are raster). Two assets still carry **uppercase**
letters and must be fixed by hand in Inkscape: the hand-edited `composite/edits/ICI_Fig1_full.svg`
(a/b/c) and the `_c` Fig 3 panels (`ICI_Fig3_atlas_c` → b, `ICI_Fig3b_separation_c` → c,
`ICI_Fig3c_doseresponse_c` → d).
| **Fig5** | `composite/ICI_Fig_predictive` | `make_fig_predictive.R` | yes |
| **Fig6** | `composite/ICI_Fig_cascade_biology` | `make_fig_cascade_biology.R` | yes |
| **Fig7** | `composite/ICI_Fig_archetype_biology` | `make_fig_archetype_biology.R` | yes |

## Supplementary figures

| Canonical | Source asset(s) | Builder | In set? |
|---|---|---|---|
| **FigS1** | `composite/edits/ICI_SuppFig_coexpression_modules` (network + module heatmap, hand-assembled) | `make_supp_coexpression.R` + Inkscape | yes |
| **FigS2** | `composite/edits/ICI_SuppFig_prisma` (study-selection flow, hand-assembled) | PRISMA flow (Inkscape) | yes |
| **FigS3** | `composite/ICI_SuppFig_scoring_flow` | `make_supp_scoring_flow.R` | yes |
| **FigS4** | `composite/ICI_Fig3b_auc` (per-cohort forest) | `make_fig3b_auc.R` | yes |
| **FigS5** | merge `ICI_SuppFig_separation_sensory` + `ICI_SuppFig_separation_execution` | `make_fig3b_separation.R` + Inkscape | **pending (merge)** |
| **FigS6** | merge `ICI_Fig_addedvalue` + `ICI_Fig_addedvalue_rates` | `make_fig_addedvalue*.R` + Inkscape | **pending (merge)** |
| **FigS7** | `supplementary/ICI_SuppFig_stability` | `make_supp_stability.R` | yes |
| **FigS8** | `composite/ICI_Fig_io_landscape` | `make_fig_io_landscape.R` | yes |
| **FigS9** | `supplementary/ICI_SuppFig_cascadescope_tool` | `make_supp_cascadescope_tool.R` | yes |

## Still to do before submission
- **Fig 3, Fig 4** — assemble editable-SVG composites in Inkscape (your step, per rule #1),
  then add `copy` lines for them to `assemble_final.sh`.
- **FigS5, FigS6** — merge the listed multi-panel components into one figure each.
- **Manuscript cross-refs** — prose cites "Supplementary Fig. 1" for the scoring-invariance
  property (lines 27, 131, 223); S1 is now co-expression, so those become **Supplementary Fig. 3**
  (scoring flow). Left for the writing pass — see below.
- **Prose fixes** (from the value-check, see `FIGURE_RECONCILIATION_v1.md` V2–V6):
  OR/SD caption 1.82 vs prose 1.81; GeparNuevo "Fig. 5c" should be 5a; Execution P 0.705→;
  decide whether to state the "91% independent" and cross-context AUCs in Results.
