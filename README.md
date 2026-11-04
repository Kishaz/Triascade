# Immune Cascade Index (ICI)

A transcriptomic, three-layer (**Sensory / Integration / Execution**) immune-checkpoint-blockade
(ICB) response biomarker built on the legacy 98-gene panel, being repackaged as a
single-patient-scorable, pre-treatment clinical tool.

## CascadeScope — the app

**[`app/`](app/) is the self-contained [CascadeScope](app/README.md) tool**: an R/Shiny app that
scores tumor RNA against the frozen ICI reference, one patient at a time, and reads it out
(patient report, 98-gene cascade view, cohort explorer). It runs locally and ships its own frozen
reference and example data — see **[`app/README.md`](app/README.md)** for install, run, and
reproducibility.

```bash
cd app && Rscript run.R      # installs deps if needed, launches CascadeScope
```

## Repository layout

```
ICI/
├── analysis/          R analysis & figure code (scoring, validation, figure builders)
│   ├── assets/icons/  Servier Medical Art / Bioicons artwork used in figures
│   ├── data/          cohort expression + clinical data  (NOT in git — see below)
│   ├── *.R            pipeline + figure scripts
│   ├── *.csv          panel definitions (ICI_98gene_panel.csv, module assignment, lock)
│   └── svg_join_text.py   makes saved SVGs hand-editable (one text object per label)
├── results/
│   ├── figures/       figure outputs — composite/ (full figs), panels/, edited/, supplementary/
│   └── tables/        result tables (CSV, incl. ICI_SuppTable_*)
├── manuscript/        CURRENT manuscript materials (study plan, panel lock, data requests)
├── review/            systematic-review / provenance tooling (separate sub-project)
└── archive/           superseded snapshots & old deliverables (Oct-29 final, Aug figures)
```

## Not under version control (local on disk only)
Set in `.gitignore`:
- `analysis/data/` — downloaded cohort data (~280 MB)
- `manuscript/`, `archive/` — manuscript drafts and superseded snapshots (kept locally)
- `review/` — systematic-review sub-project (its own git repo)
- `**/.venv/`, `*.pdf`, `*.zip`, `*.pptx`, `*.gmt` — environments and large binaries

## Reproducing the figures
From `analysis/`:
```bash
Rscript make_fig1_full.R      # → results/figures/composite/ICI_Fig1_full.{png,svg}
Rscript make_fig2_combined.R  # → results/figures/composite/ICI_Fig2.{png,svg}
Rscript make_fig3_atlas.R     # → results/figures/composite/ICI_Fig3_atlas.{png,svg}   (Fig 3a: anatomy atlas)
Rscript run_pancancer_extra.R # → analysis/pancancer_validation.rds   (16-cohort ICB validation; per-sample 3-layer scores)
Rscript run_pooled_auc.R      # → results/tables/pooled_integration_auc.csv, pancancer_cohort_table.csv
Rscript make_fig3b_auc.R      # → results/figures/composite/ICI_Fig3b_auc.{png,svg}  (forest, summary/supp)
Rscript make_fig3b_separation.R # → Fig 3b (Integration) + SUPPLEMENTARY Sensory/Execution separation figures
Rscript make_supp_tables.R    # → results/tables/ICI_SuppTable_layer_auc_by_cohort.csv (all 3 layers)
```

**Fig 3b** is the responder-vs-non-responder violin separation for the **Integration** layer
(`results/figures/composite/ICI_Fig3b_separation.{png,svg}`). The **same plot for the Sensory and
Execution layers** is written to `results/figures/supplementary/` — these show Integration
separates responders best (Sensory/Execution AUCs are markedly lower), the evidence that
Integration is the discriminating layer. All three share one cohort order (Integration-AUC sorted).

**Fig 3a** is a hand-tuned atlas: the committed SVG is the author's Inkscape-edited copy
(`results/figures/composite/edits/` holds the backup). Re-running `make_fig3_atlas.R`
reproduces it *except* a manual vertical-framing tweak that lives only in the edited copy —
so regenerate onto a new name, don't overwrite the edited deliverable.
`save_fig()` (in `analysis/ici_theme.R`) writes a 600-dpi PNG plus an editable SVG whose text is
post-processed by `svg_join_text.py` (sentences stay one object; no fixed `textLength`).

## Attribution
Organ illustrations in Figure 3 are **Servier Medical Art (smart.servier.com), CC BY 4.0** —
this credit must appear in the Figure 3 caption.
