# Immune Cascade Index (ICI)

A transcriptomic, three-layer (**Sensory / Integration / Execution**) immune-checkpoint-blockade
(ICB) response biomarker built on the legacy 98-gene panel, being repackaged as a
single-patient-scorable, pre-treatment clinical tool.

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
│   ├── figures/       figure outputs — composite/ (full figs), panels/, edited/
│   └── tables/        result tables (CSV)
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
Rscript run_pancancer_extra.R # → analysis/pancancer_validation.rds   (15-cohort ICB validation)
Rscript run_pooled_auc.R      # → results/tables/pooled_integration_auc.csv, pancancer_cohort_table.csv
Rscript make_fig3b_auc.R      # → results/figures/composite/ICI_Fig3b_auc.{png,svg}  (Fig 3b: AUC forest)
```

**Fig 3a** is a hand-tuned atlas: the committed SVG is the author's Inkscape-edited copy
(`results/figures/composite/edits/` holds the backup). Re-running `make_fig3_atlas.R`
reproduces it *except* a manual vertical-framing tweak that lives only in the edited copy —
so regenerate onto a new name, don't overwrite the edited deliverable.
`save_fig()` (in `analysis/ici_theme.R`) writes a 600-dpi PNG plus an editable SVG whose text is
post-processed by `svg_join_text.py` (sentences stay one object; no fixed `textLength`).

## Attribution
Organ illustrations in Figure 3 are **Servier Medical Art (smart.servier.com), CC BY 4.0** —
this credit must appear in the Figure 3 caption.
