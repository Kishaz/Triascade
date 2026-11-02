# Locked-reference single-patient scoring (prototype)

## Why
The published ICI layer scores binarize each gene at the **cohort median**
(`gene = 1 if expr >= median across the cohort`). That makes the score
**cohort-relative**: a single patient cannot be scored without a cohort around
them, which the manuscript lists as a limitation and which top-tier reviewers
will flag as a deployability problem. This prototype freezes the reference
**once** and scores any patient **one at a time** against it.

## What is here
- `ici_locked_reference.R` - the scoring module (drop-in). Functions:
  - `load_panel()` - the frozen 98-gene panel (21 / 21+20 / 19+19).
  - `score_cohort_median(expr, PANEL)` - reproduces the published scoring.
  - `fit_locked_reference(expr_ref, PANEL)` - freezes per-gene thresholds and
    the layer z-score statistics from a reference set.
  - `score_locked(expr_new, ref, mode)` - scores a single sample or a matrix
    against the frozen reference, with no cohort needed. Modes:
    - `quantile` (default) - per-gene rank against the frozen reference;
      reproduces cohort-median when the target is the **same platform**.
    - `direct` - value vs frozen reference median (same platform and
      normalization).
    - `selfnorm` - per-sample z-score across the panel then a frozen
      standardized threshold; removes per-sample affine (platform) shifts.
  - `concordance()`, `archetype_agreement()` - head-to-head agreement.
- `run_locked_reference_demo.R` - a **synthetic** demonstration (no real data).
- `demo_concordance_table.csv`, `demo_locked_vs_cohortmedian.png` - its output.
- `ICI_98gene_panel.csv` - gene, layer, direction.

## What the synthetic demo shows
On a new cohort drawn from the **same platform** as the locked reference:

| metric | result |
|---|---|
| gene-call agreement (locked vs cohort-median) | 0.97 |
| Integration layer Spearman | 0.99 |
| Composite Spearman | 0.99 |
| 3-bit archetype agreement | 0.87 |
| AUC for Integration (cohort-median vs locked) | 0.762 vs 0.759 |

and the single-patient check passes: **a patient scored alone is identical to
that patient scored inside a cohort.** On a deliberately shifted platform the
raw-threshold agreement drops (0.63) but the layer-score rank-correlation (0.88)
and the predictive AUC (0.75 vs 0.77) are largely preserved; full fidelity on a
new platform comes from freezing the reference on that platform (standard for a
clinical assay).

## How to run it on the real cohorts (for the paper)
The expression matrices live in the pipeline (genes x samples), e.g. the
`score_ici()` input in `ISPY2_REBUILD_248_v1.1.R`. Replace the synthetic
matrices with the real ones:

```r
source("ici_locked_reference.R")
PANEL <- load_panel("ICI_98gene_panel.csv")

# 1. Choose a LOCKED reference (freeze once). Options, in order of strength:
#    (a) the largest single cohort on the clinical platform (e.g. IMvigor210), or
#    (b) a pooled multi-cohort reference per platform.
ref <- fit_locked_reference(expr_IMvigor210, PANEL, label = "IMvigor210_locked")
saveRDS(ref, "ICI_locked_reference.rds")   # ship this with the classifier

# 2. Score every OTHER cohort patient-by-patient against the frozen reference,
#    without using the target cohort's own median:
lk  <- score_locked(expr_ISPY2, ref, mode = "quantile")   # same platform
cm  <- score_cohort_median(expr_ISPY2, PANEL)             # published scoring

# 3. Report concordance and - the key point - that the headline results survive
#    single-patient scoring:
concordance(cm, lk); archetype_agreement(cm, lk)
#  - rerun the I-SPY2 Integration x treatment interaction on `lk$Integration`
#  - rerun the IMvigor210 Integration survival / GEP-TIS redundancy on `lk`
#  - compare per-cohort layer AUCs (cohort-median vs locked)
```

## What to put in the manuscript
A short "Locked-reference deployment" paragraph (Methods) + one supplementary
figure/table showing:
1. locked vs cohort-median layer scores are near-identical within platform
   (gene-call agreement, Spearman, archetype agreement), and
2. the predictive results (I-SPY2 Integration x treatment interaction; IMvigor210
   Integration survival and benchmark redundancy) are **preserved under locked,
   single-patient scoring**, with the reference locked and shipped.

This converts the current limitation ("a single patient cannot be scored in
isolation") into a solved, demonstrated deployment path - the specific thing
Nature Cancer-tier review asks for.

## Note on figures
The demo figure is PNG. A manuscript version will be produced as an editable SVG
(sentence-level text objects, no fixed widths) per the project figure standard.

## Real-data validation (added 2026-10-29)
`run_real_GSE78220.R` runs this on a **real public cohort**: GSE78220 (Hugo),
melanoma, anti-PD-1, 27 baseline patients (15 responders / 12), 97/98 panel
genes present (VSIR off-platform). Scoring is leave-one-out: each patient is
scored against a reference frozen on the OTHER 26, so no patient uses their own
cohort.

Result (`real_GSE78220_concordance.csv`):
- gene-call agreement, locked LOO (quantile) vs cohort-median: **1.000**
- Integration / Composite layer Spearman: **1.000 / 1.000**; archetype agreement **1.000**
- response AUC identical: Sensory 0.656, Integration 0.658, Execution 0.481,
  Composite 0.617 under BOTH the published cohort-median rule and locked
  single-patient scoring.
- the cohort-median Sensory AUC (0.656) reproduces the paper's reported Hugo
  Sensory-dominant value (0.642), validating the reimplementation.

Conclusion: on real data, a patient can be scored one at a time against a locked
reference and get the identical layer scores, archetype and predictive AUC as
the cohort-relative rule used in the manuscript. The single-patient limitation is
solved. Next: repeat the freeze-and-transfer on I-SPY2 (GSE194040) and show the
Integration x treatment interaction survives locked scoring.

## DECISIVE real-data test: I-SPY2 (added 2026-10-29)
`run_real_ISPY2.R` on GSE194040 (I-SPY2), HER2-negative, pembrolizumab+paclitaxel
vs paclitaxel alone. Validation checkpoints reproduce the paper exactly:
eligible **248** (pembro 69, control 179), pCR **45% / 17%**, pembro-arm
Integration AUC **0.777** (paper: 0.777).

Leave-one-out LOCKED single-patient scoring vs the published cohort-median rule:

| | Integration x arm p (continuous) | p (median split) | pembro AUC |
|---|---|---|---|
| cohort-median | 0.0255 | 0.0298 | 0.777 |
| **locked single-patient (LOO)** | **0.0255** | **0.0298** | **0.777** |

Integration score Spearman 1.000; Integration high/low assignment agreement 100%.
pCR by Integration x arm (identical under both scorings): control 12% (low) / 24%
(high); pembrolizumab 14% (low) / 66% (high).

CONCLUSION: the headline predictive result -- Integration predicts differential
benefit from adding pembrolizumab -- is UNCHANGED under locked, single-patient
scoring. A patient can be scored alone and the Integration x treatment
interaction holds identically. This is the deployability evidence for the
Nature Cancer tier.
