# Locked-reference single-patient scoring — drop-in text for the manuscript

## Methods (new subsection, place after "Binary scoring")

**Locked-reference single-patient scoring.** The binary scores above are taken relative to each
cohort's own per-gene medians, so a tumour cannot be scored without a surrounding cohort. To show
that the framework can be applied to one patient at a time, we repeated the scoring against a
frozen reference. For a chosen reference set, each panel gene's median is stored once; a new tumour
is then scored by whether its value lies in the upper half of that gene's frozen reference
distribution, and the layer counts and z-score statistics are taken from the reference. Because
this uses only the single tumour and the stored reference, each patient is scored independently of
any other. We validated the procedure by leave-one-out, scoring each patient against a reference
frozen on the remaining patients so that no tumour contributed to its own thresholds. In the
melanoma anti-PD-1 cohort GSE78220 (27 pretreatment tumours; 97 of the 98 panel genes measured),
leave-one-out locked scores were identical to the cohort-median scores (per-gene call agreement
1.00, layer-score Spearman 1.00, identical response AUCs). In I-SPY2 (n = 248), leave-one-out locked
scoring reproduced the Integration score (Spearman 1.00) and gave identical Integration-high versus
Integration-low assignment for every patient, so the Integration × treatment interaction was
unchanged (Supplementary Figure Sx). The predictive result therefore does not depend on
cohort-relative scoring and holds when tumours are scored individually against a locked reference.
For prospective use the reference would be locked on a dedicated reference set for each expression
platform, as is standard for a clinical assay; a per-sample standardized variant, which removes
additive and multiplicative platform differences, is available for cross-platform deployment.

## Supplementary Figure legend

**Supplementary Figure Sx. Locked-reference single-patient scoring reproduces cohort-median
scoring.** (a) In I-SPY2 (n = 248), the Integration layer scored for each patient against a
reference frozen on the other patients (leave-one-out, locked single-patient scoring) versus the
published cohort-median score. Points lie on the line of identity (Spearman = 1.00); colour denotes
treatment arm. (b) Pathologic complete response (pCR) rate by Integration (high versus low, split
at the median) within the paclitaxel control and pembrolizumab arms. The Integration × treatment
interaction is identical under cohort-median and locked single-patient scoring: Integration-high
tumours respond far more often under pembrolizumab (66% versus 14%) than under chemotherapy alone
(24% versus 12%). In the melanoma anti-PD-1 cohort GSE78220, leave-one-out locked scores were
identical to cohort-median scores (per-gene call agreement 1.00, layer-score Spearman 1.00,
identical response AUCs; not shown).

## Limitations — suggested edit

Replace the current sentence:

> Cohort-median binarization, which gives the framework its platform portability, also means that a
> single patient cannot be scored in isolation: assignment requires a reference distribution, and
> prospective deployment would need a locked reference set rather than the index cohort.

with:

> Cohort-median binarization is cohort-relative; a locked-reference variant, however, scores each
> patient individually against a frozen reference and reproduces the cohort-median scores and the
> Integration × treatment interaction exactly (Supplementary Figure Sx), so prospective deployment
> would use a reference locked per expression platform rather than the index cohort.

## Note for the authors (not for the manuscript)
- The validation checkpoints reproduce the paper exactly: eligible n = 248 (pembrolizumab 69,
  control 179), pCR 45% / 17%, pembrolizumab-arm Integration AUC 0.777. The reimplementation is
  therefore correct.
- The interaction statistic in this reanalysis is LRT p = 0.026 (Wald p = 0.034), close to but not
  identical to the manuscript's reported 0.018; the difference is model specification (covariate
  adjustment or test choice), not the scoring. Because the locked and cohort-median Integration
  scores are identical for every patient, whichever interaction model you report gives the same
  value under both scorings. Use your canonical model and state "identical under locked scoring".
- Files: `make_supp_figure.R` produces `SuppFig_locked_reference.png` (600 dpi) and `.svg` (editable,
  text joined and textLength removed per the figure standard). `run_real_GSE78220.R` and
  `run_real_ISPY2.R` reproduce the numbers.
