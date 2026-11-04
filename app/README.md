# ICI scorer — Shiny app

Score tumour RNA-seq against a locked **Immune Cascade Index** reference, one
patient at a time (platform-agnostic, cohort-independent), and read it out at
patient or cohort level. Runs locally — your data never leaves the machine.

## Run

```r
shiny::runApp("app")   # from the project root
```

Requires: `shiny`, `bslib`, `plotly`, `DT`, `ggplot2` (+ `data.table` to rebuild assets).

## Two views

- **Patient report** (clinician): a plain-language call, a 3-layer percentile
  gauge with the rate-limiting layer flagged, the limiting-layer biological
  rationale (investigational), a version stamp, and a research-use-only footer.
- **Cohort explorer** (trial / biomarker teams): distributions of the three
  layer scores and a per-patient beeswarm (hover for detail) + table + CSV.

## Input

RNA-seq expression — **TPM, log-TPM or counts** (linear data is auto-log2'd),
from whole-transcriptome or FFPE-compatible targeted RNA-seq. Provide a
**gene × sample** CSV: first column = gene symbol, remaining columns = samples.
A single patient is simply a one-sample file. Coverage of the 98 panel genes is
reported and flagged if low.

## Future-proof by design

References live in **`app/references/*.rds`**, each self-describing (version,
gene count, reference cohort). The app lists them in the **Panel / reference**
selector; drop in a new one — an expanded gene panel (v2…), or a reference built
on a newer platform — and it appears automatically, with **no code change**.
Every report is stamped with the panel version used.

## Files

- `app.R` — UI + server (two tabs, reference registry).
- `ici_locked_reference.R` — scoring engine (frozen-reference single-patient scoring).
- `score_cohort.R` — preprocessing + `score_uploaded()`.
- `references/ICI_v1.rds` — frozen locked reference (per-gene distributions + panel + metadata).
- `example_cohort.csv` — example upload (IMvigor210 subset).
- `make_app_assets.R` — rebuild the reference(s) + example (run once from project root).

## To extend later

More input adapters (Excel, counts matrices, other RNA assays), richer gene-ID
mapping (alias/Ensembl), a batch/QC panel, PDF export of the patient report, and
the per-patient cascade-profile drill-down.
