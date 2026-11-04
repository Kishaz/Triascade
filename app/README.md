# ICI scorer — Shiny app (MVP)

Score any cohort's gene-expression table against the locked 98-gene **Immune
Cascade Index** reference, one patient at a time (platform-agnostic, cohort-
independent), and explore the three layer scores at cohort and patient level.

## Run

```r
# from the project root
shiny::runApp("app")
```

Requires: `shiny`, `bslib`, `plotly`, `DT`, `ggplot2` (and `data.table` for asset
building). Runs locally — your data never leaves the machine.

## Use

1. **Load example cohort** (ships with the app) or **Browse** to a CSV.
2. Input format: a **gene × sample** table, first column = gene symbol, remaining
   columns = samples (TPM / log-TPM / counts — linear data is auto-log2'd).
3. Read the **coverage** QC (how many of the 98 panel genes were found), the
   **cohort distribution** of the 3 layer scores, and the **per-patient**
   beeswarm (hover a point for that patient's scores, floor layer and archetype).
4. **Download scores** as CSV.

## Files

- `app.R` — the Shiny app (UI + server).
- `ici_locked_reference.R` — scoring engine (frozen-reference single-patient scoring).
- `score_cohort.R` — preprocessing + scoring wrapper (`score_uploaded()`).
- `ici_reference.rds` — the frozen locked reference (per-gene distributions + panel).
- `example_cohort.csv` — example upload (IMvigor210 subset).
- `make_app_assets.R` — rebuilds `ici_reference.rds` + `example_cohort.csv` (run once from project root).

## MVP scope / to extend later

Single CSV format, one locked reference, three views. Planned: more input formats
(Excel, counts matrices, NanoString), richer gene-ID mapping, batch/QC panel,
per-patient cascade-profile drill-down, and the biological "floor explanation".
