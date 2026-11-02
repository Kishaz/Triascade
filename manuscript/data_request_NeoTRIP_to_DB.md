# NeoTRIP clinical labels — request to DB

**To:** DB
**From:** SMM
**Re:** Two de-identified columns for the NeoTRIP transcriptomic cohort (GSE319641)

Hi DB,

For the Immune Cascade Index paper I need to run the treatment-interaction analysis on NeoTRIP,
exactly as we did for I-SPY2. The public GEO record (GSE319641) has the expression matrix and
patient IDs but **not** the treatment arm or the pathologic response, so the analysis is blocked
without the trial clinical file you hold locally
(`NeoTRIP_baseline_D1C2_pheno_and_clinical.txt`).

**What I need:** a single de-identified table keyed by the GEO patient ID, with just two columns:

| column | values | notes |
|---|---|---|
| `patient_id` | `TNBCxxxxx` | the same ID used in the GEO expression columns (baseline samples) |
| `arm` | `chemo` or `chemo+atezolizumab` | the randomized assignment |
| `pCR` | `1` (pCR) or `0` (residual disease) | pathologic complete response |

No ages, dates, sites, names, or any other field — those are not needed and should be left out so
the file is fully de-identified. Baseline samples only (the `_Baseline` columns); the `_D1C2`
on-treatment samples are not used for this analysis.

**Why just these two:** the analysis is a treatment-interaction test — whether the Integration
score predicts who benefits specifically from adding atezolizumab — so I need the arm and the
outcome per patient, matched to the expression we already have. Everything else is already in GEO.

A plain `.csv` or `.txt` is perfect. Once I have it the run takes minutes (the scoring is already
built and validated on NeoTRIP expression across the 241 baseline tumors).

Thanks,
SMM
