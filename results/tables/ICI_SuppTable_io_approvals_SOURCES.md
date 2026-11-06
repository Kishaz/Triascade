# Supplementary Table SX — approved checkpoint-inhibitor regimens: data sources

Companion to `ICI_SuppTable_io_approvals_detailed.csv` (per-indication detail) and
`ICI_io_regimen_counts_by_cancer.csv` (collapsed regimen tally used in the figure).
Compiled 2026-11-06.

## Data-source statement (drop-in for Methods)

FDA-approved immune checkpoint-inhibitor regimens were enumerated from primary FDA
prescribing information via Drugs@FDA (openFDA API) and the FDA Oncology Center of
Excellence approval notifications, and cross-checked against the Cancer Research
Institute immuno-oncology approval landscape. European approvals were verified against
the European Medicines Agency (EMA) European Public Assessment Reports (EPARs). Each
indication is linked to its pivotal registrational trial. Indication-level counts were
collapsed to distinct drug + combination-partner regimens (monotherapy, +chemotherapy,
+CTLA-4, +LAG-3, or +targeted agent), pooling line and treatment setting. Withdrawn
approvals were excluded from active counts. Counts reflect status as of 6 Nov 2026.

## Regulatory references (per drug)

FDA (all): Drugs@FDA — https://www.accessdata.fda.gov/scripts/cder/daf/ (search by brand)
and FDA OCE approval notifications.

EMA European Public Assessment Reports:

| Drug (INN) | Brand | EMA EPAR |
|---|---|---|
| pembrolizumab | Keytruda | https://www.ema.europa.eu/en/medicines/human/EPAR/keytruda |
| nivolumab | Opdivo | https://www.ema.europa.eu/en/medicines/human/EPAR/opdivo |
| nivolumab + relatlimab | Opdualag | https://www.ema.europa.eu/en/medicines/human/EPAR/opdualag |
| atezolizumab | Tecentriq | https://www.ema.europa.eu/en/medicines/human/EPAR/tecentriq |
| durvalumab | Imfinzi | https://www.ema.europa.eu/en/medicines/human/EPAR/imfinzi |
| cemiplimab | Libtayo | https://www.ema.europa.eu/en/medicines/human/EPAR/libtayo |
| avelumab | Bavencio | https://www.ema.europa.eu/en/medicines/human/EPAR/bavencio |
| ipilimumab | Yervoy | https://www.ema.europa.eu/en/medicines/human/EPAR/yervoy |
| tislelizumab | Tevimbra | https://www.ema.europa.eu/en/medicines/human/EPAR/tevimbra |
| toripalimab | Loqtorzi | https://www.ema.europa.eu/en/medicines/human/EPAR/loqtorzi |
| dostarlimab | Jemperli | https://www.ema.europa.eu/en/medicines/human/EPAR/jemperli |
| retifanlimab | Zynyz | https://www.ema.europa.eu/en/medicines/human/EPAR/zynyz-0 |
| cosibelimab | Unloxcyt | No EMA authorisation (US only) |

Corroborating secondary source: Vaddepally RK, et al. Review of Indications of FDA-Approved
Immune Checkpoint Inhibitors per NCCN Guidelines. *Cancers* 2020;12(3):738.
doi:10.3390/cancers12030738.

## Reading the table

- `fda_status` / `ema_status`: year of approval, `withdrawn`, or `no` (not approved by that agency).
- `pivotal_trial`: registrational trial(s) supporting the indication.
- `ref`: number(s) of the pivotal trial(s) in the supplementary reference list
  `ICI_SuppTable_io_approvals_references.md` (two numbers for a combined-trial row;
  blank for the four recent approvals without a trial publication).
- `on_atlas`: whether the organ is one of the 16 sites shown in the anatomical atlas / figure;
  off-atlas rows (head & neck, cervix, anal, nasopharynx, lymphoma, pleura, soft tissue,
  tumor-agnostic) are tabulated here but shown only as a footnote tally in the figure.

## Reference list (supplement-scoped)

The 137 pivotal-trial references are a self-contained reference list for this supplementary
table (`ICI_SuppTable_io_approvals_references.md`), not part of the main-manuscript
Mendeley-managed bibliography. Each is numbered in order of first appearance in the table and
keyed to a row through the `ref` column. Reference metadata was resolved by DOI via Crossref
(author lists from PubMed where the trial-group byline required it) and never typed; the
build is reproducible from `docs/references.csv` through `docs/tools/fetch_references.py` →
`docs/tools/build_supp_reflist.py`. One row's PMID was found to point to the article's erratum
rather than the article and was corrected during the fetch.

## Caveats (state in the figure/table legend)

1. **Tissue-agnostic approvals.** The tumor-agnostic rows (pembrolizumab MSI-H/dMMR and
   TMB-high; dostarlimab dMMR) apply across solid tumors. Sites counted as "0" have no
   **histology-specific** approval but a biomarker-selected minority of their patients can
   still receive an approved checkpoint inhibitor.
2. **"Skin"** aggregates all cutaneous malignancies (melanoma, Merkel cell, cutaneous SCC,
   basal cell).
3. **Withdrawn approvals** retained in the detail table for transparency but excluded from
   active counts: atezolizumab urothelial 2L and 1L cisplatin-ineligible (FDA); durvalumab
   urothelial 2L (FDA); atezolizumab metastatic TNBC (FDA only — the EU metastatic TNBC
   indication remains active; the 2021 EU withdrawal was of a separate early-stage application).

## Rows left without a pivotal trial (recent approvals — verify before submission)

- pembrolizumab + belzutifan, RCC adjuvant (FDA 2026)
- dostarlimab monotherapy, rectal MSI-H / LARC (FDA 2026)
- pembrolizumab + enfortumab vedotin, MIBC perioperative (FDA 2025, EMA 2026)
- nivolumab + brentuximab vedotin, classical Hodgkin (EMA 2025)
