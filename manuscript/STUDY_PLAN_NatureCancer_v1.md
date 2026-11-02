# ICI — Study Plan v1 (Nature Cancer / Cancer Cell scope)

**Date:** 2026-10-30 · **Authors:** SMM, DB · status: proposed program

## 1. The defining idea

Every established ICB response biomarker is a **single scalar** (TMB, PD-L1 IHC, the
T-cell–inflamed GEP/TIS, CYT, TIDE, COMPASS). They measure *how inflamed* a tumor is and
stop there. None localise *where in the anti-tumor response a given patient breaks*.

The Immune Cascade Index decomposes anti-tumor immunity into three mechanistic layers and
makes the predictive claim concrete:

> Anti-tumor immunity is a cascade — the **Sensory** layer detects the tumor, **Integration** weighs
> the activating-vs-suppressive decision, **Execution** carries out cytotoxic killing. ICB does not act
> on the whole cascade; it releases the brake within **Integration**. So response is predictable not from
> how inflamed a tumor is, but from whether Integration is the rate-limiting, *releasable* step: sensing
> and execution present, integration suppressed.

This reframes the work from "another signature" to **a per-patient map of why immunotherapy
fails** — and the I-SPY2 data support it: the predictive axis is **Integration × treatment**,
not Integration alone, not Sensory, not Execution.

## 2. Why it is unique and impactful

1. **Mechanistic, not black-box** — named biology; a reviewer argues with the mechanism, not a gene list.
2. **Predictive, not prognostic** — the interaction is measured in *randomized* ICB-vs-control trials.
3. **Per-patient failure localization** — which layer is limiting → rational therapy selection.
4. **Deployable** — locked-reference single-patient scoring, validated identical across three platforms.

## 3. The three layers (hybrid 6-module → 3-layer derivation)

The systematic review yields **6 data-driven modules** (program-blind Louvain on gene
co-occurrence; enrichment-labeled). They **nest** into the three cascade layers. The two
"resistance" modules are **not** a separate axis — they are the **suppressive arm of Integration**
(empirically confirmed: removing them from the negative arm destroys the predictive signal).

| Layer | Modules (from the review) | Direction | Biology |
|---|---|---|---|
| **SENSORY** | IFN / antigen-presentation | + | antigen & danger sensing, IFN induction, MHC |
| **INTEGRATION (+)** | chemokine / trafficking | + | T-cell recruitment, costimulation, Th1 |
| **INTEGRATION (−)** | Treg/checkpoint + stromal/ECM + angiogenic/CAF-TGFβ | − | regulatory & stromal immunosuppression (the releasable brake) |
| **EXECUTION** | cytotoxic effector | + | granzyme/perforin cytolytic killing |

Integration score = (activating integration) − (suppressive integration). This is why ICB,
which removes a coinhibitory/suppressive brake, converts Integration-suppressed tumors into responders.

**Derivation discipline (non-negotiable for the target journals):** layer membership and gene
direction are fixed from biology + module structure *before* outcome is seen, and locked. We do
**not** tune gene membership to maximize significance in any one cohort. Validation is
out-of-sample across the pooled cohorts (§5).

## 4. Verification to date (I-SPY2, n = 248; de-risking run 2026-10-30)

Scoring I-SPY2 with the **derived** panel (108 genes, never trained on outcome):

- Predictive signal is **Integration × treatment** — the three layers are read separately, **no
  composite**; Integration-high pembro benefit **+42 pts** vs +12 low.
- Resistance-in-suppressive-arm confirmed: Integration p 0.069 (included) → 0.455 (excluded).
- Derived Integration alone p = 0.069 (vs legacy curated panel 0.026; ρ = 0.67). Gene
  composition of the Integration layer is the open modeling choice — to be locked on biology and
  tested on the pool, **not** tuned on I-SPY2.
- Deployability: locked single-patient scoring = cohort-median scoring, exactly, on I-SPY2
  (array), GSE78220 (FPKM), NeoTRIP (FFPE-TPM) — three platforms.

## 5. Phased program

**Phase 1 — Lock the derivation (foundation).** Formalise 6 modules → 3 layers with enrichment
evidence per module; fix gene membership + direction; produce the locked derived panel
(gene · module · layer · direction · recurrence · grade · provenance). Deliverable: locked panel
CSV + derivation figure + Methods text. *Status: structure in hand (ICI_PANEL_v1_locked.csv);
needs principled finalisation + enrichment write-up.*

**Phase 2 — Score definition.** Per-layer locked-reference single-patient score — the three layers
read **separately (no composite)** — plus a "limiting-layer" call per patient. *Status: implemented
and platform-validated.*

**Phase 3 — Pooled randomized validation (the main result).** Assemble every obtainable
randomized ICB-vs-control transcriptomic trial; test Integration × treatment
across them; meta-analyze the interaction; show the signal is predictive (interaction), not
prognostic (control-arm AUC near 0.5). Per-patient failure map. *The biggest lever on journal tier.*

**Phase 4 — Prediction & benchmark.** Locked classifier; held-out / independent validation;
head-to-head vs TMB, PD-L1, GEP/TIS, COMPASS; stated decision rule + intended-use population.

**Phase 5 (stretch) — Mechanism / actionability.** Spatial / scRNA confirmation that layer calls
reflect real tissue states; layer-matched therapy hypotheses (integration-limited → ICB;
execution-limited → other; sensing-limited → other).

## 6. Cohort pool target (Phase 3)

Randomized ICB-vs-control, with tumor transcriptomes + response:

| Cohort | Tumor / setting | Arms | Accession | Access |
|---|---|---|---|---|
| I-SPY2 | breast, neoadj | pembro vs paclitaxel | GSE194040 | **public — in hand** |
| NeoTRIP | TNBC, neoadj | atezo+chemo vs chemo | GSE319641 | expr public; **arm+pCR from DB** |
| GeparNuevo | TNBC, neoadj | durvalumab vs placebo | GSE173839 (Karn) | public — to fetch |
| IMmotion150 | mRCC | atezo±bev vs sunitinib | EGA (McDermott 2018) | controlled — request |
| POPLAR/OAK | NSCLC | atezo vs docetaxel | supp / EGA | partial — assess |

Breast-heavy; RCC + NSCLC add tumor-type breadth (important for Nature Cancer). Single-arm ICB
cohorts (IMvigor210, Gide, Riaz, GSE78220) used only for score generalization, **never** for the
predictive interaction claim.

## 7. Immediate next actions

1. Inventory + download the public cohorts (GeparNuevo first; assess NSCLC).
2. Request from DB: NeoTRIP de-identified **arm + pCR** keyed by `TNBCxxxxx` (two columns).
3. Finalise Phase 1 locked panel + enrichment write-up.
4. Pre-register the analysis (layers, directions, primary = Integration × treatment,
   meta-analysis plan) **before** touching the pooled cohorts.
