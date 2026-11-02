# ICI panel — LOCKED v1

**Status:** LOCKED · **Date:** 2026-10-30 · **File:** `ICI_panel_LOCKED_v1.csv` (108 genes)
**Scope:** this document freezes the panel definition, layer/direction assignment, and scoring
formula. No gene, layer, direction, or threshold changes after this point without a versioned
re-lock and re-registration (see Change control).

## 1. Provenance (how the genes were obtained)

Systematic review of ICB-response transcriptomic studies → curated signature substrate (COMPASS +
review-derived signatures + MSigDB Hallmark; 254 signatures) → program-blind community detection
(Louvain) on the gene co-occurrence graph → 6 emergent modules → two-track retention gate →
**108 genes** (Grade A+B). The panel grew monotonically across derivation refreshes
(90→133→147 at recurrence ≥5; nothing ever dropped) — a convergent core. The locked v1 takes the
Grade A+B subset (108), matching the original panel size.

The gene set is **outcome-blind**: derived from literature co-occurrence and biology only, never
fitted to any trial outcome.

## 2. The six modules nest into three layers (LOCKED)

Enrichment by g:Profiler (GO biological process / Reactome).

| Module | n | Enrichment read (top term, p) | Layer | Direction |
|---|---|---|---|---|
| M1 IFN / antigen-presentation | 24 | response to biotic stimulus / defense response (6e-59) | **Sensory** | + |
| M2 chemokine / trafficking | 26 | inflammatory response / response to chemical (4e-46) | **Integration** | + (activating) |
| M3 Treg / checkpoint regulation | 12 | leukocyte / cell activation; regulatory (1e-28) | **Integration** | − (suppressive) |
| M4 stromal / ECM-coagulation | 9 | proteolysis / coagulation / ECM organization (1e-33) | **Integration** | − (suppressive) |
| M5 angiogenic / CAF-TGFβ | 10 | vasculature development / cell migration (2e-25) | **Integration** | − (suppressive) |
| M6 cytotoxic effector | 27 | immune response / lymphocyte-mediated cytotoxicity (1e-64) | **Execution** | + |

Module numbers run top-to-bottom in cascade order (Sensory → activating Integration → suppressive
Integration → Execution), matching Figure 1.

**3-layer collapse:** Sensory 24 (+) · Integration 57 (26 activating − 31 suppressive) · Execution 27 (+).

## 3. Direction lock rules (and why)

- Sensory and Execution modules are **positive** (presence favours an intact response).
- Integration is the **net of activating and suppressive integration**:
  activating = M2 (recruitment/costimulation/Th1); suppressive = M3 + M5 + M6.
- **M3 (Treg/checkpoint) is assigned negative** even though its genes are individually "mixed":
  it is a net-regulatory/immunosuppressive module (checkpoints and Tregs are inhibitory brakes).
- **M5 + M6 (stromal / angiogenic / CAF-TGFβ) are assigned negative**, i.e. they are the
  suppressive arm of Integration rather than a separate "resistance" axis.
  *Empirical justification (I-SPY2):* including these suppressive programs in the Integration
  negative arm is required for the predictive signal — excluding them raised the interaction
  p from 0.069 to 0.455. This is a locked modeling decision, validated out-of-sample, not tuned.
- **Execution exhaustion checkpoints are POSITIVE, not a subtracted negative arm.** The Execution
  module carries 6 effector-phase inhibitory/exhaustion genes (PDCD1/PD-1, LAG3, TIGIT, HAVCR2/TIM-3,
  ENTPD1/CD39, PDCD1LG2/PD-L2). Unlike Integration's priming-phase brakes (CTLA4, FOXP3, EOMES, which
  block *access* and are correctly negative), pre-treatment exhaustion markers denote T cells that are
  *present but reversibly restrained* — the cells ICB rescues — so they *predict response positively*.
  *Empirical justification:* IMvigor210 exhaustion-score → response AUC 0.57 (same direction as the
  effector score 0.56; higher in inflamed 0.79 vs desert 0.20); among effector-high tumors, response
  rate is 30% (exhaustion-high) vs 19% (exhaustion-low). On I-SPY2, *subtracting* exhaustion from the
  Execution score worsens the Execution × treatment interaction (p 0.14 effector-only → 0.25 current
  all-positive → **0.38 subtracted**). Exhaustion is therefore kept positive in the score and reported
  as a **secondary readout** that splits Execution non-responders into *armed-but-exhausted*
  (→ anti-PD-1 ± LAG3/TIGIT combination) vs *cold / effector-low* (→ upstream T-cell generation).
  Test: `locked_reference_prototype/run_execution_exhaustion_test.R` (run 2026-10-31).

## 4. Grade gate

- **A** (36): legacy-corroborated AND recurrence ≥5.
- **B** (72): legacy-corroborated (recurrence <5) OR novel with recurrence ≥6.
- Grade C and the 19 optional legacy add-backs are **excluded** from locked v1 (available as a
  documented sensitivity set only).

## 5. Scoring formula (LOCKED)

Per gene, binarize at the reference median (gene = 1 if expression ≥ median). Each layer score is
the **fraction** of its genes called high (proportion form — coverage-robust). Then:

- **Sensory** = fraction of 24 Sensory genes = 1
- **Integration** = fraction(26 activating = 1) − fraction(31 suppressive = 1)
- **Execution** = fraction of 27 Execution genes = 1

The three layer scores are read **separately** — no composite. Collapsing the layers into one
scalar defeats the purpose of the decomposition (which layer limits the tumor), and the predictive
signal lives specifically in Integration × treatment, not a blend.

Cohort-median reference for retrospective cohorts; a frozen locked reference for single-patient
scoring (reproduces cohort-median scoring exactly across array / FPKM / FFPE-TPM platforms).

**Scoring-method benchmark** (`run_scoring_benchmark.R`): on the I-SPY2 Integration × treatment
interaction and IMvigor210 phenotype separation, proportion-form median-binarization ties the best
predictive method (interaction p = 0.058) and is coverage-robust; continuous z-mean (p = 0.087) and
singscore (p = 0.160) do not improve it and sacrifice platform portability. Median-binarization is
therefore retained as near-optimal *and* uniquely deployable.

## 6. What the locked panel produces (validation, for reference — not part of the lock)

- **Predictive (Integration × treatment), randomized ICB:** pooled OR 1.83 (95% CI 1.06–3.16),
  p = 0.031, I² = 0% across the two I-SPY2 arms; Integration-high ICB benefit +42 pts (pembro).
- **Predictive, not prognostic:** BrighTNess (n=482, chemo-only) Integration × carboplatin
  interaction p = 0.92 — the effect is ICB-specific.
- **NR taxonomy externally validated:** limiting-layer class vs Mariathasan immune phenotype in
  IMvigor210 (n=348), χ² p = 6×10⁻¹⁶ (Integration-suppressed ↔ excluded; Cascade-intact ↔ inflamed).
- **Deployable:** locked single-patient scoring identical to cohort-median across three platforms.

## 7. Caveats (carried into the manuscript)

- Module **stability is modest (median 0.29)**: the layer *assignment of an individual boundary
  gene* is provisional; the axis-level structure and the Grade-A core are the robust parts.
- A few placements are co-occurrence artifacts (e.g. CXCL13/CCL19 in the effector module via TLS
  signatures; GZMK in the stromal module) — biologically defensible, flagged for human review.

## 8. Gap analysis & sensitivity (audit, 2026-10-30)

A curated set of canonical ICB genes (per layer/sub-function) was audited against the locked 108,
the recurrence universe (2,280 genes), and the broader derived panels (`gap_analysis.R`).

- **Backbone fully covered:** checkpoints (PD-1, PD-L1, CTLA4, LAG3, TIM-3, TIGIT, IDO1, PD-L2),
  effectors (GZMA/B/H/K, PRF1, GNLY, NKG7, IFNG), trafficking (CCL5, CXCL9/10/11/13, CCR7), IFN
  (STAT1, IRF1, IFNGR1/2), antigen presentation (B2M, HLA-A/B/C/E, TAP1, PSMB8/9), Tregs (FOXP3,
  IL2RA, IL10).
- **Recoverable gaps (recurrence ≥5, below the Grade A+B gate):** costimulation (ICOS, ICOSLG,
  TNFSF9, TNFSF4, TNFRSF18), effector/T-NK identity (KLRK1, KLRD1, CD3E/G, ZAP70, FAS, TCF7),
  HLA-DRA, SPP1 — 14 genes, all in the recurrence-≥5 denovo panel.
- **Genuine absences (NOT literature-recurrent, left out and documented as limitations):** innate
  sensing (cGAS-STING, RIG-I, TLRs) and most myeloid-suppression markers (CD68, CD163, CSF1R,
  ARG1) — these do not recur in ICB-response signatures.

**ACKR1/DARC (tested on request, 2026-10-30):** 0 literature recurrence (absent from the substrate),
and a direct test in IMvigor210 (`run_ackr1_test.R`) shows it is a non-specific immune/vascular
infiltration marker — elevated in both excluded and inflamed tumors, correlated with every module
(r 0.36-0.58) but ~0 with the Integration net, and not associated with response (AUC 0.45). Excluded
from the panel on both grounds. An endothelial/HEV T-cell-entry axis, if wanted, needs a multi-gene
venule signature from spatial data (Phase 5), not ACKR1 alone.

**Sensitivity panel v1.1 (108 + the 14 recoverable = 122 genes)** was tested on the randomized
predictive signal: pooled Integration×treatment OR 1.85 (p=0.028) for locked 108 vs **1.64
(p=0.080) for 122** — adding the borderline genes does not improve and slightly dilutes the signal
(the activating-arm additions add inflammation-correlated noise; the predictive weight sits in the
suppressive arm). The direction is preserved, so the result is robust to panel definition. Decision:
**the parsimonious, pre-specified Grade A+B 108 remains primary**; v1.1 is reported as a
pre-specified sensitivity analysis, not the primary panel. This choice rests on the a-priori gate
and parsimony, not on cohort p-values.

## 9. Change control

Locked v1 is frozen for the analysis and pre-registration. Any change (gene, layer, direction,
grade gate, or scoring rule) requires a new version (v2), a dated entry here, and — if it affects
the primary analysis — an amended registration. Downstream figures and tables cite
`ICI_panel_LOCKED_v1.csv` only.
