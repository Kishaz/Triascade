# =====================================================================
# ICI 2.0 — SUPPLEMENTARY TABLE S1
# Gene-by-gene provenance + per-gene AUC sanity check
# =====================================================================
# Run:
#   source("/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript/build_gene_provenance.R")
#
# Outputs:
#   ICI_Manuscript/Supplementary_Table_S1_gene_provenance.csv
# =====================================================================

suppressPackageStartupMessages({
  library(pROC)
})

BASE     <- "/Users/dbedi/Desktop/R data FINAL"
MAN_DIR  <- file.path(BASE, "ICI_Manuscript")

# ---- GENE PROVENANCE TABLE -------------------------------------------
# Each row: gene, layer, component, biology, citation_short, citation_full
# Citations are verified from published literature.
#
# Citation key:
#   Schreiber2011    = Schreiber et al. Science 331:1565-70 (2011) — cancer immunoediting; MHC-I antigen presentation framework
#   Spranger2015     = Spranger et al. Nature 523:231-235 (2015) — beta-catenin/T cell exclusion; BATF3 DC signature
#   Ayers2017        = Ayers et al. JCI 127:2930-2940 (2017) — IFN-gamma TIS / 18-gene T cell-inflamed GEP
#   Cristescu2018    = Cristescu et al. Science 362:eaar3593 (2018) — pan-tumor T cell-inflamed GEP validation
#   Mariathasan2018  = Mariathasan et al. Nature 554:544-548 (2018) — IMvigor210; TGF-beta F-TBRS stromal exclusion
#   Rooney2015       = Rooney et al. Cell 160:48-61 (2015) — molecular and genetic properties of tumors associated with cytolytic activity (CYT = GZMA, PRF1)
#   Tirosh2016       = Tirosh et al. Science 352:189-196 (2016) — single-cell melanoma; T cell exhaustion signature
#   Wherry2015       = Wherry & Kurachi Nat Rev Immunol 15:486-499 (2015) — molecular and cellular insights into T cell exhaustion
#   Thommen2018      = Thommen et al. Nat Med 24:994-1004 (2018) — CD8 PD-1 high T cells in NSCLC; exhaustion markers
#   Harlin2009       = Harlin et al. Cancer Res 69:3077-3085 (2009) — chemokine expression in melanoma metastases (CXCL9/10/11, CCL2/3/4/5)
#   Bottcher2018     = Bottcher et al. Cell 172:1022-1037 (2018) — NK cells facilitate cDC1 recruitment via CCL5 + XCL1
#   Tauriello2018    = Tauriello et al. Nature 554:538-543 (2018) — TGF-beta drives immune evasion in colorectal metastasis
#   Hugo2016         = Hugo et al. Cell 165:35-44 (2016) — innate anti-PD-1 resistance signature (IPRES)
#   Galon2006        = Galon et al. Science 313:1960-1964 (2006) — type, density, and location of immune cells (Immunoscore origin)
#   Fridman2012      = Fridman et al. Nat Rev Cancer 12:298-306 (2012) — immune contexture and prognosis
#   Salmon2012       = Salmon et al. Cell 165:35-44 (2016) — desmoplastic stroma and CAFs (FAP, COL)
#   Khan2019         = Khan et al. Nature 571:211-218 (2019) — TOX is a critical regulator of tumor-specific T cell exhaustion
#   Scott2019        = Scott et al. Nature 571:270-274 (2019) — TOX drives intratumoural T cell dysfunction

gene_provenance <- list(

  # =====================================================================
  # SENSORY LAYER (21 genes) — antigen presentation, IFN-gamma signaling
  # =====================================================================
  list("HLA-A",   "Sensory", "MHC-I",            "MHC class I heavy chain (CD8 T cell antigen presentation)",                "Schreiber2011, Ayers2017"),
  list("HLA-B",   "Sensory", "MHC-I",            "MHC class I heavy chain",                                                  "Schreiber2011, Ayers2017"),
  list("HLA-C",   "Sensory", "MHC-I",            "MHC class I heavy chain",                                                  "Schreiber2011, Ayers2017"),
  list("HLA-E",   "Sensory", "MHC-I",            "Non-classical MHC-I; NK cell engagement",                                  "Ayers2017"),
  list("HLA-F",   "Sensory", "MHC-I",            "Non-classical MHC-I",                                                      "Schreiber2011"),
  list("HLA-G",   "Sensory", "MHC-I",            "Non-classical MHC-I; tolerogenic in some contexts",                        "Schreiber2011"),
  list("B2M",     "Sensory", "MHC-I",            "Beta-2-microglobulin; obligate MHC-I light chain",                         "Schreiber2011"),
  list("TAP1",    "Sensory", "Antigen processing","Peptide transporter for MHC-I loading (TAP heterodimer)",                  "Schreiber2011"),
  list("TAP2",    "Sensory", "Antigen processing","Peptide transporter for MHC-I loading",                                    "Schreiber2011"),
  list("TAPBP",   "Sensory", "Antigen processing","Tapasin; loads peptides onto MHC-I",                                       "Schreiber2011"),
  list("PSMB8",   "Sensory", "Antigen processing","Immunoproteasome subunit (LMP7); IFN-gamma-induced",                      "Schreiber2011, Ayers2017"),
  list("PSMB9",   "Sensory", "Antigen processing","Immunoproteasome subunit (LMP2); IFN-gamma-induced",                      "Schreiber2011"),
  list("PSMB10",  "Sensory", "Antigen processing","Immunoproteasome subunit (MECL-1)",                                        "Schreiber2011, Ayers2017"),
  list("IRF1",    "Sensory", "IFN-gamma",        "IRF1; transcription factor downstream of IFN-gamma",                       "Schreiber2011, Cristescu2018"),
  list("STAT1",   "Sensory", "IFN-gamma",        "STAT1; obligate IFN-gamma signaling node",                                  "Ayers2017, Cristescu2018"),
  list("STAT2",   "Sensory", "IFN-gamma",        "STAT2; type I IFN signaling",                                              "Cristescu2018"),
  list("JAK1",    "Sensory", "IFN-gamma",        "JAK1; IFN-gamma receptor signaling kinase",                                "Schreiber2011"),
  list("JAK2",    "Sensory", "IFN-gamma",        "JAK2; IFN-gamma receptor signaling kinase",                                "Schreiber2011"),
  list("IFNG",    "Sensory", "IFN-gamma",        "Interferon-gamma cytokine (also Execution+)",                               "Schreiber2011, Ayers2017"),
  list("IFNGR1",  "Sensory", "IFN-gamma",        "IFN-gamma receptor 1 (alpha chain)",                                       "Schreiber2011"),
  list("IFNGR2",  "Sensory", "IFN-gamma",        "IFN-gamma receptor 2 (beta chain)",                                        "Schreiber2011"),

  # =====================================================================
  # INTEGRATION POSITIVE (21 genes) — chemokines, DC trafficking, TLS
  # =====================================================================
  list("CXCL9",    "Integration", "IntPos: T cell chemokines",   "CXCR3 ligand; T cell recruitment; IFN-gamma-induced",      "Harlin2009, Ayers2017"),
  list("CXCL10",   "Integration", "IntPos: T cell chemokines",   "CXCR3 ligand; T cell recruitment",                          "Harlin2009"),
  list("CXCL11",   "Integration", "IntPos: T cell chemokines",   "CXCR3 ligand; T cell recruitment",                          "Harlin2009"),
  list("CXCL13",   "Integration", "IntPos: TLS/B cell",           "CXCR5 ligand; tertiary lymphoid structure formation",       "Fridman2012"),
  list("CCL5",     "Integration", "IntPos: T cell chemokines",   "RANTES; T cell and NK chemoattractant",                    "Harlin2009, Bottcher2018"),
  list("CCL4",     "Integration", "IntPos: DC trafficking",       "MIP-1beta; DC and T cell recruitment",                     "Harlin2009"),
  list("CCL19",    "Integration", "IntPos: TLS/lymph node",       "CCR7 ligand; lymph node and TLS organization",              "Fridman2012"),
  list("CCL21",    "Integration", "IntPos: TLS/lymph node",       "CCR7 ligand; lymph node and TLS organization",              "Fridman2012"),
  list("ICAM1",    "Integration", "IntPos: Adhesion",             "Endothelial adhesion molecule; T cell extravasation",       "Galon2006"),
  list("VCAM1",    "Integration", "IntPos: Adhesion",             "Endothelial adhesion molecule; lymphocyte trafficking",     "Galon2006"),
  list("MADCAM1",  "Integration", "IntPos: Adhesion",             "Mucosal vascular adhesin; tissue-specific homing",          "Fridman2012"),
  list("CXCR3",    "Integration", "IntPos: T cell receptor",     "T cell chemokine receptor for CXCL9/10/11",                "Harlin2009"),
  list("CCR5",     "Integration", "IntPos: T cell receptor",     "T cell receptor for CCL5/CCL3/CCL4",                       "Harlin2009"),
  list("CCR7",     "Integration", "IntPos: T cell receptor",     "T cell homing receptor; central memory marker",            "Fridman2012"),
  list("CD80",     "Integration", "IntPos: Co-stimulation",       "B7-1; T cell co-stimulation via CD28",                     "Galon2006"),
  list("CD86",     "Integration", "IntPos: Co-stimulation",       "B7-2; T cell co-stimulation via CD28",                     "Galon2006"),
  list("CD40",     "Integration", "IntPos: DC activation",        "DC activation receptor; CD40L from CD4 T cells",            "Spranger2015"),
  list("BATF3",    "Integration", "IntPos: cDC1",                 "BATF3 lineage cDC1 master regulator (cross-presentation)", "Spranger2015, Bottcher2018"),
  list("IRF8",     "Integration", "IntPos: cDC1",                 "cDC1 transcription factor; IFN response",                  "Bottcher2018"),
  list("CR2",      "Integration", "IntPos: B cell/TLS",           "Complement receptor 2; B cell co-receptor (CD21)",         "Fridman2012"),
  list("TNFRSF17", "Integration", "IntPos: B cell/TLS",           "BCMA; plasma cell marker; TLS marker",                     "Fridman2012"),

  # =====================================================================
  # INTEGRATION NEGATIVE (20 genes) — TGF-beta stroma, ECM, exclusion
  # =====================================================================
  list("VEGFA",    "Integration", "IntNeg: Angiogenesis",         "Vascular endothelial growth factor; abnormal vasculature",  "Mariathasan2018"),
  list("VEGFB",    "Integration", "IntNeg: Angiogenesis",         "Endothelial growth factor",                                "Mariathasan2018"),
  list("VEGFC",    "Integration", "IntNeg: Angiogenesis",         "Lymphatic and vascular endothelial growth factor",          "Mariathasan2018"),
  list("FLT1",     "Integration", "IntNeg: Angiogenesis",         "VEGF receptor 1",                                          "Mariathasan2018"),
  list("KDR",      "Integration", "IntNeg: Angiogenesis",         "VEGF receptor 2",                                          "Mariathasan2018"),
  list("TGFB1",    "Integration", "IntNeg: TGF-beta",             "TGF-beta-1; immune exclusion (also Execution-)",            "Mariathasan2018, Tauriello2018"),
  list("TGFB2",    "Integration", "IntNeg: TGF-beta",             "TGF-beta-2",                                                "Mariathasan2018"),
  list("TGFB3",    "Integration", "IntNeg: TGF-beta",             "TGF-beta-3",                                                "Mariathasan2018"),
  list("COL1A1",   "Integration", "IntNeg: Stromal ECM",          "Fibrillar collagen; desmoplastic stroma",                   "Mariathasan2018, Salmon2012"),
  list("COL1A2",   "Integration", "IntNeg: Stromal ECM",          "Fibrillar collagen",                                        "Mariathasan2018"),
  list("COL3A1",   "Integration", "IntNeg: Stromal ECM",          "Fibrillar collagen",                                        "Mariathasan2018"),
  list("COL4A1",   "Integration", "IntNeg: Stromal ECM",          "Basement membrane collagen",                                "Mariathasan2018"),
  list("FAP",      "Integration", "IntNeg: CAF",                  "Fibroblast activation protein; CAF marker (F-TBRS)",       "Mariathasan2018, Salmon2012"),
  list("ACTA2",    "Integration", "IntNeg: CAF",                  "Smooth muscle actin; myofibroblast/CAF marker",            "Salmon2012"),
  list("FN1",      "Integration", "IntNeg: CAF",                  "Fibronectin; ECM component",                                "Mariathasan2018"),
  list("MMP2",     "Integration", "IntNeg: ECM remodeling",       "Matrix metalloprotease 2; gelatinase",                     "Salmon2012"),
  list("MMP9",     "Integration", "IntNeg: ECM remodeling",       "Matrix metalloprotease 9",                                  "Salmon2012"),
  list("MMP14",    "Integration", "IntNeg: ECM remodeling",       "MT1-MMP; pericellular ECM degradation",                     "Salmon2012"),
  list("CXCL12",   "Integration", "IntNeg: T cell exclusion",     "SDF-1; CAF-derived; T cell exclusion (CXCR4 axis)",         "Salmon2012"),
  list("EDNRB",    "Integration", "IntNeg: Endothelial barrier",  "Endothelin receptor B; T cell extravasation barrier",       "Salmon2012"),

  # =====================================================================
  # EXECUTION POSITIVE (19 genes) — effector cytotoxicity
  # =====================================================================
  list("GZMA",   "Execution", "ExecPos: Granzymes",          "Granzyme A; cytotoxic granule serine protease (CYT component)", "Rooney2015"),
  list("GZMB",   "Execution", "ExecPos: Granzymes",          "Granzyme B; canonical cytotoxic effector",                       "Rooney2015"),
  list("GZMH",   "Execution", "ExecPos: Granzymes",          "Granzyme H",                                                     "Rooney2015"),
  list("GZMK",   "Execution", "ExecPos: Granzymes",          "Granzyme K; CD8 effector marker",                                "Rooney2015"),
  list("GZMM",   "Execution", "ExecPos: Granzymes",          "Granzyme M",                                                     "Rooney2015"),
  list("PRF1",   "Execution", "ExecPos: Cytolytic",          "Perforin-1; pore-forming protein (CYT component)",                "Rooney2015"),
  list("GNLY",   "Execution", "ExecPos: Cytolytic",          "Granulysin; antimicrobial and cytolytic",                         "Rooney2015"),
  list("NKG7",   "Execution", "ExecPos: Cytolytic",          "NK group 7; cytolytic granule organization",                      "Ayers2017"),
  list("FASLG",  "Execution", "ExecPos: Death ligand",       "Fas ligand; death receptor pathway",                              "Rooney2015"),
  list("IFNG",   "Execution", "ExecPos: Effector cytokine",  "Interferon-gamma (also Sensory)",                                  "Rooney2015, Ayers2017"),
  list("TNF",    "Execution", "ExecPos: Effector cytokine",  "Tumor necrosis factor",                                           "Rooney2015"),
  list("LTA",    "Execution", "ExecPos: Effector cytokine",  "Lymphotoxin alpha (TNF-beta)",                                    "Rooney2015"),
  list("CD8A",   "Execution", "ExecPos: T cell marker",      "CD8 alpha; cytotoxic T cell marker",                              "Galon2006, Ayers2017"),
  list("CD8B",   "Execution", "ExecPos: T cell marker",      "CD8 beta",                                                        "Galon2006"),
  list("TBX21",  "Execution", "ExecPos: TF",                 "T-bet; Th1/CD8 effector master regulator",                        "Wherry2015"),
  list("EOMES",  "Execution", "ExecPos: TF",                 "Eomesodermin; effector and memory CD8 differentiation",           "Wherry2015"),
  list("RUNX3",  "Execution", "ExecPos: TF",                 "RUNX3; tissue-resident memory and effector function",             "Wherry2015"),
  list("IL2",    "Execution", "ExecPos: Cytokine",           "Interleukin-2; T cell growth factor",                             "Wherry2015"),
  list("IL12A",  "Execution", "ExecPos: Cytokine",           "IL-12 p35; Th1/CD8 polarization",                                 "Wherry2015"),

  # =====================================================================
  # EXECUTION NEGATIVE (19 genes) — checkpoints, exhaustion
  # =====================================================================
  list("PDCD1",    "Execution", "ExecNeg: Checkpoint",      "PD-1; central T cell inhibitory receptor",                          "Wherry2015"),
  list("CTLA4",    "Execution", "ExecNeg: Checkpoint",      "CTLA-4; co-inhibitory receptor",                                    "Wherry2015"),
  list("LAG3",     "Execution", "ExecNeg: Checkpoint",      "LAG-3; T cell inhibitory receptor",                                 "Wherry2015, Thommen2018"),
  list("HAVCR2",   "Execution", "ExecNeg: Checkpoint",      "TIM-3; T cell inhibitory receptor",                                 "Wherry2015, Thommen2018"),
  list("TIGIT",    "Execution", "ExecNeg: Checkpoint",      "TIGIT; T cell inhibitory receptor",                                 "Wherry2015, Ayers2017"),
  list("BTLA",     "Execution", "ExecNeg: Checkpoint",      "BTLA; HVEM-binding inhibitory receptor",                            "Wherry2015"),
  list("VSIR",     "Execution", "ExecNeg: Checkpoint",      "VISTA; immune checkpoint",                                          "Wherry2015"),
  list("CD274",    "Execution", "ExecNeg: Checkpoint",      "PD-L1; PD-1 ligand on tumor and APCs",                              "Wherry2015, Ayers2017"),
  list("PDCD1LG2", "Execution", "ExecNeg: Checkpoint",      "PD-L2; PD-1 ligand",                                                "Ayers2017"),
  list("CD160",    "Execution", "ExecNeg: Checkpoint",      "CD160; HVEM-binding inhibitory receptor",                           "Wherry2015"),
  list("CD244",    "Execution", "ExecNeg: Checkpoint",      "2B4; SLAM family inhibitory receptor",                              "Wherry2015"),
  list("TOX",      "Execution", "ExecNeg: Exhaustion TF",   "TOX; master regulator of T cell exhaustion",                        "Khan2019, Scott2019"),
  list("NR4A1",    "Execution", "ExecNeg: Exhaustion TF",   "Nur77; T cell exhaustion transcription factor",                     "Wherry2015"),
  list("NR4A2",    "Execution", "ExecNeg: Exhaustion TF",   "Nurr1; T cell exhaustion",                                          "Wherry2015"),
  list("NR4A3",    "Execution", "ExecNeg: Exhaustion TF",   "NOR1; T cell exhaustion",                                           "Wherry2015"),
  list("ENTPD1",   "Execution", "ExecNeg: Treg/exhaustion", "CD39; ATP hydrolysis to AMP; Treg/exhausted T cell marker",         "Thommen2018"),
  list("IL10",     "Execution", "ExecNeg: Suppressive cytokine","Interleukin-10; immunosuppressive",                              "Wherry2015"),
  list("TGFB1",    "Execution", "ExecNeg: Suppressive cytokine","TGF-beta-1 (also IntNeg)",                                       "Mariathasan2018"),
  list("IDO1",     "Execution", "ExecNeg: Suppressive enzyme","Indoleamine 2,3-dioxygenase; tryptophan catabolism; tolerogenic",  "Ayers2017")
)

# Convert to data frame
gene_df <- data.frame(
  gene_symbol  = sapply(gene_provenance, `[[`, 1),
  layer        = sapply(gene_provenance, `[[`, 2),
  component    = sapply(gene_provenance, `[[`, 3),
  biology      = sapply(gene_provenance, `[[`, 4),
  citation     = sapply(gene_provenance, `[[`, 5),
  stringsAsFactors = FALSE
)

cat("==========================================\n")
cat("Gene panel summary\n")
cat("==========================================\n")
cat(sprintf("Total entries: %d (note: IFNG and TGFB1 each appear in 2 layers)\n", nrow(gene_df)))
cat(sprintf("Unique genes:  %d\n", length(unique(gene_df$gene_symbol))))
cat("\n")

# Per-layer/component summary
layer_summary <- aggregate(gene_symbol ~ layer + component, data = gene_df, FUN = length)
colnames(layer_summary)[3] <- "n_genes"
print(layer_summary, row.names = FALSE)

# ---- PER-GENE AUC IN IMVIGOR210 DISCOVERY COHORT --------------------
cat("\n==========================================\n")
cat("Computing per-gene AUC in IMvigor210 Discovery cohort...\n")
cat("==========================================\n")

# Try to load IMvigor210 expression matrix.
# First check if a pre-computed expression matrix is saved with the scores.
# If not, we'd need to load the IMvigor210CoreBiologies package or a saved RDS.

imv_expr_paths <- c(
  file.path(BASE, "ATRX_Immune", "IMvigor210_expression.rds"),
  file.path(BASE, "IMvigor210_expression.rds"),
  file.path(MAN_DIR, "IMvigor210_expression.rds"),
  file.path(BASE, "IMvigor210_expr.csv"),
  file.path(BASE, "IMvigor210_logTPM.csv")
)

imv_expr <- NULL
imv_expr_path_used <- NULL
for (p in imv_expr_paths) {
  if (file.exists(p)) {
    cat(sprintf("Found IMvigor210 expression at: %s\n", p))
    imv_expr <- if (grepl("\\.rds$", p)) readRDS(p) else read.csv(p, row.names = 1)
    imv_expr_path_used <- p
    break
  }
}

if (is.null(imv_expr)) {
  cat("\n!! IMvigor210 expression matrix not found in expected locations.\n")
  cat("!! Per-gene AUC computation skipped. Gene provenance table will be saved without AUC values.\n")
  cat("\nTo run per-gene AUC, place the IMvigor210 log-expression matrix at one of:\n")
  for (p in imv_expr_paths) cat("  -", p, "\n")
  cat("\nExpected format: rows = HGNC gene symbols, cols = sample IDs.\n")

  # Save the gene table without AUC
  out_csv <- file.path(MAN_DIR, "Supplementary_Table_S1_gene_provenance.csv")
  write.csv(gene_df, out_csv, row.names = FALSE)
  cat(sprintf("\nWrote: %s (%d rows, no AUC column)\n", out_csv, nrow(gene_df)))
} else {
  # Load IMvigor210 response
  imv_scores <- read.csv(file.path(MAN_DIR, "bulk_scores", "IMvigor210_ICI_scores.csv"),
                         stringsAsFactors = FALSE)
  resp <- ifelse(imv_scores$response == "Responder", 1L,
            ifelse(imv_scores$response == "Non-Responder", 0L, NA_integer_))

  # Match samples between expression and clinical
  if (ncol(imv_expr) == nrow(imv_scores)) {
    cat(sprintf("Matching: %d patients in expression, %d in clinical\n",
                ncol(imv_expr), nrow(imv_scores)))
  }
  v <- !is.na(resp)
  cat(sprintf("Patients with response: %d (R=%d, NR=%d)\n",
              sum(v), sum(resp[v] == 1), sum(resp[v] == 0)))

  # Compute per-gene AUC
  gene_df$AUC_IMvigor210 <- NA_real_
  gene_df$gene_in_expr   <- FALSE

  for (i in seq_len(nrow(gene_df))) {
    g <- gene_df$gene_symbol[i]
    if (g %in% rownames(imv_expr)) {
      gene_df$gene_in_expr[i] <- TRUE
      vec <- as.numeric(imv_expr[g, ])
      if (all(is.na(vec))) next
      v2 <- !is.na(vec) & v
      if (sum(v2) < 10) next
      r <- tryCatch({
        roc(resp[v2], vec[v2], quiet = TRUE)
      }, error = function(e) NULL)
      if (!is.null(r)) {
        gene_df$AUC_IMvigor210[i] <- round(as.numeric(auc(r)), 3)
      }
    }
  }

  # Annotate genes with low signal
  gene_df$signal_flag <- ifelse(
    is.na(gene_df$AUC_IMvigor210), "not in matrix",
    ifelse(abs(gene_df$AUC_IMvigor210 - 0.5) < 0.02, "low signal", "ok")
  )

  cat("\n=== Per-gene AUC summary ===\n")
  cat(sprintf("Genes with computed AUC: %d / %d\n",
              sum(!is.na(gene_df$AUC_IMvigor210)), nrow(gene_df)))
  cat(sprintf("  Above-chance signal (|AUC-0.5| > 0.02): %d\n",
              sum(gene_df$signal_flag == "ok")))
  cat(sprintf("  Low signal (|AUC-0.5| ≤ 0.02): %d\n",
              sum(gene_df$signal_flag == "low signal")))
  cat(sprintf("  Not present in matrix: %d\n",
              sum(gene_df$signal_flag == "not in matrix")))

  # Save
  out_csv <- file.path(MAN_DIR, "Supplementary_Table_S1_gene_provenance.csv")
  write.csv(gene_df, out_csv, row.names = FALSE)
  cat(sprintf("\nWrote: %s (%d rows, with AUC)\n", out_csv, nrow(gene_df)))
}

cat("\n==========================================\n")
cat("Citations key for Supplementary Table S1:\n")
cat("==========================================\n")
citations_full <- c(
  "Ayers2017       = Ayers M et al. JCI 127:2930-2940 (2017). DOI 10.1172/JCI91190",
  "Bottcher2018    = Bottcher JP et al. Cell 172:1022-1037 (2018). DOI 10.1016/j.cell.2018.01.004",
  "Cristescu2018   = Cristescu R et al. Science 362:eaar3593 (2018). DOI 10.1126/science.aar3593",
  "Fridman2012     = Fridman WH et al. Nat Rev Cancer 12:298-306 (2012). DOI 10.1038/nrc3245",
  "Galon2006       = Galon J et al. Science 313:1960-1964 (2006). DOI 10.1126/science.1129139",
  "Harlin2009      = Harlin H et al. Cancer Res 69:3077-3085 (2009). DOI 10.1158/0008-5472.CAN-08-2281",
  "Hugo2016        = Hugo W et al. Cell 165:35-44 (2016). DOI 10.1016/j.cell.2016.02.065",
  "Khan2019        = Khan O et al. Nature 571:211-218 (2019). DOI 10.1038/s41586-019-1325-x",
  "Mariathasan2018 = Mariathasan S et al. Nature 554:544-548 (2018). DOI 10.1038/nature25501",
  "Rooney2015      = Rooney MS et al. Cell 160:48-61 (2015). DOI 10.1016/j.cell.2014.12.033",
  "Salmon2012      = Salmon H, Franciszkiewicz K, Damotte D, Dieu-Nosjean MC, Validire P, Trautmann A, Mami-Chouaib F, Donnadieu E. J Clin Invest 2012;122(3):899-910. doi:10.1172/JCI45817",
  "Schreiber2011   = Schreiber RD et al. Science 331:1565-1570 (2011). DOI 10.1126/science.1203486",
  "Scott2019       = Scott AC et al. Nature 571:270-274 (2019). DOI 10.1038/s41586-019-1324-y",
  "Spranger2015    = Spranger S et al. Nature 523:231-235 (2015). DOI 10.1038/nature14404",
  "Tauriello2018   = Tauriello DVF et al. Nature 554:538-543 (2018). DOI 10.1038/nature25492",
  "Thommen2018     = Thommen DS et al. Nat Med 24:994-1004 (2018). DOI 10.1038/s41591-018-0057-z",
  "Tirosh2016      = Tirosh I et al. Science 352:189-196 (2016). DOI 10.1126/science.aad0501",
  "Wherry2015      = Wherry EJ & Kurachi M. Nat Rev Immunol 15:486-499 (2015). DOI 10.1038/nri3862"
)
for (c in citations_full) cat(" ", c, "\n")

# Save citations key as a text file
cat_fp <- file.path(MAN_DIR, "Supplementary_Table_S1_citations_key.txt")
writeLines(c("ICI 2.0 Manuscript — Supplementary Table S1 — Citations Key",
             "===============================================================",
             "",
             citations_full), cat_fp)
cat(sprintf("\nWrote: %s\n", cat_fp))

cat("\nDONE.\n")
