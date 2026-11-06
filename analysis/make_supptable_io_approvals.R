# =============================================================================
# Supplementary detail table: FDA + EMA approved checkpoint-inhibitor regimens
# per cancer, with brand, pivotal trial, and per-drug regulatory references.
# Enriches results/tables/ICI_io_approval_landscape.csv.
# Regulatory basis: Drugs@FDA / openFDA + FDA OCE approval notifications (FDA);
# EMA European Public Assessment Reports (EU). Compiled 2026-11-06.
# =============================================================================
setwd("/Users/smwambu1/Desktop/Projects/ICI/analysis")
D <- read.csv("../results/tables/ICI_io_approval_landscape.csv", stringsAsFactors=FALSE)

ATLAS16 <- c("Brain","Thyroid","Esophagus","Lung","Gastric","Liver","Pancreas","Kidney",
             "Colon","Bladder","Skin","Breast","Ovary","Uterus","Prostate","Testis")

# ---- per-drug brand + regulatory reference anchors (verifiable) ----
EPAR <- "https://www.ema.europa.eu/en/medicines/human/EPAR/"
dref <- read.csv(text="
drug|brand|epar_slug
pembrolizumab|Keytruda|keytruda
nivolumab|Opdivo|opdivo
atezolizumab|Tecentriq|tecentriq
durvalumab|Imfinzi|imfinzi
cemiplimab|Libtayo|libtayo
avelumab|Bavencio|bavencio
ipilimumab|Yervoy|yervoy
nivolumab+relatlimab|Opdualag|opdualag
tislelizumab|Tevimbra|tevimbra
toripalimab|Loqtorzi|loqtorzi
dostarlimab|Jemperli|jemperli
retifanlimab|Zynyz|zynyz-0
cosibelimab|Unloxcyt|NA
", sep="|", stringsAsFactors=FALSE, strip.white=TRUE)
# relatlimab appears in the data as the Opdualag combo under drug 'nivolumab'; brand handled per-row below
dref$ema_epar <- ifelse(dref$epar_slug=="NA", "No EMA authorisation (US only)", paste0(EPAR, dref$epar_slug))
FDA_SRC <- "Drugs@FDA (accessdata.fda.gov/scripts/cder/daf) + FDA OCE approval notifications"

brand_of <- function(drug, regimen){
  b <- dref$brand[match(drug, dref$drug)]
  b[is.na(b)] <- ""
  b[grepl("relatlimab|Opdualag", regimen, ignore.case=TRUE)] <- "Opdualag"
  b
}
epar_of <- function(drug, regimen){
  e <- dref$ema_epar[match(drug, dref$drug)]
  e[is.na(e)] <- ""
  e[grepl("relatlimab|Opdualag", regimen, ignore.case=TRUE)] <- paste0(EPAR,"opdualag")
  e
}

# ---- pivotal trial per approval (drug, cancer_type, regimen) ----
# Confidently sourced for atlas organs + well-established off-atlas; left "" where uncertain.
T <- function(d,c,r,t) data.frame(drug=d, cancer_type=c, regimen=r, trial=t, stringsAsFactors=FALSE)
trial <- do.call(rbind, list(
  # Skin
  T("pembrolizumab","Melanoma","monotherapy (met)","KEYNOTE-006/-002"),
  T("pembrolizumab","Melanoma","adjuvant","KEYNOTE-054/-716"),
  T("pembrolizumab","Merkel cell","monotherapy","KEYNOTE-017"),
  T("pembrolizumab","cSCC","monotherapy","KEYNOTE-629"),
  T("nivolumab","Melanoma","monotherapy (met)","CheckMate-037/-066"),
  T("nivolumab","Melanoma","+ipilimumab","CheckMate-067"),
  T("nivolumab","Melanoma","adjuvant","CheckMate-238"),
  T("nivolumab","Melanoma","+relatlimab (Opdualag)","RELATIVITY-047"),
  T("atezolizumab","Melanoma","+cobimetinib+vemurafenib","IMspire150"),
  T("ipilimumab","Melanoma","monotherapy (met)","MDX010-20"),
  T("ipilimumab","Melanoma","adjuvant","EORTC 18071"),
  T("cemiplimab","cSCC","monotherapy","EMPOWER-CSCC-1 (Study 1540)"),
  T("cemiplimab","cSCC","adjuvant","C-POST"),
  T("cemiplimab","BCC","monotherapy","Study 1620"),
  T("avelumab","Merkel cell","monotherapy","JAVELIN Merkel 200"),
  T("retifanlimab","Merkel cell","monotherapy","POD1UM-201"),
  T("cosibelimab","cSCC","monotherapy","Study CK-301-101"),
  # Lung
  T("pembrolizumab","NSCLC","+pemetrexed+platinum (1L nonsq)","KEYNOTE-189"),
  T("pembrolizumab","NSCLC","+carbo/paclitaxel (1L sq)","KEYNOTE-407"),
  T("pembrolizumab","NSCLC","monotherapy (1L)","KEYNOTE-024/-042"),
  T("pembrolizumab","NSCLC","monotherapy (2L)","KEYNOTE-010"),
  T("pembrolizumab","NSCLC","+platinum perioperative","KEYNOTE-671"),
  T("pembrolizumab","NSCLC","monotherapy adjuvant","KEYNOTE-091 (PEARLS)"),
  T("nivolumab","NSCLC","+chemo neoadjuvant","CheckMate-816"),
  T("nivolumab","NSCLC","+chemo perioperative","CheckMate-77T"),
  T("nivolumab","NSCLC","+ipilimumab (1L)","CheckMate-227"),
  T("nivolumab","NSCLC","+ipilimumab+chemo (1L)","CheckMate-9LA"),
  T("nivolumab","NSCLC","monotherapy (2L)","CheckMate-017/-057"),
  T("atezolizumab","NSCLC","monotherapy adjuvant","IMpower010"),
  T("atezolizumab","NSCLC","monotherapy (1L high PD-L1)","IMpower110"),
  T("atezolizumab","NSCLC","+bev+paclitaxel+carbo (1L)","IMpower150"),
  T("atezolizumab","NSCLC","+nab-paclitaxel+carbo (1L)","IMpower130"),
  T("atezolizumab","NSCLC","monotherapy (2L)","OAK/POPLAR"),
  T("atezolizumab","NSCLC","monotherapy (1L platinum-ineligible)","IPSOS"),
  T("atezolizumab","SCLC","+carbo+etoposide (1L)","IMpower133"),
  T("atezolizumab","SCLC","+lurbinectedin maintenance","IMforte"),
  T("durvalumab","NSCLC","monotherapy consolidation (PACIFIC)","PACIFIC"),
  T("durvalumab","NSCLC","+chemo perioperative (AEGEAN)","AEGEAN"),
  T("durvalumab","NSCLC","+tremelimumab+chemo (POSEIDON)","POSEIDON"),
  T("durvalumab","SCLC","+etoposide+platinum (CASPIAN)","CASPIAN"),
  T("durvalumab","SCLC","monotherapy limited-stage (ADRIATIC)","ADRIATIC"),
  T("cemiplimab","NSCLC","monotherapy (1L PD-L1>=50)","EMPOWER-Lung 1"),
  T("cemiplimab","NSCLC","+chemo (1L)","EMPOWER-Lung 3"),
  T("tislelizumab","NSCLC","+chemo (1L sq)","RATIONALE-307"),
  T("tislelizumab","NSCLC","+chemo (1L non-sq)","RATIONALE-304"),
  T("tislelizumab","NSCLC","monotherapy (2L)","RATIONALE-303"),
  T("tislelizumab","NSCLC","perioperative","RATIONALE-315"),
  T("tislelizumab","SCLC","+etoposide+platinum (1L)","RATIONALE-312"),
  # Gastric
  T("pembrolizumab","Gastric/GEJ","+trastuzumab+chemo (HER2+)","KEYNOTE-811"),
  T("pembrolizumab","Gastric/GEJ","+chemo (HER2-)","KEYNOTE-859"),
  T("nivolumab","Gastric/GEJ","+chemo (1L)","CheckMate-649"),
  T("durvalumab","Gastric/GEJ","+FLOT perioperative (MATTERHORN)","MATTERHORN"),
  T("tislelizumab","Gastric/GEJ","+chemo (1L)","RATIONALE-305"),
  # Esophagus
  T("pembrolizumab","Esophageal/GEJ","+platinum+FU","KEYNOTE-590"),
  T("pembrolizumab","Esophageal","monotherapy (squamous 2L)","KEYNOTE-181"),
  T("nivolumab","Esophageal/GEJ","monotherapy adjuvant","CheckMate-577"),
  T("nivolumab","Esophageal SCC","+chemo (1L)","CheckMate-648"),
  T("nivolumab","Esophageal SCC","+ipilimumab (1L)","CheckMate-648"),
  T("nivolumab","Esophageal SCC","monotherapy (2L)","ATTRACTION-3"),
  T("toripalimab","Esophageal SCC","+chemo (1L)","JUPITER-06"),
  T("tislelizumab","Esophageal SCC","monotherapy (2L)","RATIONALE-302"),
  T("tislelizumab","Esophageal SCC","+chemo (1L)","RATIONALE-306"),
  # Bladder
  T("pembrolizumab","Urothelial","+enfortumab vedotin (1L)","EV-302 (KEYNOTE-A39)"),
  T("pembrolizumab","Urothelial","monotherapy (LA/met)","KEYNOTE-045/-052"),
  T("pembrolizumab","Bladder NMIBC","monotherapy (BCG-unresponsive)","KEYNOTE-057"),
  T("nivolumab","Urothelial","monotherapy adjuvant","CheckMate-274"),
  T("nivolumab","Urothelial","+cisplatin+gemcitabine (1L)","CheckMate-901"),
  T("nivolumab","Urothelial","monotherapy (2L)","CheckMate-275"),
  T("atezolizumab","Bladder MIBC","monotherapy adjuvant ctDNA+","IMvigor011"),
  T("atezolizumab","Urothelial","monotherapy (2L)","IMvigor210 cohort 2"),
  T("atezolizumab","Urothelial","monotherapy (1L cisplatin-ineligible)","IMvigor210 cohort 1"),
  T("durvalumab","Bladder MIBC","+gem+cis perioperative (NIAGARA)","NIAGARA"),
  T("durvalumab","Bladder NMIBC","+BCG (POTOMAC)","POTOMAC"),
  T("durvalumab","Urothelial","monotherapy (2L)","Study 1108"),
  T("avelumab","Urothelial","monotherapy 1L maintenance","JAVELIN Bladder 100"),
  T("avelumab","Urothelial","monotherapy (2L)","JAVELIN Solid Tumor"),
  # Kidney
  T("pembrolizumab","RCC","+axitinib (1L)","KEYNOTE-426"),
  T("pembrolizumab","RCC","+lenvatinib (1L)","CLEAR (KEYNOTE-581)"),
  T("pembrolizumab","RCC","monotherapy adjuvant","KEYNOTE-564"),
  T("nivolumab","RCC","+ipilimumab (1L)","CheckMate-214"),
  T("nivolumab","RCC","+cabozantinib (1L)","CheckMate-9ER"),
  T("nivolumab","RCC","monotherapy (2L)","CheckMate-025"),
  T("avelumab","RCC","+axitinib (1L)","JAVELIN Renal 101"),
  # Liver
  T("pembrolizumab","HCC","monotherapy (2L)","KEYNOTE-224/-240"),
  T("pembrolizumab","Biliary tract","+gem+cis (1L)","KEYNOTE-966"),
  T("nivolumab","HCC","+ipilimumab (1L)","CheckMate-9DW"),
  T("nivolumab","HCC","+ipilimumab (2L post-sorafenib)","CheckMate-040"),
  T("atezolizumab","HCC","+bevacizumab (1L)","IMbrave150"),
  T("durvalumab","Biliary tract","+gem+cis (TOPAZ-1)","TOPAZ-1"),
  T("durvalumab","HCC","+tremelimumab (STRIDE)","HIMALAYA"),
  T("durvalumab","HCC","monotherapy (1L)","HIMALAYA"),
  # Colon
  T("pembrolizumab","Colorectal MSI-H","monotherapy (incl 1L)","KEYNOTE-177/-164"),
  T("nivolumab","Colorectal MSI-H","+ipilimumab (1L)","CheckMate-8HW"),
  T("nivolumab","Colorectal MSI-H","+ipilimumab (2L)","CheckMate-142"),
  T("nivolumab","Colorectal MSI-H","monotherapy","CheckMate-142"),
  T("atezolizumab","Colon dMMR","+FOLFOX adjuvant (Stage III)","ATOMIC"),
  # Breast
  T("pembrolizumab","TNBC","+chemo neoadjuvant (early)","KEYNOTE-522"),
  T("pembrolizumab","TNBC","+sacituzumab govitecan (1L)","ASCENT-04 (KEYNOTE-D19)"),
  T("pembrolizumab","TNBC","+chemo (met)","KEYNOTE-355"),
  T("atezolizumab","TNBC","+nab-paclitaxel (1L)","IMpassion130"),
  # Uterus
  T("pembrolizumab","Endometrial","+carbo+paclitaxel","KEYNOTE-868 (NRG-GY018)"),
  T("pembrolizumab","Endometrial","+lenvatinib (pMMR)","KEYNOTE-775 (Study 309)"),
  T("pembrolizumab","Endometrial","monotherapy (MSI-H)","KEYNOTE-158"),
  T("durvalumab","Endometrial","+chemo->mono (dMMR)","DUO-E"),
  T("durvalumab","Endometrial","+chemo->+olaparib (pMMR)","DUO-E"),
  T("dostarlimab","Endometrial","monotherapy (dMMR)","GARNET"),
  T("dostarlimab","Endometrial","+chemo->mono (all-comers)","RUBY"),
  # Ovary
  T("pembrolizumab","Ovarian","+paclitaxel+-bevacizumab","KEYNOTE-B96 (ENGOT-ov65)"),
  # Head and neck
  T("pembrolizumab","HNSCC","neoadjuvant->adjuvant","KEYNOTE-689"),
  T("pembrolizumab","HNSCC","+platinum+FU (1L)","KEYNOTE-048"),
  T("pembrolizumab","HNSCC","monotherapy (1L)","KEYNOTE-048"),
  T("pembrolizumab","HNSCC","monotherapy (2L)","KEYNOTE-040/-012"),
  T("nivolumab","HNSCC","monotherapy (2L)","CheckMate-141"),
  # Cervix
  T("pembrolizumab","Cervical","+chemoradiotherapy","KEYNOTE-A18"),
  T("pembrolizumab","Cervical","+chemo+-bev","KEYNOTE-826"),
  T("pembrolizumab","Cervical","monotherapy","KEYNOTE-158"),
  T("cemiplimab","Cervical","monotherapy","EMPOWER-Cervical 1"),
  # Anal
  T("retifanlimab","Anal SCC","+chemo (1L)","POD1UM-303"),
  T("retifanlimab","Anal SCC","monotherapy (2L)","POD1UM-202"),
  # Nasopharynx
  T("toripalimab","Nasopharyngeal","+chemo (1L)","JUPITER-02"),
  T("toripalimab","Nasopharyngeal","monotherapy","POLARIS-02"),
  T("tislelizumab","Nasopharyngeal","+chemo (1L)","RATIONALE-309"),
  # Lymphoma
  T("pembrolizumab","Classical Hodgkin","monotherapy","KEYNOTE-204/-087"),
  T("pembrolizumab","PMBCL","monotherapy","KEYNOTE-170"),
  T("nivolumab","Classical Hodgkin","+AVD (1L)","SWOG S1826"),
  T("nivolumab","Classical Hodgkin","monotherapy","CheckMate-205"),
  # Pleura
  T("pembrolizumab","Mesothelioma","+pemetrexed+platinum","IND.227/KEYNOTE-483"),
  T("nivolumab","Mesothelioma","+ipilimumab","CheckMate-743"),
  # Soft tissue
  T("atezolizumab","Alveolar soft part sarcoma","monotherapy","Study ML39345"),
  # Tumor-agnostic
  T("pembrolizumab","MSI-H/dMMR solid tumors","monotherapy","KEYNOTE-158/-164"),
  T("pembrolizumab","TMB-high solid tumors","monotherapy","KEYNOTE-158"),
  T("dostarlimab","dMMR solid tumors","monotherapy","GARNET")
))

M <- merge(D, trial, by=c("drug","cancer_type","regimen"), all.x=TRUE)
M$trial[is.na(M$trial)] <- ""

# ---- assemble the detail table ----
out <- data.frame(
  on_atlas   = ifelse(M$organ %in% ATLAS16, "yes", "no"),
  organ      = M$organ,
  cancer_type= M$cancer_type,
  drug       = M$drug,
  brand      = brand_of(M$drug, M$regimen),
  regimen    = M$regimen,
  pivotal_trial = M$trial,
  fda_status = M$fda,
  ema_status = M$ema,
  fda_source = FDA_SRC,
  ema_epar   = epar_of(M$drug, M$regimen),
  stringsAsFactors=FALSE)

ord_org <- c(ATLAS16, sort(setdiff(unique(out$organ), ATLAS16)))
out <- out[order(match(out$organ, ord_org), out$drug, out$cancer_type), ]

write.csv(out, "../results/tables/ICI_SuppTable_io_approvals_detailed.csv", row.names=FALSE)

n_unmatched <- sum(out$pivotal_trial=="")
cat("rows:", nrow(out), " | pivotal trial assigned:", nrow(out)-n_unmatched,
    " | blank (uncertain):", n_unmatched, "\n")
if (n_unmatched>0){
  cat("\n--- rows with no pivotal trial (left blank, verify before submission) ---\n")
  print(out[out$pivotal_trial=="", c("organ","drug","cancer_type","regimen","fda_status","ema_status")], row.names=FALSE)
}
cat("\nwrote ../results/tables/ICI_SuppTable_io_approvals_detailed.csv\n")
