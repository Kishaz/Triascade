# =============================================================================
SCRIPT_VERSION <- "ICI_patient_landscape  v2.5  2026-10-29"
# v2.5: analysis names as in your scripts (Rose BLCA, Gide-Pembro, Gide-IpiPembro, Hugo SKCM,
#       Riaz SKCM, Braun ccRCC ...), the same names as Figure 2B. Nothing else changed.
# v2.4: row labels "ICB target"/"ICB drug" renamed "Target"/"Treatment", since two
#       analyses did not receive checkpoint blockade (I-SPY2 paclitaxel only; Friedrich,
#       a bispecific T-cell engager). Nothing else changed.
# v2.3: treatment labels corrected from the source studies: Friedrich myeloma
#       (GSE217245) received a BCMA×CD3 bispecific T-cell engager, not PD-(L)1;
#       EAC (GSE165252, PERFECT) received atezolizumab with chemoradiotherapy.
# v2.2: High = at or above the median (>=), the rule your archetype, survival and
#       gene-level scripts use, so Figure 2A matches the archetype numbers exactly.
# v2.1: the two NSCLC studies drawn and binarized separately (13 analyses).
# Patient landscape across all cohort x treatment-arm splits (Figure 2A)
# Same design as v1. Changes:
#   - paper cohorts only: I-SPY2 at n=248 (246 with Integration), Huang out
#   - I-SPY2 control arm: ICB target/drug = None (paclitaxel only)
#   - Gide Ipi+Pembro: target PD-1 + CTLA-4, drug Ipi+Pembro
#   - NSCLC drawn as its two studies, each binarized at its own median (13 analyses);
#     GSE126044 labelled nivolumab / PD-1
#   - Tier row removed
#   - one-line column labels so small cohorts do not overlap
#   - saved to Oct 16/figures
# =============================================================================
need <- c("ggplot2","dplyr","tidyr")
miss <- need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")
suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })
cat("\n", SCRIPT_VERSION, "\n", sep = "")

MAN  <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"
GRID <- file.path(MAN, "MASTER_GRID")
OUT  <- file.path(MAN, "Oct 16", "figures")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- cohort metadata --------------------------------------------------------
META <- data.frame(
  key    = c("imvigor","ispy2","braun","rose","gide","hugo","riaz","nsclc","eac","friedrich"),
  label  = c("IMvigor210","I-SPY2","Braun","Rose","Gide","Hugo","Riaz","NSCLC","EAC","Friedrich"),
  cancer = c("Urothelial","Breast","Renal","Urothelial","Melanoma","Melanoma","Melanoma",
             "NSCLC","Esophageal","Myeloma"),
  target = c("PD-L1","PD-1","PD-1","PD-1/PD-L1","PD-1","PD-1","PD-1",
             "PD-1/PD-L1","PD-1/PD-L1","PD-1/PD-L1"),
  drug   = c("Atezolizumab","Pembrolizumab","Nivolumab","PD-(L)1","Pembrolizumab",
             "Pembrolizumab","Nivolumab","PD-(L)1","PD-(L)1","PD-(L)1"),
  stringsAsFactors = FALSE)

# paper cohort files (I-SPY2 = n=248 rebuild, 2 controls without CXCL9 removed)
FILES <- c(
  imvigor   = "ICI_binary_imvigor_scores.csv",
  rose      = "ICI_binary_rose_scores.csv",
  braun     = "ICI_binary_braun_scores.csv",
  ispy2     = "ICI_binary_ispy2_scores_246complete.csv",
  gide      = "ICI_binary_gide_scores.csv",
  hugo      = "ICI_binary_hugo_scores.csv",
  riaz      = "ICI_binary_riaz_scores.csv",
  nsclc     = "ICI_binary_nsclc_scores.csv",
  eac       = "ICI_binary_eac_scores.csv",
  friedrich = "ICI_binary_friedrich_scores.csv")

ID_CANDS   <- c("sample","sample_id","RNA_ID","vst_id","pid_str","pid","ID","patient_id","patient","Sample_ID")
RESP_CANDS <- c("response_bin","pcr_binary","responder","response","Response","BOR","binaryResponse","resp")

read_cohort <- function(key) {
  f <- file.path(GRID, FILES[[key]])
  x <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
  id <- intersect(ID_CANDS, names(x))[1]; if (is.na(id)) id <- names(x)[1]
  x$.id <- as.character(x[[id]])
  rc <- intersect(RESP_CANDS, names(x))[1]
  if (is.na(rc)) stop("no response column in ", basename(f))
  r <- x[[rc]]
  if (is.character(r) || is.factor(r)) {
    r <- toupper(trimws(as.character(r)))
    r <- ifelse(r %in% c("R","CR","PR","RESPONDER","YES","TRUE","1"), 1L,
         ifelse(r %in% c("NR","SD","PD","NON-RESPONDER","NO","FALSE","0"), 0L, NA_integer_)) }
  x$.resp <- as.integer(r)
  lay <- c("Sensory","Integration","Execution")
  armcol <- if (key %in% c("ispy2","gide")) intersect(c("arm","Arm","treatment"), names(x))[1] else
            if (key == "nsclc") intersect(c("cohort","study","dataset"), names(x))[1] else NA
  if (key == "nsclc" && is.na(armcol)) stop("NSCLC file has no study column (cohort/study) - paste this to Claude")
  x$.arm <- if (!is.na(armcol)) as.character(x[[armcol]]) else key       # NSCLC: one panel per study
  out <- x[, c(".id",".resp",".arm", lay)]; out$key <- key
  out <- out[!is.na(out$.resp) & complete.cases(out[, lay]), ]
  cat(sprintf("  %-10s n=%4d  id=%-10s resp=%-12s arms=%s\n", key, nrow(out), id, rc,
              paste(unique(out$.arm), collapse = "/")))
  out
}

cat("=== reading cohort score files ===\n")
all <- do.call(rbind, lapply(names(FILES), read_cohort))

# ---- binarize layers at cohort-arm median -----------------------------------
all <- all %>% group_by(key, .arm) %>%
  mutate(across(c(Sensory, Integration, Execution),
                ~ as.integer(.x >= median(.x, na.rm = TRUE)), .names = "{.col}_hi")) %>%
  ungroup()
all <- left_join(all, META, by = "key")

# arm-specific treatment
ctl <- all$key == "ispy2" & grepl("control", all$.arm, ignore.case = TRUE)
all$target[ctl] <- "None"; all$drug[ctl] <- "Paclitaxel only"
ip  <- all$key == "gide" & grepl("ipi", all$.arm, ignore.case = TRUE)
all$target[ip] <- "PD-1 + CTLA-4"; all$drug[ip] <- "Ipi+Pembro"
n126 <- all$key == "nsclc" & grepl("126044", all$.arm)                  # GSE126044 = nivolumab
all$target[n126] <- "PD-1"; all$drug[n126] <- "Nivolumab"
fr <- all$key == "friedrich"; all$target[fr] <- "BCMA×CD3"; all$drug[fr] <- "Bispecific T-cell engager"
ea <- all$key == "eac"; all$target[ea] <- "PD-L1"; all$drug[ea] <- "Atezolizumab + CRT"

# one-line panel labels, fixed order (by cancer type)
armlab <- ifelse(all$.arm == all$key, "",
          ifelse(ctl, " control", ifelse(ip, " Ipi+Pembro",
          ifelse(all$key == "gide", " Pembro", ifelse(all$key == "ispy2", " pembro", paste0(" ", all$.arm))))))
all$panel0 <- paste0(all$label, armlab)
NAMES <- c("Rose" = "Rose BLCA", "Gide Pembro" = "Gide-Pembro", "Gide Ipi+Pembro" = "Gide-IpiPembro",
           "Hugo" = "Hugo SKCM", "Riaz" = "Riaz SKCM", "Braun" = "Braun ccRCC")
all$panel0 <- ifelse(all$panel0 %in% names(NAMES), NAMES[all$panel0], all$panel0)
ORDER <- c("IMvigor210","Rose BLCA","I-SPY2 control","I-SPY2 pembro","Gide-Pembro","Gide-IpiPembro",
           "Hugo SKCM","Riaz SKCM","NSCLC GSE126044","NSCLC GSE135222","EAC","Braun ccRCC","Friedrich")
ORDER <- c(ORDER, setdiff(unique(all$panel0), ORDER))
nn <- table(all$panel0)
all$panel <- factor(paste0(all$panel0, " (n=", nn[all$panel0], ")"),
                    levels = paste0(ORDER, " (n=", nn[ORDER], ")")[ORDER %in% names(nn)])

all <- all %>% arrange(panel, .resp, Sensory_hi, Integration_hi, Execution_hi) %>%
  group_by(panel) %>% mutate(x = row_number()) %>% ungroup()

# ---- long format ------------------------------------------------------------
rows <- list(
  data.frame(row = "Cancer",      val = all$cancer),
  data.frame(row = "Target",      val = all$target),
  data.frame(row = "Treatment",   val = all$drug),
  data.frame(row = "Sensory",     val = ifelse(all$Sensory_hi == 1, "High", "Low")),
  data.frame(row = "Integration", val = ifelse(all$Integration_hi == 1, "High", "Low")),
  data.frame(row = "Execution",   val = ifelse(all$Execution_hi == 1, "High", "Low")),
  data.frame(row = "Response",    val = ifelse(all$.resp == 1, "Responder", "Non-responder")))
L <- do.call(rbind, lapply(rows, function(r) cbind(r, x = all$x, panel = all$panel)))
L$row <- factor(L$row, levels = rev(c("Cancer","Target","Treatment",
                                      "Sensory","Integration","Execution","Response")))
L$panel <- factor(L$panel, levels = levels(all$panel))

PAL <- c(
  Urothelial="#4C72B0", Breast="#DD8452", Renal="#55A868", Melanoma="#C44E52",
  NSCLC="#8172B3", Esophageal="#937860", Myeloma="#DA8BC3",
  `PD-L1`="#7B9FD4", `PD-1`="#5FA37C", `PD-1 + CTLA-4`="#7E5EA8", `PD-1/PD-L1`="#A8C4E0", None="#D9D9D9",
  Atezolizumab="#E1A3C8", Nivolumab="#D2B48C", Pembrolizumab="#F0C419",
  `Ipi+Pembro`="#E08A5F", `PD-(L)1`="#B9C4CC", `Paclitaxel only`="#D9D9D9",
  `BCMA×CD3`="#B07AA1", `Bispecific T-cell engager`="#8E6C8A", `Atezolizumab + CRT`="#C46FA6",
  High="#2C6FA6", Low="#E8E8E8",
  Responder="#2E7D4F", `Non-responder`="#E8A33D")
BRK <- c("Urothelial","Breast","Melanoma","NSCLC","Esophageal","Renal","Myeloma",
         "PD-L1","PD-1","PD-1 + CTLA-4","PD-1/PD-L1","BCMA×CD3","None",
         "Atezolizumab","Atezolizumab + CRT","Pembrolizumab","Nivolumab","Ipi+Pembro","PD-(L)1",
         "Bispecific T-cell engager","Paclitaxel only",
         "High","Low","Responder","Non-responder")

g <- ggplot(L, aes(x = x, y = row, fill = val)) +
  geom_tile() +
  facet_grid(. ~ panel, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = PAL, breaks = BRK, na.value = "white", name = NULL) +
  labs(title = sprintf("Patient landscape across %d cohort × treatment-arm analyses (n = %s)",
                       nlevels(all$panel), format(nrow(all), big.mark = ",")),
       subtitle = "Each column is one patient. High = at or above the median of its analysis.",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_blank(), panel.grid = element_blank(),
        panel.spacing = unit(1.5, "pt"),
        strip.text.x = element_text(angle = 90, hjust = 0, vjust = 0.5, size = 6.5),
        strip.clip = "off",
        legend.position = "right", legend.key.size = unit(8, "pt"),
        legend.text = element_text(size = 7),
        plot.margin = margin(10, 10, 10, 10))

ggsave(file.path(OUT, "Figure2A_landscape_v2_5.png"), g, width = 15, height = 5.4, dpi = 300)
ggsave(file.path(OUT, "Figure2A_landscape_v2_5.pdf"), g, width = 15, height = 5.4)
print(g)

cat("\n=== coverage by cohort ===\n")
print(as.data.frame(all %>% count(panel, cancer, target, drug)), row.names = FALSE)
cat(sprintf("\nTOTAL: %d patients, %d analyses\nwritten -> %s\n", nrow(all), nlevels(all$panel),
            file.path(OUT, "Figure2A_landscape_v2_5.png")))
