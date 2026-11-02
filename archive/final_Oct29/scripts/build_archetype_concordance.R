# =====================================================================
# ICI 2.0 — Archetype assignment (binary) + continuous concordance
# =====================================================================
# Mirrors prior framework's threshold-based pattern definition.
#   Coordinated         = S+ I+ E+
#   Sensory-Blind       = S- I+ E+
#   Spatially-Excluded  = S+ I- E+
#   Execution-Failed    = S+ I+ E-
#   Multi-Layer Ataxia  = S- I- E-
#   Mixed               = (any other pattern)
#
# Computes archetype under both:
#   (a) Continuous ssGSEA scores median-split per cohort
#   (b) Binary cohort-median binarization scores median-split per cohort
# Compares concordance.
#
# Run:
#   source("/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript/build_archetype_concordance.R")
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(gridExtra)
})

BASE     <- "/Users/dbedi/Desktop/R data FINAL"
MAN_DIR  <- file.path(BASE, "ICI_Manuscript")
GRID_DIR <- file.path(MAN_DIR, "MASTER_GRID")
OUT_DIR  <- file.path(GRID_DIR, "Figures_v2", "archetypes")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("==========================================\n")
cat("Archetype assignment (binary + continuous concordance)\n")
cat("==========================================\n\n")

# ---- Archetype classifier (matches prior framework) ---------------
classify_archetype <- function(S, I, E) {
  S_high <- S >= median(S, na.rm = TRUE)
  I_high <- I >= median(I, na.rm = TRUE)
  E_high <- E >= median(E, na.rm = TRUE)
  dplyr::case_when(
     S_high &  I_high &  E_high  ~ "Coordinated",
    !S_high &  I_high &  E_high  ~ "Sensory-Blind",
     S_high & !I_high &  E_high  ~ "Spatially-Excluded",
     S_high &  I_high & !E_high  ~ "Execution-Failed",
    !S_high & !I_high & !E_high  ~ "Multi-Layer Ataxia",
    TRUE                          ~ "Mixed"
  )
}

ARCH_ORDER <- c("Coordinated", "Spatially-Excluded", "Sensory-Blind",
                "Execution-Failed", "Multi-Layer Ataxia", "Mixed")
ARCH_COLORS <- c(
  "Coordinated"        = "#1B7837",
  "Spatially-Excluded" = "#A23B72",
  "Sensory-Blind"      = "#2E86AB",
  "Execution-Failed"   = "#F18F01",
  "Multi-Layer Ataxia" = "#999999",
  "Mixed"              = "#D5D5D5"
)

# ---- Cohort registry (binary scores) ------------------------------
cohort_specs <- list(
  "IMvigor210"      = list(file = "ICI_binary_imvigor_scores.csv", tumor = "Urothelial",
                           cont_csv = file.path(MAN_DIR, "bulk_scores", "IMvigor210_ICI_scores.csv")),
  "I-SPY2 pembro"   = list(file = "ICI_binary_ispy2_scores.csv", tumor = "Breast",
                           filter = function(d) d[d$arm == "pembro", ],
                           cont_csv = NULL),
  "I-SPY2 control"  = list(file = "ICI_binary_ispy2_scores.csv", tumor = "Breast",
                           filter = function(d) d[d$arm == "control", ],
                           cont_csv = NULL),
  "Braun ccRCC"     = list(file = "ICI_binary_braun_scores.csv", tumor = "Renal",
                           cont_csv = file.path(MAN_DIR, "bulk_scores",
                                                 "Braun_ccRCC_ICI_scores.csv")),
  "Rose BLCA"       = list(file = "ICI_binary_rose_scores.csv", tumor = "Urothelial"),
  "Gide-Pembro"     = list(file = "ICI_binary_gide_scores.csv", tumor = "Melanoma",
                           filter = function(d) d[d$arm == "Pembro", ]),
  "Gide-IpiPembro"  = list(file = "ICI_binary_gide_scores.csv", tumor = "Melanoma",
                           filter = function(d) d[d$arm == "Ipi+Pembro", ]),
  "Hugo SKCM"       = list(file = "ICI_binary_hugo_scores.csv", tumor = "Melanoma"),
  "Riaz SKCM"       = list(file = "ICI_binary_riaz_scores.csv", tumor = "Melanoma"),
  "NSCLC"           = list(file = "ICI_binary_nsclc_scores.csv", tumor = "NSCLC",
                           cont_csv = file.path(MAN_DIR, "bulk_scores",
                                                 "NSCLC_ICI_scores.csv")),
  "EAC"             = list(file = "ICI_binary_eac_scores.csv", tumor = "Esophageal"),
  "Friedrich"       = list(file = "ICI_binary_friedrich_scores.csv", tumor = "Multiple Myeloma")
)

# ---- Build archetype table per cohort -----------------------------
all_archetypes <- list()

for (cn in names(cohort_specs)) {
  spec <- cohort_specs[[cn]]

  # BINARY assignment
  d_bin <- read.csv(file.path(GRID_DIR, spec$file), stringsAsFactors = FALSE)
  if (!is.null(spec$filter)) d_bin <- spec$filter(d_bin)

  resp_col <- intersect(c("response_bin", "pcr_binary"), colnames(d_bin))[1]
  if (is.na(resp_col)) next
  d_bin$resp <- d_bin[[resp_col]]
  d_bin <- d_bin[!is.na(d_bin$resp), ]

  d_bin$archetype_binary <- classify_archetype(d_bin$Sensory,
                                                d_bin$Integration,
                                                d_bin$Execution)
  d_bin$cohort <- cn
  d_bin$tumor  <- spec$tumor

  # Find a sample ID column for matching to continuous later
  id_col <- intersect(c("sample_id", "vst_id", "sample", "patient_id", "RNA_ID"),
                      colnames(d_bin))[1]
  if (!is.na(id_col)) d_bin$id <- d_bin[[id_col]] else d_bin$id <- seq_len(nrow(d_bin))

  cat(sprintf("[%s, n=%d] binary archetype distribution:\n",
              cn, nrow(d_bin)))
  print(table(d_bin$archetype_binary))

  all_archetypes[[cn]] <- data.frame(
    cohort = cn, tumor = spec$tumor,
    id = d_bin$id, resp = d_bin$resp,
    archetype_binary = d_bin$archetype_binary,
    Sensory_bin = d_bin$Sensory,
    Integration_bin = d_bin$Integration,
    Execution_bin = d_bin$Execution,
    stringsAsFactors = FALSE
  )
}

archetype_df <- do.call(rbind, all_archetypes)
archetype_df$archetype_binary <- factor(archetype_df$archetype_binary,
                                         levels = ARCH_ORDER)

cat(sprintf("\n\nTotal patients across 12 cohorts: %d\n", nrow(archetype_df)))
cat("Overall binary archetype distribution:\n")
print(table(archetype_df$archetype_binary))

# ---- Compute archetype-level response rates -----------------------
archetype_summary <- archetype_df %>%
  dplyr::group_by(archetype_binary) %>%
  dplyr::summarise(
    n            = dplyr::n(),
    n_responders = sum(resp, na.rm = TRUE),
    response_rate = mean(resp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(response_rate))

cat("\nArchetype response rates (across all 12 cohorts):\n")
print(archetype_summary)

write.csv(archetype_df, file.path(GRID_DIR, "archetype_assignments_binary.csv"),
          row.names = FALSE)
write.csv(archetype_summary, file.path(GRID_DIR, "archetype_summary_binary.csv"),
          row.names = FALSE)

# ---- Plot 1: Response rate by archetype ---------------------------
p1 <- ggplot(archetype_summary,
             aes(x = reorder(archetype_binary, response_rate),
                 y = response_rate, fill = archetype_binary)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(aes(label = sprintf("%.0f%%\n(%d/%d)",
                                100 * response_rate,
                                n_responders, n)),
            vjust = -0.3, size = 3) +
  scale_fill_manual(values = ARCH_COLORS, guide = "none") +
  scale_y_continuous(labels = scales::percent_format(),
                     limits = c(0, max(archetype_summary$response_rate) * 1.25),
                     expand = c(0, 0)) +
  coord_flip() +
  labs(title = "Response rate by binary archetype",
       subtitle = "All 12 cohorts pooled (n = 971)",
       x = NULL, y = "Response rate") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "archetype_response_rate_binary.pdf"), p1,
       width = 8, height = 5)
ggsave(file.path(OUT_DIR, "archetype_response_rate_binary.png"), p1,
       width = 8, height = 5, dpi = 200)

# ---- Plot 2: Archetype distribution by tumor type ----------------
arch_by_tumor <- archetype_df %>%
  dplyr::group_by(tumor, archetype_binary) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(tumor) %>%
  dplyr::mutate(pct = n / sum(n)) %>%
  dplyr::ungroup()

tumor_order <- archetype_df %>%
  dplyr::group_by(tumor) %>%
  dplyr::summarise(n_total = dplyr::n(), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(n_total)) %>%
  dplyr::pull(tumor)

arch_by_tumor$tumor <- factor(arch_by_tumor$tumor, levels = rev(tumor_order))
arch_by_tumor$archetype_binary <- factor(arch_by_tumor$archetype_binary,
                                          levels = ARCH_ORDER)

p2 <- ggplot(arch_by_tumor,
             aes(x = pct, y = tumor, fill = archetype_binary)) +
  geom_col(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = ARCH_COLORS, name = "Archetype") +
  scale_x_continuous(labels = scales::percent_format(),
                     expand = c(0, 0)) +
  labs(title = "Archetype composition by tumor type (binary)",
       subtitle = "Each tumor type's distinct cascade signature",
       x = "% of patients", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(face = "bold"),
        legend.position = "right",
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "archetype_by_tumor_binary.pdf"), p2,
       width = 10, height = 5.5)
ggsave(file.path(OUT_DIR, "archetype_by_tumor_binary.png"), p2,
       width = 10, height = 5.5, dpi = 200)

# ---- Plot 3: Archetype × response heatmap by cohort --------------
arch_cohort_resp <- archetype_df %>%
  dplyr::group_by(cohort, archetype_binary) %>%
  dplyr::summarise(
    n = dplyr::n(),
    response_rate = mean(resp, na.rm = TRUE),
    .groups = "drop"
  )

# Cohort order from manuscript
cohort_order_lvl <- c("IMvigor210", "I-SPY2 pembro", "I-SPY2 control", "Braun ccRCC",
                       "Rose BLCA", "Gide-Pembro", "Gide-IpiPembro",
                       "Hugo SKCM", "Riaz SKCM", "NSCLC", "EAC", "Friedrich")
arch_cohort_resp$cohort <- factor(arch_cohort_resp$cohort,
                                   levels = rev(cohort_order_lvl))
arch_cohort_resp$archetype_binary <- factor(arch_cohort_resp$archetype_binary,
                                             levels = ARCH_ORDER)

p3 <- ggplot(arch_cohort_resp,
             aes(x = archetype_binary, y = cohort, fill = response_rate)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%\n(n=%d)",
                                100 * response_rate, n)),
            size = 2.5, color = "black") +
  scale_fill_gradient2(low = "#3B7AB8", mid = "white", high = "#C8443A",
                       midpoint = 0.3, limits = c(0, 1),
                       labels = scales::percent_format(),
                       name = "Response\nrate") +
  scale_x_discrete(position = "top") +
  labs(title = "Archetype × cohort response heatmap (binary)",
       subtitle = "Response rate (and n) per archetype × cohort cell",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(panel.grid = element_blank(),
        axis.text.x.top = element_text(angle = 35, hjust = 0, face = "bold"),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "archetype_x_cohort_response_binary.pdf"), p3,
       width = 11, height = 7)
ggsave(file.path(OUT_DIR, "archetype_x_cohort_response_binary.png"), p3,
       width = 11, height = 7, dpi = 200)

# ---- Continuous-vs-binary concordance (where continuous available)
cat("\n==========================================\n")
cat("CONTINUOUS vs BINARY archetype concordance\n")
cat("==========================================\n")

concordance_results <- list()

for (cn in names(cohort_specs)) {
  spec <- cohort_specs[[cn]]
  if (is.null(spec$cont_csv) || !file.exists(spec$cont_csv)) next

  cont <- read.csv(spec$cont_csv, stringsAsFactors = FALSE, check.names = FALSE)

  # Find S/I/E continuous columns and an ID
  s_col <- intersect(c("Layer1_Sensory", "Sensory"), colnames(cont))[1]
  i_col <- intersect(c("Layer2_Integration", "Integration"), colnames(cont))[1]
  e_col <- intersect(c("Layer3_Execution", "Execution"), colnames(cont))[1]
  if (any(is.na(c(s_col, i_col, e_col)))) next

  cont_id <- intersect(c("sample_id", "vst_id", "sample", "patient_id",
                          "RNA_ID", "expr_id"), colnames(cont))[1]
  if (is.na(cont_id)) next

  cont$archetype_continuous <- classify_archetype(
    cont[[s_col]], cont[[i_col]], cont[[e_col]]
  )
  cont$id <- cont[[cont_id]]

  # Pull binary archetype labels
  bin_subset <- archetype_df[archetype_df$cohort == cn, ]
  merged <- merge(bin_subset, cont[, c("id", "archetype_continuous")],
                  by = "id", all.x = FALSE)

  if (nrow(merged) < 5) next

  # Concordance: same archetype label between continuous and binary?
  match_pct <- mean(merged$archetype_binary == merged$archetype_continuous,
                    na.rm = TRUE) * 100
  cat(sprintf("[%s, n_match=%d] continuous-binary archetype agreement: %.1f%%\n",
              cn, nrow(merged), match_pct))

  concordance_results[[cn]] <- merged
}

# Combined concordance
if (length(concordance_results) > 0) {
  conc_df <- dplyr::bind_rows(concordance_results)

  # Confusion matrix
  conf_tbl <- table(continuous = conc_df$archetype_continuous,
                    binary     = conc_df$archetype_binary)
  cat("\nConfusion matrix (rows = continuous, cols = binary):\n")
  print(conf_tbl)

  conf_df <- as.data.frame(conf_tbl)
  conf_df$continuous <- factor(conf_df$continuous, levels = ARCH_ORDER)
  conf_df$binary     <- factor(conf_df$binary,     levels = ARCH_ORDER)
  conf_df$pct <- conf_df$Freq /
    ave(conf_df$Freq, conf_df$continuous, FUN = sum)

  p4 <- ggplot(conf_df, aes(x = binary, y = continuous, fill = pct)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(Freq > 0,
                                 sprintf("%d\n(%.0f%%)", Freq, 100 * pct),
                                 "")),
              size = 2.5, color = "black") +
    scale_fill_gradient(low = "white", high = "#3B7AB8",
                        labels = scales::percent_format(),
                        name = "% of\ncontinuous\nclass") +
    labs(title = "Continuous vs binary archetype concordance",
         subtitle = sprintf("%d patients across cohorts with both score types; overall agreement: %.1f%%",
                            nrow(conc_df),
                            mean(conc_df$archetype_continuous == conc_df$archetype_binary,
                                 na.rm = TRUE) * 100),
         x = "Binary archetype",
         y = "Continuous archetype") +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 30, hjust = 1, face = "bold"),
          axis.text.y = element_text(face = "bold"),
          plot.title = element_text(face = "bold"))

  ggsave(file.path(OUT_DIR, "archetype_continuous_vs_binary_confusion.pdf"), p4,
         width = 8, height = 7)
  ggsave(file.path(OUT_DIR, "archetype_continuous_vs_binary_confusion.png"), p4,
         width = 8, height = 7, dpi = 200)

  write.csv(conc_df, file.path(GRID_DIR, "archetype_concordance.csv"),
            row.names = FALSE)
}

cat("\n==========================================\n")
cat("DONE. Files in: ", OUT_DIR, "\n")
list.files(OUT_DIR, pattern = "\\.pdf$") |> sapply(function(f) cat("  -", f, "\n"))
