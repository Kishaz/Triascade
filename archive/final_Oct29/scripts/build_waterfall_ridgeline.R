# =====================================================================
# ICI 2.0 — Cross-cohort waterfall and ridgeline (binary)
# =====================================================================
# Visualizes:
#   1. Per-cohort R - NR delta (binary composite score) — waterfall
#   2. Binary composite score density per cohort — ridgelines
#
# Run:
#   source("/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript/build_waterfall_ridgeline.R")
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  if (!requireNamespace("ggridges", quietly = TRUE)) {
    install.packages("ggridges", repos = "https://cloud.r-project.org")
  }
  library(ggridges)
})

BASE     <- "/Users/dbedi/Desktop/R data FINAL"
MAN_DIR  <- file.path(BASE, "ICI_Manuscript")
GRID_DIR <- file.path(MAN_DIR, "MASTER_GRID")
OUT_DIR  <- file.path(GRID_DIR, "Figures_v2", "waterfall_ridgeline")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("==========================================\n")
cat("Cross-cohort waterfall + ridgeline (binary)\n")
cat("==========================================\n\n")

# ---- Cohort registry ---------------------------------------------
cohort_specs <- list(
  "IMvigor210"      = list(file = "ICI_binary_imvigor_scores.csv", tumor = "Urothelial"),
  "I-SPY2 pembro"   = list(file = "ICI_binary_ispy2_scores.csv", tumor = "Breast",
                           filter = function(d) d[d$arm == "pembro", ]),
  "I-SPY2 control"  = list(file = "ICI_binary_ispy2_scores.csv", tumor = "Breast",
                           filter = function(d) d[d$arm == "control", ]),
  "Braun ccRCC"     = list(file = "ICI_binary_braun_scores.csv", tumor = "Renal"),
  "Rose BLCA"       = list(file = "ICI_binary_rose_scores.csv", tumor = "Urothelial"),
  "Gide-Pembro"     = list(file = "ICI_binary_gide_scores.csv", tumor = "Melanoma",
                           filter = function(d) d[d$arm == "Pembro", ]),
  "Gide-IpiPembro"  = list(file = "ICI_binary_gide_scores.csv", tumor = "Melanoma",
                           filter = function(d) d[d$arm == "Ipi+Pembro", ]),
  "Hugo SKCM"       = list(file = "ICI_binary_hugo_scores.csv", tumor = "Melanoma"),
  "Riaz SKCM"       = list(file = "ICI_binary_riaz_scores.csv", tumor = "Melanoma"),
  "NSCLC"           = list(file = "ICI_binary_nsclc_scores.csv", tumor = "NSCLC"),
  "EAC"             = list(file = "ICI_binary_eac_scores.csv", tumor = "Esophageal"),
  "Friedrich"       = list(file = "ICI_binary_friedrich_scores.csv", tumor = "Multiple Myeloma")
)

# ---- Load + assemble pooled long-format binary scores -----------
all_long <- list()
delta_summary <- list()

for (cn in names(cohort_specs)) {
  spec <- cohort_specs[[cn]]
  d <- read.csv(file.path(GRID_DIR, spec$file), stringsAsFactors = FALSE)
  if (!is.null(spec$filter)) d <- spec$filter(d)
  resp_col <- intersect(c("response_bin", "pcr_binary"), colnames(d))[1]
  if (is.na(resp_col)) next
  d$resp <- d[[resp_col]]
  d <- d[!is.na(d$resp) & d$resp %in% c(0, 1), ]
  if (nrow(d) < 5) next

  long <- data.frame(
    cohort = cn,
    tumor  = spec$tumor,
    response = factor(ifelse(d$resp == 1, "R", "NR"), levels = c("NR", "R")),
    Sensory     = d$Sensory,
    Integration = d$Integration,
    Execution   = d$Execution,
    Composite   = d$composite,
    stringsAsFactors = FALSE
  )
  all_long[[cn]] <- long

  # Compute per-cohort R-NR delta for each layer
  for (lay in c("Sensory", "Integration", "Execution", "Composite")) {
    v <- long[[lay]]
    delta <- mean(v[long$response == "R"], na.rm = TRUE) -
             mean(v[long$response == "NR"], na.rm = TRUE)
    sd_pool <- sd(v, na.rm = TRUE)
    cohens_d <- if (sd_pool > 0) delta / sd_pool else NA
    p_val <- tryCatch(wilcox.test(v ~ long$response)$p.value, error = function(e) NA)
    delta_summary[[length(delta_summary) + 1]] <- data.frame(
      cohort = cn, tumor = spec$tumor, layer = lay,
      n = nrow(long), n_R = sum(long$response == "R"), n_NR = sum(long$response == "NR"),
      delta = delta, cohens_d = cohens_d, p = p_val,
      stringsAsFactors = FALSE
    )
  }
  cat(sprintf("[%s] n=%d (R=%d, NR=%d) ✓\n", cn, nrow(long),
              sum(long$response == "R"), sum(long$response == "NR")))
}

pooled <- do.call(rbind, all_long)
delta_df <- do.call(rbind, delta_summary)

cat(sprintf("\nTotal pooled: %d patients across %d cohorts\n",
            nrow(pooled), length(unique(pooled$cohort))))

# Save underlying data
write.csv(delta_df, file.path(GRID_DIR, "waterfall_delta_summary.csv"), row.names = FALSE)

# ============================================================
# PANEL: Waterfall — per-cohort R-NR delta per layer
# ============================================================

cohort_order <- delta_df %>%
  filter(layer == "Composite") %>%
  arrange(delta) %>%
  pull(cohort)

delta_df$cohort <- factor(delta_df$cohort, levels = cohort_order)
delta_df$layer <- factor(delta_df$layer,
                          levels = c("Sensory", "Integration", "Execution", "Composite"))
delta_df$sig <- delta_df$p < 0.05

LAYER_COLORS <- c(
  "Sensory"     = "#2E86AB",
  "Integration" = "#A23B72",
  "Execution"   = "#F18F01",
  "Composite"   = "#3B7AB8"
)

p_water <- ggplot(delta_df, aes(x = cohort, y = delta, fill = layer)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_col(aes(alpha = sig), color = "black", linewidth = 0.3) +
  scale_alpha_manual(values = c(`TRUE` = 1.0, `FALSE` = 0.4),
                     guide = "none") +
  scale_fill_manual(values = LAYER_COLORS, guide = "none") +
  facet_wrap(~ layer, ncol = 2, scales = "free_x") +
  coord_flip() +
  geom_text(aes(label = sprintf("%s%s",
                                ifelse(p < 0.001, "***",
                                       ifelse(p < 0.01, "**",
                                              ifelse(p < 0.05, "*", ""))),
                                ifelse(is.na(p), "", ""))),
            hjust = ifelse(delta_df$delta > 0, -0.2, 1.2), size = 3) +
  labs(title = "Per-cohort R−NR delta in binary layer scores",
       subtitle = "Solid = p < 0.05; faded = NS. Cohorts ordered by Composite delta.",
       y = "Mean(R) − Mean(NR), binary score (raw delta)",
       x = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey92", color = NA),
        axis.text.y = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "waterfall_RminusNR_binary.pdf"), p_water,
       width = 12, height = 9)
ggsave(file.path(OUT_DIR, "waterfall_RminusNR_binary.png"), p_water,
       width = 12, height = 9, dpi = 200)
cat("\nWrote: waterfall_RminusNR_binary.pdf\n")

# ============================================================
# PANEL: Cohen's d effect size waterfall (standardized)
# ============================================================
p_cohen <- ggplot(delta_df, aes(x = cohort, y = cohens_d, fill = layer)) +
  geom_hline(yintercept = 0, color = "grey50", linewidth = 0.4) +
  geom_col(aes(alpha = sig), color = "black", linewidth = 0.3) +
  scale_alpha_manual(values = c(`TRUE` = 1.0, `FALSE` = 0.4), guide = "none") +
  scale_fill_manual(values = LAYER_COLORS, guide = "none") +
  facet_wrap(~ layer, ncol = 2, scales = "free_x") +
  coord_flip() +
  labs(title = "Per-cohort R−NR effect size (Cohen's d, binary scoring)",
       subtitle = "Standardized effect size; |d| > 0.8 = large effect (dashed lines)",
       y = "Cohen's d (mean(R) − mean(NR) / sd_pooled)",
       x = NULL) +
  geom_hline(yintercept = c(0.8, -0.8), linetype = "dashed", color = "grey60",
             linewidth = 0.3) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey92", color = NA),
        axis.text.y = element_text(face = "bold", size = 9),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "waterfall_cohens_d_binary.pdf"), p_cohen,
       width = 12, height = 9)
ggsave(file.path(OUT_DIR, "waterfall_cohens_d_binary.png"), p_cohen,
       width = 12, height = 9, dpi = 200)
cat("Wrote: waterfall_cohens_d_binary.pdf\n")

# ============================================================
# PANEL: Ridgelines — binary Composite distribution per cohort
# ============================================================
ridge_df <- pooled |>
  pivot_longer(cols = c(Sensory, Integration, Execution, Composite),
               names_to = "layer", values_to = "score") |>
  mutate(layer = factor(layer, levels = c("Sensory", "Integration", "Execution", "Composite")))

# Within-cohort z-score for layer comparability across cohorts
ridge_df <- ridge_df |>
  group_by(cohort, layer) |>
  mutate(score_z = scale(score)[, 1]) |>
  ungroup()

# Order cohorts by sample size
cohort_n_order <- pooled |>
  group_by(cohort) |>
  summarise(n = n(), .groups = "drop") |>
  arrange(n) |>
  pull(cohort)
ridge_df$cohort <- factor(ridge_df$cohort, levels = cohort_n_order)

p_ridge <- ggplot(ridge_df, aes(x = score_z, y = cohort, fill = response,
                                 color = response)) +
  geom_density_ridges(alpha = 0.5, scale = 1.2, linewidth = 0.4,
                      rel_min_height = 0.005) +
  facet_wrap(~ layer, ncol = 2) +
  scale_fill_manual(values = c("NR" = "steelblue", "R" = "darkred"),
                    name = "Response") +
  scale_color_manual(values = c("NR" = "steelblue", "R" = "darkred"),
                     guide = "none") +
  labs(title = "Layer score distributions by cohort and response (binary, z-scored within cohort)",
       subtitle = "R distributions shifted right = layer is response-discriminating",
       x = "Binary layer score (z-scored within cohort)",
       y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.background = element_rect(fill = "grey92", color = NA),
        axis.text.y = element_text(face = "bold", size = 9),
        legend.position = "top",
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "ridgeline_distributions_binary.pdf"), p_ridge,
       width = 13, height = 11)
ggsave(file.path(OUT_DIR, "ridgeline_distributions_binary.png"), p_ridge,
       width = 13, height = 11, dpi = 200)
cat("Wrote: ridgeline_distributions_binary.pdf\n")

# Summary stats
cat("\n==========================================\n")
cat("PER-COHORT EFFECT SIZE SUMMARY (Composite layer)\n")
cat("==========================================\n")
print(delta_df |> filter(layer == "Composite") |>
        select(cohort, n, n_R, delta, cohens_d, p) |>
        arrange(desc(cohens_d)))

cat("\n==========================================\n")
cat("DONE. Files in:\n", OUT_DIR, "\n")
list.files(OUT_DIR, pattern = "\\.pdf$") |>
  sapply(function(f) cat("  -", f, "\n"))
