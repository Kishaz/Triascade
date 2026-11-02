# =====================================================================
# ICI 2.0 — scRNA cell-state decomposition (Fig 6, binary)
# =====================================================================
# Single-cell T-cell decomposition of Execution layer in melanoma.
#   2,182 T cells, 15 patients (5 NR / 10 R), 3 cell-state clusters
#
# Validates that:
#   1. Execution layer score correlates with T-cell activation phenotype
#   2. Execution-positive cell fractions differ between R vs NR
#   3. Cell-state clusters show distinct layer dominance
#
# Run:
#   source("/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript/build_scRNA_celltype.R")
# =====================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(grid)
  library(gridExtra)
})

BASE     <- "/Users/dbedi/Desktop/R data FINAL"
MAN_DIR  <- file.path(BASE, "ICI_Manuscript")
GRID_DIR <- file.path(MAN_DIR, "MASTER_GRID")
OUT_DIR  <- file.path(GRID_DIR, "Figures_v2", "scRNA")

cat("==========================================\n")
cat("scRNA T-cell-state Execution decomposition\n")
cat("==========================================\n\n")

LAYER_COLORS <- c(
  "Sensory"     = "#2E86AB",
  "Integration" = "#A23B72",
  "Execution"   = "#F18F01"
)

# ---- Load Seurat object ------------------------------------------
cat("Loading Seurat object...\n")
sc <- readRDS(file.path(BASE, "CereNet_Validation_Seurat_TumorOnly.rds"))
cat(sprintf("  %d cells × %d genes, %d patients, %d clusters\n",
            ncol(sc), nrow(sc),
            length(unique(sc$patient_id)),
            length(unique(sc$cluster))))

# ---- Define ICI panels (matched to bulk) ------------------------
SENSORY <- c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G","B2M",
             "TAP1","TAP2","TAPBP","PSMB8","PSMB9","PSMB10",
             "IRF1","STAT1","STAT2","JAK1","JAK2","IFNG","IFNGR1","IFNGR2")
INTPOS  <- c("CXCL9","CXCL10","CXCL11","CXCL13","CCL5","CCL4","CCL19","CCL21",
             "ICAM1","VCAM1","MADCAM1","CXCR3","CCR5","CCR7",
             "CD80","CD86","CD40","BATF3","IRF8","CR2","TNFRSF17")
INTNEG  <- c("VEGFA","VEGFB","VEGFC","FLT1","KDR","TGFB1","TGFB2","TGFB3",
             "COL1A1","COL1A2","COL3A1","COL4A1","FAP","ACTA2","FN1",
             "MMP2","MMP9","MMP14","CXCL12","EDNRB")
EXECPOS <- c("GZMA","GZMB","GZMH","GZMK","GZMM","PRF1","GNLY","NKG7","FASLG",
             "IFNG","TNF","LTA","CD8A","CD8B","TBX21","EOMES","RUNX3","IL2","IL12A")
EXECNEG <- c("PDCD1","CTLA4","LAG3","HAVCR2","TIGIT","BTLA","VSIR","CD274",
             "PDCD1LG2","CD160","CD244","TOX","NR4A1","NR4A2","NR4A3",
             "ENTPD1","IL10","TGFB1","IDO1")

# ---- Extract expression matrix -----------------------------------
# Get default layer (counts or data depending on Seurat version)
expr <- tryCatch(GetAssayData(sc, layer = "data"),
                 error = function(e) GetAssayData(sc, slot = "data"))
expr <- as.matrix(expr)
cat(sprintf("Expression matrix: %d genes × %d cells\n", nrow(expr), ncol(expr)))
cat(sprintf("Range: %.2f to %.2f (log-normalized)\n",
            min(expr, na.rm = TRUE), max(expr, na.rm = TRUE)))

# Coverage check
all_panel <- unique(c(SENSORY, INTPOS, INTNEG, EXECPOS, EXECNEG))
in_panel <- intersect(all_panel, rownames(expr))
cat(sprintf("Panel coverage: %d / %d genes\n", length(in_panel), length(all_panel)))

# ---- Cell-level binarization at cell-population median ----------
binarize_gene <- function(g, m) {
  if (!g %in% rownames(m)) return(NULL)
  x <- as.numeric(m[g, ])
  thresh <- median(x, na.rm = TRUE)
  as.integer(x >= thresh)
}

sum_binary <- function(gs, m) {
  total <- rep(0L, ncol(m)); hits <- 0
  for (g in gs) {
    b <- binarize_gene(g, m)
    if (!is.null(b)) { total <- total + b; hits <- hits + 1 }
  }
  list(scores = total, n_used = hits)
}

cat("Binarizing genes at cell-population medians...\n")
S_res     <- sum_binary(SENSORY, expr)
IntPos_r  <- sum_binary(INTPOS,  expr)
IntNeg_r  <- sum_binary(INTNEG,  expr)
ExecPos_r <- sum_binary(EXECPOS, expr)
ExecNeg_r <- sum_binary(EXECNEG, expr)
cat(sprintf("  S=%d/21, IntP=%d/21, IntN=%d/20, ExecP=%d/19, ExecN=%d/19\n",
            S_res$n_used, IntPos_r$n_used, IntNeg_r$n_used,
            ExecPos_r$n_used, ExecNeg_r$n_used))

# Per-cell layer scores
sc$Sensory_bin     <- S_res$scores
sc$Integration_bin <- IntPos_r$scores - IntNeg_r$scores
sc$Execution_bin   <- ExecPos_r$scores - ExecNeg_r$scores

# Z-score per layer for composite
z <- function(x) {
  s <- sd(x, na.rm = TRUE); if (s == 0) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / s
}
sc$composite_bin <- (z(sc$Sensory_bin) + z(sc$Integration_bin) + z(sc$Execution_bin)) / 3

# ---- Validate Execution layer correlates with prior activation -
cor_act <- cor(sc$Execution_bin, sc$activation, method = "spearman")
cor_exh <- cor(sc$Execution_bin, sc$exhaustion, method = "spearman")
cor_net <- cor(sc$Execution_bin, sc$net_activation, method = "spearman")
cat("\nExecution layer (binary) vs prior T-cell phenotype scores:\n")
cat(sprintf("  vs activation:     ρ = %.3f\n", cor_act))
cat(sprintf("  vs exhaustion:     ρ = %.3f\n", cor_exh))
cat(sprintf("  vs net_activation: ρ = %.3f\n", cor_net))

# ---- Cell-state cluster summary ---------------------------------
cluster_summary <- as.data.frame(sc@meta.data) |>
  group_by(cluster) |>
  summarise(
    n_cells       = n(),
    Sensory_mean      = mean(Sensory_bin, na.rm = TRUE),
    Integration_mean  = mean(Integration_bin, na.rm = TRUE),
    Execution_mean    = mean(Execution_bin, na.rm = TRUE),
    activation_mean   = mean(activation, na.rm = TRUE),
    exhaustion_mean   = mean(exhaustion, na.rm = TRUE),
    .groups = "drop"
  )
cat("\nPer-cluster layer score means:\n")
print(cluster_summary)

# ---- Patient-level aggregation: % Execution-positive cells -----
patient_summary <- as.data.frame(sc@meta.data) |>
  group_by(patient_id, response) |>
  summarise(
    n_cells = n(),
    Sensory_pct_high     = mean(Sensory_bin >= median(sc$Sensory_bin, na.rm = TRUE)) * 100,
    Integration_pct_high = mean(Integration_bin >= median(sc$Integration_bin, na.rm = TRUE)) * 100,
    Execution_pct_high   = mean(Execution_bin >= median(sc$Execution_bin, na.rm = TRUE)) * 100,
    Execution_mean       = mean(Execution_bin, na.rm = TRUE),
    .groups = "drop"
  )
cat("\nPer-patient summary:\n")
print(patient_summary)

# Wilcoxon R vs NR
for (col in c("Sensory_pct_high", "Integration_pct_high", "Execution_pct_high",
              "Execution_mean")) {
  v_R  <- patient_summary[[col]][patient_summary$response == "R"]
  v_NR <- patient_summary[[col]][patient_summary$response == "NR"]
  if (length(v_R) > 0 && length(v_NR) > 0 &&
      sd(v_R, na.rm = TRUE) > 0 && sd(v_NR, na.rm = TRUE) > 0) {
    p <- tryCatch(wilcox.test(patient_summary[[col]] ~ patient_summary$response)$p.value,
                  error = function(e) NA)
    cat(sprintf("  %s: R median = %.1f, NR median = %.1f, p = %.3f\n",
                col, median(v_R, na.rm = TRUE), median(v_NR, na.rm = TRUE), p))
  }
}

# =================================================================
# PANEL A — Execution layer score correlates with T-cell activation
# =================================================================
df_cells <- data.frame(
  cell        = colnames(expr),
  Execution   = sc$Execution_bin,
  activation  = sc$activation,
  exhaustion  = sc$exhaustion,
  net_act     = sc$net_activation,
  cluster     = factor(sc$cluster),
  patient     = sc$patient_id,
  response    = sc$response,
  stringsAsFactors = FALSE
)

p_A <- ggplot(df_cells, aes(x = activation, y = Execution)) +
  geom_jitter(aes(color = response), width = 0.05, alpha = 0.4, size = 0.6) +
  geom_smooth(method = "lm", se = TRUE, color = "darkred", linewidth = 0.7,
              alpha = 0.2) +
  scale_color_manual(values = c("R" = "darkred", "NR" = "steelblue"),
                     name = NULL) +
  annotate("text", x = max(df_cells$activation, na.rm = TRUE) * 0.9,
           y = min(df_cells$Execution, na.rm = TRUE) + 1,
           label = sprintf("Spearman ρ = %.2f", cor_act),
           hjust = 1, fontface = "bold", size = 3.5) +
  labs(title = "A. Execution layer (binary) tracks T-cell activation phenotype",
       subtitle = sprintf("%d single T cells, melanoma (Sade-Feldman)", ncol(expr)),
       x = "T-cell activation score (curated phenotype)",
       y = "Execution layer score (ICI binary)") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold"),
        legend.position = "top")

ggsave(file.path(OUT_DIR, "Fig6A_Execution_vs_activation.pdf"), p_A,
       width = 7.5, height = 6.5)
ggsave(file.path(OUT_DIR, "Fig6A_Execution_vs_activation.png"), p_A,
       width = 7.5, height = 6.5, dpi = 200)
cat("\nWrote: Fig6A_Execution_vs_activation.pdf\n")

# =================================================================
# PANEL B — Cell-state cluster signature heatmap
# =================================================================
cluster_long <- cluster_summary |>
  pivot_longer(cols = c(Sensory_mean, Integration_mean, Execution_mean,
                        activation_mean, exhaustion_mean),
               names_to = "Metric", values_to = "value") |>
  mutate(Metric = factor(Metric,
                          levels = c("Sensory_mean", "Integration_mean", "Execution_mean",
                                     "activation_mean", "exhaustion_mean"),
                          labels = c("Sensory (binary)", "Integration (binary)",
                                     "Execution (binary)", "Activation (curated)",
                                     "Exhaustion (curated)")),
         cluster_label = sprintf("Cluster %s\n(n=%d)", cluster, n_cells))

# Z-score each metric across clusters for visualization
cluster_long_z <- cluster_long |>
  group_by(Metric) |>
  mutate(value_z = scale(value)[, 1]) |>
  ungroup()

p_B <- ggplot(cluster_long_z, aes(x = cluster_label, y = Metric, fill = value_z)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", value)), size = 3) +
  scale_fill_gradient2(low = "#3B7AB8", mid = "white", high = "#C8443A",
                       midpoint = 0,
                       name = "z-score\n(across\nclusters)") +
  labs(title = "B. T-cell cluster signatures (binary ICI layers vs curated phenotypes)",
       subtitle = "Cell-state clusters show distinct layer dominance",
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank(),
        axis.text.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold"))

ggsave(file.path(OUT_DIR, "Fig6B_cluster_signatures.pdf"), p_B,
       width = 8, height = 5)
ggsave(file.path(OUT_DIR, "Fig6B_cluster_signatures.png"), p_B,
       width = 8, height = 5, dpi = 200)
cat("Wrote: Fig6B_cluster_signatures.pdf\n")

# =================================================================
# PANEL C — Patient-level Execution-high % by R/NR
# =================================================================
p_C <- ggplot(patient_summary,
              aes(x = response, y = Execution_pct_high, fill = response)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 2) +
  geom_text(aes(label = patient_id), nudge_x = 0.15, size = 2.3, color = "grey40") +
  scale_fill_manual(values = c("R" = "darkred", "NR" = "steelblue"),
                    guide = "none") +
  labs(title = sprintf("C. Patient-level Execution-high cell fraction by response\n(n_R=%d, n_NR=%d patients)",
                       sum(patient_summary$response == "R"),
                       sum(patient_summary$response == "NR")),
       subtitle = "Each point = 1 patient; y-axis = % of cells above population median Execution score",
       x = "Response", y = "% Execution-positive cells") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = 11))

ggsave(file.path(OUT_DIR, "Fig6C_patient_ExecutionPct.pdf"), p_C,
       width = 7, height = 6)
ggsave(file.path(OUT_DIR, "Fig6C_patient_ExecutionPct.png"), p_C,
       width = 7, height = 6, dpi = 200)
cat("Wrote: Fig6C_patient_ExecutionPct.pdf\n")

# Save underlying data
write.csv(patient_summary, file.path(GRID_DIR, "scRNA_celltype_patient_summary.csv"),
          row.names = FALSE)
write.csv(cluster_summary, file.path(GRID_DIR, "scRNA_celltype_cluster_summary.csv"),
          row.names = FALSE)
write.csv(df_cells, file.path(GRID_DIR, "scRNA_celltype_per_cell.csv"),
          row.names = FALSE)

cat("\n==========================================\n")
cat("DONE. Files in:", OUT_DIR, "\n")
list.files(OUT_DIR, pattern = "Fig6[ABC].*\\.pdf$") |>
  sapply(function(f) cat("  -", f, "\n"))
