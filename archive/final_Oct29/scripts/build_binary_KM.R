# =====================================================================
# ICI 2.0 — Multi-cohort BINARY KM survival curves
# =====================================================================
# Cohorts with survival data accessible:
#   IMvigor210 (OS, n=298)
#   Braun ccRCC (OS, n=181)
#   NSCLC (PFS, n=43)
#   EAC (OS, n=32)
#
# For each: median-split KM by Sensory, Integration, Composite (binary).
# Plus: pooled cross-cohort forest of HRs.
# Plus: IMvigor inflamed-only panel (Integration vs CD8A vs TIS vs Teff).
#
# Run:
#   source("/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript/build_binary_KM.R")
# =====================================================================

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(grid)
  library(gridExtra)
  library(dplyr)
})

BASE     <- "/Users/dbedi/Desktop/R data FINAL"
MAN_DIR  <- file.path(BASE, "ICI_Manuscript")
GRID_DIR <- file.path(MAN_DIR, "MASTER_GRID")
OUT_DIR  <- file.path(GRID_DIR, "Figures_v2", "KM_panels")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cat("==========================================\n")
cat("Multi-cohort BINARY KM survival\n")
cat("Output: ", OUT_DIR, "\n")
cat("==========================================\n\n")

LAYER_COLORS <- c(
  "Sensory"     = "#2E86AB",
  "Integration" = "#A23B72",
  "Execution"   = "#F18F01",
  "Composite"   = "#3B7AB8"
)

# ---- KM helper ----------------------------------------------------
make_km <- function(d, layer_col, time_col, event_col, layer_name,
                    cohort_name, time_label = "Time (months)") {
  d <- d[!is.na(d[[layer_col]]) & !is.na(d[[time_col]]) & !is.na(d[[event_col]]), ]
  if (nrow(d) < 10) return(NULL)
  med <- median(d[[layer_col]], na.rm = TRUE)
  d$grp <- factor(ifelse(d[[layer_col]] >= med, "High", "Low"), levels = c("Low", "High"))
  if (length(unique(d$grp)) < 2) return(NULL)

  # Standardize time/event into well-known column names so ggsurvplot doesn't
  # choke on dynamic formulas
  d$.time  <- d[[time_col]]
  d$.event <- d[[event_col]]

  fit <- survfit(Surv(.time, .event) ~ grp, data = d)
  cox <- coxph(Surv(.time, .event) ~ grp, data = d)

  s <- summary(cox)
  hr <- s$coefficients[1, "exp(coef)"]
  ci_lo <- s$conf.int[1, "lower .95"]
  ci_hi <- s$conf.int[1, "upper .95"]
  p <- s$coefficients[1, "Pr(>|z|)"]

  n_high <- sum(d$grp == "High")
  n_low  <- sum(d$grp == "Low")

  hr_label <- sprintf("HR = %.2f (%.2f–%.2f)\np = %s",
                      hr, ci_lo, ci_hi,
                      ifelse(p < 0.001, sprintf("%.0e", p), sprintf("%.3f", p)))

  list(
    fit = fit, d = d, cox = cox,
    hr = hr, ci_lo = ci_lo, ci_hi = ci_hi, p = p,
    n_high = n_high, n_low = n_low,
    label = hr_label,
    layer = layer_name,
    cohort = cohort_name,
    time_label = time_label,
    time_col = ".time"
  )
}

plot_km <- function(km_result) {
  if (is.null(km_result)) return(NULL)
  p <- ggsurvplot(km_result$fit,
                  data = km_result$d,
                  pval = FALSE,
                  conf.int = FALSE,
                  legend.labs = c(sprintf("Low (n=%d)", km_result$n_low),
                                  sprintf("High (n=%d)", km_result$n_high)),
                  legend.title = "",
                  palette = c("steelblue", "darkred"),
                  risk.table = FALSE,
                  ggtheme = theme_minimal(base_size = 10) +
                    theme(panel.grid.minor = element_blank(),
                          legend.background = element_rect(fill = "white", color = "grey80"),
                          plot.title = element_text(face = "bold", size = 11)),
                  title = sprintf("%s (binary, %s)", km_result$layer, km_result$cohort),
                  xlab = km_result$time_label,
                  ylab = "Overall Survival")
  tmax <- max(km_result$d$.time, na.rm = TRUE)
  p_plot <- p$plot +
    annotate("text", x = tmax * 0.55, y = 0.15,
             label = km_result$label, size = 3, hjust = 0)
  p_plot
}

# Storage for HR forest data
hr_results <- list()

# =====================================================
# IMvigor210
# =====================================================
cat("[IMvigor210]\n")
imv_clin_path <- file.path(MAN_DIR, "SIGNIOS_TMB_IMvigor210.csv")
imv_bin_path  <- file.path(GRID_DIR, "ICI_binary_imvigor_scores.csv")

imv_clin <- read.csv(imv_clin_path, stringsAsFactors = FALSE, check.names = FALSE)
imv_bin  <- read.csv(imv_bin_path, stringsAsFactors = FALSE)
cat("  Clinical cols:", paste(colnames(imv_clin), collapse = ", "), "\n")
cat("  Binary cols:", paste(colnames(imv_bin), collapse = ", "), "\n")

# Determine ID column in clinical
imv_id_col <- intersect(c("sample_id", "Sample_ID", "sample", "ID"),
                        colnames(imv_clin))[1]
if (is.na(imv_id_col)) {
  # First col is usually IDs
  imv_id_col <- colnames(imv_clin)[1]
}
imv_clin$sample <- imv_clin[[imv_id_col]]

imv_merged <- merge(imv_bin, imv_clin[, c("sample", "os_time", "os_event")],
                    by = "sample", all.x = FALSE)
cat(sprintf("  Merged: %d patients with OS\n", nrow(imv_merged)))

# Convert time to months if needed (some IMvigor data is in months already, some days)
# Check range
if (max(imv_merged$os_time, na.rm = TRUE) > 100) {
  cat("  os_time appears to be in days; converting to months\n")
  imv_merged$os_time <- imv_merged$os_time / 30.44
}

km_imv_S <- make_km(imv_merged, "Sensory",     "os_time", "os_event",
                    "Sensory (SEE)", "IMvigor210")
km_imv_I <- make_km(imv_merged, "Integration", "os_time", "os_event",
                    "Integration (REACH)", "IMvigor210")
km_imv_E <- make_km(imv_merged, "Execution",   "os_time", "os_event",
                    "Execution (KILL)", "IMvigor210")
km_imv_C <- make_km(imv_merged, "composite",   "os_time", "os_event",
                    "Composite (ICI)", "IMvigor210")

for (k in list(km_imv_S, km_imv_I, km_imv_E, km_imv_C)) {
  if (!is.null(k)) {
    hr_results[[length(hr_results) + 1]] <- data.frame(
      cohort = k$cohort, layer = k$layer,
      hr = k$hr, ci_lo = k$ci_lo, ci_hi = k$ci_hi, p = k$p,
      n_high = k$n_high, n_low = k$n_low,
      stringsAsFactors = FALSE
    )
  }
}

# Save IMvigor 4-panel
imv_panels <- list()
for (k in list(km_imv_S, km_imv_I, km_imv_E, km_imv_C)) {
  if (!is.null(k)) imv_panels[[length(imv_panels) + 1]] <- plot_km(k)
}
if (length(imv_panels) > 0) {
  g <- arrangeGrob(grobs = imv_panels, ncol = 2,
                   top = textGrob("IMvigor210 — Binary KM by ICI layers",
                                  gp = gpar(fontsize = 12, fontface = "bold")))
  ggsave(file.path(OUT_DIR, "KM_IMvigor210_binary.pdf"), g, width = 11, height = 9)
  ggsave(file.path(OUT_DIR, "KM_IMvigor210_binary.png"), g, width = 11, height = 9, dpi = 200)
  cat(sprintf("  Wrote: KM_IMvigor210_binary.pdf (%d panels)\n", length(imv_panels)))
}

# =====================================================
# Braun ccRCC
# =====================================================
cat("\n[Braun ccRCC]\n")
braun_full_path <- file.path(MAN_DIR, "Braun_ICI_FULL_results", "Braun_FULL_ICI_scores_nivolumab.csv")
if (!file.exists(braun_full_path)) {
  braun_full_path <- file.path(MAN_DIR, "Braun_ICI_FULL_results", "Braun_FULL_ICI_scores_all.csv")
}

braun_clin <- read.csv(braun_full_path, stringsAsFactors = FALSE, check.names = FALSE)
braun_bin  <- read.csv(file.path(GRID_DIR, "ICI_binary_braun_scores.csv"),
                       stringsAsFactors = FALSE)
cat("  Clinical cols:", paste(head(colnames(braun_clin), 30), collapse = ", "), "...\n")
cat("  Binary cols:", paste(colnames(braun_bin), collapse = ", "), "\n")

# Find an ID column that matches between binary scores and clinical
braun_bin_id <- "RNA_ID"  # confirmed earlier
candidates <- intersect(c("RNA_ID", "ID", "Sample_ID", "patient_id", "sample"),
                        colnames(braun_clin))
if (length(candidates) == 0) {
  cat("  No matching ID column. First clinical col:", colnames(braun_clin)[1], "\n")
} else {
  braun_clin_id <- candidates[1]
  braun_clin$RNA_ID <- braun_clin[[braun_clin_id]]

  # Take only relevant cols
  if ("os_time" %in% colnames(braun_clin) && "os_event" %in% colnames(braun_clin)) {
    keep_cols <- c("RNA_ID", "os_time", "os_event")
    if ("pfs_time" %in% colnames(braun_clin)) keep_cols <- c(keep_cols, "pfs_time", "pfs_event")
    braun_merged <- merge(braun_bin, braun_clin[, keep_cols],
                          by = "RNA_ID", all.x = FALSE)
    cat(sprintf("  Merged: %d patients with OS\n", nrow(braun_merged)))

    # Convert time if needed
    if (max(braun_merged$os_time, na.rm = TRUE) > 200) {
      cat("  os_time appears to be in days; converting to months\n")
      braun_merged$os_time <- braun_merged$os_time / 30.44
    }

    km_braun_S <- make_km(braun_merged, "Sensory",     "os_time", "os_event",
                          "Sensory (SEE)", "Braun ccRCC")
    km_braun_I <- make_km(braun_merged, "Integration", "os_time", "os_event",
                          "Integration (REACH)", "Braun ccRCC")
    km_braun_E <- make_km(braun_merged, "Execution",   "os_time", "os_event",
                          "Execution (KILL)", "Braun ccRCC")
    km_braun_C <- make_km(braun_merged, "composite",   "os_time", "os_event",
                          "Composite (ICI)", "Braun ccRCC")

    for (k in list(km_braun_S, km_braun_I, km_braun_E, km_braun_C)) {
      if (!is.null(k)) {
        hr_results[[length(hr_results) + 1]] <- data.frame(
          cohort = k$cohort, layer = k$layer,
          hr = k$hr, ci_lo = k$ci_lo, ci_hi = k$ci_hi, p = k$p,
          n_high = k$n_high, n_low = k$n_low,
          stringsAsFactors = FALSE
        )
      }
    }

    braun_panels <- list()
    for (k in list(km_braun_S, km_braun_I, km_braun_E, km_braun_C)) {
      if (!is.null(k)) braun_panels[[length(braun_panels) + 1]] <- plot_km(k)
    }
    if (length(braun_panels) > 0) {
      g <- arrangeGrob(grobs = braun_panels, ncol = 2,
                       top = textGrob("Braun ccRCC — Binary KM by ICI layers (negative control)",
                                      gp = gpar(fontsize = 12, fontface = "bold")))
      ggsave(file.path(OUT_DIR, "KM_Braun_binary.pdf"), g, width = 11, height = 9)
      ggsave(file.path(OUT_DIR, "KM_Braun_binary.png"), g, width = 11, height = 9, dpi = 200)
      cat(sprintf("  Wrote: KM_Braun_binary.pdf (%d panels)\n", length(braun_panels)))
    }
  } else {
    cat("  os_time/os_event columns not found in Braun clinical\n")
  }
}

# =====================================================
# NSCLC (PFS)
# =====================================================
cat("\n[NSCLC]\n")
nsclc_clin <- read.csv(file.path(BASE, "datasets", "NSCLC_ICI_results",
                                  "NSCLC_combined_ICI_scores.csv"),
                       stringsAsFactors = FALSE)
nsclc_bin <- read.csv(file.path(GRID_DIR, "ICI_binary_nsclc_scores.csv"),
                      stringsAsFactors = FALSE)
cat("  Clinical cols:", paste(colnames(nsclc_clin), collapse = ", "), "\n")
cat("  Binary cols:", paste(colnames(nsclc_bin), collapse = ", "), "\n")

# Try to match
nsclc_id_clin <- intersect(c("sample_id", "expr_id", "patient_id", "ID"),
                           colnames(nsclc_clin))[1]
nsclc_id_bin  <- intersect(c("sample_id", "expr_id", "patient_id", "ID"),
                           colnames(nsclc_bin))[1]

if (!is.na(nsclc_id_clin) && !is.na(nsclc_id_bin)) {
  if ("pfs_time" %in% colnames(nsclc_clin) && "pfs_event" %in% colnames(nsclc_clin)) {
    nsclc_bin$.id <- nsclc_bin[[nsclc_id_bin]]
    nsclc_clin$.id <- nsclc_clin[[nsclc_id_clin]]
    nsclc_merged <- merge(nsclc_bin, nsclc_clin[, c(".id", "pfs_time", "pfs_event")],
                          by = ".id", all.x = FALSE)
    cat(sprintf("  Merged: %d patients with PFS\n", nrow(nsclc_merged)))

    if (max(nsclc_merged$pfs_time, na.rm = TRUE) > 100) {
      nsclc_merged$pfs_time <- nsclc_merged$pfs_time / 30.44
    }

    km_nsclc_S <- make_km(nsclc_merged, "Sensory",     "pfs_time", "pfs_event",
                          "Sensory (SEE)", "NSCLC", "PFS Time (months)")
    km_nsclc_I <- make_km(nsclc_merged, "Integration", "pfs_time", "pfs_event",
                          "Integration (REACH)", "NSCLC", "PFS Time (months)")
    km_nsclc_E <- make_km(nsclc_merged, "Execution",   "pfs_time", "pfs_event",
                          "Execution (KILL)", "NSCLC", "PFS Time (months)")
    km_nsclc_C <- make_km(nsclc_merged, "composite",   "pfs_time", "pfs_event",
                          "Composite (ICI)", "NSCLC", "PFS Time (months)")

    for (k in list(km_nsclc_S, km_nsclc_I, km_nsclc_E, km_nsclc_C)) {
      if (!is.null(k)) {
        hr_results[[length(hr_results) + 1]] <- data.frame(
          cohort = k$cohort, layer = k$layer,
          hr = k$hr, ci_lo = k$ci_lo, ci_hi = k$ci_hi, p = k$p,
          n_high = k$n_high, n_low = k$n_low,
          stringsAsFactors = FALSE
        )
      }
    }

    nsclc_panels <- list()
    for (k in list(km_nsclc_S, km_nsclc_I, km_nsclc_E, km_nsclc_C)) {
      if (!is.null(k)) nsclc_panels[[length(nsclc_panels) + 1]] <- plot_km(k)
    }
    if (length(nsclc_panels) > 0) {
      g <- arrangeGrob(grobs = nsclc_panels, ncol = 2,
                       top = textGrob("NSCLC — Binary KM by ICI layers (PFS)",
                                      gp = gpar(fontsize = 12, fontface = "bold")))
      ggsave(file.path(OUT_DIR, "KM_NSCLC_binary.pdf"), g, width = 11, height = 9)
      ggsave(file.path(OUT_DIR, "KM_NSCLC_binary.png"), g, width = 11, height = 9, dpi = 200)
      cat(sprintf("  Wrote: KM_NSCLC_binary.pdf (%d panels)\n", length(nsclc_panels)))
    }
  } else {
    cat("  pfs_time/pfs_event not found\n")
  }
} else {
  cat("  ID column matching failed\n")
}

# =====================================================
# EAC
# =====================================================
cat("\n[EAC]\n")
eac_path_candidates <- c(
  file.path(BASE, "Esophageal_cere_predictions.csv")
)
eac_clin <- NULL
for (p in eac_path_candidates) {
  if (file.exists(p)) {
    eac_clin <- read.csv(p, stringsAsFactors = FALSE, check.names = FALSE)
    cat(sprintf("  Loaded: %s (%d cols)\n", basename(p), ncol(eac_clin)))
    break
  }
}

if (!is.null(eac_clin)) {
  eac_bin <- read.csv(file.path(GRID_DIR, "ICI_binary_eac_scores.csv"),
                      stringsAsFactors = FALSE)
  cat("  Binary cols:", paste(colnames(eac_bin), collapse = ", "), "\n")

  # Find ID column
  eac_id_cands <- intersect(c("vst_id", "sample_id", "patient_id", "patient", "ID", "Sample_ID"),
                            colnames(eac_clin))
  if (length(eac_id_cands) > 0 && "os_time" %in% colnames(eac_clin)
      && "os_event" %in% colnames(eac_clin)) {
    eac_id <- eac_id_cands[1]
    eac_clin$vst_id <- eac_clin[[eac_id]]
    eac_merged <- merge(eac_bin, eac_clin[, c("vst_id", "os_time", "os_event")],
                        by = "vst_id", all.x = FALSE)
    cat(sprintf("  Merged: %d patients with OS\n", nrow(eac_merged)))

    if (nrow(eac_merged) >= 10) {
      if (max(eac_merged$os_time, na.rm = TRUE) > 200) {
        eac_merged$os_time <- eac_merged$os_time / 30.44
      }
      km_eac_S <- make_km(eac_merged, "Sensory",     "os_time", "os_event",
                          "Sensory (SEE)", "EAC")
      km_eac_I <- make_km(eac_merged, "Integration", "os_time", "os_event",
                          "Integration (REACH)", "EAC")
      km_eac_E <- make_km(eac_merged, "Execution",   "os_time", "os_event",
                          "Execution (KILL)", "EAC")
      km_eac_C <- make_km(eac_merged, "composite",   "os_time", "os_event",
                          "Composite (ICI)", "EAC")

      for (k in list(km_eac_S, km_eac_I, km_eac_E, km_eac_C)) {
        if (!is.null(k)) {
          hr_results[[length(hr_results) + 1]] <- data.frame(
            cohort = k$cohort, layer = k$layer,
            hr = k$hr, ci_lo = k$ci_lo, ci_hi = k$ci_hi, p = k$p,
            n_high = k$n_high, n_low = k$n_low,
            stringsAsFactors = FALSE
          )
        }
      }

      eac_panels <- list()
      for (k in list(km_eac_S, km_eac_I, km_eac_E, km_eac_C)) {
        if (!is.null(k)) eac_panels[[length(eac_panels) + 1]] <- plot_km(k)
      }
      if (length(eac_panels) > 0) {
        g <- arrangeGrob(grobs = eac_panels, ncol = 2,
                         top = textGrob("EAC — Binary KM by ICI layers",
                                        gp = gpar(fontsize = 12, fontface = "bold")))
        ggsave(file.path(OUT_DIR, "KM_EAC_binary.pdf"), g, width = 11, height = 9)
        ggsave(file.path(OUT_DIR, "KM_EAC_binary.png"), g, width = 11, height = 9, dpi = 200)
        cat(sprintf("  Wrote: KM_EAC_binary.pdf (%d panels)\n", length(eac_panels)))
      }
    } else {
      cat("  Insufficient merged rows for KM\n")
    }
  } else {
    cat("  ID col or OS cols missing in EAC clinical\n")
  }
} else {
  cat("  EAC clinical file not found\n")
}

# =====================================================
# Cross-cohort HR forest
# =====================================================
cat("\n[Cross-cohort HR forest]\n")
if (length(hr_results) > 0) {
  hr_df <- do.call(rbind, hr_results)
  hr_df$layer <- factor(hr_df$layer,
                         levels = c("Sensory (SEE)", "Integration (REACH)",
                                    "Execution (KILL)", "Composite (ICI)"))
  hr_df$cohort <- factor(hr_df$cohort,
                          levels = c("IMvigor210", "Braun ccRCC", "NSCLC", "EAC"))
  hr_df$lbl <- sprintf("%s (n=%d)", hr_df$cohort, hr_df$n_high + hr_df$n_low)
  hr_df$sig <- hr_df$p < 0.05

  write.csv(hr_df, file.path(GRID_DIR, "HR_summary_binary.csv"), row.names = FALSE)

  p_forest <- ggplot(hr_df, aes(x = hr, y = cohort, color = sig)) +
    geom_vline(xintercept = 1, color = "grey60", linetype = "dashed") +
    geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi), height = 0.2, linewidth = 0.5) +
    geom_point(aes(size = n_high + n_low), shape = 22, fill = NA, stroke = 1.5) +
    facet_wrap(~ layer, ncol = 2) +
    scale_x_log10(limits = c(0.3, 5),
                  breaks = c(0.5, 1, 2, 5)) +
    scale_color_manual(values = c(`TRUE` = "darkred", `FALSE` = "grey50"),
                       labels = c("p ≥ 0.05", "p < 0.05"),
                       name = NULL) +
    scale_size_continuous(range = c(2, 5), guide = "none") +
    geom_text(aes(label = sprintf("HR=%.2f, p=%s",
                                  hr,
                                  ifelse(p < 0.001, sprintf("%.0e", p),
                                         sprintf("%.3f", p)))),
              x = log10(3), hjust = 0, size = 2.6, color = "black") +
    labs(title = "Cross-cohort HR forest (binary, High vs Low layer score)",
         subtitle = "High score = above cohort median; OS or PFS Cox regression",
         x = "Hazard Ratio (log scale)",
         y = NULL) +
    theme_minimal(base_size = 11) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(face = "bold"),
          strip.background = element_rect(fill = "grey92", color = NA),
          legend.position = "bottom")

  ggsave(file.path(OUT_DIR, "KM_HR_forest_crosscohort.pdf"), p_forest,
         width = 11, height = 8)
  ggsave(file.path(OUT_DIR, "KM_HR_forest_crosscohort.png"), p_forest,
         width = 11, height = 8, dpi = 200)
  cat("  Wrote: KM_HR_forest_crosscohort.pdf\n")

  cat("\nHR summary:\n")
  print(hr_df[, c("cohort", "layer", "hr", "ci_lo", "ci_hi", "p", "n_high", "n_low")],
        row.names = FALSE)
}

cat("\n==========================================\n")
cat("DONE. Files in: ", OUT_DIR, "\n")
files_made <- list.files(OUT_DIR, pattern = "\\.pdf$")
for (f in files_made) cat("  -", f, "\n")
