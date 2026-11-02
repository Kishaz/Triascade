# =============================================================================
# 05_tumor_purity_sensitivity.R   (v4 — handles integer sample-ID columns)
# Only change from v3: coerce scored$sample to character before join,
# so integer-typed pid columns (ispy2) don't blow up the inner_join with
# "Can't join <integer> with <character>".
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
})

DATA_DIR <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript/MASTER_GRID/"
PURITY_FILE <- file.path(DATA_DIR, "per_cohort_tumor_purity.csv")

if (!file.exists(PURITY_FILE)) {
  stop(paste0("\n\nMissing input: ", PURITY_FILE,
              "\nRun 08_extract_tumor_purity_ESTIMATE.R first.\n"))
}

purity <- read.csv(PURITY_FILE, stringsAsFactors = FALSE)

# Prefer string-form IDs first (pid_str before pid)
SAMPLE_COL_CANDIDATES <- c(
  "sample", "Sample", "sample_id", "SampleID", "Sample_ID", "sample_ID",
  "RNA_ID", "rna_id", "vst_id",
  "patient", "Patient", "patient_id", "PatientID", "Patient_ID", "patientID",
  "pid_str", "pid", "PID",
  "ID", "id", "subject_id", "Subject_ID"
)

find_sample_col <- function(df, candidates = SAMPLE_COL_CANDIDATES) {
  hit <- intersect(candidates, colnames(df))
  if (length(hit) > 0) return(hit[1])
  NA_character_
}

purity_sample_col <- find_sample_col(purity)
if (is.na(purity_sample_col)) {
  stop("\nCould not identify sample-ID column in purity file. ",
       "Columns: ", paste(colnames(purity), collapse=", "), "\n")
}
if (purity_sample_col != "sample") {
  names(purity)[names(purity) == purity_sample_col] <- "sample"
}
purity$sample <- as.character(purity$sample)
purity$sample_stripped <- sub("^X", "", purity$sample)

try_join <- function(scored, purity) {
  j1 <- inner_join(scored, purity[, c("sample", "purity", "ESTIMATE_score",
                                       "stromal_score", "immune_score")],
                   by = "sample")
  if (nrow(j1) >= nrow(scored) * 0.5) return(j1)
  purity2 <- purity
  purity2$sample <- purity2$sample_stripped
  j2 <- inner_join(scored, purity2[, c("sample", "purity", "ESTIMATE_score",
                                        "stromal_score", "immune_score")],
                   by = "sample")
  if (nrow(j2) > nrow(j1)) return(j2)
  j1
}

scored_files <- list.files(DATA_DIR,
                            pattern = "ICI_binary_.*_scores\\.csv$",
                            full.names = TRUE)
cat(sprintf("Found %d scored cohort files\n", length(scored_files)))

results <- list()
skipped <- character()

for (f in scored_files) {
  cohort_name <- gsub("ICI_binary_|_scores\\.csv", "", basename(f))
  scored <- tryCatch(read.csv(f, stringsAsFactors = FALSE),
                     error = function(e) NULL)
  if (is.null(scored)) { skipped <- c(skipped, cohort_name); next }

  sample_col <- find_sample_col(scored)
  if (is.na(sample_col)) {
    cat(sprintf("\n[%s] no sample-ID column among: %s -- skipping\n",
                cohort_name, paste(colnames(scored), collapse=", ")))
    skipped <- c(skipped, cohort_name); next
  }
  if (sample_col != "sample") {
    names(scored)[names(scored) == sample_col] <- "sample"
  }
  # === v4 fix: ensure character type before join ============================
  scored$sample <- as.character(scored$sample)

  resp_col_candidates <- c("response_bin", "response", "responder",
                            "pcr_binary", "pCR", "pcr")
  resp_col <- intersect(resp_col_candidates, colnames(scored))
  if (length(resp_col) == 0) {
    cat(sprintf("\n[%s] no response column -- skipping\n", cohort_name))
    skipped <- c(skipped, cohort_name); next
  }
  if (resp_col[1] != "response_bin") {
    names(scored)[names(scored) == resp_col[1]] <- "response_bin"
  }

  joined <- try_join(scored, purity)

  if (nrow(joined) < 8) {
    cat(sprintf("\n[%s] sample_col='%s'; %d matched -- skipping\n",
                cohort_name, sample_col, nrow(joined)))
    skipped <- c(skipped, cohort_name); next
  }
  cat(sprintf("\n=== %s (sample_col='%s'; n=%d matched, R=%d) ===\n",
              cohort_name, sample_col, nrow(joined), sum(joined$response_bin)))

  qs <- unique(quantile(joined$purity, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE))
  if (length(qs) >= 3) {
    joined$purity_tertile <- cut(joined$purity, breaks = qs,
                                  labels = c("low", "mid", "high")[seq_len(length(qs) - 1)],
                                  include.lowest = TRUE)
  } else {
    joined$purity_tertile <- factor("all")
  }
  joined$purity_z <- as.numeric(scale(joined$purity))

  layer_candidates <- intersect(c("Sensory", "Integration", "Execution",
                                   "composite", "Composite", "ICI_composite"),
                                 colnames(joined))

  for (lyr in layer_candidates) {
    auc_all <- tryCatch(
      as.numeric(suppressMessages(pROC::roc(joined$response_bin, joined[[lyr]])$auc)),
      error = function(e) NA_real_)
    if (!is.na(auc_all)) {
      cat(sprintf("  %-14s ALL  (n=%d): AUC = %.3f\n", lyr, nrow(joined), auc_all))
      results[[length(results) + 1]] <- data.frame(
        cohort = cohort_name, layer = lyr, tertile = "ALL",
        n = nrow(joined), n_R = sum(joined$response_bin), AUC = auc_all
      )
    }
    for (tt in levels(joined$purity_tertile)) {
      sub <- joined[joined$purity_tertile == tt, ]
      if (length(unique(sub$response_bin)) < 2 || nrow(sub) < 5) next
      auc_val <- tryCatch(
        as.numeric(suppressMessages(pROC::roc(sub$response_bin, sub[[lyr]])$auc)),
        error = function(e) NA_real_)
      if (!is.na(auc_val)) {
        cat(sprintf("  %-14s %s tertile (n=%d, R=%d): AUC = %.3f\n",
                    lyr, tt, nrow(sub), sum(sub$response_bin), auc_val))
        results[[length(results) + 1]] <- data.frame(
          cohort = cohort_name, layer = lyr, tertile = tt,
          n = nrow(sub), n_R = sum(sub$response_bin), AUC = auc_val
        )
      }
    }
    fmla <- as.formula(sprintf("response_bin ~ scale(%s) + purity_z", lyr))
    fit <- tryCatch(glm(fmla, family = binomial(), data = joined),
                    error = function(e) NULL)
    if (!is.null(fit)) {
      or_layer <- exp(coef(fit)[2])
      p_layer  <- summary(fit)$coefficients[2, "Pr(>|z|)"]
      or_pur   <- exp(coef(fit)[3])
      p_pur    <- summary(fit)$coefficients[3, "Pr(>|z|)"]
      cat(sprintf("  %-14s purity-adjusted: layer OR/SD = %.2f (p = %.3f), purity OR/SD = %.2f (p = %.3f)\n",
                  lyr, or_layer, p_layer, or_pur, p_pur))
    }
  }
}

if (length(results) > 0) {
  out_table <- do.call(rbind, results)
  out_file <- file.path(DATA_DIR, "tumor_purity_sensitivity_results.csv")
  write.csv(out_table, out_file, row.names = FALSE)
  cat(sprintf("\nSaved %d rows to:\n  %s\n", nrow(out_table), out_file))
}

if (length(skipped) > 0) {
  cat(sprintf("\nSkipped: %s\n", paste(skipped, collapse=", ")))
}

cat("\n---- sessionInfo ----\n")
print(sessionInfo())
