# Score an uploaded expression matrix (genes x samples) against the frozen ICI
# locked reference. Single-patient, cohort-independent (selfnorm). Returns per-
# patient layer z-scores + floor / archetype, plus coverage QC.
score_uploaded <- function(E, assets){
  ref <- assets$ref; PG <- assets$panel_genes
  E <- as.matrix(E); storage.mode(E) <- "numeric"
  rownames(E) <- toupper(rownames(E)); E <- E[!duplicated(rownames(E)), , drop=FALSE]
  logged <- FALSE
  if (suppressWarnings(quantile(E, 0.99, na.rm=TRUE)) > 50) { E <- log2(E + 1); logged <- TRUE }  # auto-log linear data
  cov_n <- length(intersect(PG, rownames(E)))
  sc <- score_locked(E, ref, mode = "selfnorm")
  zf <- function(x, L){ m <- ref$zstats[[L]]["mean"]; s <- ref$zstats[[L]]["sd"]
                        if (is.na(s) || s == 0) x*0 else (x - m)/s }
  Z <- data.frame(sample = as.character(sc$sample),
                  Sensory     = as.numeric(zf(sc$Sensory,     "Sensory")),
                  Integration = as.numeric(zf(sc$Integration, "Integration")),
                  Execution   = as.numeric(zf(sc$Execution,   "Execution")),
                  stringsAsFactors = FALSE)
  lay <- c("Sensory","Integration","Execution"); mat <- as.matrix(Z[, lay])
  wmin <- lay[max.col(-mat, ties.method="first")]; mz <- apply(mat, 1, min)
  Z$floor_layer <- ifelse(mz >= -0.5, "none", wmin)
  Z$archetype   <- ifelse(mz >= -0.5, "intact",
                     ifelse(wmin=="Sensory","cold", ifelse(wmin=="Integration","blocked","no-kill")))
  list(scores = Z, coverage = cov_n, n_panel = length(PG), n_samples = ncol(E), logged = logged)
}
