# =====================================================================
# FIGURE 2 — the layers are separable, and different patients are
# limited at different steps.
#
# A  per-cohort correlation between layer pairs (are the three layers
#    three names for the same thing?)
# B  which layer is lowest, per patient, per cohort
# C  example patient cascade profiles
#
# Reads every MASTER_GRID score file. Nothing is hardcoded per cohort.
# Output: ICI_Manuscript/Figures/Fig2_layer_separability.pdf (+ .png)
#         ICI_Manuscript/Figures/Fig2_source_numbers.csv
# =====================================================================
SCRIPT_VERSION <- "fig2_layer_separability v1.0  2026-08-21"
cat(SCRIPT_VERSION, "\n\n")

suppressPackageStartupMessages(library(ggplot2))

MAN    <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"
GRID   <- file.path(MAN, "MASTER_GRID")
OUTDIR <- file.path(MAN, "Figures")
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

MIN_N <- 20   # cohorts below this are computed but drawn faintly

pick <- function(df, opts) {
  hit <- opts[tolower(opts) %in% tolower(names(df))]
  if (!length(hit)) return(NA_character_)
  names(df)[match(tolower(hit[1]), tolower(names(df)))]
}
OPTS <- list(Sensory     = c("Sensory","Score_S","Layer1","S_z"),
             Integration = c("Integration","Score_I","Layer2","I_z"),
             Execution   = c("Execution","Score_E","Layer3","E_z"))

files <- list.files(GRID, pattern = "^ICI_binary_.*_scores\\.csv$", full.names = TRUE)
cat("score files:", length(files), "\n")

cors <- list(); lim <- list(); zall <- list()

for (f in files) {
  co <- sub("^ICI_binary_(.*)_scores\\.csv$", "\\1", basename(f))
  d  <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(d)) { cat(sprintf("  %-12s unreadable, skipped\n", co)); next }
  cn <- vapply(OPTS, function(o) pick(d, o), "")
  if (any(is.na(cn))) { cat(sprintf("  %-12s missing layer columns, skipped\n", co)); next }

  M <- as.matrix(d[, cn]); colnames(M) <- names(OPTS)
  M <- M[stats::complete.cases(M), , drop = FALSE]
  if (nrow(M) < 10) { cat(sprintf("  %-12s n=%d too small, skipped\n", co, nrow(M))); next }

  C <- stats::cor(M)
  cors[[co]] <- data.frame(
    cohort = co, n = nrow(M),
    pair = c("Sensory \u2013 Integration", "Integration \u2013 Execution",
             "Sensory \u2013 Execution"),
    r = c(C[1,2], C[2,3], C[1,3]), stringsAsFactors = FALSE)

  Z <- scale(M)
  who <- colnames(Z)[apply(Z, 1, which.min)]
  lim[[co]] <- data.frame(cohort = co, n = nrow(M),
                          layer = names(OPTS),
                          prop = as.numeric(table(factor(who, levels = names(OPTS)))) / nrow(M),
                          stringsAsFactors = FALSE)
  if (co == "ispy2") zall$ispy2 <- as.data.frame(Z)
  cat(sprintf("  %-12s n=%-4d S-I %.2f  I-E %.2f  S-E %.2f\n",
              co, nrow(M), C[1,2], C[2,3], C[1,3]))
}

cors <- do.call(rbind, cors); lim <- do.call(rbind, lim)
rownames(cors) <- NULL; rownames(lim) <- NULL
write.csv(merge(cors, aggregate(prop ~ cohort + layer, lim, sum),
                by = "cohort", all = TRUE),
          file.path(OUTDIR, "Fig2_source_numbers.csv"), row.names = FALSE)

ord <- cors$cohort[order(-cors$n)]
ord <- ord[!duplicated(ord)]
cors$cohort <- factor(cors$cohort, levels = rev(ord))
lim$cohort  <- factor(lim$cohort,  levels = rev(ord))
cors$small  <- cors$n < MIN_N
cors$pair   <- factor(cors$pair, levels = c("Sensory \u2013 Integration",
                                            "Integration \u2013 Execution",
                                            "Sensory \u2013 Execution"))
lim$layer   <- factor(lim$layer, levels = c("Sensory","Integration","Execution"))

PAL3 <- c("#1D9E75", "#378ADD", "#D4537E")

pA <- ggplot(cors, aes(r, cohort, colour = pair, alpha = !small, size = n)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
  geom_point() +
  scale_colour_manual(values = PAL3, name = NULL) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
  scale_size(range = c(1.1, 3.4), guide = "none") +
  scale_x_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Pearson r between layer scores", y = NULL,
       subtitle = "A   Layer pairs are correlated but not interchangeable") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top", panel.grid.minor = element_blank(),
        plot.subtitle = element_text(face = "bold", size = 9))

pB <- ggplot(lim, aes(prop, cohort, fill = layer)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.4) +
  scale_fill_manual(values = PAL3, name = NULL) +
  scale_x_continuous(labels = function(v) paste0(v * 100, "%"), expand = c(0, 0)) +
  labs(x = "Patients whose lowest z-scored layer is\u2026", y = NULL,
       subtitle = "B   Different patients are limited at different steps") +
  theme_minimal(base_size = 9) +
  theme(legend.position = "top", panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.subtitle = element_text(face = "bold", size = 9))

pC <- NULL
if (!is.null(zall$ispy2)) {
  Z <- zall$ispy2
  Z$lim <- colnames(Z)[apply(as.matrix(Z[, 1:3]), 1, which.min)]
  ex <- do.call(rbind, lapply(c("Sensory","Integration","Execution"), function(L) {
    s <- Z[Z$lim == L, ]
    s[order(s[[L]])[1], , drop = FALSE]
  }))
  ex$id <- paste(ex$lim, "limited")
  long <- data.frame(
    id    = rep(ex$id, each = 3),
    layer = factor(rep(c("Sensory","Integration","Execution"), nrow(ex)),
                   levels = c("Sensory","Integration","Execution")),
    z     = as.numeric(t(as.matrix(ex[, 1:3]))))
  long$id <- factor(long$id, levels = ex$id)
  pC <- ggplot(long, aes(layer, z, group = id)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey80") +
    geom_line(colour = "grey55", linewidth = 0.5) +
    geom_point(aes(colour = layer), size = 2.4) +
    scale_colour_manual(values = PAL3, guide = "none") +
    facet_wrap(~id, nrow = 1) +
    labs(x = NULL, y = "layer z-score",
         subtitle = "C   Example cascade profiles (I-SPY2)") +
    theme_minimal(base_size = 9) +
    theme(panel.grid.minor = element_blank(),
          strip.text = element_text(size = 8),
          plot.subtitle = element_text(face = "bold", size = 9))
}

fp <- file.path(OUTDIR, "Fig2_layer_separability")
combine <- function(a, b, c) {
  if (requireNamespace("patchwork", quietly = TRUE)) {
    library(patchwork)
    if (is.null(c)) (a | b) else (a | b) / c + plot_layout(heights = c(2.3, 1))
  } else if (requireNamespace("gridExtra", quietly = TRUE)) {
    if (is.null(c)) gridExtra::arrangeGrob(a, b, ncol = 2)
    else gridExtra::arrangeGrob(gridExtra::arrangeGrob(a, b, ncol = 2), c,
                                ncol = 1, heights = c(2.3, 1))
  } else NULL
}
g <- combine(pA, pB, pC)
if (is.null(g)) {
  cat("\npatchwork/gridExtra not installed - writing panels separately\n")
  ggsave(paste0(fp, "_A.pdf"), pA, width = 5.2, height = 4.4, device = cairo_pdf)
  ggsave(paste0(fp, "_B.pdf"), pB, width = 5.2, height = 4.4, device = cairo_pdf)
  if (!is.null(pC)) ggsave(paste0(fp, "_C.pdf"), pC, width = 7.5, height = 2.4,
                           device = cairo_pdf)
} else {
  ggsave(paste0(fp, ".pdf"), g, width = 10.4, height = 7.4, device = cairo_pdf)
  ggsave(paste0(fp, ".png"), g, width = 10.4, height = 7.4, dpi = 400)
  cat("\nWrote:\n ", paste0(fp, ".pdf"), "\n ", paste0(fp, ".png"), "\n")
}
cat(" ", file.path(OUTDIR, "Fig2_source_numbers.csv"), "\n")
