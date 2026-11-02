# =====================================================================
SCRIPT_VERSION <- "ICI_FIGS_v12  v1.1  2026-10-16"
# Builds the two figures the v12 manuscript cites but that do not exist:
#   Figure 3  I-SPY2, n = 248 (from ICI_binary_ispy2_scores_248.csv)
#     A  AUC by arm per layer, DeLong CI, interaction p
#     B  OR per SD by arm per layer
#     C  pCR by Integration tertile, by arm
#     D  Integration OR by arm within TNBC and HR+/HER2-
#     E  pCR by lowest layer, by arm (shape null)
#     F  sequential gating in the pembrolizumab arm
#   Figure 4  NeoTRIP (GSE319641) baseline, same scoring as I-SPY2
#     A  AUC by arm   B  OR per SD by arm
#     C  pCR by Integration tertile   D  layer x arm interaction ORs
# Every panel's numbers are written to a source CSV next to the figure.
# Panel 'high' = score >= cohort median, the same rule used for binarizing.
# =====================================================================
suppressPackageStartupMessages({ library(data.table); library(pROC) })

BASE      <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"
ISPY      <- file.path(BASE, "MASTER_GRID/ICI_binary_ispy2_scores_248.csv")
NEO_DIR   <- path.expand("~/datasets/GSE319641")
NEO_EXPR  <- file.path(NEO_DIR, "GSE319641_NeoTRIP_baseline_D1C2_TPM_ComBat_all_samples.txt.gz")
NEO_PHENO <- file.path(NEO_DIR, "Early_GEP_NeoTRIP/data/NeoTRIP_baseline_D1C2_pheno_and_clinical.txt")
OUT       <- file.path(BASE, "Figures_v12")
dir.create(OUT, showWarnings = FALSE)
cat("\n", SCRIPT_VERSION, "\n", strrep("=", 62), "\n", sep = "")

COL <- c(control = "#8C8C8C", drug = "#1F6FB4")
fmtp <- function(p) ifelse(is.na(p), "NA", ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
                                                  formatC(p, format = "f", digits = 3)))
wilson <- function(x, n) {
  if (n == 0) return(c(NA, NA, NA)); z <- 1.96; p <- x / n; d <- 1 + z^2 / n
  cc <- (p + z^2 / (2 * n)) / d; h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / d
  c(p, cc - h, cc + h)
}
tag <- function(t) mtext(t, side = 3, adj = -0.12, font = 2, line = 0.6, cex = 1.1)

# ---------------- statistics: one row per layer -----------------------
stats_layer <- function(d, L) {
  d <- d[!is.na(d[[L]]), ]
  d$X <- as.numeric(scale(d[[L]])); d$A <- factor(d$arm, levels = c("control", "drug"))
  m  <- glm(y ~ X * A, binomial, d); m0 <- glm(y ~ X + A, binomial, d)
  ci <- confint.default(m)["X:Adrug", ]
  res <- data.frame(layer = L,
                    int_OR = unname(exp(coef(m)["X:Adrug"])), int_lo = unname(exp(ci[1])), int_hi = unname(exp(ci[2])),
                    p_wald = coef(summary(m))["X:Adrug", 4],
                    p_lrt  = anova(m0, m, test = "LRT")[2, "Pr(>Chi)"])
  for (a in c("control", "drug")) {
    s <- d[d$arm == a, ]; mm <- glm(y ~ X, binomial, s); cc <- confint.default(mm)["X", ]
    r <- roc(s$y, s[[L]], levels = c(0, 1), direction = "<", quiet = TRUE)
    ca <- as.numeric(ci.auc(r, method = "delong"))
    res[[paste0("OR_", a)]]   <- unname(exp(coef(mm)["X"]))
    res[[paste0("ORlo_", a)]] <- unname(exp(cc[1])); res[[paste0("ORhi_", a)]] <- unname(exp(cc[2]))
    res[[paste0("AUC_", a)]]  <- ca[2]; res[[paste0("AUClo_", a)]] <- ca[1]; res[[paste0("AUChi_", a)]] <- ca[3]
    res[[paste0("n_", a)]]    <- nrow(s); res[[paste0("ev_", a)]] <- sum(s$y)
  }
  res
}

# ---------------- panels ---------------------------------------------
panel_auc <- function(S, drug, t) {
  k <- nrow(S); x <- seq_len(k)
  plot(NA, xlim = c(0.5, k + 0.5), ylim = c(0.25, 1), xaxt = "n", xlab = "", ylab = "AUC for pCR (95% CI)", las = 1)
  abline(h = 0.5, lty = 2, col = "grey60")
  for (a in c("control", "drug")) {
    xx <- x + ifelse(a == "control", -0.15, 0.15)
    segments(xx, S[[paste0("AUClo_", a)]], xx, S[[paste0("AUChi_", a)]], col = COL[a], lwd = 2)
    points(xx, S[[paste0("AUC_", a)]], pch = 19, col = COL[a], cex = 1.3)
  }
  axis(1, at = x, labels = S$layer, cex.axis = 0.85)
  text(x, 0.28, paste0("interaction\np = ", fmtp(S$p_wald)), cex = 0.7)
  legend("topleft", c("Chemotherapy alone", drug), col = COL, pch = 19, bty = "n", cex = 0.8)
  tag(t)
}

panel_or <- function(S, drug, t, xlab = "OR for pCR per SD (95% CI)") {
  k <- nrow(S); y <- rev(seq_len(k)); XMAX <- 40
  plot(NA, log = "x", xlim = c(0.2, XMAX * 3), ylim = c(0.5, k + 0.5), yaxt = "n", ylab = "", xlab = xlab)
  abline(v = 1, lty = 2, col = "grey60")
  for (a in c("control", "drug")) {
    yy <- y + ifelse(a == "control", -0.15, 0.15)
    segments(pmax(S[[paste0("ORlo_", a)]], 0.2), yy, pmin(S[[paste0("ORhi_", a)]], XMAX), yy, col = COL[a], lwd = 2)
    points(S[[paste0("OR_", a)]], yy, pch = 19, col = COL[a])
  }
  axis(2, at = y, labels = S$layer, las = 1, cex.axis = 0.75)
  text(XMAX * 2.9, y, paste0("p = ", fmtp(S$p_wald)), cex = 0.72, adj = 1)
  legend("bottomleft", c("Chemotherapy alone", drug), col = COL, pch = 19, bty = "n", cex = 0.75)
  tag(t)
}

bars_by_arm <- function(tab, drug, t, xlab, ptxt) {   # tab: columns grp, arm, x, n
  grps <- unique(tab$grp)
  M  <- sapply(grps, function(g) sapply(c("control", "drug"), function(a) { r <- tab[tab$grp == g & tab$arm == a, ]; wilson(r$x, r$n) }), simplify = "array")
  P  <- 100 * M[1, , ]; LO <- 100 * M[2, , ]; HI <- 100 * M[3, , ]
  bp <- barplot(P, beside = TRUE, col = COL, border = NA, ylim = c(0, 110), las = 1, ylab = "pCR (%)", xlab = xlab, names.arg = grps)
  arrows(bp, LO, bp, HI, angle = 90, code = 3, length = 0.03, col = "grey30")
  nn <- sapply(grps, function(g) sapply(c("control", "drug"), function(a) tab$n[tab$grp == g & tab$arm == a]))
  text(bp, 3, paste0("n=", nn), cex = 0.6, col = "white")
  legend("topleft", c("Chemotherapy alone", drug), fill = COL, border = NA, bty = "n", cex = 0.75)
  mtext(ptxt, side = 3, line = -1.2, adj = 0.98, cex = 0.72)
  tag(t)
}

tertile_tab <- function(d, L = "Integration") {
  d <- d[!is.na(d[[L]]), ]
  q <- quantile(d[[L]], c(1/3, 2/3)); d$T <- cut(d[[L]], c(-Inf, q, Inf), labels = c("Low", "Mid", "High"))
  tab <- do.call(rbind, lapply(levels(d$T), function(g) do.call(rbind, lapply(c("control", "drug"), function(a) {
    s <- d[d$T == g & d$arm == a, ]; data.frame(grp = g, arm = a, x = sum(s$y), n = nrow(s)) }))))
  d$Tn <- as.integer(d$T); d$A <- factor(d$arm, levels = c("control", "drug"))
  p <- coef(summary(glm(y ~ Tn * A, binomial, d)))["Tn:Adrug", 4]
  list(tab = tab, p = p)
}

# =====================================================================
# FIGURE 3 — I-SPY2, n = 248
# =====================================================================
cat("\nFIGURE 3: I-SPY2\n")
d <- fread(ISPY, data.table = FALSE)
comp <- intersect(names(d), c("ICI_composite", "composite", "Composite"))[1]
d$Composite <- d[[comp]]
d$y <- as.integer(d$pcr_binary); d$arm <- ifelse(d$arm == "pembro", "drug", "control")
cat("patients:", nrow(d), "| arms:", paste(names(table(d$arm)), table(d$arm), collapse = " "), "\n")
LAY <- c("Sensory", "Integration", "Execution", "Composite")

S3 <- do.call(rbind, lapply(LAY, function(L) stats_layer(d, L)))
cat("\nlayer x arm (expect Integration p_wald ~0.025):\n")
print(S3[, c("layer", "OR_control", "OR_drug", "p_wald", "p_lrt", "AUC_control", "AUC_drug")], digits = 3)

T3 <- tertile_tab(d)

rec <- do.call(rbind, lapply(c("TNBC", "HRpos"), function(r) {
  s <- stats_layer(d[d$receptor == r, ], "Integration")
  s$layer <- sprintf("%s (n=%d)", ifelse(r == "TNBC", "TNBC", "HR+/HER2-"), s$n_control + s$n_drug); s }))

z <- scale(d[, c("Sensory", "Integration", "Execution")])
ok <- complete.cases(z)
d$lowest <- NA; d$lowest[ok] <- c("Sensory", "Integration", "Execution")[apply(z[ok, ], 1, which.min)]
dl <- d[ok, ]; dl$A <- factor(dl$arm, levels = c("control", "drug")); dl$lowest <- factor(dl$lowest)
p_low <- anova(glm(y ~ lowest + A, binomial, dl), glm(y ~ lowest * A, binomial, dl), test = "LRT")[2, "Pr(>Chi)"]
low_tab <- do.call(rbind, lapply(c("Sensory", "Integration", "Execution"), function(g) do.call(rbind, lapply(c("control", "drug"), function(a) {
  s <- dl[dl$lowest == g & dl$arm == a, ]; data.frame(grp = g, arm = a, x = sum(s$y), n = nrow(s)) }))))

hi <- function(L) d[[L]] >= median(d[[L]], na.rm = TRUE)
pb <- d$arm == "drug"
gates <- list("All"             = pb,
              "Sensory low"     = pb & !hi("Sensory"),
              "Sensory high"    = pb & hi("Sensory"),
              "+ Integration high" = pb & hi("Sensory") & hi("Integration"),
              "+ Execution high"   = pb & hi("Sensory") & hi("Integration") & hi("Execution"))
gate_tab <- do.call(rbind, lapply(names(gates), function(g) {
  k <- which(gates[[g]] %in% TRUE); w <- wilson(sum(d$y[k]), length(k))
  data.frame(gate = g, n = length(k), pcr = sum(d$y[k]), pct = 100 * w[1], lo = 100 * w[2], hi = 100 * w[3]) }))
cat("\ngating, pembrolizumab arm (defined post hoc):\n"); print(gate_tab, digits = 3)

draw3 <- function() {
  par(mfrow = c(2, 3), mar = c(4.2, 6, 2.6, 1), mgp = c(2.4, 0.7, 0))
  panel_auc(S3, "Paclitaxel + pembrolizumab", "A")
  panel_or(S3, "Paclitaxel + pembrolizumab", "B")
  bars_by_arm(T3$tab, "Paclitaxel + pembrolizumab", "C", "Integration tertile", paste0("tertile x arm p = ", fmtp(T3$p)))
  panel_or(rec, "Paclitaxel + pembrolizumab", "D", "Integration OR per SD (95% CI)")
  bars_by_arm(low_tab, "Paclitaxel + pembrolizumab", "E", "Lowest layer", paste0("lowest layer x arm p = ", fmtp(p_low)))
  par(mar = c(7, 4.5, 2.6, 1))
  bp <- barplot(gate_tab$pct, col = COL["drug"], border = NA, ylim = c(0, 110), ylab = "pCR (%)",
                names.arg = gate_tab$gate, las = 2, cex.names = 0.75)
  arrows(bp, gate_tab$lo, bp, gate_tab$hi, angle = 90, code = 3, length = 0.03, col = "grey30")
  text(bp, 4, paste0("n=", gate_tab$n), cex = 0.65, col = "white")
  mtext("Pembrolizumab arm; gates defined post hoc", side = 3, line = -1.2, adj = 0.98, cex = 0.65)
  tag("F")
}
pdf(file.path(OUT, "Figure3_ISPY2_n248.pdf"), width = 13, height = 8); draw3(); invisible(dev.off())
png(file.path(OUT, "Figure3_ISPY2_n248.png"), width = 13, height = 8, units = "in", res = 300); draw3(); invisible(dev.off())
fwrite(S3, file.path(OUT, "Figure3AB_source.csv")); fwrite(T3$tab, file.path(OUT, "Figure3C_source.csv"))
fwrite(rec, file.path(OUT, "Figure3D_source.csv")); fwrite(low_tab, file.path(OUT, "Figure3E_source.csv"))
fwrite(gate_tab, file.path(OUT, "Figure3F_source.csv"))
cat("Figure 3 written. tertile p =", fmtp(T3$p), "| lowest-layer p =", fmtp(p_low), "\n")

# =====================================================================
# FIGURE 4 — NeoTRIP baseline, frozen panel, same scoring
# =====================================================================
cat("\n", strrep("=", 62), "\nFIGURE 4: NeoTRIP\n", sep = "")
panel_src <- list.files(BASE, pattern = "^build_scRNA_celltype\\.R$", recursive = TRUE, full.names = TRUE)[1]
pe <- new.env(); lines <- readLines(panel_src, warn = FALSE)
keep <- grep("^\\s*(SENSORY|INTPOS|INTNEG|EXECPOS|EXECNEG)\\s*<-", lines); blk <- c()
for (k in keep) { j <- k
  repeat { blk <- c(blk, lines[j]); s <- paste(lines[k:j], collapse = "")
    if (lengths(regmatches(s, gregexpr("\\(", s))) == lengths(regmatches(s, gregexpr("\\)", s)))) break
    j <- j + 1; if (j > k + 25) stop("panel parse") } }
eval(parse(text = blk), envir = pe)
PANEL <- mget(c("SENSORY", "INTPOS", "INTNEG", "EXECPOS", "EXECNEG"), envir = pe)
stopifnot(length(unique(unlist(PANEL))) == 98L)

ex <- fread(NEO_EXPR, data.table = FALSE)
# v1.1: keep only sample columns; drop any annotation columns (e.g. a second ID column)
smp <- grepl("^TNBC", names(ex))
cat("expression columns:", ncol(ex), "| sample columns:", sum(smp),
    "| dropped:", paste(head(names(ex)[!smp], 5), collapse = ", "), "\n")
g  <- as.character(ex[[1]])
ex <- as.matrix(as.data.frame(lapply(ex[, smp, drop = FALSE], function(v) suppressWarnings(as.numeric(v)))))
keep_g <- !is.na(g) & g != ""; ex <- ex[keep_g, , drop = FALSE]; g <- g[keep_g]
cat("non-numeric cells set to NA:", sum(is.na(ex)), "\n")
ex <- log2(ex + 1)                                                    # file is LINEAR TPM
ex <- rowsum(ex, g, na.rm = TRUE) / as.vector(table(g)[sort(unique(g))])  # collapse duplicate symbols
ph <- fread(NEO_PHENO, data.table = FALSE)
ph <- ph[grepl("base", ph$Timepoint, ignore.case = TRUE) & ph$Sample_ID %in% colnames(ex), ]
cat("baseline samples matched:", nrow(ph), " (expect 241)\n")
ph$arm <- ifelse(grepl("A", sub("^CT", "", ph$Arm)), "drug", "control")
cat("arm mapping - check it:\n"); print(table(ph$Arm, ph$arm))
stopifnot(length(unique(ph$arm)) == 2)
ycol <- if ("pCR.num" %in% names(ph)) "pCR.num" else "pCR"
ph$y <- as.integer(as.character(ph[[ycol]]) %in% c("1", "pCR", "yes", "Yes"))

m <- ex[, ph$Sample_ID]
bin <- function(gs) { gg <- intersect(gs, rownames(m)); t(apply(m[gg, , drop = FALSE], 1, function(v) as.integer(v >= median(v, na.rm = TRUE)))) }
ph$Sensory     <- colSums(bin(PANEL$SENSORY))
ph$Integration <- colSums(bin(PANEL$INTPOS)) - colSums(bin(PANEL$INTNEG))
ph$Execution   <- colSums(bin(PANEL$EXECPOS)) - colSums(bin(PANEL$EXECNEG))
ph$Composite   <- rowMeans(scale(ph[, c("Sensory", "Integration", "Execution")]))
cat("panel genes off-platform:", paste(setdiff(unlist(PANEL), rownames(m)), collapse = ", "), "\n")

S4 <- do.call(rbind, lapply(LAY, function(L) stats_layer(ph, L)))
cat("\nlayer x arm (expect Integration interaction OR ~1.76, LRT p ~0.040):\n")
print(S4[, c("layer", "OR_control", "OR_drug", "int_OR", "p_wald", "p_lrt", "AUC_control", "AUC_drug")], digits = 3)
T4 <- tertile_tab(ph)

draw4 <- function() {
  par(mfrow = c(2, 2), mar = c(4.2, 6, 2.6, 1), mgp = c(2.4, 0.7, 0))
  panel_auc(S4, "Chemotherapy + atezolizumab", "A")
  panel_or(S4, "Chemotherapy + atezolizumab", "B")
  bars_by_arm(T4$tab, "Chemotherapy + atezolizumab", "C", "Integration tertile", paste0("tertile x arm p = ", fmtp(T4$p)))
  k <- nrow(S4); y <- rev(seq_len(k))
  plot(NA, log = "x", xlim = c(0.3, 8), ylim = c(0.5, k + 0.5), yaxt = "n", ylab = "", xlab = "Layer x arm interaction OR per SD (95% CI)")
  abline(v = 1, lty = 2, col = "grey60")
  segments(S4$int_lo, y, S4$int_hi, y, lwd = 2, col = COL["drug"]); points(S4$int_OR, y, pch = 19, col = COL["drug"])
  axis(2, at = y, labels = S4$layer, las = 1, cex.axis = 0.8)
  text(7.8, y, paste0("LRT p = ", fmtp(S4$p_lrt)), cex = 0.72, adj = 1)
  tag("D")
}
pdf(file.path(OUT, "Figure4_NeoTRIP.pdf"), width = 9.5, height = 8); draw4(); invisible(dev.off())
png(file.path(OUT, "Figure4_NeoTRIP.png"), width = 9.5, height = 8, units = "in", res = 300); draw4(); invisible(dev.off())
fwrite(S4, file.path(OUT, "Figure4ABD_source.csv")); fwrite(T4$tab, file.path(OUT, "Figure4C_source.csv"))
fwrite(ph[, c("Sample_ID", "Arm", "arm", "y", LAY)], file.path(OUT, "NeoTRIP_baseline_scores.csv"))
cat("Figure 4 written.\n\nAll outputs in", OUT, "\n")
