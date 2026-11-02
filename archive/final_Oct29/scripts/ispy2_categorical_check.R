# =====================================================================
# Does ANY categorical encoding of the layer profile stratify
# pembrolizumab benefit in I-SPY2?
#
# The continuous Integration x arm interaction (p=0.008) is NOT in
# question here. This asks a narrower thing: does knowing WHICH layer
# is limiting, or the SHAPE of the profile, add anything?
#
# Every test is a layer-profile variable x arm interaction, by LRT.
# Nine tests -> Bonferroni 0.0056. Read the pCR cell rates first.
#
# Output: ICI_Manuscript/ISPY2/CATEGORICAL/ispy2_categorical_interactions.csv
# =====================================================================
SCRIPT_VERSION <- "ispy2_categorical_check v1.0  2026-08-21"
cat(SCRIPT_VERSION, "\n\n")

MAN <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"
F   <- file.path(MAN, "MASTER_GRID", "ICI_binary_ispy2_scores.csv")
OUT <- file.path(MAN, "ISPY2", "CATEGORICAL")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

d <- read.csv(F, stringsAsFactors = FALSE)
d$y   <- as.integer(d[["pcr_binary"]])
d$arm <- factor(d$arm, levels = c("control", "pembro"))
d <- d[!is.na(d$y) & !is.na(d$arm), ]

L <- c("Sensory", "Integration", "Execution")
Z <- scale(as.matrix(d[, L]));  colnames(Z) <- L
cat(sprintf("n=%d  (control %d, pembro %d)\n\n",
            nrow(d), sum(d$arm == "control"), sum(d$arm == "pembro")))

# ---- candidate encodings --------------------------------------------
enc <- list()

enc$min_layer <- factor(L[apply(Z, 1, which.min)], levels = c("Execution","Sensory","Integration"))

enc$int_tertile <- cut(Z[, "Integration"], breaks = quantile(Z[, "Integration"], c(0, 1/3, 2/3, 1)),
                       include.lowest = TRUE, labels = c("low","mid","high"))
enc$int_median  <- factor(ifelse(Z[, "Integration"] >= median(Z[, "Integration"]), "high", "low"),
                          levels = c("low","high"))

# profile SHAPE, independent of overall level
enc$S_minus_E   <- Z[, "Sensory"] - Z[, "Execution"]
enc$I_minus_mean<- Z[, "Integration"] - rowMeans(Z)
enc$spread      <- apply(Z, 1, max) - apply(Z, 1, min)

# is the minimum layer DEEP, i.e. a real bottleneck rather than a tie?
mn <- apply(Z, 1, min); md <- apply(Z, 1, function(v) sort(v)[2])
enc$gap         <- md - mn
enc$deep_bottle <- factor(ifelse(enc$gap >= median(enc$gap) & mn < 0, "deep", "shallow"),
                          levels = c("shallow","deep"))
# min layer, but only among patients with a genuinely deep minimum
dl <- as.character(enc$min_layer); dl[enc$deep_bottle == "shallow"] <- NA
enc$min_layer_deep <- factor(dl, levels = c("Execution","Sensory","Integration"))

# ---- test each --------------------------------------------------------
res <- list()
for (nm in names(enc)) {
  v  <- enc[[nm]]
  ok <- !is.na(v)
  dd <- d[ok, ]; dd$v <- v[ok]
  if (is.factor(dd$v)) dd$v <- droplevels(dd$v)
  if (is.factor(dd$v) && nlevels(dd$v) < 2) next
  m0 <- glm(y ~ v + arm, data = dd, family = binomial)
  m1 <- glm(y ~ v * arm, data = dd, family = binomial)
  a  <- anova(m0, m1, test = "Chisq")
  res[[nm]] <- data.frame(
    encoding = nm, type = if (is.factor(dd$v)) "categorical" else "continuous",
    n = nrow(dd), df = a$Df[2], deviance = round(a$Deviance[2], 3),
    interaction_p = signif(a$`Pr(>Chi)`[2], 3), stringsAsFactors = FALSE)
  cat(sprintf("%-16s %-11s n=%-4d df=%d  dev=%6.3f  p=%.4f\n",
              nm, res[[nm]]$type, nrow(dd), a$Df[2], a$Deviance[2], a$`Pr(>Chi)`[2]))
}
res <- do.call(rbind, res); rownames(res) <- NULL

# ---- reference: the continuous Integration result ---------------------
d$Iz <- Z[, "Integration"]
mc0 <- glm(y ~ Iz + arm, data = d, family = binomial)
mc1 <- glm(y ~ Iz * arm, data = d, family = binomial)
ac  <- anova(mc0, mc1, test = "Chisq")
cat(sprintf("\nREFERENCE  Integration continuous x arm:  p=%.4f\n", ac$`Pr(>Chi)`[2]))

cat("\n---- pCR rates, deep-bottleneck patients only ----\n")
sub <- d[!is.na(enc$min_layer_deep), ]; sub$v <- enc$min_layer_deep[!is.na(enc$min_layer_deep)]
print(aggregate(y ~ arm + v, sub,
      function(z) sprintf("%.0f%% (%d/%d)", 100*mean(z), sum(z), length(z))), row.names = FALSE)

res <- rbind(res, data.frame(encoding = "Integration_continuous", type = "continuous",
             n = nrow(d), df = ac$Df[2], deviance = round(ac$Deviance[2], 3),
             interaction_p = signif(ac$`Pr(>Chi)`[2], 3)))
res$bonferroni_0.0056 <- ifelse(res$interaction_p < 0.0056, "clears", "")
write.csv(res, file.path(OUT, "ispy2_categorical_interactions.csv"), row.names = FALSE)
cat("\nWrote", file.path(OUT, "ispy2_categorical_interactions.csv"), "\n")
