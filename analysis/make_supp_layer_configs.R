# =============================================================================
# Joint layer-state configurations (reviewer comment C110, Fig 6b): observed vs
# expected-under-independence frequencies for the eight Sensory/Integration/
# Execution high-low combinations, with standardized (Pearson) residuals, so the
# specific configurations driving the significant chi-square can be identified.
# "Forward-concordant" = states whose set of high layers is a prefix of the
# ordered cascade Sensory -> Integration -> Execution (nested: none, S, S+I,
# S+I+E); all other states are discordant. Reuses ICI_cascade_patterns.csv.
# =============================================================================
D <- read.csv("../results/tables/ICI_cascade_patterns.csv", stringsAsFactors = FALSE)
D$pattern <- sprintf("%03d", as.integer(D$pattern))         # keep SIE bit order
D$std_residual <- round((D$observed - D$expected) / sqrt(D$expected), 2)
D$expected <- round(D$expected, 1)
D <- D[order(-D$std_residual), c("pattern","type","observed","expected","std_residual")]
names(D)[1] <- "SIE_pattern"
out <- "../results/tables/ICI_SuppTable_layer_configs.csv"
write.csv(D, out, row.names = FALSE)
cat("wrote", out, "  (chi-square total =", round(sum((D$observed-D$expected)^2/D$expected),1), ")\n\n")
print(D, row.names = FALSE)
cat("\nEnriched (obs > exp): forward-concordant 111 and 000; depleted: discordant states.\n")
