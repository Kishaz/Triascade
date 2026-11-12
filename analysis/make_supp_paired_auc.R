# =============================================================================
# Pooled paired AUC comparison: ICI Integration vs each established biomarker
# (reviewer comment C50 / R17). The per-cohort paired DeLong differences
# (ICI_benchmark_paired_delong.csv) are meta-analyzed per comparator so we can
# distinguish Integration's NUMERICAL advantage (higher AUC in most cohorts)
# from a statistically significant paired difference. Random-effects (logit-free,
# on the AUC-difference scale), with the per-cohort SE backed out of the DeLong P.
# =============================================================================
D <- read.csv("../results/tables/ICI_benchmark_paired_delong.csv", stringsAsFactors = FALSE)
D$se <- abs(D$dAUC) / qnorm(1 - D$p/2)                      # DeLong SE from two-sided P
D$se[!is.finite(D$se) | D$se == 0] <- NA

dl_pool <- function(d){                                     # DerSimonian-Laird on dAUC
  d <- d[is.finite(d$dAUC) & is.finite(d$se) & d$se > 0, ]
  y <- d$dAUC; s2 <- d$se^2; w <- 1/s2
  yf <- sum(w*y)/sum(w); Q <- sum(w*(y-yf)^2); k <- nrow(d)
  C <- sum(w) - sum(w^2)/sum(w); tau2 <- max(0,(Q-(k-1))/C)
  wr <- 1/(s2+tau2); yr <- sum(wr*y)/sum(wr); se <- sqrt(1/sum(wr))
  data.frame(pooled_dAUC = yr, lo = yr-1.96*se, hi = yr+1.96*se,
             p = 2*pnorm(-abs(yr/se)),
             n_cohorts = k, n_higher = sum(d$dAUC > 0), n_sig = sum(d$p < 0.05),
             row.names = NULL) }

COMP <- c(GEP="T-cell-inflamed GEP", IFNG="IFN-gamma signature",
          CYT="Cytolytic activity", CD274="PD-L1 mRNA (CD274)")
TAB <- do.call(rbind, lapply(names(COMP), function(c){
  r <- dl_pool(D[D$comp == c, ]); cbind(comparator = COMP[[c]], r) }))
num <- c("pooled_dAUC","lo","hi"); TAB[num] <- round(TAB[num], 3); TAB$p <- signif(TAB$p, 2)
out <- "../results/tables/ICI_SuppTable_paired_auc_pooled.csv"
write.csv(TAB, out, row.names = FALSE)
cat("wrote", out, "\n\n")
cat("Integration vs comparator, pooled paired AUC difference (random-effects):\n")
print(TAB[, c("comparator","pooled_dAUC","lo","hi","p","n_higher","n_sig","n_cohorts")], row.names = FALSE)
cat("\nn_higher = cohorts where Integration AUC is higher; n_sig = cohorts with paired DeLong P<0.05\n")
