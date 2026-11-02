# =====================================================================
SCRIPT_VERSION <- "ICI_ISPY2_LRT  v1.0  2026-10-28"
# Every I-SPY2 interaction p-value by likelihood-ratio test, at n = 248
# (246 with Integration), so the paper uses one test throughout (NeoTRIP,
# PRINCE and Table S4 already use it). Uses YOUR scripts:
#   ispy2_categorical_check.R  -> Supplementary Table S4 (all encodings), rerun on n = 248
#   fig2_ispy2_predictive.R    -> layer x arm tests, ORs by arm, receptor subgroups,
#                                 tertiles (its own LRT function), rerun on n = 248
# and adds, with the same likelihood-ratio test: Composite, the median split,
# the subgroup interaction OR, the three-way test, and the old 189-patient subset.
# Prints each text value next to its likelihood-ratio value.
# Outputs: Oct 16/source_data/rerun_split_v1/ispy2_LRT; new Table S4 in
# supplementary_tables (old one moved to supplementary_tables/_removed).
# MASTER_GRID files are not changed.
# =====================================================================
suppressPackageStartupMessages(library(data.table))
BASE <- Sys.getenv("ICI_BASE", "/Users/dbedi/Desktop/R data FINAL")
MAN  <- file.path(BASE, "ICI_Manuscript"); GRID <- file.path(MAN, "MASTER_GRID")
OCT  <- file.path(MAN, "Oct 16"); ST <- file.path(OCT, "supplementary_tables")
W    <- file.path(OCT, "source_data", "rerun_split_v1", "ispy2_LRT"); dir.create(W, showWarnings = FALSE, recursive = TRUE)
cat("\n", SCRIPT_VERSION, "\n", sep = "")
N248 <- file.path(GRID, "ICI_binary_ispy2_scores_246complete.csv"); OLD189 <- file.path(GRID, "ICI_binary_ispy2_scores.csv")
FC <- file.path(MAN, "ispy2_categorical_check.R"); FF <- file.path(MAN, "Oct4_ICI", "scripts", "fig2_ispy2_predictive.R")
for (f in c(N248, OLD189, FC, FF)) if (!file.exists(f)) stop("not found: ", f)

run_txt <- function(file, subs, tag) {
  txt <- paste(readLines(file, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (BASE != "/Users/dbedi/Desktop/R data FINAL")                       # test runs only
    subs <- c(list(c('"/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"', deparse(MAN), 1)), subs)
  for (s in subs) {
    k <- lengths(regmatches(txt, gregexpr(s[1], txt, fixed = TRUE)))
    if (k != as.integer(s[3])) stop(sprintf("%s: expected %s match(es) of [%s], found %d - paste this to Claude", tag, s[3], s[1], k))
    txt <- gsub(s[1], s[2], txt, fixed = TRUE) }
  log <- file.path(W, paste0(tag, ".log")); zz <- file(log, open = "wt"); sink(zz); sink(zz, type = "message")
  e <- new.env()
  res <- tryCatch({ eval(parse(text = txt, keep.source = FALSE, encoding = "UTF-8"), envir = e); NULL }, error = function(err) err)
  sink(type = "message"); sink(); close(zz)
  if (!is.null(res)) { cat("ERROR in", tag, ":", conditionMessage(res), "\n"); cat(tail(readLines(log, warn = FALSE), 15), sep = "\n"); stop("paste this to Claude") }
  e }
lrt <- function(m0, m1) anova(m0, m1, test = "Chisq")$`Pr(>Chi)`[2]
fmt <- function(p) ifelse(p < 0.001, formatC(p, format = "e", digits = 1), sprintf("%.3f", p))

# ---- 1. Table S4: your categorical script on n = 248 --------------------------
cat1 <- run_txt(FC, list(
  c('F   <- file.path(MAN, "MASTER_GRID", "ICI_binary_ispy2_scores.csv")', paste0("F   <- ", deparse(N248)), 1),
  c('OUT <- file.path(MAN, "ISPY2", "CATEGORICAL")', paste0("OUT <- ", deparse(W)), 1)), "table_S4")
S4 <- cat1$res
cat(sprintf("\n==== Table S4 (likelihood-ratio), n = %d ====\n", max(S4$n))); print(S4[, c("encoding", "type", "n", "df", "interaction_p")], row.names = FALSE)
dir.create(file.path(ST, "_removed"), showWarnings = FALSE)
for (f in list.files(ST, pattern = "^Supp_Table_S4", full.names = TRUE)) file.rename(f, file.path(ST, "_removed", basename(f)))
file.copy(file.path(W, "ispy2_categorical_interactions.csv"), file.path(ST, "Supp_Table_S4_ispy2_interactions_LRT_n248.csv"), overwrite = TRUE)

# ---- 2. your Figure 3 script on n = 248 ------------------------------------------
f3 <- run_txt(FF, list(
  c('SCOR <- file.path(MAN,"MASTER_GRID","ICI_binary_ispy2_scores.csv")', paste0("SCOR <- ", deparse(N248)), 1),
  c('OUT  <- file.path(MAN,"Aug_1_ICI_FINAL_PAPER","Figures")', paste0("OUT  <- ", deparse(W)), 1)), "figure3_script")
d <- f3$d; inter <- f3$inter
cat(sprintf("\nI-SPY2: n = %d (control %d, pembro %d); receptor subtype found for %d\n", nrow(d),
            sum(d$arm == "control"), sum(d$arm == "pembro"), sum(!is.na(d$subtype))))

# ---- 3. the text's I-SPY2 tests, likelihood-ratio -------------------------------
cc <- intersect(c("composite", "Composite", "ICI_composite"), names(d))[1]
d$Composite_z <- as.numeric(scale(d[[cc]]))
row <- function(what, now, r) data.frame(test = what, text_now = now, OR_control = round(r$or_c, 2),
                                         OR_pembro = round(r$or_p, 2), LRT_p = fmt(r$p), stringsAsFactors = FALSE)
Tt <- rbind(row("Integration x arm, per SD", "OR 3.16 vs 1.31; p = 0.025", inter(d, "Integration_z")),
            row("Composite x arm, per SD",   "OR 4.00 vs 1.61; p = 0.032", inter(d, "Composite_z")),
            row("Sensory x arm, per SD",     "OR 3.10 vs 1.54; p = 0.068", inter(d, "Sensory_z")),
            row("Execution x arm, per SD",   "OR 2.24 vs 1.94; p = 0.71",  inter(d, "Execution_z")))
d$Int_high <- factor(ifelse(d$Integration >= median(d$Integration), "high", "low"), levels = c("low", "high"))
Tt <- rbind(Tt, row("Integration median split x arm", "p = 0.026", inter(d, "Int_high")))
Tt <- rbind(Tt, data.frame(test = "Integration tertiles x arm", text_now = "p = 0.018", OR_control = NA, OR_pembro = NA,
                           LRT_p = fmt(f3$src$tertile_p)))
cat("\n==== I-SPY2 interaction tests (text now -> likelihood-ratio) ====\n"); print(Tt, row.names = FALSE)
D <- f3$src$panel_d; cat("\npCR by Integration tertile (text: pembro 11/46/68%, control 14/16/20%):\n")
print(data.frame(tertile = D$tert, arm = D$arm, pCR = sprintf("%.0f%%", 100 * D$rate), n = D$n), row.names = FALSE)

# receptor subgroups: interaction OR with profile-likelihood CI, and the three-way test
if (!"subtype" %in% names(d)) cat("\nreceptor subtype not found (series matrices missing in ISPY2/) - subgroup tests skipped\n") else {
cat("\n==== receptor subgroups (text: TNBC n = 112, OR 4.45 (1.32-14.98), p = 0.016; HR+ n = 134, OR 1.22, p = 0.71; three-way p = 0.105) ====\n")
for (st in c("TNBC", "HR+/HER2-")) {
  s <- d[!is.na(d$subtype) & d$subtype == st, ]; s$Iz <- as.numeric(scale(s$Integration))
  m1 <- glm(y ~ Iz * arm, s, family = binomial); m0 <- glm(y ~ Iz + arm, s, family = binomial)
  k <- grep(":", names(coef(m1))); ci <- suppressMessages(confint(m1))[k, ]
  cat(sprintf("  %-10s n = %3d  interaction OR %.2f (95%% CI %.2f-%.2f)  LRT p = %s\n", st, nrow(s),
              exp(coef(m1)[k]), exp(ci[1]), exp(ci[2]), fmt(lrt(m0, m1)))) }
h <- d[!is.na(d$subtype) & d$subtype %in% c("TNBC", "HR+/HER2-"), ]
h$X <- as.numeric(scale(h$Integration)); h$S <- factor(h$subtype, levels = c("TNBC", "HR+/HER2-"))
cat(sprintf("  three-way Integration x arm x receptor (n = %d): LRT p = %s\n", nrow(h),
            fmt(lrt(glm(y ~ X * arm + X * S + arm * S, h, family = binomial), glm(y ~ X * arm * S, h, family = binomial)))))
}

# the earlier 189-patient subset ("binarized p = 0.008")
o <- read.csv(OLD189, stringsAsFactors = FALSE); o$y <- as.integer(o$pcr_binary)
o$arm <- factor(o$arm, levels = c("control", "pembro")); o <- o[!is.na(o$y) & !is.na(o$arm) & !is.na(o$Integration), ]
o$Ih <- as.integer(o$Integration >= median(o$Integration)); o$Iz <- as.numeric(scale(o$Integration))
cat(sprintf("\n==== earlier 189-patient subset (n = %d; text: binarized p = 0.008) ====\n  median split LRT p = %s | per SD LRT p = %s\n", nrow(o),
            fmt(lrt(glm(y ~ Ih + arm, o, family = binomial), glm(y ~ Ih * arm, o, family = binomial))),
            fmt(lrt(glm(y ~ Iz + arm, o, family = binomial), glm(y ~ Iz * arm, o, family = binomial)))))

fwrite(Tt, file.path(W, "ispy2_LRT_text_values.csv"))
cat("\nnew Table S4: supplementary_tables/Supp_Table_S4_ispy2_interactions_LRT_n248.csv (old one in supplementary_tables/_removed)\n")
cat("Figure 3 script output (n = 248) for reference:", W, "\nDONE - paste everything above to Claude.\n")
