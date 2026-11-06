# =============================================================================
# Figure 3, opening panel: "from a single transcriptome to a response prediction."
# A generic flow schematic (illustration) paired with Fig 3's OWN pooled real data
# (reinforcement), de-duplicated so the randomized interaction stays in Fig 5 and
# the archetype breakdown stays in Fig 7.
#   left  : transcriptome -> score vs frozen reference -> S/I/E layer scores ->
#           Integration score -> predicted R / NR
#   right : pooled responder-vs-non-responder Integration (within-cohort z) across
#           the 11 validated ICB cohorts, annotated with the random-effects pooled
#           AUC (NOT a pooled-ROC AUC, which cohort offsets would distort).
# Numbers read from pancancer_validation.rds + ICI_benchmark_pooled_auc.csv; never typed.
# House style via ici_theme.R; all text black. Build to an editable SVG for Inkscape
# assembly into the Fig 3 composite (intended panel a).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
setwd("/Users/smwambu1/Desktop/Projects/ICI/analysis")
source("ici_theme.R"); .pt <- 2.845276
col <- pal_layer
cS <- unname(col["Sense"]); cI <- unname(col["Integration"]); cE <- unname(col["Execution"])

# ---- real data: pooled within-cohort-z Integration, R vs NR (11 validated cohorts) ----
R <- readRDS("pancancer_validation.rds")
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
dd <- do.call(rbind, lapply(WORKS, function(k){ s <- R[[k]]$scores
  z <- (s$integ_single - mean(s$integ_single)) / sd(s$integ_single)
  data.frame(z = z, resp = factor(ifelse(s$y == 1, "R", "NR"), c("NR", "R"))) }))
ic <- read.csv("../results/tables/ICI_benchmark_pooled_auc.csv", check.names = FALSE)
ic <- ic[ic$biom == "ICI", ]
auclab <- sprintf("Pooled AUC %.2f (%.2f–%.2f)", ic$auc, ic$lo, ic$hi)
nlab   <- sprintf("%d cohorts · %d patients", ic$k, nrow(dd))

# ---- flow schematic (left) ----
txt  <- function(x, y, lab, size = 9, face = 1, colour = INK, hjust = 0.5)
  annotate("text", x = x, y = y, label = lab, family = PUB_FONT, size = size/.pt, fontface = face, colour = colour, hjust = hjust)
rect <- function(xmin, xmax, ymin, ymax, fill, colour, lwd = 0.6)
  annotate("rect", xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill, colour = colour, linewidth = lwd)
arrw <- function(x, xend, y) annotate("segment", x = x, xend = xend, y = y, yend = y,
  colour = INK, linewidth = 0.7, arrow = arrow(length = unit(1.4, "mm"), type = "closed"))
seg  <- function(x, xend, y, colour) annotate("segment", x = x, xend = xend, y = y, yend = y, colour = colour, linewidth = 0.5)

flow <- list(
  annotate("text", x = 0, y = 36, label = "a", family = PUB_FONT, size = PT_PANEL_TAG/.pt, fontface = 2, hjust = 0),
  # 1 transcriptome
  rect(2, 33, 15, 33, CREAM, MUTED), txt(17.5, 27, "bulk", 10), txt(17.5, 21.5, "transcriptome", 9.2),
  arrw(33, 41, 24),
  # 2 score vs frozen reference (the single-sample step)
  rect(41, 77, 14, 34, "white", MUTED, 0.5), txt(59, 29, "score vs", 9.5, face = 3),
  txt(59, 24, "frozen reference", 9.5, face = 3), txt(59, 19, "(single sample)", 7.8, face = 3, colour = MUTED),
  arrw(77, 85, 24),
  # 3 S/I/E layer scores (mini bars)
  rect(85, 120, 13, 35, "white", MUTED, 0.5), txt(102.5, 31.5, "layer scores", 9.5, face = 3),
  txt(88, 26, "S", 10, face = 2, colour = cS), txt(88, 22, "I", 10, face = 2, colour = cI), txt(88, 18, "E", 10, face = 2, colour = cE),
  seg(91, 116, 26, cS),   annotate("point", x = 108, y = 26, size = 1.9, colour = cS),
  seg(91, 116, 22, cI),   annotate("point", x = 100, y = 22, size = 1.9, colour = cI),
  seg(91, 116, 18, cE),   annotate("point", x = 112, y = 18, size = 1.9, colour = cE),
  arrw(120, 128, 24),
  # 4 Integration (emphasized)
  rect(128, 160, 15, 33, "white", cI, 0.9), txt(144, 26.5, "Integration", 10, face = 2, colour = cI),
  txt(144, 21, "score", 9.5, face = 2, colour = cI),
  arrw(160, 168, 24),
  # 5 predicted R / NR
  txt(181, 27, "predicted", 9.5, hjust = 0.5), txt(181, 21, "R / NR", 11, face = 2, colour = cI, hjust = 0.5)
)
pFlow <- ggplot() + flow +
  coord_fixed(ratio = 1, xlim = c(0, 197), ylim = c(12, 37), expand = FALSE, clip = "off") +
  theme_void() + theme(plot.margin = margin(2, 1, 2, 2, "mm"))

# ---- pooled real data (right) ----
pData <- ggplot(dd, aes(resp, z, fill = resp)) +
  geom_violin(colour = NA, width = 0.92, alpha = 0.85) +
  geom_boxplot(width = 0.15, fill = "white", colour = INK, linewidth = 0.4, outlier.shape = NA) +
  scale_fill_manual(values = c(NR = NEUTRAL, R = cI), guide = "none") +
  scale_y_continuous("Integration (z within cohort)", breaks = seq(-2, 2, 1)) +
  scale_x_discrete(NULL, labels = c(NR = "NR", R = "R")) +
  labs(subtitle = paste0(auclab, "\n", nlab)) +
  theme_publication(grid = "y") +
  theme(axis.text.x = element_text(size = PT_BASE + 1, colour = INK, face = 2),
        axis.text.y = element_text(size = PT_BASE, colour = INK),
        axis.title.y = element_text(size = PT_AXIS_TITLE, colour = INK),
        plot.subtitle = element_text(size = PT_BASE - 0.5, colour = INK, hjust = 0.5, margin = margin(b = 2)),
        plot.margin = margin(2, 3, 2, 3, "mm"))

p <- pFlow + pData + plot_layout(widths = c(2.15, 1))
FIG_W_DOUBLE <- 320                       # lift save_fig's double-column width cap for this wide panel
save_fig(p, "ICI_Fig3_flow", width_mm = 234, height_mm = 58, dir = "../results/figures/composite")
cat(sprintf("done - ICI_Fig3_flow: pooled %d patients (R=%d, NR=%d), %s\n",
            nrow(dd), sum(dd$resp == "R"), sum(dd$resp == "NR"), auclab))
