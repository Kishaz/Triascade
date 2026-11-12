# =============================================================================
# Figure 2 (combined): top row = 2A overview (LEFT) + 2B math stacked VERTICALLY
# (right); bottom row = 2C validation, spanning the width below them.  All drawn on
# ONE shared canvas so the coordinate scale (0.935 mm/unit) is identical everywhere.
#   2A occupies local x 0-186 (sourced, unchanged).
#   2B occupies local x OFF..OFF+88, two steps stacked: binarize (top) / sum (bottom).
#   2C is sourced as a layer list and shifted DOWN below the A/B row (y-offset).
# House style and helpers (txt/rect, Arial, warm palette) come from 2A via source().
# RULE: all text BLACK (INK) except layer-colored formulas; no em-dashes; editable SVG.
# =============================================================================
FIG2_NOSAVE <- TRUE                     # source 2A for pA without re-saving it
source("make_fig2A_overview.R")         # -> pA, and helpers txt()/rect()/col/.pt ...

OFF <- 188                              # x-offset of the 2B column (tight gutter to 2A)
tb  <- function(x, y, ...) txt(x + OFF, y, ...)

# ---- geometry: three stacked steps in the 2B column -------------------------
axy <- 120; thx <- 40                     # number-line y ; frozen-threshold x
pts2 <- data.frame(x = c(16,23,30,37, 46,53,60,67) + OFF, y = axy,
                   fill = c(rep("white",4), rep(NEUTRAL,4)))

lab_bin <- "b[g]==bgroup('{', atop(1*','~~'if'~~z[g]>=t[g], 0*','~~'otherwise'), '.')"
S_lab <- "bold(S)==sum(b[g], g %in% 'Sensory')"
I_lab <- "bold(I)==sum(b[g], g %in% 'activating') - sum(b[g], g %in% 'suppressive')"
E_lab <- "bold(E)==sum(b[g], g %in% 'effector') - sum(b[g], g %in% 'exhaustion')"
std_lab <- "bold(L)*'*' == frac(bold(L) - mu[L], sigma[L])"

b2 <- list(
  # panel tag
  annotate("text", x = OFF, y = 148, label = "b", family = PUB_FONT, size = PT_PANEL_TAG/.pt, fontface = 2, hjust = 0),

  # ========== STEP 1 : flag each gene vs the frozen reference ==========
  annotate("point", x = 5+OFF, y = 144, size = 5.2, colour = INK),
  tb(5, 144, "1", size = 8, face = 2, colour = "white"),
  tb(11, 145,   "Flag each gene against", size = 9.5, face = 2, hjust = 0),
  tb(11, 140.3, "the frozen reference",   size = 9.5, face = 2, hjust = 0),

  # number line + frozen-reference threshold
  annotate("segment", x = 10+OFF, xend = 78+OFF, y = axy, yend = axy, colour = INK, linewidth = 0.4,
           arrow = arrow(length = unit(1.4,"mm"), type = "closed")),
  annotate("text", x = 80+OFF, y = axy, label = "standardized~expression", angle = 90,
           family = PUB_FONT, size = 7.2*FSCALE/.pt, colour = INK, hjust = 0.5, parse = TRUE),
  annotate("segment", x = thx+OFF, xend = thx+OFF, y = axy-5, yend = axy+5, colour = ACCENT, linewidth = 0.5, linetype = "22"),
  tb(thx, axy+8,   "t[g]", size = 8.5, colour = ACCENT, parse = TRUE),
  tb(thx, axy+11.5,"frozen reference", size = 7.3),
  geom_point(data = pts2, aes(x = x, y = y), shape = 21, size = 3.0, fill = pts2$fill, colour = INK, stroke = 0.5),
  tb(23, axy-6, "z[g] < t[g]  (b[g]==0)", size = 7.2, parse = TRUE),
  tb(57, axy-6, "z[g] >= t[g]  (b[g]==1)", size = 7.2, parse = TRUE),

  # indicator formula
  tb(40, 102, lab_bin, size = 10, parse = TRUE),

  # legend
  tb(2, 95, "z[g]*' = within-sample standardized expression'", size = 7.4, hjust = 0, parse = TRUE),
  tb(2, 91, "t[g]*' = frozen gene-specific reference threshold'", size = 7.4, hjust = 0, parse = TRUE),

  # divider
  annotate("segment", x = 2+OFF, xend = 86+OFF, y = 85, yend = 85, colour = "grey80", linewidth = 0.4),

  # ========== STEP 2 : count within each layer ==========
  annotate("point", x = 5+OFF, y = 79.5, size = 5.2, colour = INK),
  tb(5, 79.5, "2", size = 8, face = 2, colour = "white"),
  tb(11, 80.5, "Count within each layer", size = 9.5, face = 2, hjust = 0),
  tb(11, 76,   "(activating − suppressive)", size = 8, hjust = 0),

  tb(3, 66, S_lab, size = 9.5, hjust = 0, parse = TRUE, colour = unname(col["Sense"])),
  tb(3, 57, I_lab, size = 9.5, hjust = 0, parse = TRUE, colour = unname(col["Integration"])),
  tb(3, 48, E_lab, size = 9.5, hjust = 0, parse = TRUE, colour = unname(col["Execution"])),

  # divider
  annotate("segment", x = 2+OFF, xend = 86+OFF, y = 41, yend = 41, colour = "grey80", linewidth = 0.4),

  # ========== STEP 3 : standardize against the frozen reference ==========
  annotate("point", x = 5+OFF, y = 35, size = 5.2, colour = INK),
  tb(5, 35, "3", size = 8, face = 2, colour = "white"),
  tb(11, 35.5, "Standardize vs the frozen reference", size = 9.0, face = 2, hjust = 0),

  tb(22, 24, std_lab, size = 11, parse = TRUE),
  tb(2, 12, "mu[L]*', '*sigma[L]*' = frozen reference layer mean, SD'", size = 7.4, hjust = 0, parse = TRUE),
  tb(2, 7.5, "bold(L)*'* in reference SD units  (0 = reference mean)'", size = 7.4, hjust = 0, parse = TRUE))

# ---- per-panel frames: A | B (taller for legibility). A gutter separates A | B.
frames <- list(
  annotate("rect", xmin =  -2, xmax = 184.5, ymin = 1.5, ymax = 150, fill = NA, colour = INK, linewidth = 0.7),  # A
  annotate("rect", xmin = 186.5, xmax = 276, ymin = 1.5, ymax = 150, fill = NA, colour = INK, linewidth = 0.7))  # B

fig2 <- pA + b2 + frames +
  coord_fixed(ratio = 1, xlim = c(-4, 278), ylim = c(-1, 152), expand = FALSE, clip = "off")

# Two-panel top row (2A + a vertical 2B column); the 2C validation panel now lives in
# Fig 3. Lift save_fig's double-column width cap for it only.
FIG_W_DOUBLE <- 320
save_fig(fig2, "ICI_Fig2", width_mm = 282*174/186, height_mm = 153*174/186, dir = "../results/figures/composite")
cat("done\n")
