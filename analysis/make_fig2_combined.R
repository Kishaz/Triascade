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

# ---- step 1: binarization number line ---------------------------------------
axy <- 102; medx <- 40
pts2 <- data.frame(x = c(16,23,30,37, 46,53,60,67) + OFF, y = axy,
                   fill = c(rep("white",4), rep(NEUTRAL,4)))

lab_bin <- "b[g]==bgroup('{', atop(1*','~~'if'~~x[g]>=m[g], 0*','~~'otherwise'), '.')"
S_lab <- "bold(S)==frac(1,21)~sum(b[g], g %in% 'Sensory')"
I_lab <- "bold(I)==frac(1,21)~sum(b[g], g %in% 'activating') - frac(1,20)~sum(b[g], g %in% 'suppressive')"
E_lab <- "bold(E)==frac(1,19)~sum(b[g], g %in% 'effector') - frac(1,19)~sum(b[g], g %in% 'exhaustion')"

b2 <- list(
  # panel tag
  annotate("text", x = OFF, y = 129, label = "B", family = PUB_FONT, size = PT_PANEL_TAG/.pt, fontface = 2, hjust = 0),

  # ========== STEP 1 : binarize ==========
  annotate("point", x = 5+OFF, y = 124, size = 5.2, colour = INK),
  tb(5, 124, "1", size = 8, face = 2, colour = "white"),
  tb(11, 125,   "Binarize each gene at its", size = 9.5, face = 2, hjust = 0),
  tb(11, 120.3, "cohort median",             size = 9.5, face = 2, hjust = 0),

  # number line + median
  annotate("segment", x = 10+OFF, xend = 78+OFF, y = axy, yend = axy, colour = INK, linewidth = 0.4,
           arrow = arrow(length = unit(1.4,"mm"), type = "closed")),
  annotate("text", x = 79+OFF, y = axy, label = "expression", angle = 90,
           family = PUB_FONT, size = 8*FSCALE/.pt, colour = INK, hjust = 0.5),
  annotate("segment", x = medx+OFF, xend = medx+OFF, y = axy-5, yend = axy+5, colour = ACCENT, linewidth = 0.5, linetype = "22"),
  tb(medx, axy+8,   "m[g]", size = 8.5, colour = ACCENT, parse = TRUE),
  tb(medx, axy+11.5,"cohort median", size = 7.3),
  geom_point(data = pts2, aes(x = x, y = y), shape = 21, size = 3.0, fill = pts2$fill, colour = INK, stroke = 0.5),
  tb(25, axy-5, "x[g] < m[g]", size = 7.3, parse = TRUE),
  tb(25, axy-9, "b[g]==0",     size = 8, face = 2, parse = TRUE),
  tb(55, axy-5, "x[g] >= m[g]", size = 7.3, parse = TRUE),
  tb(55, axy-9, "b[g]==1",      size = 8, face = 2, parse = TRUE),

  # indicator formula
  tb(40, 85, lab_bin, size = 10, parse = TRUE),

  # legend
  tb(2, 76, "x[g]*' = gene expression in the tumor'", size = 8, hjust = 0, parse = TRUE),
  tb(2, 71.5, "m[g]*' = cohort median'",     size = 8, hjust = 0, parse = TRUE),
  tb(2, 67, "b[g]*' = 1 (high) or 0 (low)'",           size = 8, hjust = 0, parse = TRUE),

  # divider between the two steps
  annotate("segment", x = 2+OFF, xend = 86+OFF, y = 61, yend = 61, colour = "grey80", linewidth = 0.4),

  # ========== STEP 2 : sum within each layer ==========
  annotate("point", x = 5+OFF, y = 56, size = 5.2, colour = INK),
  tb(5, 56, "2", size = 8, face = 2, colour = "white"),
  tb(11, 56, "Combine within each layer", size = 9.5, face = 2, hjust = 0),

  tb(3, 46, S_lab, size = 10, hjust = 0, parse = TRUE, colour = unname(col["Sense"])),
  tb(3, 31, I_lab, size = 10, hjust = 0, parse = TRUE, colour = unname(col["Integration"])),
  tb(3, 16, E_lab, size = 10, hjust = 0, parse = TRUE, colour = unname(col["Execution"])),

  tb(3, 7.5, "S %in% group('[',list(0,1),']')*';   '*I*','~E %in% group('[',list(-1,+1),']')",
      size = 8, hjust = 0, parse = TRUE))

# ---- Panel C: validation, sourced as layers and shifted BELOW the A/B row -------
# pA and b2 are already materialized, so it is safe for 2C to redefine the helpers.
FIG2C_XOFF     <- 0                      # left-align C with A
FIG2C_YOFF     <- -95                    # drop C below the A/B row
FIG2C_ASLAYERS <- TRUE                   # expose C_layers; do not save 2C on its own
source("make_fig2C_validation.R")        # -> C_layers

# ---- per-panel frames, so the three panels never read as one continuous figure --
# Each box encloses one panel's content with padding; a gutter separates A | B and
# a clear band separates the top row from C.
frames <- list(
  annotate("rect", xmin =  -2, xmax = 184.5, ymin =   1.5, ymax = 131.5, fill = NA, colour = INK, linewidth = 0.7),  # A
  annotate("rect", xmin = 186.5, xmax = 276, ymin =   1.5, ymax = 131.5, fill = NA, colour = INK, linewidth = 0.7),  # B
  annotate("rect", xmin =  -2, xmax = 276,   ymin = -70,   ymax =  -4,   fill = NA, colour = INK, linewidth = 0.7))  # C

fig2 <- pA + b2 + C_layers + frames +
  coord_fixed(ratio = 1, xlim = c(-4, 278), ylim = c(-72, 133.5), expand = FALSE, clip = "off")

# This combined figure is wider and taller than one journal page on purpose
# (2A + a vertical 2B column, with 2C below); lift save_fig's double-column width
# cap for it only. Scale is kept identical to 2A standalone: 186 units -> 174 mm.
FIG_W_DOUBLE <- 320
save_fig(fig2, "ICI_Fig2", width_mm = 282*174/186, height_mm = 205.5*174/186, dir = "../results/figures/composite")
cat("done\n")
