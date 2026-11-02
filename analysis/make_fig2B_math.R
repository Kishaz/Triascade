# =============================================================================
# Figure 2B - HOW the layer scores are calculated (LOCKED v1 formula).
#   Step 1  binarize each panel gene at its locked reference median:
#             b_g = 1 if x_g >= m_g, else 0            (indicator function)
#   Step 2  each layer score = normalized sum of its b_g (LEGACY 98-gene panel):
#             S = (1/21) * sum_{g in Sensory}     b_g
#             I = (1/21) * sum_{g in activating}  b_g  -  (1/20) * sum_{g in suppressive} b_g
#             E = (1/19) * sum_{g in effector}    b_g  -  (1/19) * sum_{g in exhaustion}  b_g
#   three scores read SEPARATELY -> no composite.  S in [0,1]; I,E in [-1,+1].
# House style: Arial, warm layer palette, geometric lines (ici_theme.R).
# Formulas are plotmath -> true Sigma / fractions, kept as editable SVG text.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R")
.pt <- 2.845276

# ---- worked example (illustrative counts; console sanity check only) --------
# legacy panel sizes: Sensory 21, Integration +21/-20, Execution +19/-19
Sv<-18; Ipv<-15; Inv<-8; Epv<-14; Env<-6
S_sc<-Sv/21; I_sc<-Ipv/21-Inv/20; E_sc<-Epv/19-Env/19

col <- pal_layer    # Sense / Integration / Execution deep tokens

# step badges
badge <- data.frame(x=c(4,96), y=c(74,74), n=c("1","2"))

# ---- binarization schematic (left) ------------------------------------------
axy <- 56; ax0 <- 10; ax1 <- 70; medx <- 40
pts <- data.frame(
  x = c(15,22,29,36, 45,52,59,66),
  hi = c(0,0,0,0, 1,1,1,1))
pts$fill <- ifelse(pts$hi==1, NEUTRAL, "white")

# ---- plotmath formula labels ------------------------------------------------
lab_bin <- "b[g]==bgroup('{', atop(1*','~~'if'~~x[g]>=m[g], 0*','~~'otherwise'), '.')"
forms <- data.frame(
  x = 100, y = c(58, 42, 26),
  col = unname(col[c("Sense","Integration","Execution")]),
  lab = c(
    "bold(S)==frac(1,21)~sum(b[g], g %in% 'Sensory')",
    "bold(I)==frac(1,21)~sum(b[g], g %in% 'activating') - frac(1,20)~sum(b[g], g %in% 'suppressive')",
    "bold(E)==frac(1,19)~sum(b[g], g %in% 'effector') - frac(1,19)~sum(b[g], g %in% 'exhaustion')"))

p <- ggplot() +
  # panel tag
  annotate("text", x = 1, y = 79, label = "B", family = PUB_FONT, size = PT_PANEL_TAG/.pt, fontface = 2, hjust = 0) +
  # step badges + headers
  geom_point(data = badge, aes(x = x, y = y), shape = 21, size = 5.4, fill = INK, colour = INK) +
  geom_text(data = badge, aes(x = x, y = y, label = n), family = PUB_FONT, size = 8/.pt, fontface = 2, colour = "white") +
  annotate("text", x = 9,  y = 74, hjust = 0, family = PUB_FONT, size = 10.5/.pt, fontface = 2, colour = INK,
           label = "Binarize each gene at its cohort median") +
  annotate("text", x = 101, y = 74, hjust = 0, family = PUB_FONT, size = 10.5/.pt, fontface = 2, colour = INK,
           label = "Combine within each layer, per gene count") +

  # --- binarization schematic ---
  annotate("segment", x = ax0, xend = ax1, y = axy, yend = axy, colour = INK, linewidth = 0.4,
           arrow = arrow(length = unit(1.4,"mm"), type = "closed")) +
  annotate("text", x = ax1+1.5, y = axy, hjust = 0, family = PUB_FONT, size = 8/.pt, colour = INK, label = "expression") +
  annotate("segment", x = medx, xend = medx, y = axy-6, yend = axy+6, colour = ACCENT, linewidth = 0.5, linetype = "22") +
  annotate("text", x = medx, y = axy+8.5, family = PUB_FONT, size = 8.5/.pt, colour = ACCENT,
           label = "m[g]", parse = TRUE) +
  annotate("text", x = medx, y = axy+12, family = PUB_FONT, size = 7.3/.pt, colour = INK, label = "cohort median") +
  geom_point(data = pts, aes(x = x, y = axy), shape = 21, size = 3.1, fill = pts$fill, colour = INK, stroke = 0.5) +
  # each point's binary call: below median -> 0 (open), at/above -> 1 (filled)
  annotate("text", x = 25, y = axy-5,   family = PUB_FONT, size = 7.3/.pt, colour = INK, parse = TRUE, label = "x[g] < m[g]") +
  annotate("text", x = 25, y = axy-9,   family = PUB_FONT, size = 8/.pt,   colour = INK, parse = TRUE, fontface = 2, label = "b[g]==0") +
  annotate("text", x = 55, y = axy-5,   family = PUB_FONT, size = 7.3/.pt, colour = INK, parse = TRUE, label = "x[g] >= m[g]") +
  annotate("text", x = 55, y = axy-9,   family = PUB_FONT, size = 8/.pt,   colour = INK, parse = TRUE, fontface = 2, label = "b[g]==1") +
  # indicator formula
  annotate("text", x = 40, y = 34, family = PUB_FONT, size = 12/.pt, colour = INK, label = lab_bin, parse = TRUE) +
  # symbol legend
  annotate("text", x = 4, y = 20, hjust = 0, family = PUB_FONT, size = 8.5/.pt, colour = INK, parse = TRUE,
           label = "x[g]*' = gene expression in the tumor'") +
  annotate("text", x = 4, y = 15.5, hjust = 0, family = PUB_FONT, size = 8.5/.pt, colour = INK, parse = TRUE,
           label = "m[g]*' = cohort median,  '*b[g]*' = 1 (high) or 0 (low)'") +

  # --- divider ---
  annotate("segment", x = 91, xend = 91, y = 10, yend = 70, colour = "grey80", linewidth = 0.4) +

  # --- layer-score formulas ---
  geom_text(data = forms, aes(x = x, y = y, label = lab), parse = TRUE, hjust = 0,
            family = PUB_FONT, size = 12/.pt, colour = forms$col) +

  # footer: ranges + no composite (right col) + worked example (full-width bottom)
  annotate("text", x = 101, y = 15.5, hjust = 0, family = PUB_FONT, size = 8.5/.pt, colour = INK, parse = TRUE,
           label = "'N = 21 / 21,20 / 19,19 genes.   '*S %in% group('[',list(0,1),']')*';   '*I*','~E %in% group('[',list(-1,+1),']')") +

  coord_fixed(ratio = 1, xlim = c(0, 185), ylim = c(9, 82), expand = FALSE, clip = "off") +
  theme_void() + theme(plot.background = element_rect(fill = "white", colour = NA),
                       plot.margin = margin(1,1,1,1,"mm"))

save_fig(p, "ICI_Fig2B_math", width_mm = 174, height_mm = 174*73/185)
cat(sprintf("S=%.3f  I=%.3f  E=%.3f\n", S_sc, I_sc, E_sc))
cat("done\n")
