# =============================================================================
# ici_theme.R - one visual system for every ICI figure.
# Adapted from scMultiome prostate_multiome R/functions/theme_publication.R
# (author: Samuel Mwamburi). Self-contained copy per project convention.
#
# RULES (carried over, verbatim intent):
#  - Author at FINAL printed size, never scale: 90 mm single col, 180 mm double,
#    170 mm max height. 7 pt stays 7 pt because nothing is resized.
#  - Helvetica throughout. One knob PT_BASE; nothing below 5 pt.
#  - NO titles or subtitles - everything lives in axis labels, legend, caption.
#  - Okabe-Ito colour-blind-safe palette. Two-group contrast = PLUM vs OCHRE,
#    never red/blue (overused + unintended good/bad valence).
#  - INK black text/axes. Save PNG 600 dpi + editable SVG (svg_join_text.py).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(svglite) })

# Cell Press / Cancer Cell figure spec: Arial font; panels labelled UPPERCASE
# bold A, B, C; column widths 85 / 114 / 174 mm (1 / 1.5 / 2 col); max height 234 mm.
FIG_W_SINGLE <- 85; FIG_W_1_5 <- 114; FIG_W_DOUBLE <- 174; FIG_H_MAX <- 234; MM_TO_IN <- 1/25.4
PUB_FONT <- "Arial"
PT_BASE <- 10
PT_AXIS_TEXT <- PT_BASE; PT_AXIS_TITLE <- PT_BASE + 1; PT_LEGEND <- PT_BASE
PT_ANNOT <- PT_BASE - 1; PT_PANEL_TAG <- PT_BASE + 2
INK <- "black"; RULE <- "grey25"

OKABE_ITO <- c(orange="#E69F00", sky_blue="#56B4E9", bluish_green="#009E73",
               yellow="#F0E442", blue="#0072B2", vermillion="#D55E00",
               purple="#CC79A7", black="#000000")
# ---- ICI warm house palette: plum / pink / sand / cream-to-white ----
# (set by the user 2026-10-30). Deep tokens read on cream; light tints are fills.
PLUM <- "#6A2C5A"; PLUM_DEEP <- "#4C1F41"
ROSE <- "#B2567E"; PINK <- "#D98FB0"; PINK_LT <- "#F1C9DB"
SAND <- "#C9A66B"; SAND_DEEP <- "#9C7A3E"
CREAM <- "#F7F1E6"; CREAM_LT <- "#FBF7F0"
ACCENT  <- PLUM_DEEP          # reference lines / thresholds
NEUTRAL <- "#C9BEB0"          # warm sand-grey for unremarkable data
MUTED   <- "#8A7F71"

# cascade-layer palette (borders + titles use the deep tokens so they read on cream)
pal_layer      <- c(Sense = SAND_DEEP, Integration = PLUM,     Execution = ROSE)
pal_layer_fill <- c(Sense = "#F4ECDB", Integration = "#EAD8E5", Execution = "#F7E1EC")
# two-group contrast (arm / high-low): plum (emphasis) vs sand
pal_two <- c(low = SAND, high = PLUM, control = SAND, treat = PLUM)
# sequential ramp cream -> plum for continuous fills
EXPR_RAMP <- c(low = CREAM, high = PLUM)
scale_expr_fill   <- function(...) ggplot2::scale_fill_gradient(low=EXPR_RAMP[["low"]], high=EXPR_RAMP[["high"]], ...)
scale_expr_colour <- function(...) ggplot2::scale_colour_gradient(low=EXPR_RAMP[["low"]], high=EXPR_RAMP[["high"]], ...)

theme_publication <- function(base_size = PT_AXIS_TEXT, grid = c("none","y","x","both")) {
  grid <- match.arg(grid)
  th <- theme_bw(base_size = base_size, base_family = PUB_FONT) +
    theme(text = element_text(colour = INK, family = PUB_FONT),
      axis.text = element_text(size = PT_AXIS_TEXT, colour = INK),
      axis.title = element_text(size = PT_AXIS_TITLE, colour = INK),
      axis.ticks = element_line(linewidth = 0.25, colour = RULE),
      axis.ticks.length = unit(1.2, "mm"),
      axis.line = element_line(linewidth = 0.3, colour = RULE),
      panel.border = element_blank(), panel.background = element_blank(),
      panel.grid.minor = element_blank(), panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "white", colour = NA),
      plot.title = element_blank(), plot.subtitle = element_blank(),
      plot.margin = margin(2,2,2,2,"mm"),
      legend.text = element_text(size = PT_LEGEND, colour = INK),
      legend.title = element_text(size = PT_LEGEND, colour = INK),
      legend.key = element_blank(), legend.key.size = unit(3.2,"mm"),
      legend.background = element_blank(), legend.margin = margin(0,0,0,0),
      strip.background = element_blank(),
      strip.text = element_text(size = PT_AXIS_TEXT, colour = INK))
  # Grid lines disabled globally (house style: no grid lines on any plot).
  # The `grid` argument is retained for call-site compatibility but is a no-op.
  th
}

tag_panels <- function(p, tags = "A") {   # Cell Press: UPPERCASE bold A, B, C
  p + patchwork::plot_annotation(tag_levels = tags) &
    theme(plot.tag = element_text(size = PT_PANEL_TAG, face = "bold",
                                  family = PUB_FONT, colour = INK, hjust = 0, vjust = 1))
}

.join_svg <- function(path) {
  tool <- file.path(dirname(sys.frame(1)$ofile %||% "."), "svg_join_text.py")
  if (!file.exists(tool)) tool <- "svg_join_text.py"
  suppressWarnings(system2("python3", c(shQuote(tool), shQuote(path)), stdout=TRUE, stderr=TRUE))
  invisible(path)
}
`%||%` <- function(a,b) if (is.null(a)) b else a

#' Save at FINAL printed size: PNG 600 dpi + editable SVG (text joined).
save_fig <- function(p, name, width_mm = FIG_W_DOUBLE, height_mm = 100,
                     dir = "../results/figures/panels", dpi = 600) {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
  stopifnot(height_mm <= FIG_H_MAX, width_mm <= FIG_W_DOUBLE)
  w <- width_mm*MM_TO_IN; h <- height_mm*MM_TO_IN
  png_dev <- if (requireNamespace("ragg", quietly=TRUE)) ragg::agg_png else "png"
  ggsave(file.path(dir, paste0(name,".png")), p, width=w, height=h, dpi=dpi, bg="white",
         device=png_dev, units="in")
  ggsave(file.path(dir, paste0(name,".svg")), p, width=w, height=h, device=svglite::svglite,
         bg="white", units="in")
  tool <- file.path("svg_join_text.py")
  if (file.exists(tool)) suppressWarnings(system2("python3", c(shQuote(tool),
         shQuote(file.path(dir, paste0(name,".svg")))), stdout=TRUE, stderr=TRUE))
  message("  FIGURE: ", name, " (", width_mm, " x ", height_mm, " mm)")
  invisible(p)
}
