# =============================================================================
# Figure 3 composite (layout preview):
#   banner: a flow (transcriptome -> scoring -> prediction + pooled separation), full width
#   row 1:  b atlas (where)       |  d dose-response (how well)
#   row 2:  c violin (separates)  |  e alluvial (why)
# Panels keep their own a/b/c/d/e tags. Raster composite for layout review; the
# EDITABLE figure is assembled in Illustrator/Inkscape from the panel SVGs.
# =============================================================================
suppressPackageStartupMessages({ library(grid); library(png) })
D  <- "../results/figures/composite/"
rd <- function(f){ a<-readPNG(paste0(D,f,".png"))
  list(g=rasterGrob(a, interpolate=TRUE, width=unit(1,"npc"), height=unit(1,"npc")), asp=dim(a)[2]/dim(a)[1]) }
fl <- rd("ICI_Fig3_flow")                      # A banner
a  <- rd("ICI_Fig3_atlas"); b <- rd("ICI_Fig3b_separation")   # B atlas, C violin
cc <- rd("ICI_Fig3c_doseresponse"); dd <- rd("ICI_Fig_archetypes")  # D dose, E alluvial

W <- 234; gap <- 4                             # per-row cell widths (b gets extra width + font)
wa <- 128; wc <- W-gap-wa                       # row 1:  B atlas | D dose-response
wb <- 146; wd <- W-gap-wb                       # row 2:  C violin (wider) | E alluvial
dim <- function(p, cw) c(cw, cw/p$asp)          # display w,h (fit by width)
bh  <- W/fl$asp                                 # banner height (full width)
dA<-dim(a,wa); dC<-dim(cc,wc); dB<-dim(b,wb); dD<-dim(dd,wd)
row1 <- max(dA[2], dC[2]); row2 <- max(dB[2], dD[2]); H <- bh + gap + row1 + gap + row2

place <- function(p, leftx, topy, d){
  pushViewport(viewport(x=unit(leftx,"mm"), y=unit(topy,"mm"),
                        width=unit(d[1],"mm"), height=unit(d[2],"mm"), just=c("left","top")))
  grid.draw(p$g); popViewport() }

png_dev <- if (requireNamespace("ragg", quietly=TRUE)) ragg::agg_png else png
MM <- 1/25.4
png_dev(paste0(D,"ICI_Fig3_composite.png"), width=W*MM, height=H*MM, units="in", res=300, background="white")
grid.newpage()
topgrid <- H - bh - gap
place(fl, 0,       H,                  c(W, bh))     # A flow banner (full width)
place(a,  0,       topgrid,            dA)           # B atlas top-left
place(cc, W-dC[1], topgrid,            dC)           # D dose-response top-right
place(b,  0,       topgrid-row1-gap,   dB)           # C violin bottom-left (wider)
place(dd, W-dD[1], topgrid-row1-gap,   dD)           # E alluvial bottom-right
invisible(dev.off())
cat(sprintf("done — ICI_Fig3_composite.png  (%.0f x %.0f mm)  banner %.0f mm\n", W, H, bh))
