# =============================================================================
# Figure 3 composite (layout preview): assembles the four panel PNGs into one
# figure.  a atlas (where) | b violin (separates) | c dose-response (how well) |
# d alluvial (why).  Panels keep their own A/B/C/D tags.  Raster composite for
# layout review; the EDITABLE figure is assembled in Illustrator from the four
# panel SVGs (this just shows the arrangement at true relative sizes).
# =============================================================================
suppressPackageStartupMessages({ library(grid); library(png) })
D <- "../results/figures/composite/"
rd  <- function(f){ a<-readPNG(paste0(D,f,".png")); list(g=rasterGrob(a, interpolate=TRUE), asp=dim(a)[2]/dim(a)[1]) }
a <- rd("ICI_Fig3_atlas"); b <- rd("ICI_Fig3b_separation")
cc<- rd("ICI_Fig3c_doseresponse"); dd<- rd("ICI_Fig_archetypes")

W   <- 180                                  # figure width (mm)
gap <- 4
hA  <- W/a$asp                               # atlas, full width
hB  <- W/b$asp                               # violin, full width
wR  <- (W-gap)/2                             # bottom row: c + d side by side
hRow<- max(wR/cc$asp, wR/dd$asp)
H   <- hA + hB + hRow + 2*gap

png_dev <- if (requireNamespace("ragg", quietly=TRUE)) ragg::agg_png else png
MM <- 1/25.4
png_dev(paste0(D,"ICI_Fig3_composite.png"), width=W*MM, height=H*MM, units="in", res=300, background="white")
grid.newpage()
pushViewport(viewport(layout=grid.layout(5, 1, heights=unit(c(hA,gap,hB,gap,hRow),"mm"))))
# a (atlas)
pushViewport(viewport(layout.pos.row=1)); grid.draw(a$g); popViewport()
# b (violin)
pushViewport(viewport(layout.pos.row=3)); grid.draw(b$g); popViewport()
# c + d
pushViewport(viewport(layout.pos.row=5, layout=grid.layout(1, 3, widths=unit(c(wR,gap,wR),"mm"))))
pushViewport(viewport(layout.pos.col=1)); grid.draw(cc$g); popViewport()
pushViewport(viewport(layout.pos.col=3)); grid.draw(dd$g); popViewport()
popViewport(2)
invisible(dev.off())
cat(sprintf("done — ICI_Fig3_composite.png  (%.0f x %.0f mm)\n", W, H))
