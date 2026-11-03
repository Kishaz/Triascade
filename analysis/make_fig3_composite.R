# =============================================================================
# Figure 3 composite (layout preview): 2x2 grid of the four panel PNGs.
#   row 1:  a atlas (where)      |  c dose-response (how well)
#   row 2:  b violin (separates) |  d alluvial (why)
# Panels keep their own A/B/C/D tags. Raster composite for layout review; the
# EDITABLE figure is assembled in Illustrator from the four panel SVGs.
# =============================================================================
suppressPackageStartupMessages({ library(grid); library(png) })
D  <- "../results/figures/composite/"
rd <- function(f){ a<-readPNG(paste0(D,f,".png"))
  list(g=rasterGrob(a, interpolate=TRUE, width=unit(1,"npc"), height=unit(1,"npc")), asp=dim(a)[2]/dim(a)[1]) }
a <- rd("ICI_Fig3_atlas"); b <- rd("ICI_Fig3b_separation")
cc<- rd("ICI_Fig3c_doseresponse"); dd<- rd("ICI_Fig_archetypes")

W   <- 180; gap <- 3
wL  <- 100; wR <- W-gap-wL                  # left col (a,b) wider; right col (c,d)
row1<- max(wL/a$asp, wR/cc$asp)
row2<- max(wL/b$asp, wR/dd$asp)
H   <- row1 + gap + row2

## fit image (preserve aspect) into a cell, top-centred
fit <- function(p, cw, ch){
  if (cw/ch > p$asp) { h<-ch; w<-h*p$asp } else { w<-cw; h<-w/p$asp }
  pushViewport(viewport(x=0.5, y=1, width=unit(w,"mm"), height=unit(h,"mm"), just=c("centre","top")))
  grid.draw(p$g); popViewport()
}

png_dev <- if (requireNamespace("ragg", quietly=TRUE)) ragg::agg_png else png
MM <- 1/25.4
png_dev(paste0(D,"ICI_Fig3_composite.png"), width=W*MM, height=H*MM, units="in", res=300, background="white")
grid.newpage()
pushViewport(viewport(layout=grid.layout(3, 3,
  heights=unit(c(row1,gap,row2),"mm"), widths=unit(c(wL,gap,wR),"mm"))))
vp <- function(r,c) viewport(layout.pos.row=r, layout.pos.col=c)
pushViewport(vp(1,1)); fit(a,  wL, row1); popViewport()
pushViewport(vp(1,3)); fit(cc, wR, row1); popViewport()
pushViewport(vp(3,1)); fit(b,  wL, row2); popViewport()
pushViewport(vp(3,3)); fit(dd, wR, row2); popViewport()
popViewport()
invisible(dev.off())
cat(sprintf("done — ICI_Fig3_composite.png  (%.0f x %.0f mm)\n", W, H))
