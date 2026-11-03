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

W <- 234; gap <- 4                            # per-row cell widths (b gets extra width + font)
wa <- 128; wc <- W-gap-wa                      # top row:  a | c
wb <- 146; wd <- W-gap-wb                      # bottom row: b (wider) | d
dim <- function(p, cw) c(cw, cw/p$asp)         # display w,h (fit by width)
dA<-dim(a,wa); dC<-dim(cc,wc); dB<-dim(b,wb); dD<-dim(dd,wd)
row1 <- max(dA[2], dC[2]); row2 <- max(dB[2], dD[2]); H <- row1 + gap + row2

place <- function(p, leftx, topy, d){
  pushViewport(viewport(x=unit(leftx,"mm"), y=unit(topy,"mm"),
                        width=unit(d[1],"mm"), height=unit(d[2],"mm"), just=c("left","top")))
  grid.draw(p$g); popViewport() }

png_dev <- if (requireNamespace("ragg", quietly=TRUE)) ragg::agg_png else png
MM <- 1/25.4
png_dev(paste0(D,"ICI_Fig3_composite.png"), width=W*MM, height=H*MM, units="in", res=300, background="white")
grid.newpage()
place(a,  0,        H,            dA)          # a top-left
place(cc, W-dC[1],  H,            dC)          # c top-right
place(b,  0,        H-row1-gap,   dB)          # b bottom-left (wider)
place(dd, W-dD[1],  H-row1-gap,   dD)          # d bottom-right
invisible(dev.off())
cat(sprintf("done — ICI_Fig3_composite.png  (%.0f x %.0f mm)  b width %d mm\n", W, H, wb))
