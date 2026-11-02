# Figure 1C (prototype) - biological cascade with real Servier cell icons
# (BioIcons, CC BY 3.0). Icons rasterised via qlmanage to assets/icons/*.png.
suppressPackageStartupMessages({ library(ggplot2); library(grid); library(png) })
source("ici_theme.R")
.pt <- 2.845276
ICO <- function(f) rasterGrob(readPNG(file.path("assets/icons", f)), interpolate=TRUE)
img <- function(f, cx, cy, s, ar=1) annotation_custom(ICO(f), xmin=cx-s/2, xmax=cx+s/2, ymin=cy-s/2*ar, ymax=cy+s/2*ar)
tumor<-ICO("tumour.svg.png"); tcell<-ICO("tcell.svg.png"); dc<-ICO("dc.svg.png"); treg<-ICO("treg.svg.png")

zc<-c(Sense=30,Integration=90,Execution=150); zw<-56; zy<-37; zh<-62
zones<-within(data.frame(layer=names(zc),x=unname(zc)),{xmin<-x-zw/2;xmax<-x+zw/2;ymin<-zy-zh/2;ymax<-zy+zh/2})
flow<-data.frame(x=c(59,119),y=zy,xend=c(61,121),yend=zy)
apop<-data.frame(x=c(169,172,167),y=c(48,43,52))       # apoptotic bodies
gran<-data.frame(x=c(145,148.5,152),y=c(43,41,42))     # cytotoxic granules

p <- ggplot() +
  geom_rect(data=zones, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, colour=layer), fill=pal_layer_fill[zones$layer], linewidth=1.0) +
  scale_colour_manual(values=pal_layer, guide="none") +
  geom_text(data=zones, aes(x=x, y=ymax-4, label=toupper(layer), colour=layer), family=PUB_FONT, size=12/.pt, fontface=2) +
  geom_segment(data=flow, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.9, colour=INK, arrow=arrow(length=unit(2.6,"mm"),type="closed")) +
  ## SENSE: tumor cell presents antigen -> dendritic cell
  img("tumour.svg.png", 18, 46, 17) +
  annotate("point", x=c(27,29.5,32), y=c(52,53,52), shape=17, size=1.6, colour=INK) +
  annotate("segment", x=27, y=47, xend=33, yend=45, linewidth=0.5, colour=INK, arrow=arrow(length=unit(1.4,"mm"),type="closed")) +
  img("dc.svg.png", 42, 44, 18) +
  annotate("text", x=18, y=34, label="tumor cell", family=PUB_FONT, size=6.5/.pt, colour=MUTED) +
  annotate("text", x=44, y=31, label="dendritic cell", family=PUB_FONT, size=6.5/.pt, colour=MUTED) +
  annotate("text", x=30, y=17, label="Antigen presentation + IFN;\nT cells recruited", family=PUB_FONT, size=8/.pt, colour=INK, lineheight=1.15) +
  annotate("text", x=30, y=9.5, label="MHC-I, STAT1, CXCL9/10", family=PUB_FONT, size=6.8/.pt, fontface=3, colour=MUTED) +
  ## INTEGRATION: T cell, activating (+) vs suppressive (-) with a Treg
  img("tcell.svg.png", 90, 41, 16) +
  annotate("segment", x=90, y=58, xend=90, yend=50, linewidth=0.85, colour=pal_layer["Integration"], arrow=arrow(length=unit(1.9,"mm"),type="closed")) +
  annotate("text", x=90, y=61, label="activating (+)", family=PUB_FONT, size=7.2/.pt, colour=pal_layer["Integration"]) +
  annotate("text", x=90, y=54.5, label="costim, Th1", family=PUB_FONT, size=6.5/.pt, fontface=3, colour=MUTED) +
  img("treg.svg.png", 90, 18, 11) +
  annotate("segment", x=90, y=24, xend=90, yend=32, linewidth=0.85, colour=MUTED, arrow=arrow(length=unit(1.9,"mm"),type="closed")) +
  annotate("text", x=90, y=10.5, label="suppressive (−): Treg, checkpoint, TGF-β, CAF", family=PUB_FONT, size=6.8/.pt, fontface=3, colour=MUTED) +
  ## EXECUTION: CD8 T cell -> granules -> dying tumor cell
  img("tcell.svg.png", 133, 42, 15) +
  annotate("text", x=133, y=32, label="CD8 T cell", family=PUB_FONT, size=6.5/.pt, colour=MUTED) +
  geom_point(data=gran, aes(x,y), size=1.4, colour=pal_layer["Execution"]) +
  img("tumour.svg.png", 164, 44, 15) +
  geom_point(data=apop, aes(x,y), size=1.1, colour=MUTED) +
  annotate("text", x=150, y=17, label="Cytotoxic killing", family=PUB_FONT, size=8/.pt, colour=INK) +
  annotate("text", x=150, y=9.5, label="GZMB, PRF1, NKG7", family=PUB_FONT, size=6.8/.pt, fontface=3, colour=MUTED) +
  annotate("text", x=1, y=71, label="C", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +
  coord_fixed(ratio=1, xlim=c(0,180), ylim=c(0,73), expand=FALSE, clip="off") +
  theme_void() + theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(1,1,1,1,"mm"))

save_fig(p, "ICI_Fig1C_icons", width_mm=174, height_mm=174*73/180)
cat("done\n")
