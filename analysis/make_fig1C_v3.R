# Figure 1C (v3) - real TME + activating/suppressive balance + two-line fate.
# Servier icons (BioIcons, CC BY 3.0). Treg recoloured plum; T-cell nucleus navy.
suppressPackageStartupMessages({ library(ggplot2); library(grid); library(png) })
source("ici_theme.R")
.pt <- 2.845276
ICO <- function(f) rasterGrob(readPNG(file.path("assets/icons", f)), interpolate=TRUE)
img <- function(f, cx, cy, s) annotation_custom(ICO(f), xmin=cx-s/2, xmax=cx+s/2, ymin=cy-s/2, ymax=cy+s/2)
S_HDR<-7.6; S_GENE<-6.8; S_CELL<-7; S_FLOW<-7.2

zc<-c(Sense=30,Integration=90,Execution=150); zw<-56; zy<-51; zh<-58
zones<-within(data.frame(layer=names(zc),x=unname(zc)),{xmin<-x-zw/2;xmax<-x+zw/2;ymin<-zy-zh/2;ymax<-zy+zh/2})
killX<-data.frame(x=c(160,170,156),y=c(62,55,52))
ct<-seq(0,1,length.out=10); chemo<-data.frame(x=34+46*ct, y=51-4.5*sin(pi*ct)+3.5*ct)  # curved chemokine trail tumor -> T cell

p <- ggplot() +
  geom_rect(data=zones, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, colour=layer), fill=pal_layer_fill[zones$layer], linewidth=1.0) +
  scale_colour_manual(values=pal_layer, guide="none") +
  geom_text(data=zones, aes(x=x, y=ymax-3.5, label=toupper(layer), colour=layer), family=PUB_FONT, size=12/.pt, fontface=2) +
  ## ===== SENSE: real TME =====
  img("tumour.svg.png",20,62,12) + img("tumour2.svg.png",31,54,11) + img("tumour3.svg.png",19,50,10) +
  img("macro.svg.png",41,45,12) + img("dc.svg.png",50,61,12) + img("fibro.svg.png",24,39,21) +
  annotate("text", x=20, y=69.5, label="tumor cells", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  annotate("text", x=24, y=33.5, label="fibroblast (CAF)", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  annotate("text", x=39, y=36.5, label="macrophage", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  annotate("text", x=52, y=53.5, label="DC", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  geom_point(data=chemo, aes(x,y), size=1.3, colour="grey45") +
  annotate("text", x=52, y=44.5, label="chemokines", family=PUB_FONT, size=6.6/.pt, fontface=3, colour=INK) +
  annotate("text", x=30, y=28, label="Antigen sensing, IFN", family=PUB_FONT, size=S_HDR/.pt, fontface=2, colour=INK) +
  annotate("text", x=30, y=24, label="CXCL9/10, STAT1, B2M, HLA-A", family=PUB_FONT, size=S_GENE/.pt, fontface=3, colour=INK) +
  ## extended curved antigen-presentation arrow DC -> centered T cell, label ON the arrow
  geom_curve(aes(x=57,y=60,xend=82,yend=56), curvature=-0.26, linewidth=0.85, colour=pal_layer["Sense"], arrow=arrow(length=unit(2.6,"mm"),type="closed")) +
  annotate("text", x=66, y=64.5, label="antigen presentation", family=PUB_FONT, size=6.6/.pt, fontface=3, colour=INK) +
  ## ===== INTEGRATION: activating (+) vs suppressive (-) balance =====
  img("tcell.svg.png",90,55,15) +
  annotate("text", x=90, y=46, label="T cell activated", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  annotate("text", x=88, y=71.5, label="Activation (+)", family=PUB_FONT, size=S_HDR/.pt, fontface=2, colour=pal_layer["Execution"]) +
  annotate("text", x=88, y=67.5, label="CCL5, CCR7, CD80, TNF", family=PUB_FONT, size=S_GENE/.pt, fontface=3, colour=INK) +
  # SOLID line: infiltration -> execution
  annotate("segment", x=98,y=57,xend=125,yend=57, linewidth=1.0, colour=pal_layer["Execution"], arrow=arrow(length=unit(2.4,"mm"),type="closed")) +
  annotate("text", x=110, y=61.5, label="infiltration", family=PUB_FONT, size=S_FLOW/.pt, colour=pal_layer["Execution"]) +
  # suppressive cells: Treg + CAF
  img("treg.svg.png",95,39,9) + img("fibro.svg.png",110,38,14) +
  annotate("text", x=95, y=33, label="Treg", family=PUB_FONT, size=5.8/.pt, colour=INK) +
  annotate("text", x=110, y=32, label="CAF", family=PUB_FONT, size=5.8/.pt, colour=INK) +
  # DOTTED line: suppressive signals act on the T cell
  annotate("segment", x=87,y=48,xend=93,yend=43, linewidth=0.8, colour=PLUM, linetype=3) +
  # Suppression inhibits infiltration (blunt bar on the solid line)
  annotate("segment", x=101,y=43,xend=106,yend=54.5, linewidth=0.8, colour=ACCENT) +
  annotate("segment", x=102.5,y=57,xend=109.5,yend=57, linewidth=1.5, colour=ACCENT) +
  annotate("text", x=90, y=27.5, label="Suppression (−)", family=PUB_FONT, size=S_HDR/.pt, fontface=2, colour=ACCENT) +
  annotate("text", x=90, y=23.5, label="FOXP3, CTLA4, TGFB1, VEGFA, FAP", family=PUB_FONT, size=S_GENE/.pt, fontface=3, colour=INK) +
  ## ===== EXECUTION: CD8 infiltrates the TME, kills tumor cells (red X) =====
  img("tcell.svg.png",129,57,12) +
  annotate("text", x=129, y=48, label="CD8 T cell", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  img("tumour.svg.png",160,62,11) + img("tumour2.svg.png",170,55,10) + img("tumour3.svg.png",156,52,9) +
  img("macro.svg.png",150,44,9) + img("dc.svg.png",137,44,9) + img("fibro.svg.png",168,42,15) +
  annotate("segment", x=135,y=57,xend=152,yend=53, linewidth=0.5, colour="#CC0000") +
  annotate("segment", x=135,y=58,xend=165,yend=61, linewidth=0.5, colour="#CC0000") +
  annotate("segment", x=135,y=56,xend=166,yend=55, linewidth=0.5, colour="#CC0000") +
  geom_point(data=killX, aes(x,y), shape=4, size=4.2, stroke=1.5, colour="#CC0000") +
  annotate("text", x=150, y=28, label="Tumor-cell killing", family=PUB_FONT, size=S_HDR/.pt, fontface=2, colour=INK) +
  annotate("text", x=150, y=24, label="GZMB, GZMA, PRF1, NKG7, CD8A", family=PUB_FONT, size=S_GENE/.pt, fontface=3, colour=INK) +
  annotate("text", x=1, y=81, label="C", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +
  coord_fixed(ratio=1, xlim=c(0,180), ylim=c(19,82), expand=FALSE, clip="off") +
  theme_void() + theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(1,1,1,1,"mm"))

save_fig(p, "ICI_Fig1C_v3", width_mm=174, height_mm=174*63/180)
cat("done\n")
