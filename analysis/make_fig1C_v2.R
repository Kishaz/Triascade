# Figure 1C (v2) - DYNAMIC biological cascade with inter-layer flow.
# TME -> tumor releases antigen -> DC binds -> wavy antigen presentation to
# T cell -> (Treg BLOCKS) -> activated T cell -> Execution -> returns to TME.
# Servier icons (BioIcons, CC BY 3.0). Treg recoloured plum to differ from T cell.
suppressPackageStartupMessages({ library(ggplot2); library(grid); library(png) })
source("ici_theme.R")
.pt <- 2.845276
ICO <- function(f) rasterGrob(readPNG(file.path("assets/icons", f)), interpolate=TRUE)
img <- function(f, cx, cy, s) annotation_custom(ICO(f), xmin=cx-s/2, xmax=cx+s/2, ymin=cy-s/2, ymax=cy+s/2)
wav <- function(x1,x2,y,amp=1.8,cyc=3,n=90){ x<-seq(x1,x2,length.out=n); data.frame(x=x,y=y+amp*sin(seq(0,cyc*2*pi,length.out=n))) }
# return arc (semi-ellipse) that runs BELOW the three zones
rt <- function(xl,xr,ytop,dip,n=120){ t<-seq(0,pi,length.out=n); data.frame(x=(xl+xr)/2+(xr-xl)/2*cos(t), y=ytop-dip*sin(t)) }

zc<-c(Sense=30,Integration=90,Execution=150); zw<-56; zy<-53; zh<-54
zones<-within(data.frame(layer=names(zc),x=unname(zc)),{xmin<-x-zw/2;xmax<-x+zw/2;ymin<-zy-zh/2;ymax<-zy+zh/2})
antigenP<-data.frame(x=c(33,36,39),y=c(60,61,60))
granP<-data.frame(x=c(147,150,153),y=c(57,55,56))
apopP<-data.frame(x=c(171,173,169),y=c(60,55,63))
retA <- rt(34,148,26,13)

p <- ggplot() +
  geom_rect(data=zones, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, colour=layer), fill=pal_layer_fill[zones$layer], linewidth=1.0) +
  scale_colour_manual(values=pal_layer, guide="none") +
  geom_text(data=zones, aes(x=x, y=ymax-3.5, label=toupper(layer), colour=layer), family=PUB_FONT, size=12/.pt, fontface=2) +
  ## ---- SENSE: TME with tumor cells; antigen -> DC ----
  annotate("rect", xmin=5,xmax=52,ymin=42,ymax=74, fill="#EADFC6", colour=SAND_DEEP, linewidth=0.4, linetype=3) +
  annotate("text", x=8.5, y=71, label="TME", family=PUB_FONT, size=7/.pt, fontface=2, colour=SAND_DEEP, hjust=0) +
  img("tumour.svg.png",18,58,13) + img("tumour.svg.png",27,49,11) +
  geom_point(data=antigenP, aes(x,y), shape=17, size=1.5, colour=INK) +
  annotate("segment", x=32,y=57,xend=39,yend=57, linewidth=0.45, colour=INK, linetype=2, arrow=arrow(length=unit(1.1,"mm"),type="closed")) +
  img("dc.svg.png",45,56,15) +
  annotate("text", x=21, y=41, label="tumor cells", family=PUB_FONT, size=6.3/.pt, colour=MUTED) +
  annotate("text", x=47, y=45.5, label="dendritic cell", family=PUB_FONT, size=6.3/.pt, colour=MUTED) +
  ## ---- wavy antigen presentation SENSE -> INTEGRATION ----
  geom_path(data=wav(53,82,57), aes(x,y), linewidth=0.8, colour=pal_layer["Sense"], arrow=arrow(length=unit(2.2,"mm"),type="closed")) +
  annotate("text", x=67, y=64, label="antigen presentation", family=PUB_FONT, size=6.6/.pt, fontface=3, colour=MUTED) +
  ## ---- INTEGRATION: T cell activated; Treg BLOCKS exit ----
  img("tcell.svg.png",90,57,15) +
  annotate("text", x=90, y=48, label="T cell activated", family=PUB_FONT, size=6.3/.pt, colour=MUTED) +
  img("treg.svg.png",108,38,10) +
  annotate("text", x=108, y=31.5, label="Treg (FoxP3+)", family=PUB_FONT, size=6.3/.pt, colour=PLUM) +
  ## ---- INTEGRATION -> EXECUTION with Treg inhibition ----
  annotate("segment", x=100,y=57,xend=127,yend=57, linewidth=0.9, colour=INK, arrow=arrow(length=unit(2.2,"mm"),type="closed")) +
  annotate("segment", x=108,y=43,xend=113,yend=54.5, linewidth=0.8, colour=ACCENT) +
  annotate("segment", x=109.5,y=57,xend=116.5,yend=57, linewidth=1.4, colour=ACCENT) +
  annotate("text", x=118.5, y=64, label="if Treg:\nno execution", family=PUB_FONT, size=6.3/.pt, fontface=3, colour=ACCENT, lineheight=1.05, hjust=0) +
  ## ---- EXECUTION: CD8 kills tumor ----
  img("tcell.svg.png",138,57,14) +
  annotate("text", x=138, y=48, label="CD8 T cell", family=PUB_FONT, size=6.3/.pt, colour=MUTED) +
  geom_point(data=granP, aes(x,y), size=1.3, colour=pal_layer["Execution"]) +
  img("tumour.svg.png",166,58,14) +
  geom_point(data=apopP, aes(x,y), size=1.0, colour=MUTED) +
  annotate("text", x=166, y=48, label="tumor killed", family=PUB_FONT, size=6.3/.pt, colour=MUTED) +
  annotate("text", x=150, y=38, label="Cytotoxic killing  (GZMB, PRF1, NKG7)", family=PUB_FONT, size=6.6/.pt, fontface=3, colour=MUTED) +
  ## ---- return arc EXECUTION -> TME (below the zones) ----
  geom_path(data=retA, aes(x,y), linewidth=0.9, colour=pal_layer["Execution"], arrow=arrow(length=unit(2.6,"mm"),type="closed")) +
  annotate("text", x=91, y=7.5, label="activated T cell returns to the TME and kills tumor cells",
           family=PUB_FONT, size=7/.pt, colour=INK) +
  annotate("text", x=1, y=81, label="C", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +
  coord_fixed(ratio=1, xlim=c(0,180), ylim=c(0,82), expand=FALSE, clip="off") +
  theme_void() + theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(1,1,1,1,"mm"))

save_fig(p, "ICI_Fig1C_v2", width_mm=174, height_mm=174*82/180)
cat("done\n")
