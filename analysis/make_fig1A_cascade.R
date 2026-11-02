# Figure 1A - what ICI represents biologically: the anti-tumor immune cascade
# drawn as three mechanistic layers (Sense -> Integration -> Execution).
# House style (Cancer Cell: Arial, warm palette, geometric connectors).
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R")
.pt <- 2.845276

circ <- function(cx,cy,r,n=72){ t<-seq(0,2*pi,length.out=n); data.frame(x=cx+r*cos(t), y=cy+r*sin(t)) }
star <- function(cx,cy,ro,ri,p=11){ a<-seq(0,2*pi,length.out=2*p+1); rr<-rep(c(ro,ri),length.out=length(a))
  data.frame(x=cx+rr*cos(a), y=cy+rr*sin(a)) }

## zones (one per layer)
zc <- c(Sense=30, Integration=90, Execution=150); zw<-56; zy<-37; zh<-64
zones <- within(data.frame(layer=names(zc), x=unname(zc)),
                {xmin<-x-zw/2; xmax<-x+zw/2; ymin<-zy-zh/2; ymax<-zy+zh/2})
# flow arrows between zones (geometric, straight)
flow <- data.frame(x=c(59,119), y=zy, xend=c(61,121), yend=zy)

## ---------- glyphs ----------
# SENSE: tumor cell (spiky) presenting antigen; recruits a T cell
tumourS <- star(22,45,8,5.2); antigen <- data.frame(x=c(19,22,25), y=c(53,54,53))  # peptide marks
tcellS  <- circ(40,30,5)
senseArr<- data.frame(x=29,y=38,xend=36,yend=32)
# INTEGRATION: central T cell with opposing activating (+) / suppressive (-) inputs
tcellI  <- circ(90,37,6.5)
# EXECUTION: CD8 T cell -> cytotoxic granules -> dying tumor cell
cd8     <- circ(135,40,5.5)
dying   <- star(162,38,8,3.8,13)
gran    <- data.frame(x=c(144,147.5,151), y=c(41,39,40.5))

p <- ggplot() +
  # zone panels + headers
  geom_rect(data=zones, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, colour=layer),
            fill=pal_layer_fill[zones$layer], linewidth=1.0) +
  geom_text(data=zones, aes(x=x, y=ymax-4.5, label=toupper(layer), colour=layer),
            family=PUB_FONT, size=12/.pt, fontface=2) +
  scale_colour_manual(values=pal_layer, guide="none") +
  geom_segment(data=flow, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.9, colour=INK,
               arrow=arrow(length=unit(2.6,"mm"), type="closed")) +
  # ---- SENSE ----
  geom_polygon(data=tumourS, aes(x,y), fill=pal_layer_fill["Execution"], colour=pal_layer["Sense"], linewidth=0.5) +
  geom_point(data=antigen, aes(x,y), shape=17, size=1.6, colour=INK) +
  geom_polygon(data=tcellS, aes(x,y), fill="white", colour=pal_layer["Sense"], linewidth=0.6) +
  annotate("text", x=40, y=30, label="T", family=PUB_FONT, size=8/.pt, fontface=2, colour=INK) +
  geom_segment(data=senseArr, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.5, colour=INK, arrow=arrow(length=unit(1.5,"mm"),type="closed")) +
  annotate("text", x=22, y=36, label="tumor cell", family=PUB_FONT, size=7/.pt, colour=MUTED) +
  annotate("text", x=30, y=20, label="Antigen presentation + IFN;\nT cells recruited", family=PUB_FONT, size=8/.pt, colour=INK, lineheight=1.15) +
  annotate("text", x=30, y=12, label="MHC-I, STAT1, CXCL9/10", family=PUB_FONT, size=7/.pt, fontface=3, colour=MUTED) +
  # ---- INTEGRATION ----
  geom_polygon(data=tcellI, aes(x,y), fill="white", colour=pal_layer["Integration"], linewidth=0.7) +
  annotate("text", x=90, y=37, label="T", family=PUB_FONT, size=9/.pt, fontface=2, colour=INK) +
  annotate("segment", x=90, y=57, xend=90, yend=45, linewidth=0.8, colour=pal_layer["Integration"], arrow=arrow(length=unit(2,"mm"),type="closed")) +
  annotate("segment", x=90, y=17, xend=90, yend=29, linewidth=0.8, colour=MUTED, arrow=arrow(length=unit(2,"mm"),type="closed")) +
  annotate("text", x=90, y=60, label="activating (+)", family=PUB_FONT, size=7.5/.pt, colour=pal_layer["Integration"]) +
  annotate("text", x=90, y=14, label="suppressive (−)", family=PUB_FONT, size=7.5/.pt, colour=MUTED) +
  annotate("text", x=90, y=52, label="costim, Th1", family=PUB_FONT, size=6.8/.pt, fontface=3, colour=MUTED) +
  annotate("text", x=90, y=21.5, label="Treg, checkpoint, TGF-β, CAF", family=PUB_FONT, size=6.8/.pt, fontface=3, colour=MUTED) +
  # ---- EXECUTION ----
  geom_polygon(data=cd8, aes(x,y), fill="white", colour=pal_layer["Execution"], linewidth=0.7) +
  annotate("text", x=135, y=40, label="CD8", family=PUB_FONT, size=7/.pt, fontface=2, colour=INK) +
  geom_point(data=gran, aes(x,y), size=1.5, colour=pal_layer["Execution"]) +
  geom_polygon(data=dying, aes(x,y), fill="white", colour=MUTED, linewidth=0.5, linetype=2) +
  annotate("text", x=162, y=27, label="tumor\ncell death", family=PUB_FONT, size=7/.pt, colour=MUTED, lineheight=1.1) +
  annotate("text", x=150, y=14, label="Cytotoxic killing", family=PUB_FONT, size=8/.pt, colour=INK) +
  annotate("text", x=150, y=8.5, label="GZMB, PRF1, NKG7", family=PUB_FONT, size=7/.pt, fontface=3, colour=MUTED) +
  annotate("text", x=1, y=71, label="A", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +
  coord_fixed(ratio=1, xlim=c(0,180), ylim=c(0,73), expand=FALSE, clip="off") +
  theme_void() + theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(1,1,1,1,"mm"))

save_fig(p, "ICI_Fig1A_cascade", width_mm=174, height_mm=174*73/180)
cat("done\n")
