# Figure 1 Panel C (scoring, gene-grid) - one tumor's panel genes are binarized
# at the cohort median; layer scores sum the calls. House style (Cancer Cell).
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R")
.pt <- 2.845276
SQ<-1.25; PITCH<-3.1   # square half-size, pitch
# illustrative binary calls (filled = gene scored 1)
sq_row <- function(y, x0, patt, layer){
  data.frame(xc=x0+(seq_along(patt)-1)*PITCH, y=y,
             fill=ifelse(patt==1, unname(pal_layer[layer]),"white"),
             border=unname(pal_layer[layer])) }
pS <- c(1,1,0,1,1,1,0,1,1,0,1,1)
pA <- c(1,1,0,1,1,0,1,1,0,1)
pN <- c(1,0,1,0,1,0,0,1,0,1)
pE <- c(1,1,1,0,1,1,1,0,1,1,0,1)
x0<-30
S  <- sq_row(30, x0, pS, "Sense")
A  <- sq_row(19, x0, pA, "Integration")
N  <- sq_row(19, x0+length(pA)*PITCH+6, pN, "Integration")
E  <- sq_row(8,  x0, pE, "Execution")
sq <- rbind(S,A,N,E)
minus_x <- x0+length(pA)*PITCH+2.5
rowend  <- max(N$xc)+3
chip <- within(data.frame(layer=c("Sense","Integration","Execution"), x=14, y=c(30,19,8), w=22, h=8.5),
               {xmin<-x-w/2;xmax<-x+w/2;ymin<-y-h/2;ymax<-y+h/2})
# geometric bus -> composite
seg_h<-data.frame(x=rowend, y=c(30,19,8), xend=rowend+5, yend=c(30,19,8))
seg_v<-data.frame(x=rowend+5, y=8, xend=rowend+5, yend=30)
arr <-data.frame(x=rowend+5, y=19, xend=rowend+13, yend=19)
comp<-within(data.frame(x=140, y=19, w=54, h=22),{xmin<-x-w/2;xmax<-x+w/2;ymin<-y-h/2;ymax<-y+h/2})

p <- ggplot() +
  annotate("text", x=7, y=40.5, hjust=0, family=PUB_FONT, size=7.2/.pt, fontface=3, colour=INK,
           label="Each panel gene is scored 1 (filled) if its expression is ≥ the cohort median, else 0 (open)") +
  # chips
  geom_rect(data=chip, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, fill=layer), colour=NA) +
  geom_text(data=chip, aes(x=x,y=y,label=layer), family=PUB_FONT, size=7.4/.pt, fontface=2, colour="white") +
  scale_fill_manual(values=pal_layer, guide="none") +
  # gene squares
  geom_rect(data=sq, aes(xmin=xc-SQ,xmax=xc+SQ,ymin=y-SQ,ymax=y+SQ), fill=sq$fill, colour=sq$border, linewidth=0.5) +
  annotate("text", x=minus_x, y=19, label="−", family=PUB_FONT, size=13/.pt, colour=INK) +
  annotate("text", x=x0, y=24, hjust=0, label="activating (+)", family=PUB_FONT, size=6/.pt, colour=INK) +
  annotate("text", x=x0+length(pA)*PITCH+6, y=24, hjust=0, label="suppressive (−)", family=PUB_FONT, size=6/.pt, colour=INK) +
  # score formulas
  annotate("text", x=rowend-1, y=33.5, hjust=1, label="= Σ high", family=PUB_FONT, size=6.4/.pt, colour=INK) +
  annotate("text", x=rowend-1, y=11.5, hjust=1, label="= Σ high", family=PUB_FONT, size=6.4/.pt, colour=INK) +
  # bus to composite
  geom_segment(data=seg_h, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.5, colour=INK) +
  geom_segment(data=seg_v, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.5, colour=INK) +
  geom_segment(data=arr, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.5, colour=INK, arrow=arrow(length=unit(1.8,"mm"),type="closed")) +
  geom_rect(data=comp, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill="#EDE4D3", colour=MUTED, linewidth=0.7) +
  annotate("text", x=140, y=22.5, label="Composite score", family=PUB_FONT, size=7.2/.pt, fontface=2, colour=INK) +
  annotate("text", x=140, y=16.5, label="mean of the three\nz-scored layer scores", family=PUB_FONT, size=6.4/.pt, colour=INK, lineheight=1.1) +
  annotate("text", x=1, y=42, label="C", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +
  coord_fixed(ratio=1, xlim=c(0,180), ylim=c(0,43), expand=FALSE, clip="off") +
  theme_void() + theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(1,1,1,1,"mm"))

save_fig(p, "ICI_Fig1C_scoring", width_mm=174, height_mm=174*43/180)
cat("done\n")
