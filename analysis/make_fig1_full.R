# Figure 1 (full): A derivation ; B modules -> 3 layers ; C biological cascade
# with Servier cell icons (BioIcons, CC BY 3.0). One canvas, house style.
suppressPackageStartupMessages({ library(ggplot2); library(grid); library(png) })
source("ici_theme.R")
.pt <- 2.845276
rectify <- function(d){ within(d,{xmin<-x-w/2; xmax<-x+w/2; ymin<-y-h/2; ymax<-y+h/2}) }
ICO <- function(f) rasterGrob(readPNG(file.path("assets/icons", f)), interpolate=TRUE)
img <- function(f, cx, cy, s) annotation_custom(ICO(f), xmin=cx-s/2, xmax=cx+s/2, ymin=cy-s/2, ymax=cy+s/2)
imgwh <- function(f, cx, cy, w, h) annotation_custom(ICO(f), xmin=cx-w/2, xmax=cx+w/2, ymin=cy-h/2, ymax=cy+h/2)
VES_AR <- 2.87          # blood-vessel icon aspect (w/h)
TEAL   <- "#1F7A72"     # checkpoint / exhaustion = ICB target (therapeutic hook)
shiftY <- function(d,off){ for(cc in intersect(c("y","yend","ymin","ymax"),names(d))) d[[cc]]<-d[[cc]]+off; d }
S_HDR<-7.6; S_GENE<-6.8; S_CELL<-7; S_FLOW<-7.2
OFF2 <- 83   # A/B (derivation) lifted above the biological cascade (C, native coords)

# ---- legacy-98 counts (data, not typed); substrate count from CORE run ----
asg <- read.csv("ICI_98_module_assignment.csv", stringsAsFactors=FALSE)
mc  <- table(factor(asg$module, paste0("M",1:6)))
cnt <- function(m) as.integer(mc[m])
L_sense <- sum(asg$layer=="Sense")
I_pos <- sum(asg$layer=="Integration" & asg$direction=="positive")
I_neg <- sum(asg$layer=="Integration" & asg$direction=="negative"); I_tot <- I_pos+I_neg
E_pos <- sum(asg$layer=="Execution" & asg$direction=="positive")
E_neg <- sum(asg$layer=="Execution" & asg$direction=="negative"); E_tot <- E_pos+E_neg
n_genes <- length(unique(asg$gene))
.subf <- "../ICI-study-redesign/provenance-review-tight/derivation/substrate_core.csv"
n_sig <- if (file.exists(.subf)) nrow(read.csv(.subf, stringsAsFactors=FALSE)) else 210L

## ===================== C: biological cascade (native, bottom) =====================
zc<-c(Sense=30,Integration=90,Execution=150); zw<-56; zy<-51; zh<-58
zones<-within(data.frame(layer=names(zc),x=unname(zc)),{xmin<-x-zw/2;xmax<-x+zw/2;ymin<-zy-zh/2;ymax<-zy+zh/2})
zones$disp<-ifelse(zones$layer=="Sense","Sensory",zones$layer)
zq<-data.frame(x=unname(zc), layer=names(zc),
  q=c("Is the tumor seen?","Can lymphocytes reach it?","Can lymphocytes kill it?"))
killX<-data.frame(x=c(160,170,156),y=c(62,55,52))
ct<-seq(0,1,length.out=10); chemo<-data.frame(x=34+46*ct, y=51-4.5*sin(pi*ct)+3.5*ct)

## ===================== A: provenance =====================
prov<-rectify(data.frame(x=25,y=c(94,74,54,34,14),w=47,h=15,
  label=c("Systematic review of\nICB transcriptomic studies",sprintf("%d curated signatures\n(substrate)", n_sig),
          "Gene co-occurrence\ngraph","Program-blind Louvain\n+ g:Profiler GO",
          sprintf("6 emergent modules\n(%d genes)", n_genes)),
  out=c(FALSE,FALSE,FALSE,FALSE,TRUE)))
prov_arr<-data.frame(x=25,y=prov$y[-5]-7.5,xend=25,yend=prov$y[-1]+7.5)
## ===================== B: modules -> layers =====================
subn <- function(m, sign, desc) sprintf("'%d genes'~~~bold('%s')~~~'%s'", cnt(m), sign, desc)
mod<-rectify(data.frame(x=90,y=seq(96,12,length.out=6),w=72,h=13.5,
  title=c("M1   IFN / antigen-presentation","M2   chemokine / trafficking","M3   stromal / ECM",
          "M4   angiogenic / CAF-TGFβ","M5   cytotoxic effector","M6   checkpoint / exhaustion"),
  sub=c(subn("M1","+","antigen processing / IFN"), subn("M2","+","chemotaxis / adhesion"),
        subn("M3","−","ECM / collagen"), subn("M4","−","angiogenesis / TGFβ"),
        subn("M5","+","lymphocyte cytotoxicity"), subn("M6","−","checkpoint / exhaustion")),
  layer=c("Sense","Integration","Integration","Integration","Execution","Execution"),lty=c(1,1,2,2,1,2)))
lay_cfg<-data.frame(layer=c("Sense","Integration","Execution"),x=157,cy=c(96,62.4,20.4),h=c(18,30,30))
lay<-rectify(data.frame(x=lay_cfg$x,y=lay_cfg$cy,w=42,h=lay_cfg$h,layer=lay_cfg$layer))
lay_lines<-list(
  Sense       = c("bold('Sensory')", sprintf("'%d genes  ('*bold('+')*')'", L_sense), "'antigen sensing / IFN'"),
  Integration = c("bold('Integration')", sprintf("'%d genes'", I_tot),
                  sprintf("'%d activating ('*bold('+')*')'", I_pos),
                  sprintf("'%d suppressive ('*bold('−')*')'", I_neg), "'access balance'"),
  Execution   = c("bold('Execution')", sprintf("'%d genes'", E_tot),
                  sprintf("'%d effector ('*bold('+')*')'", E_pos),
                  sprintf("'%d exhaustion ('*bold('−')*')'", E_neg), "'killing balance'"))
PITCH<-4.8
layer_txt<-do.call(rbind,lapply(seq_len(nrow(lay_cfg)),function(i){Lr<-lay_cfg$layer[i];ln<-lay_lines[[Lr]];n<-length(ln)
  first<-lay_cfg$cy[i]+(n-1)*PITCH/2
  data.frame(x=lay_cfg$x[i],y=first-(seq_len(n)-1)*PITCH,label=ln,
             sz=ifelse(seq_len(n)==1,13,9.5),
             col=ifelse(seq_len(n)==1,pal_layer[Lr],INK))}))
mod_xr<-126;bus_x<-131;box_xl<-136
int_y<-mod$y[mod$layer=="Integration"]; exec_y<-mod$y[mod$layer=="Execution"]
sense_arr<-data.frame(x=mod_xr,y=mod$y[mod$layer=="Sense"],xend=box_xl,yend=mod$y[mod$layer=="Sense"])
seg_hi<-data.frame(x=mod_xr,y=int_y,xend=bus_x,yend=int_y)
seg_vi<-data.frame(x=bus_x,y=min(int_y),xend=bus_x,yend=max(int_y))
arr_i<-data.frame(x=bus_x,y=lay_cfg$cy[2],xend=box_xl,yend=lay_cfg$cy[2])
seg_he<-data.frame(x=mod_xr,y=exec_y,xend=bus_x,yend=exec_y)
seg_ve<-data.frame(x=bus_x,y=min(exec_y),xend=bus_x,yend=max(exec_y))
arr_e<-data.frame(x=bus_x,y=lay_cfg$cy[3],xend=box_xl,yend=lay_cfg$cy[3])
for(nm in c("prov","prov_arr","mod","lay","layer_txt","sense_arr","seg_hi","seg_vi","arr_i","seg_he","seg_ve","arr_e")) assign(nm, shiftY(get(nm),OFF2))

# scattered signalling dots (Integration + Execution), edited-style density (C native coords)
dots<-data.frame(x=c(80,84,100,105,92,99,   150,158,166,173,153),
                 y=c(50,62,61,51,46,43,      60,46,58,50,52))

p <- ggplot() +
  ## ---- C: cascade ----
  geom_rect(data=zones, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, colour=layer), fill=pal_layer_fill[zones$layer], linewidth=1.0) +
  scale_colour_manual(values=pal_layer, guide="none") +
  geom_text(data=zones, aes(x=x, y=ymax-3.5, label=toupper(disp), colour=layer), family=PUB_FONT, size=12/.pt, fontface=2) +
  geom_text(data=zq, aes(x=x, y=72.5, label=q, colour=layer), family=PUB_FONT, size=9.3/.pt, fontface=3) +
  img("tumour.svg.png",20,62,12) + img("tumour2.svg.png",31,54,11) + img("tumour3.svg.png",19,50,10) +
  img("tumour.svg.png",28,46,8.5) + img("dc.svg.png",33,62,8) +
  img("macro.svg.png",41,45,12) + img("dc.svg.png",50,61,12) +
  annotate("text", x=20, y=69.5, label="tumor cells", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  annotate("text", x=39, y=36.5, label="macrophage", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  annotate("text", x=52, y=53.5, label="DC", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  geom_point(data=chemo, aes(x,y), size=1.3, colour="grey45") +
  geom_point(data=dots, aes(x,y), size=1.2, colour="grey55") +
  annotate("text", x=52, y=44.5, label="chemokines", family=PUB_FONT, size=6.6/.pt, fontface=3, colour=INK) +
  annotate("text", x=30, y=28, label="Antigen sensing, IFN", family=PUB_FONT, size=S_HDR/.pt, fontface=2, colour=INK) +
  annotate("text", x=30, y=25, label="CXCL9/10, STAT1, B2M, HLA-A", family=PUB_FONT, size=S_GENE/.pt, fontface=3, colour=INK) +
  imgwh("vessel.png", 71, 59.5, 30, 30/VES_AR) +
  img("tcell.svg.png", 73, 59, 7.5) +
  annotate("text", x=72, y=51.5, label="recruitment", family=PUB_FONT, size=S_FLOW/.pt, colour=INK) +
  img("tcell.svg.png",90,55,15) +
  annotate("text", x=90, y=46, label="activated T / NK", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  annotate("text", x=88, y=69.0, label="Activation (+)", family=PUB_FONT, size=S_HDR/.pt, fontface=2, colour=pal_layer["Execution"]) +
  annotate("text", x=88, y=65.2, label="CCL5, CCR7, CD80, TNF", family=PUB_FONT, size=S_GENE/.pt, fontface=3, colour=INK) +
  imgwh("vessel.png", 112, 57, 30, 30/VES_AR) +
  img("tcell.svg.png", 114, 56.5, 7.5) +
  annotate("text", x=112, y=64, label="infiltration", family=PUB_FONT, size=S_FLOW/.pt, colour=INK) +
  img("stroma.svg.png",95,39,10) + img("fibro.svg.png",110,38,14) +
  annotate("text", x=95, y=33, label="stroma", family=PUB_FONT, size=5.8/.pt, colour=INK) +
  annotate("text", x=110, y=32, label="CAF", family=PUB_FONT, size=5.8/.pt, colour=INK) +
  annotate("segment", x=87,y=48,xend=93,yend=43, linewidth=0.8, colour=PLUM, linetype=3) +
  annotate("segment", x=101,y=43,xend=106,yend=54.5, linewidth=0.8, colour=ACCENT) +
  annotate("segment", x=102.5,y=57,xend=109.5,yend=57, linewidth=1.5, colour=ACCENT) +
  annotate("text", x=90, y=27.5, label="Suppression (−)", family=PUB_FONT, size=S_HDR/.pt, fontface=2, colour=ACCENT) +
  annotate("text", x=90, y=24.5, label="TGFB1, VEGFA, FAP, COL1A1, MMP9", family=PUB_FONT, size=S_GENE/.pt, fontface=3, colour=INK) +
  img("tcell.svg.png",128,57,12) +
  img("nk.svg.png",142,57,10) +
  annotate("text", x=128, y=48.5, label="CD8 T", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  annotate("text", x=142, y=48.5, label="NK", family=PUB_FONT, size=S_CELL/.pt, colour=INK) +
  annotate("text", x=135, y=68, label="ICB target", family=PUB_FONT, size=6.4/.pt, fontface=2, colour=TEAL) +
  annotate("segment", x=135, y=66.2, xend=135, yend=62, linewidth=0.6, colour=TEAL) +
  img("tumour.svg.png",160,62,11) + img("tumour2.svg.png",170,55,10) + img("tumour3.svg.png",156,52,9) +
  img("macro.svg.png",174,61,8) + img("macro.svg.png",172,47,8) +
  img("treg.svg.png",135,43,9) + annotate("text", x=135, y=38, label="Treg", family=PUB_FONT, size=5.8/.pt, colour=INK) +
  annotate("segment", x=134,y=57,xend=152,yend=52, linewidth=0.5, colour="#CC0000") +
  annotate("segment", x=134,y=58,xend=165,yend=61, linewidth=0.5, colour="#CC0000") +
  annotate("segment", x=147,y=56,xend=167,yend=55, linewidth=0.6, colour="#CC0000", linetype=3) +
  annotate("segment", x=147,y=55,xend=159,yend=61, linewidth=0.6, colour="#CC0000", linetype=3) +
  geom_point(data=killX, aes(x,y), shape=4, size=4.2, stroke=1.5, colour="#CC0000") +
  annotate("text", x=150, y=34.0, label="Effector (+)", family=PUB_FONT, size=S_HDR/.pt, fontface=2, colour=pal_layer["Execution"]) +
  annotate("text", x=150, y=30.6, label="GZMB, PRF1, IFNG, NKG7, GNLY", family=PUB_FONT, size=S_GENE/.pt, fontface=3, colour=INK) +
  annotate("text", x=150, y=27.2, label="Checkpoint / exhaustion", family=PUB_FONT, size=S_HDR/.pt, fontface=2, colour=TEAL) +
  annotate("text", x=150, y=23.8, label="PD-1, CTLA4, LAG3, TIGIT, TIM-3", family=PUB_FONT, size=S_GENE/.pt, fontface=3, colour=INK) +
  ## ---- A: provenance ----
  geom_rect(data=prov, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill=CREAM, colour=ifelse(prov$out,PLUM,MUTED), linewidth=ifelse(prov$out,1.1,0.6)) +
  geom_text(data=prov, aes(x,y,label=label), family=PUB_FONT, size=ifelse(prov$out,10.5,9.5)/.pt, fontface=ifelse(prov$out,2,1), colour=INK, lineheight=1.25) +
  geom_segment(data=prov_arr, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK, arrow=arrow(length=unit(1.8,"mm"),type="closed")) +
  ## ---- B: modules -> layers ----
  geom_segment(data=sense_arr, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK, arrow=arrow(length=unit(1.9,"mm"),type="closed")) +
  geom_segment(data=seg_hi, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK) +
  geom_segment(data=seg_vi, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK) +
  geom_segment(data=arr_i, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK, arrow=arrow(length=unit(1.9,"mm"),type="closed")) +
  geom_segment(data=seg_he, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK) +
  geom_segment(data=seg_ve, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK) +
  geom_segment(data=arr_e, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK, arrow=arrow(length=unit(1.9,"mm"),type="closed")) +
  geom_rect(data=mod, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, colour=layer, linetype=I(lty)), fill=pal_layer_fill[mod$layer], linewidth=0.8) +
  geom_text(data=mod, aes(x=x, y=y+2.8, label=title), family=PUB_FONT, size=11/.pt, fontface=2, colour=INK) +
  geom_text(data=mod, aes(x=x, y=y-2.9, label=sub), family=PUB_FONT, size=9/.pt, colour=INK, parse=TRUE) +
  geom_rect(data=lay, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, colour=layer), fill=pal_layer_fill[lay$layer], linewidth=1.1) +
  geom_text(data=layer_txt, aes(x=x,y=y,label=label), family=PUB_FONT, size=layer_txt$sz/.pt, colour=layer_txt$col, parse=TRUE) +
  ## tags
  annotate("text", x=1,  y=OFF2+106, label="A", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +
  annotate("text", x=52, y=OFF2+106, label="B", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +
  annotate("text", x=1,  y=84,        label="C", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +
  coord_fixed(ratio=1, xlim=c(0,180), ylim=c(19,OFF2+108), expand=FALSE, clip="off") +
  theme_void() + theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(1,1,1,1,"mm"))

H <- (OFF2+108)-19
save_fig(p, "ICI_Fig1_full", width_mm=174, height_mm=174*H/180, dir="../results/figures/composite")
cat(sprintf("done (canvas height %d units)\n", H))
