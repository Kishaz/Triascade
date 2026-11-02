# Figure 1 (derivation): systematic review -> 6 co-occurrence modules (M1-M6,
# cascade order) -> 3 layers (S/I/E), populated by the OFFICIAL 98-gene panel.
# Module / layer counts are read from ICI_98_module_assignment.csv (gen_98_modules.R);
# nothing is typed.  House style ici_theme.R (Arial, UPPERCASE panel letters, warm
# palette, geometric square-corner connectors).
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R")
.pt <- 2.845276
rectify <- function(d){ within(d,{xmin<-x-w/2; xmax<-x+w/2; ymin<-y-h/2; ymax<-y+h/2}) }
S_MODT<-11; S_SUB<-9; S_PROV<-9.5; S_OUT<-10.5; S_LAYT<-13; S_LAYS<-9.5
PITCH <- 4.8   # uniform line pitch for layer-box text

# ---- counts from the 98-gene module assignment (data, not typed) ------------
asg <- read.csv("ICI_98_module_assignment.csv", stringsAsFactors=FALSE)
mc  <- table(factor(asg$module, paste0("M",1:6)))
cnt <- function(m) as.integer(mc[m])
L_sense  <- sum(asg$layer=="Sense")
I_pos <- sum(asg$layer=="Integration" & asg$direction=="positive")
I_neg <- sum(asg$layer=="Integration" & asg$direction=="negative"); I_tot <- I_pos+I_neg
E_pos <- sum(asg$layer=="Execution"  & asg$direction=="positive")
E_neg <- sum(asg$layer=="Execution"  & asg$direction=="negative"); E_tot <- E_pos+E_neg
n_genes <- length(unique(asg$gene))
# substrate signature count: CORE run that produced the 6 emergent modules (DERIVATION_REPORT.md)
.subf <- "../ICI-study-redesign/provenance-review-tight/derivation/substrate_core.csv"
n_sig <- if (file.exists(.subf)) nrow(read.csv(.subf, stringsAsFactors=FALSE)) else 210L

## ---- Panel A: provenance, ending in the connected output box ----
prov <- rectify(data.frame(x=25, y=c(94,74,54,34,14), w=47, h=15,
  label=c("Systematic review of\nICB transcriptomic studies",sprintf("%d curated signatures\n(substrate)", n_sig),
          "Gene co-occurrence\ngraph","Program-blind Louvain\n+ g:Profiler GO",
          sprintf("6 emergent modules\n(%d genes)", n_genes)),
  out=c(FALSE,FALSE,FALSE,FALSE,TRUE)))
prov_arr <- data.frame(x=25, y=prov$y[-5]-7.5, xend=25, yend=prov$y[-1]+7.5)

## ---- Panel B: six modules in cascade order (M1..M6), 98-gene counts ----
subn <- function(m, sign, desc) sprintf("'%d genes'~~~bold('%s')~~~'%s'", cnt(m), sign, desc)
mod <- rectify(data.frame(x=90, y=seq(96,12,length.out=6), w=72, h=13.5,
  title=c("M1   IFN / antigen-presentation","M2   chemokine / trafficking","M3   stromal / ECM",
          "M4   angiogenic / CAF-TGFβ","M5   cytotoxic effector","M6   checkpoint / exhaustion"),
  sub=c(subn("M1","+","antigen processing / IFN"), subn("M2","+","chemotaxis / adhesion"),
        subn("M3","−","ECM / collagen"),       subn("M4","−","angiogenesis / TGFβ"),
        subn("M5","+","lymphocyte cytotoxicity"),   subn("M6","−","checkpoint / exhaustion")),
  layer=c("Sense","Integration","Integration","Integration","Execution","Execution"),
  lty=c(1,1,2,2,1,2)))

## ---- three layers: box + text on ONE uniform pitch ----
lay_cfg <- data.frame(layer=c("Sense","Integration","Execution"), x=157, cy=c(96,62.4,20.4), h=c(18,30,30))
lay <- rectify(data.frame(x=lay_cfg$x, y=lay_cfg$cy, w=42, h=lay_cfg$h, layer=lay_cfg$layer))
# labels are plotmath so only the +/- direction sign is bold; titles baked bold()
lay_lines <- list(
  Sense       = c("bold('Sensory')", sprintf("'%d genes  ('*bold('+')*')'", L_sense), "'antigen sensing / IFN'"),
  Integration = c("bold('Integration')", sprintf("'%d genes'", I_tot),
                  sprintf("'%d activating ('*bold('+')*')'", I_pos),
                  sprintf("'%d suppressive ('*bold('−')*')'", I_neg), "'access balance'"),
  Execution   = c("bold('Execution')", sprintf("'%d genes'", E_tot),
                  sprintf("'%d effector ('*bold('+')*')'", E_pos),
                  sprintf("'%d exhaustion ('*bold('−')*')'", E_neg), "'killing balance'"))
layer_txt <- do.call(rbind, lapply(seq_len(nrow(lay_cfg)), function(i){
  Lr<-lay_cfg$layer[i]; ln<-lay_lines[[Lr]]; n<-length(ln)
  first <- lay_cfg$cy[i] + (n-1)*PITCH/2
  data.frame(x=lay_cfg$x[i], y=first-(seq_len(n)-1)*PITCH, label=ln,
             sz=ifelse(seq_len(n)==1,S_LAYT,S_LAYS),
             col=ifelse(seq_len(n)==1, pal_layer[Lr], INK)) }))

## ---- GEOMETRIC connectors: Sense straight; Integration + Execution each a bus ----
mod_xr<-126; bus_x<-131; box_xl<-136
int_y  <- mod$y[mod$layer=="Integration"]
exec_y <- mod$y[mod$layer=="Execution"]
sense_arr <- data.frame(x=mod_xr, y=mod$y[mod$layer=="Sense"], xend=box_xl, yend=mod$y[mod$layer=="Sense"])
seg_hi <- data.frame(x=mod_xr, y=int_y,  xend=bus_x, yend=int_y)
seg_vi <- data.frame(x=bus_x, y=min(int_y), xend=bus_x, yend=max(int_y))
arr_i  <- data.frame(x=bus_x, y=lay_cfg$cy[2], xend=box_xl, yend=lay_cfg$cy[2])
seg_he <- data.frame(x=mod_xr, y=exec_y, xend=bus_x, yend=exec_y)
seg_ve <- data.frame(x=bus_x, y=min(exec_y), xend=bus_x, yend=max(exec_y))
arr_e  <- data.frame(x=bus_x, y=lay_cfg$cy[3], xend=box_xl, yend=lay_cfg$cy[3])

p <- ggplot() +
  geom_rect(data=prov, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),
            fill=CREAM, colour=ifelse(prov$out, PLUM, MUTED), linewidth=ifelse(prov$out,1.1,0.6)) +
  geom_text(data=prov, aes(x,y,label=label), family=PUB_FONT,
            size=ifelse(prov$out,S_OUT,S_PROV)/.pt, fontface=ifelse(prov$out,2,1), colour=INK, lineheight=1.25) +
  geom_segment(data=prov_arr, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK, arrow=arrow(length=unit(1.8,"mm"),type="closed")) +
  geom_segment(data=sense_arr, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK, arrow=arrow(length=unit(1.9,"mm"),type="closed")) +
  geom_segment(data=seg_hi, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK) +
  geom_segment(data=seg_vi, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK) +
  geom_segment(data=arr_i,  aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK, arrow=arrow(length=unit(1.9,"mm"),type="closed")) +
  geom_segment(data=seg_he, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK) +
  geom_segment(data=seg_ve, aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK) +
  geom_segment(data=arr_e,  aes(x=x,y=y,xend=xend,yend=yend), linewidth=0.55, colour=INK, arrow=arrow(length=unit(1.9,"mm"),type="closed")) +
  geom_rect(data=mod, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, colour=layer, linetype=I(lty)),
            fill=pal_layer_fill[mod$layer], linewidth=0.8) +
  geom_text(data=mod, aes(x=x, y=y+2.8, label=title), family=PUB_FONT, size=S_MODT/.pt, fontface=2, colour=INK) +
  geom_text(data=mod, aes(x=x, y=y-2.9, label=sub),   family=PUB_FONT, size=S_SUB/.pt, colour=INK, parse=TRUE) +
  geom_rect(data=lay, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax, colour=layer),
            fill=pal_layer_fill[lay$layer], linewidth=1.1) +
  geom_text(data=layer_txt, aes(x=x,y=y,label=label), family=PUB_FONT,
            size=layer_txt$sz/.pt, colour=layer_txt$col, parse=TRUE) +
  scale_colour_manual(values=pal_layer, guide="none") +
  annotate("text", x=2,  y=108, label="A", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +
  annotate("text", x=52, y=108, label="B", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +
  coord_fixed(ratio=1, xlim=c(0,180), ylim=c(0,110), expand=FALSE, clip="off") +
  theme_void() + theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(1,1,1,1,"mm"))

save_fig(p, "ICI_Fig1_derivation", width_mm=174, height_mm=174*110/180)
cat("done\n")
