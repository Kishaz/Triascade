# =============================================================================
# SUPPLEMENTARY: how each patient's gene expression is DISSOCIATED and scored
# into the three ICI layer scores. Unlike a learned-weight network, the ICI is a
# DETERMINISTIC, SIGNED cascade: each panel gene maps to a fixed sub-module and
# layer, and a layer = activating (+) minus suppressive (-). Flow, left -> right:
#   gene z  ->  indicator (>= locked ref median)  ->  signed sub-module
#   ->  the 3 layer scores (floor ringed).
# 2x2: one clean intact RESPONDER (no floor) + three floored NON-responders
# (cold / blocked / no-kill). All 98 gene dots shown; ~20 marker genes labelled.
# House style; genes grouped BY module so edges read cleanly; editable SVG.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("ici_theme.R"); source("ici_locked_reference.R"); .pt <- 2.845276
PANEL <- load_panel("ICI_98gene_panel.csv")
obj <- readRDS("data/imvigor210/imvigor210_clean.rds"); E <- obj$E; pd <- obj$pd
rownames(E) <- toupper(rownames(E))
y <- ifelse(pd$binaryResponse=="CR/PR",1L, ifelse(pd$binaryResponse=="SD/PD",0L,NA))
Zc <- function(v){ s<-sd(v,na.rm=TRUE); if(is.na(s)||s==0) v*0 else (v-mean(v,na.rm=TRUE))/s }

## ---- cascade wiring --------------------------------------------------------
layer_of <- c(SENSORY="Sensory", INTPOS="Integration", INTNEG="Integration",
              EXECPOS="Execution", EXECNEG="Execution")
sign_of  <- c(SENSORY=1, INTPOS=1, INTNEG=-1, EXECPOS=1, EXECNEG=-1)
smlab    <- c(SENSORY="Sensory", INTPOS="Trafficking", INTNEG="Stroma/Treg",
              EXECPOS="Effector", EXECNEG="Exhaustion")
SMORD <- c("SENSORY","INTPOS","INTNEG","EXECPOS","EXECNEG")
LORD  <- c("Sensory","Integration","Execution")
LABEL <- c("HLA-A","B2M","STAT1","TAP1","CXCL9","CXCL10","CXCR3","CCL5","TGFB1","FAP",
           "VEGFA","COL1A1","GZMB","PRF1","CD8A","NKG7","PDCD1","LAG3","TOX","IDO1")
slots <- do.call(rbind, lapply(SMORD, function(sm){
  g <- toupper(PANEL[[sm]]); data.frame(submod=sm, gene=g, layer=layer_of[sm],
       sign=sign_of[sm], present=g %in% rownames(E), row.names=NULL) }))
GAP <- 1.4; slots$row <- NA_real_; yy <- 0
for(sm in SMORD){ idx <- which(slots$submod==sm)
  slots$row[idx] <- yy + seq_along(idx); yy <- max(slots$row[idx]) + GAP }
slots$y <- max(slots$row) - slots$row + 1
sm_y <- tapply(slots$y, slots$submod, mean)[SMORD]
ly_y <- sapply(LORD, function(L) mean(sm_y[names(layer_of)[layer_of==L]]))
x_gene<-0; x_ind<-0.85; x_sm<-2.6; x_lay<-4.3
zmat <- t(apply(E,1,Zc)); med <- apply(E,1,median,na.rm=TRUE)
cm <- score_cohort_median(E, PANEL)
LZall <- cbind(Sensory=Zc(cm$Sensory), Integration=Zc(cm$Integration), Execution=Zc(cm$Execution))
rownames(LZall) <- cm$sample

draw_patient <- function(samp, pid, sub){
  d <- slots; d$z <- NA_real_; d$on <- NA_integer_; pg <- d$gene[d$present]
  d$z[d$present] <- zmat[pg, samp]; d$on[d$present] <- as.integer(E[pg,samp] >= med[pg])
  d$lab <- ifelse(d$gene %in% LABEL, d$gene, NA)
  LZ <- LZall[samp,]
  onv <- function(sm) mean(d$on[d$submod==sm & d$present], na.rm=TRUE)   # fraction of genes "on"
  raw <- c(Sensory=onv("SENSORY"), Integration=onv("INTPOS")-onv("INTNEG"),
           Execution=onv("EXECPOS")-onv("EXECNEG"))                       # normalised layer scores
  smv <- sapply(SMORD, function(sm) mean(d$z[d$submod==sm & d$present], na.rm=TRUE))
  SN <- data.frame(submod=SMORD, y=sm_y, lab=paste0(smlab[SMORD], c(" +"," +"," −"," +"," −")[match(SMORD,SMORD)]),
                   v=smv[SMORD], sign=sign_of[SMORD], ly=ly_y[layer_of[SMORD]])
  LN <- data.frame(layer=LORD, y=ly_y, z=LZ[LORD], raw=raw[LORD])
  LN$scorelab <- sprintf("%s  %.2f", LN$layer, LN$raw)
  dd <- d[d$present,]
  ggplot() +
    geom_segment(data=dd, aes(x=x_ind,y=y,xend=x_sm,yend=sm_y[submod], colour=z, alpha=on), linewidth=0.15) +
    geom_segment(data=SN, aes(x=x_sm,y=y,xend=x_lay,yend=ly, colour=sign*2.5), linewidth=0.45, alpha=0.9) +
    geom_point(data=dd, aes(x_gene,y,fill=z), shape=21, size=0.75, stroke=0.1, colour="white") +
    geom_point(data=dd, aes(x_ind,y,fill=on), shape=22, size=0.7, stroke=0.1, colour=MUTED) +
    ggrepel::geom_text_repel(data=subset(dd,!is.na(lab)), aes(x_gene, y, label=lab),
              family=PUB_FONT, size=4.4/.pt, colour=INK, direction="y", hjust=1, nudge_x=-0.55,
              segment.size=0.12, segment.colour=MUTED, min.segment.length=0, box.padding=0.02,
              point.padding=0, max.overlaps=Inf, xlim=c(NA, x_gene-0.1)) +
    geom_point(data=SN, aes(x_sm,y,fill=v), shape=21, size=2.8, stroke=0.25, colour="white") +
    geom_text(data=SN, aes(x_sm, y+2.4, label=lab), hjust=0.5, vjust=0, family=PUB_FONT, size=5/.pt, colour=INK, fontface=2) +
    geom_point(data=LN, aes(x_lay,y,fill=z), shape=21, size=5.2, stroke=0.25, colour="white") +
    geom_text(data=LN, aes(x_lay,y,label=substr(layer,1,1)), family=PUB_FONT, size=5/.pt, fontface=2,
              colour=ifelse(abs(LN$z)>0.8,"white",INK)) +
    geom_text(data=LN, aes(x_lay+0.38,y,label=scorelab), hjust=0, family=PUB_FONT, size=5/.pt, fontface=2, colour=INK) +
    scale_fill_gradient2(low=SAND_DEEP, mid="#F4ECE0", high=PLUM, midpoint=0, limits=c(-1.8,1.8), oob=scales::squish, guide="none") +
    scale_colour_gradient2(low=SAND_DEEP, mid="#DCCFBE", high=PLUM, midpoint=0, limits=c(-1.8,1.8), oob=scales::squish, guide="none") +
    scale_alpha_continuous(range=c(0.22,0.9), guide="none") +
    scale_x_continuous(limits=c(-1.7, 6.4)) + scale_y_continuous(expand=expansion(add=c(1.5,3))) +
    labs(title=sprintf("%s · %s", pid, sub)) +
    theme_void() +
    theme(plot.title=element_text(size=PT_BASE-1, face=2, family=PUB_FONT, colour=INK, hjust=0.5),
          plot.margin=margin(1,2,1,2,"mm"))
}

## ---- pick the four patients ------------------------------------------------
arch <- apply(LZall,1,function(r){ w<-LORD[which.min(r)]; if(min(r)>=-0.5) "intact" else
  c(Sensory="cold",Integration="blocked",Execution="no-kill")[w] })
ys <- setNames(y, colnames(E))
pickNR <- function(a){ i<-which(arch==a & ys==0); j<-i[which.max(apply(LZall[i,,drop=FALSE],1,function(r) sort(r)[2]-min(r)))]; names(arch)[j] }
pickR_intact <- { i<-which(arch=="intact" & ys==1); names(arch)[i[which.max(apply(LZall[i,,drop=FALSE],1,min))]] }
P <- c(R=pickR_intact, cold=pickNR("cold"), blocked=pickNR("blocked"), nokill=pickNR("no-kill"))

pR  <- draw_patient(P["R"],      "R1", "intact responder")
pC  <- draw_patient(P["cold"],   "N1", "cold non-responder")
pB  <- draw_patient(P["blocked"],"N2", "blocked non-responder")
pK  <- draw_patient(P["nokill"], "N3", "no-kill non-responder")

## ---- shared legend ---------------------------------------------------------
leg <- ggplot() + xlim(0,10) + ylim(0.6,2.4) + theme_void() +
  # z gradient swatches
  annotate("text", x=0, y=2, label="Gene / layer z:", hjust=0, family=PUB_FONT, size=(PT_BASE-1.5)/.pt, fontface=2, colour=INK) +
  annotate("point", x=1.9+(0:4)*0.33, y=2, shape=21, size=2.2, colour="white",
           fill=colorRampPalette(c(SAND_DEEP,"#F4ECE0",PLUM))(5)) +
  annotate("text", x=1.75, y=2, label="low", hjust=1, family=PUB_FONT, size=(PT_BASE-2)/.pt, colour=INK) +
  annotate("text", x=3.7, y=2, label="high", hjust=0, family=PUB_FONT, size=(PT_BASE-2)/.pt, colour=INK) +
  # indicator
  annotate("point", x=5.3, y=2, shape=22, size=1.8, fill=PLUM, colour=MUTED, stroke=0.2) +
  annotate("text", x=5.5, y=2, label="≥ ref median (on)", hjust=0, family=PUB_FONT, size=(PT_BASE-2)/.pt, colour=INK) +
  annotate("point", x=0.1, y=1, shape=22, size=1.8, fill="#FBF7F0", colour=MUTED, stroke=0.2) +
  annotate("text", x=0.3, y=1, label="< ref median (off)", hjust=0, family=PUB_FONT, size=(PT_BASE-2)/.pt, colour=INK) +
  # edge sign
  annotate("segment", x=2.9, xend=3.4, y=1, yend=1, colour=PLUM, linewidth=0.8) +
  annotate("text", x=3.5, y=1, label="activating (+)", hjust=0, family=PUB_FONT, size=(PT_BASE-2)/.pt, colour=INK) +
  annotate("segment", x=5.3, xend=5.8, y=1, yend=1, colour=SAND_DEEP, linewidth=0.8) +
  annotate("text", x=5.9, y=1, label="suppressive (−)", hjust=0, family=PUB_FONT, size=(PT_BASE-2)/.pt, colour=INK) +
  theme(plot.margin=margin(1,2,1,2,"mm"))

tag_th <- theme(plot.tag=element_text(size=PT_PANEL_TAG, face="bold", family=PUB_FONT, colour=INK, hjust=0, vjust=1))
pR<-pR+labs(tag="a")+tag_th; pC<-pC+labs(tag="b")+tag_th; pB<-pB+labs(tag="c")+tag_th; pK<-pK+labs(tag="d")+tag_th
grid <- (wrap_elements(full=pR) + wrap_elements(full=pC)) / (wrap_elements(full=pB) + wrap_elements(full=pK))
fig  <- grid / leg + plot_layout(heights=c(1,1,0.14)) & theme(plot.margin=margin(1,1,1,1,"mm"))
save_fig(fig, "ICI_SuppFig_scoring_flow", width_mm=174, height_mm=224, dir="../results/figures/composite")
cat("patients:\n"); for(n in names(P)) cat(sprintf("  %-7s %s\n", n, P[n]))
cat("done — 2x2 scoring-flow (1 responder + 3 floored NRs)\n")
