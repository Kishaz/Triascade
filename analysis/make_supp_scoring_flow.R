# =============================================================================
# SUPPLEMENTARY: how one patient's gene expression is DISSOCIATED and scored into
# the three ICI layer scores. Unlike a learned-weight network, the ICI is a
# DETERMINISTIC, SIGNED cascade: each panel gene maps to a fixed sub-module and
# layer, and a layer = activating (+) minus suppressive (-). Flow, left -> right:
#   gene expression (within-cohort z)  ->  gene indicator (>= locked ref median)
#   ->  signed sub-module  ->  the 3 layer scores  ->  floor / archetype call.
# House style; genes grouped BY module so edges read cleanly; floor marked.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R"); source("ici_locked_reference.R"); .pt <- 2.845276
PANEL <- load_panel("ICI_98gene_panel.csv")
obj <- readRDS("data/imvigor210/imvigor210_clean.rds"); E <- obj$E; rownames(E) <- toupper(rownames(E))
Zc <- function(v){ s<-sd(v,na.rm=TRUE); if(is.na(s)||s==0) v*0 else (v-mean(v,na.rm=TRUE))/s }

## ---- fixed gene -> submodule -> layer map (the cascade wiring) --------------
layer_of <- c(SENSORY="Sensory", INTPOS="Integration", INTNEG="Integration",
              EXECPOS="Execution", EXECNEG="Execution")
sign_of  <- c(SENSORY=1, INTPOS=1, INTNEG=-1, EXECPOS=1, EXECNEG=-1)
smlab    <- c(SENSORY="Sensory (+)", INTPOS="Trafficking (+)", INTNEG="Stroma/Treg (−)",
              EXECPOS="Effector (+)", EXECNEG="Exhaustion (−)")
SMORD <- c("SENSORY","INTPOS","INTNEG","EXECPOS","EXECNEG")
slots <- do.call(rbind, lapply(SMORD, function(sm){
  g <- toupper(PANEL[[sm]]); data.frame(submod=sm, gene=g,
       layer=layer_of[sm], sign=sign_of[sm], present=g %in% rownames(E), row.names=NULL) }))

## ---- vertical layout: genes stacked by module, with gaps between modules ----
GAP <- 2.2
slots$row <- NA_real_; yy <- 0
for(sm in SMORD){ idx <- which(slots$submod==sm)
  slots$row[idx] <- yy + seq_along(idx); yy <- max(slots$row[idx]) + GAP }
slots$y <- max(slots$row) - slots$row + 1                  # top = first module
sm_y <- tapply(slots$y, slots$submod, mean)[SMORD]
LORD <- c("Sensory","Integration","Execution")
ly_y <- sapply(LORD, function(L) mean(sm_y[names(layer_of)[layer_of==L]]))

## ---- columns (x positions) -------------------------------------------------
x_gene <- 0; x_ind <- 1.5; x_sm <- 3.4; x_lay <- 5.2; x_call <- 6.3

draw_patient <- function(samp, resp, phen, pid){
  med <- apply(E,1,median,na.rm=TRUE)
  zmat <- t(apply(E,1,Zc))
  d <- slots; d$z <- NA_real_; d$on <- NA_integer_
  pg <- d$gene[d$present]
  d$z[d$present]  <- zmat[pg, samp]
  d$on[d$present] <- as.integer(E[pg, samp] >= med[pg])
  # layer scores for THIS patient (cohort-median scoring, z within cohort) + floor
  cm <- score_cohort_median(E, PANEL)
  LZ <- c(Sensory=Zc(cm$Sensory)[cm$sample==samp], Integration=Zc(cm$Integration)[cm$sample==samp],
          Execution=Zc(cm$Execution)[cm$sample==samp])
  floorL <- names(LZ)[which.min(LZ)]; isfloor <- min(LZ) < -0.5
  # submodule node value = signed fraction "on"
  smv <- sapply(SMORD, function(sm){ mean(d$z[d$submod==sm & d$present], na.rm=TRUE) })  # expression magnitude

  LN <- data.frame(layer=LORD, y=ly_y, z=LZ[LORD], floor=(LORD==floorL & isfloor))
  LN$namelab <- ifelse(LN$floor, paste0(LN$layer, "  ← floor"), LN$layer)
  SN <- data.frame(submod=SMORD, y=sm_y, lab=smlab[SMORD], v=smv[SMORD],
                   sign=sign_of[SMORD], ly=ly_y[layer_of[SMORD]])
  EG <- within(d[d$present,], { x1<-x_ind; y1<-y; x2<-x_sm; y2<-sm_y[submod] })

  ggplot() +
    # edges: gene -> submodule (colour by expression z, emphasise ON genes)
    geom_segment(data=EG, aes(x=x1,y=y1,xend=x2,yend=y2, colour=z, alpha=on),
                 linewidth=0.22) +
    # edges: submodule -> layer (colour by sign)
    geom_segment(data=SN, aes(x=x_sm, y=y, xend=x_lay, yend=ly, colour=sign*2),
                 linewidth=0.5, alpha=0.9) +
    # gene expression dots + labels
    geom_point(data=d[d$present,], aes(x_gene, y, fill=z), shape=21, size=1.5, stroke=0.15, colour="white") +
    geom_text(data=d[d$present,], aes(x_gene-0.12, y, label=gene), hjust=1,
              family=PUB_FONT, size=3.2/.pt, colour=INK) +
    # gene indicator (filled=on)
    geom_point(data=d[d$present,], aes(x_ind, y, fill=on), shape=22, size=1.3, stroke=0.15, colour=MUTED) +
    # submodule nodes + labels
    geom_point(data=SN, aes(x_sm, y, fill=v), shape=21, size=4.4, stroke=0.3, colour="white") +
    geom_text(data=SN, aes(x_sm+0.18, y, label=lab), hjust=0, family=PUB_FONT,
              size=6/.pt, colour=INK, fontface=2) +
    # layer nodes + labels + floor ring
    geom_point(data=LN, aes(x_lay, y, fill=z), shape=21, size=8, stroke=0.3, colour="white") +
    geom_point(data=subset(LN,floor), aes(x_lay, y), shape=21, size=9.2, stroke=1.1, colour=INK, fill=NA) +
    geom_text(data=LN, aes(x_lay, y, label=substr(layer,1,1)), family=PUB_FONT,
              size=7/.pt, fontface=2, colour=ifelse(abs(LN$z)>0.8,"white",INK)) +
    geom_text(data=subset(LN,!floor), aes(x_lay+0.42, y, label=namelab), hjust=0, family=PUB_FONT,
              size=7/.pt, fontface=2, colour=INK) +
    geom_text(data=subset(LN,floor), aes(x_lay+0.42, y, label=namelab), hjust=0, family=PUB_FONT,
              size=7/.pt, fontface=2, colour=PLUM) +
    { if(!isfloor) annotate("text", x=x_lay+0.42, y=min(ly_y)-4, hjust=0, family=PUB_FONT,
             size=6.5/.pt, fontface=2, colour=PLUM, label="no floor (intact)") } +
    scale_fill_gradient2(low=SAND_DEEP, mid="#FBF7F0", high=PLUM, midpoint=0, guide="none") +
    scale_colour_gradient2(low=SAND_DEEP, mid="#E8DFD3", high=PLUM, midpoint=0, guide="none") +
    scale_alpha_continuous(range=c(0.12,0.7), guide="none") +
    annotate("text", x=c(x_gene, x_ind, x_sm, x_lay), y=max(slots$y)+5,
             label=c("Gene (z)","Indicator","Sub-module","Layer"),
             family=PUB_FONT, size=(PT_BASE-2)/.pt, fontface=2, colour=INK, lineheight=0.92) +
    scale_x_continuous(limits=c(-1.4, 7.3)) +
    scale_y_continuous(expand=expansion(add=c(2,7))) +
    labs(title=sprintf("%s  ·  %s, %s", pid, phen, resp)) +
    theme_void() +
    theme(plot.title=element_text(size=PT_BASE, face=2, family=PUB_FONT, colour=INK, hjust=0.5),
          plot.margin=margin(2,3,2,3,"mm"))
}

p <- draw_patient("SAM18bc1078bc15", "non-responder", "immune-cold", pid="P1")
save_fig(p, "ICI_SuppFig_scoring_flow", width_mm=118, height_mm=220, dir="../results/figures/composite")
cat("done — single-patient scoring-flow prototype\n")
