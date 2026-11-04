# =============================================================================
# Per-patient cascade profiling: NAME the block, then EXPLAIN it (IMvigor210, n=348).
# Each tumour is scored on the 3 layers; the weakest (floor, z < -0.5) is its
# rate-limiting step / archetype. We then decompose each floor into the panel's
# own signed sub-modules to pinpoint the lesion.
#   (a) exemplar patients, one per archetype, S/I/E profile with the floor marked
#   (b) floor decomposition: archetype x 5 signed sub-modules (median z)
# Reads data/imvigor210/imvigor210_clean.rds. House style; editable SVG.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("ici_theme.R"); source("ici_locked_reference.R"); .pt <- 2.845276
PANEL <- load_panel("ICI_98gene_panel.csv")
obj <- readRDS("data/imvigor210/imvigor210_clean.rds"); E <- obj$E; pd <- obj$pd
rownames(E) <- toupper(rownames(E))
y <- ifelse(pd$binaryResponse=="CR/PR",1L, ifelse(pd$binaryResponse=="SD/PD",0L,NA))
Z  <- function(v){ s<-sd(v,na.rm=TRUE); if(is.na(s)||s==0) v*0 else (v-mean(v,na.rm=TRUE))/s }
sub<- function(gs){ g<-intersect(toupper(gs),rownames(E)); colMeans(t(apply(E[g,,drop=FALSE],1,Z))) }

cm <- score_cohort_median(E, PANEL); cm <- cm[match(colnames(E),cm$sample),]
S <- Z(cm$Sensory); I <- Z(cm$Integration); X <- Z(cm$Execution)
ti <- cbind(Sensory=S, Integration=I, Execution=X)
wmin <- colnames(ti)[apply(ti,1,which.min)]; mz <- apply(ti,1,min)
arch <- factor(ifelse(mz>=-0.5,"intact",ifelse(wmin=="Sensory","cold",
         ifelse(wmin=="Integration","blocked","no-kill"))), c("cold","blocked","no-kill","intact"))
ALAB <- c(cold="Cold\n(Sensory floor)", blocked="Blocked\n(Integration floor)",
          `no-kill`="No-kill\n(Execution floor)", intact="Intact\n(no floor)")
AWHY <- c(cold="luminal/UroA immune-cold", blocked="TGF-β / stromal exclusion",
          `no-kill`="effector deficit", intact="inflamed yet resistant")

## ---- panel a: exemplar patient per archetype -------------------------------
# pick the clearest exemplar: a non-responder whose floor is unambiguous
pick <- function(a){ idx <- which(arch==a & y==0 & !is.na(y))
  prof <- ti[idx,,drop=FALSE]
  floorv <- apply(prof,1,min); nextv <- apply(prof,1,function(r) sort(r)[2])
  idx[ which.max((nextv-floorv) - 0.15*abs(floorv-median(mz))) ] }       # clear gap, typical depth
ex <- sapply(levels(arch), pick)
EA <- do.call(rbind, lapply(levels(arch), function(a){ i<-ex[a]
  data.frame(arch=a, layer=c("Sensory","Integration","Execution"), z=ti[i,],
             is_floor=(c("Sensory","Integration","Execution")==wmin[i]) & a!="intact", row.names=NULL) }))
EA$arch  <- factor(ALAB[EA$arch], levels=unname(ALAB))
EA$layer <- factor(EA$layer, c("Sensory","Integration","Execution"))

FL <- EA[EA$is_floor, ]; FL$yl <- pmin(FL$z,0) - 0.30
pa <- ggplot(EA, aes(layer, z, fill=z>=0)) +
  geom_hline(yintercept=0, colour=RULE, linewidth=0.3) +
  geom_col(width=0.68, colour="white", linewidth=0.3) +
  geom_hline(yintercept=-0.5, linetype="22", colour=MUTED, linewidth=0.3) +
  geom_text(data=FL, aes(layer, yl, label="floor"), inherit.aes=FALSE, fontface=3,
            family=PUB_FONT, size=(PT_BASE-2.5)/.pt, colour=INK) +
  facet_wrap(~arch, nrow=1) +
  scale_fill_manual(values=c(`TRUE`=PLUM, `FALSE`=SAND_DEEP), guide="none") +
  scale_y_continuous("Layer score (within-cohort z)", limits=c(-1.85,1.85), breaks=seq(-1.5,1.5,1.5)) +
  scale_x_discrete(NULL, labels=c("S","I","E")) +
  theme_publication(grid="y") +
  theme(axis.text=element_text(size=PT_BASE, colour=INK),
        axis.title.y=element_text(size=PT_AXIS_TITLE, colour=INK, margin=margin(r=1)),
        strip.text=element_text(size=PT_BASE-0.5, face=2, colour=INK, lineheight=0.95),
        panel.spacing=unit(3,"mm"), plot.margin=margin(2,3,1,2,"mm"))

## ---- panel b: floor decomposition (archetype x signed sub-modules) ---------
SUBM <- list(`Sensory (+)`=PANEL$SENSORY, `Trafficking (+)`=PANEL$INTPOS,
             `Stroma/Treg (−)`=PANEL$INTNEG, `Effector (+)`=PANEL$EXECPOS,
             `Exhaustion (−)`=PANEL$EXECNEG)
SM <- sapply(SUBM, sub)                                    # n x 5
DM <- do.call(rbind, lapply(levels(arch), function(a)
        data.frame(arch=a, module=colnames(SM), z=apply(SM[arch==a,,drop=FALSE],2,median), row.names=NULL)))
DM$arch   <- factor(DM$arch, levels=rev(levels(arch)))      # cold at bottom -> intact top
DM$module <- factor(DM$module, levels=names(SUBM))
lim <- max(abs(DM$z))
pb <- ggplot(DM, aes(module, arch, fill=z)) +
  geom_tile(colour="white", linewidth=1.1) +
  geom_text(aes(label=sprintf("%+.2f", z), colour=abs(z)>0.45*lim), size=(PT_BASE-1.5)/.pt, family=PUB_FONT) +
  scale_fill_gradient2(low=SAND_DEEP, mid="#FBF7F0", high=PLUM, midpoint=0, limits=c(-lim,lim),
                       name="median z", guide=guide_colourbar(barheight=unit(22,"mm"), barwidth=unit(2.6,"mm"))) +
  scale_colour_manual(values=c(`TRUE`="white",`FALSE`=INK), guide="none") +
  scale_x_discrete(NULL, position="top") + scale_y_discrete(NULL) +
  coord_equal() +
  theme_publication() +
  theme(axis.text.x=element_text(size=PT_BASE-0.5, colour=INK, angle=35, hjust=0),
        axis.text.y=element_text(size=PT_BASE, colour=INK, face=2),
        axis.ticks=element_blank(), panel.border=element_blank(),
        legend.title=element_text(size=PT_BASE-1, colour=INK),
        plot.margin=margin(2,3,2,2,"mm"))

## ---- panel c: floor validates against the Mariathasan immune phenotype -----
ph <- pd[["Immune phenotype"]]; keep <- !is.na(ph)
PC <- as.data.frame(prop.table(table(arch[keep], ph[keep]), 1))
names(PC) <- c("arch","phenotype","frac")
PC$arch <- factor(PC$arch, levels=levels(arch))
PC$phenotype <- factor(PC$phenotype, levels=c("desert","excluded","inflamed"))
chip <- chisq.test(table(arch[keep], ph[keep]))$p.value
pc <- ggplot(PC, aes(arch, frac, fill=phenotype)) +
  geom_col(width=0.74, colour="white", linewidth=0.3) +
  annotate("text", x=0.5, y=1.07, hjust=0, vjust=0, family=PUB_FONT, fontface=3,
           size=(PT_BASE-1.5)/.pt, colour=INK, label=sprintf("χ² p = %.0e", chip)) +
  scale_fill_manual(values=c(desert=SAND, excluded=ROSE, inflamed=PLUM), name=NULL) +
  scale_y_continuous("fraction of tumours", limits=c(0,1.12), breaks=seq(0,1,0.5), expand=expansion(mult=0)) +
  scale_x_discrete(NULL) +
  theme_publication(grid="y") +
  theme(axis.text.x=element_text(size=PT_BASE-0.5, colour=INK, angle=20, hjust=1),
        axis.text.y=element_text(size=PT_BASE, colour=INK),
        axis.title.y=element_text(size=PT_AXIS_TITLE, colour=INK, margin=margin(r=1)),
        legend.position="top", legend.justification="left", legend.box.spacing=unit(1,"mm"),
        legend.text=element_text(size=PT_BASE-1, colour=INK), plot.margin=margin(2,2,2,2,"mm"))

## ---- panel d: flagship ORTHOGONAL driver — TGF-beta/stromal exclusion -------
# genes outside the 98-gene panel; explains the Integration ("blocked") floor
PG <- unique(toupper(unlist(PANEL)))
tgfb_g <- setdiff(toupper(c("ACTA2","ACTG2","ADAM12","CNN1","COL4A1","COL5A1","COL5A2","COL6A3",
  "COMP","CTGF","FAP","FBN1","FN1","LRRC15","POSTN","PDGFRB","TAGLN","TGFB1","TGFBR2","VCAN",
  "THBS2","SPARC")), PG)
DD <- data.frame(arch=arch, tgfb=sub(tgfb_g))
pbl <- wilcox.test(DD$tgfb[arch=="blocked"], DD$tgfb[arch!="blocked"])$p.value
pd_ <- ggplot(DD, aes(arch, tgfb, fill=arch=="blocked")) +
  geom_hline(yintercept=0, colour=RULE, linewidth=0.3) +
  geom_boxplot(outlier.size=0.4, linewidth=0.3, width=0.62, colour=INK) +
  annotate("text", x=2, y=2.5, label=sprintf("blocked vs rest\np = %.0e", pbl), hjust=0.5, vjust=1,
           family=PUB_FONT, fontface=3, size=(PT_BASE-2)/.pt, colour=INK, lineheight=0.95) +
  scale_fill_manual(values=c(`TRUE`=PLUM, `FALSE`="#EFE7DD"), guide="none") +
  scale_y_continuous("TGF-β / stroma (orthogonal z)", limits=c(-1.6,2.7), breaks=seq(-1,2,1)) +
  scale_x_discrete(NULL) +
  theme_publication(grid="y") +
  theme(axis.text.x=element_text(size=PT_BASE-0.5, colour=INK, angle=20, hjust=1),
        axis.text.y=element_text(size=PT_BASE, colour=INK),
        axis.title.y=element_text(size=PT_AXIS_TITLE-0.5, colour=INK, margin=margin(r=1)),
        plot.margin=margin(2,3,2,2,"mm"))

## ---- panel e: no-kill REFINEMENT — two actionable flavours ------------------
nk <- arch=="no-kill"
effM <- SM[,"Effector (+)"]; exhM <- SM[,"Exhaustion (−)"]
tex_g <- setdiff(c("LAYN","TNFRSF9","TOX","CXCL13"), PG); tex <- sub(tex_g)
flav <- ifelse(effM[nk] < median(effM), "effector-absent", "primed-but-exhausted")
NK <- rbind(
  data.frame(flavour=flav, feature="Effector (+)",     z=effM[nk]),
  data.frame(flavour=flav, feature="Exhaustion (−)", z=exhM[nk]),   # higher = more exhaustion expressed
  data.frame(flavour=flav, feature="Term-exhaustion\n(LAYN/TNFRSF9)", z=tex[nk]))
NK$feature <- factor(NK$feature, levels=c("Effector (+)","Exhaustion (−)","Term-exhaustion\n(LAYN/TNFRSF9)"))
AGG <- aggregate(z~flavour+feature, NK, median)
ptex <- wilcox.test(tex[nk]~flav)$p.value
nlab <- sprintf("%s (n=%d)", c("effector-absent","primed-but-exhausted"), c(sum(flav=="effector-absent"),sum(flav=="primed-but-exhausted")))
AGG$flavour <- factor(AGG$flavour, c("effector-absent","primed-but-exhausted"), nlab)
pe <- ggplot(AGG, aes(feature, z, fill=flavour)) +
  geom_hline(yintercept=0, colour=RULE, linewidth=0.3) +
  geom_col(position=position_dodge(width=0.72), width=0.64, colour="white", linewidth=0.3) +
  annotate("text", x=3, y=0.95, label=sprintf("p = %.0e", ptex), hjust=0.5, vjust=1,
           family=PUB_FONT, fontface=3, size=(PT_BASE-2.5)/.pt, colour=INK) +
  scale_fill_manual(values=setNames(c(SAND, PLUM), nlab), name=NULL) +
  scale_y_continuous("median z", limits=c(-0.7,1.0), breaks=seq(-0.5,1,0.5)) +
  scale_x_discrete(NULL) +
  theme_publication(grid="y") +
  theme(axis.text.x=element_text(size=PT_BASE-1.5, colour=INK, lineheight=0.9),
        axis.text.y=element_text(size=PT_BASE, colour=INK),
        axis.title.y=element_text(size=PT_AXIS_TITLE-0.5, colour=INK, margin=margin(r=1)),
        legend.position="top", legend.justification="left", legend.box.spacing=unit(0.5,"mm"),
        legend.text=element_text(size=PT_BASE-1.5, colour=INK), plot.margin=margin(2,2,2,2,"mm"))

tag_th <- theme(plot.tag=element_text(size=PT_PANEL_TAG, face="bold", family=PUB_FONT, colour=INK, hjust=0, vjust=1))
pa <- pa + labs(tag="a") + tag_th; pb <- pb + labs(tag="b") + tag_th; pc <- pc + labs(tag="c") + tag_th
pd_ <- pd_ + labs(tag="d") + tag_th; pe <- pe + labs(tag="e") + tag_th
row2 <- wrap_elements(full=pb) + wrap_elements(full=pc) + plot_layout(widths=c(1.5, 1))
row3 <- wrap_elements(full=pd_) + wrap_elements(full=pe) + plot_layout(widths=c(1, 1.25))
fig <- wrap_elements(full=pa) / row2 / row3 + plot_layout(heights=c(0.9, 1.15, 1.0)) &
  theme(plot.margin=margin(2,2,2,2,"mm"))
save_fig(fig, "ICI_Fig_archetype_biology", width_mm=174, height_mm=200, dir="../results/figures/composite")
cat("exemplars (sample | archetype):\n"); for(a in levels(arch)) cat(sprintf("  %-8s %s\n", a, colnames(E)[ex[a]]))
cat("done — panels a (exemplars) + b (floor decomposition)\n")
