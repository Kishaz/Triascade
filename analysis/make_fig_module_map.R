# =============================================================================
# SUPPLEMENTARY: module map / reference of the locked 98-gene ICI panel.
# The 6 co-occurrence modules (M1-M6) nested into the 3-layer cascade and signed
# (activating + / suppressive -), each listing its member genes and count.
# Reads ICI_98_module_assignment.csv (aligned to the 98 panel). House style.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R"); .pt <- 2.845276
asg <- read.csv("ICI_98_module_assignment.csv", stringsAsFactors=FALSE)
asg$gene <- toupper(asg$gene)
asg$Layer <- c(Sense="Sensory", Integration="Integration", Execution="Execution")[asg$layer]

MODCOL <- c(M1="#C9A66B", M2="#6A2C5A", M3="#8A6D30", M4="#E79FC4", M5="#C44E6E", M6="#2B1C33")
LIGHT  <- c(M1=TRUE, M2=FALSE, M3=FALSE, M4=TRUE, M5=FALSE, M6=FALSE)   # black vs white header text
PROC   <- c(M1="IFN / antigen-pres.", M2="chemokine / traffic.", M3="stromal / ECM",
            M4="angiogenic / CAF", M5="cytotoxic effector", M6="checkpoint / exhaust.")
layer_of <- c(M1="Sensory", M2="Integration", M3="Integration", M4="Integration",
              M5="Execution", M6="Execution")
sign_of  <- c(M1="+", M2="+", M3="−", M4="−", M5="+", M6="−")
col_of   <- c(Sensory=1, Integration=2, Execution=3)
ord_in_col <- list(Sensory=c("M1"), Integration=c("M2","M3","M4"), Execution=c("M5","M6"))
LAYCOL <- c(Sensory=SAND_DEEP, Integration=PLUM, Execution=ROSE)

PER <- 3; LH <- 1.0; HDR <- 2.3; GAP <- 0.7; CW <- 2.0; XSP <- 2.2
wrapgenes <- function(gs){ gs<-sort(gs); n<-length(gs)
  paste(sapply(split(gs, ceiling(seq_along(gs)/PER)), paste, collapse=", "), collapse="\n") }

rects <- list(); htxt <- list(); btxt <- list(); lhead <- list()
for(L in names(ord_in_col)){
  cx <- col_of[L]*XSP; y <- 0
  lhead[[L]] <- data.frame(x=cx, y=1.6, lab=L, col=LAYCOL[L])
  for(m in ord_in_col[[L]]){
    gs <- asg$gene[asg$module==m]; nl <- ceiling(length(gs)/PER)
    bodyh <- nl*LH + 0.5; toth <- HDR + bodyh
    top <- y; bot <- y - toth
    rects[[length(rects)+1]] <- data.frame(xmin=cx-CW/2, xmax=cx+CW/2, ymin=bot, ymax=top,
        hymin=top-HDR, module=m, fill=MODCOL[m])
    htxt[[length(htxt)+1]] <- data.frame(x=cx-CW/2+0.08, y=top-0.75,
        lab1=sprintf("%s  %s", m, PROC[m]), lab2=sprintf("%s  ·  %s  ·  %d genes",
        c("activating","suppressive")[(sign_of[m]=="−")+1], sign_of[m], length(gs)),
        light=LIGHT[m])
    btxt[[length(btxt)+1]] <- data.frame(x=cx-CW/2+0.08, y=top-HDR-0.25, lab=wrapgenes(gs))
    y <- bot - GAP
  }
}
R <- do.call(rbind, rects); H <- do.call(rbind, htxt); B <- do.call(rbind, btxt); LHd <- do.call(rbind, lhead)

p <- ggplot() +
  # layer headers + accent underline
  geom_text(data=LHd, aes(x, y, label=lab, colour=I(col)), family=PUB_FONT, fontface=2, size=(PT_BASE+2)/.pt) +
  geom_segment(data=LHd, aes(x=x-CW/2, xend=x+CW/2, y=0.9, yend=0.9, colour=I(col)), linewidth=1.1) +
  # card body + header
  geom_rect(data=R, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), fill="white", colour="grey80", linewidth=0.3) +
  geom_rect(data=R, aes(xmin=xmin, xmax=xmax, ymin=hymin, ymax=ymax, fill=I(fill)), colour=NA) +
  geom_text(data=H, aes(x, y, label=lab1), hjust=0, vjust=1, family=PUB_FONT, fontface=2,
            size=(PT_BASE-1.5)/.pt, colour=ifelse(H$light, INK, "white")) +
  geom_text(data=H, aes(x, y-0.95, label=lab2), hjust=0, vjust=1, family=PUB_FONT, fontface=3,
            size=(PT_BASE-3)/.pt, colour=ifelse(H$light, INK, "white")) +
  geom_text(data=B, aes(x, y, label=lab), hjust=0, vjust=1, family=PUB_FONT,
            size=(PT_BASE-2.5)/.pt, colour=INK, lineheight=0.95) +
  coord_cartesian(clip="off") +
  scale_y_continuous(expand=expansion(add=c(1,1.5))) +
  scale_x_continuous(expand=expansion(add=c(0.4,0.4))) +
  theme_void() + theme(plot.margin=margin(3,3,3,3,"mm"))

save_fig(p, "ICI_SuppFig_module_map", width_mm=174, height_mm=150, dir="../results/figures/composite")
cat("done - module map. genes/module:\n"); print(table(asg$module))
