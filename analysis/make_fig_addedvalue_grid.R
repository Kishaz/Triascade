# =============================================================================
# Added-value as a 2x2 grid (teaching version, headline marker = T-cell-inflamed GEP).
# Every validated patient is split two ways at the within-cohort median:
#   columns = GEP (established marker) Low / High
#   rows    = ICI Integration score   Low / High
# Each cell = ICB response rate (%) for patients in that box, with n.
# Margins show each split's OWN effect: the ICI (row) gap is large, the GEP (column)
# gap is small -> ICI adds information GEP misses. Reads pancancer_validation.rds.
# House style; single canvas; editable SVG (no textLength).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(grid) })
source("ici_theme.R"); .pt <- 2.845276
R <- readRDS("pancancer_validation.rds")
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
hi <- function(v) ifelse(v > median(v), "High", "Low")
D <- do.call(rbind, lapply(WORKS, function(k){ s<-R[[k]]$scores
  data.frame(y=s$y, ICI=hi(s$integ_single), GEP=hi(s$GEP)) }))
base <- 100*mean(D$y)
rate <- function(sub){ v<-D$y[sub]; c(pct=100*mean(v), n=length(v)) }

## four cells + margins -------------------------------------------------------
cell <- expand.grid(GEP=c("Low","High"), ICI=c("Low","High"), stringsAsFactors=FALSE)
cell <- cbind(cell, t(mapply(function(g,i) rate(D$GEP==g & D$ICI==i), cell$GEP, cell$ICI)))
cell$cx <- ifelse(cell$GEP=="Low",1,2); cell$cy <- ifelse(cell$ICI=="Low",1,2)
mICI <- sapply(c("Low","High"), function(i) rate(D$ICI==i))      # row margins (ICI effect)
mGEP <- sapply(c("Low","High"), function(g) rate(D$GEP==g))      # col margins (GEP effect)

p <- ggplot(cell, aes(cx, cy)) +
  geom_tile(aes(fill=pct), colour="white", linewidth=2.2, width=0.96, height=0.96) +
  geom_text(aes(label=sprintf("%.0f%%", pct), colour=pct>32), vjust=-0.1,
            family=PUB_FONT, size=(PT_BASE+9)/.pt, fontface=2) +
  geom_text(aes(label=sprintf("n = %d", n), colour=pct>32), vjust=2.9,
            family=PUB_FONT, size=(PT_BASE-1)/.pt) +
  scale_fill_gradient(low="#F3EAF0", high=PLUM, limits=c(12,50), oob=scales::squish, guide="none") +
  scale_colour_manual(values=c(`TRUE`="white",`FALSE`=INK), guide="none") +
  ## axis labels
  annotate("text", x=1.5, y=2.74, label="T-cell-inflamed GEP  (established marker)",
           family=PUB_FONT, size=(PT_BASE+1)/.pt, fontface=2, colour=INK) +
  annotate("text", x=1, y=2.60, label="Low", family=PUB_FONT, size=PT_BASE/.pt, colour=MUTED) +
  annotate("text", x=2, y=2.60, label="High", family=PUB_FONT, size=PT_BASE/.pt, colour=MUTED) +
  annotate("text", x=0.26, y=1.5, label="ICI Integration score", angle=90,
           family=PUB_FONT, size=(PT_BASE+1)/.pt, fontface=2, colour=INK) +
  annotate("text", x=0.40, y=1, label="Low",  angle=90, family=PUB_FONT, size=PT_BASE/.pt, colour=MUTED) +
  annotate("text", x=0.40, y=2, label="High", angle=90, family=PUB_FONT, size=PT_BASE/.pt, colour=MUTED) +
  ## key comparison 1: down the GEP-high column, ICI low -> high (big jump)
  annotate("segment", x=2.60, xend=2.60, y=1.2, yend=1.8,
           arrow=arrow(length=unit(2.4,"mm"), type="closed"), linewidth=0.9, colour=PLUM) +
  ## key comparison 2: across the ICI-low row, GEP low -> high (flat)
  annotate("segment", x=1.2, xend=1.8, y=0.44, yend=0.44,
           arrow=arrow(length=unit(2.0,"mm"), type="closed"), linewidth=0.7, colour=SAND_DEEP) +
  coord_fixed(clip="off", xlim=c(0.2,4.5), ylim=c(0.35,2.82)) +
  theme_void(base_family=PUB_FONT) +
  theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(3,3,3,3,"mm"))

## ---- plain-words read-out to the right (conditional = the real added value) --
LX <- 2.74
p <- p +
  annotate("text", x=LX, y=2.60, hjust=0, family=PUB_FONT, size=PT_BASE/.pt, fontface=2, colour=INK,
           label="How to read it") +
  annotate("point", x=LX+0.03, y=2.28, size=4, colour=PLUM) +
  annotate("text", x=LX+0.12, y=2.28, hjust=0, lineheight=0.95, family=PUB_FONT,
           size=(PT_BASE-0.5)/.pt, colour=INK,
           label="Down a column (same GEP, ICI differs):") +
  annotate("text", x=LX+0.12, y=2.07, hjust=0, family=PUB_FONT, size=(PT_BASE-0.5)/.pt, colour=PLUM, fontface=2,
           label="GEP-high goes 21% → 46% as ICI rises") +
  annotate("text", x=LX+0.12, y=1.90, hjust=0, family=PUB_FONT, size=(PT_BASE-1)/.pt, colour=MUTED,
           label="→ ICI finds responders GEP alone misses") +
  annotate("point", x=LX+0.03, y=1.46, size=4, colour=SAND_DEEP) +
  annotate("text", x=LX+0.12, y=1.46, hjust=0, family=PUB_FONT, size=(PT_BASE-0.5)/.pt, colour=INK,
           label="Across a row (same ICI, GEP differs):") +
  annotate("text", x=LX+0.12, y=1.25, hjust=0, family=PUB_FONT, size=(PT_BASE-0.5)/.pt, colour=SAND_DEEP, fontface=2,
           label="ICI-low barely moves  18% → 21%") +
  annotate("text", x=LX+0.12, y=1.08, hjust=0, family=PUB_FONT, size=(PT_BASE-1)/.pt, colour=MUTED,
           label="→ GEP adds little once ICI is known") +
  annotate("text", x=LX, y=0.62, hjust=0, lineheight=0.95, family=PUB_FONT, size=(PT_BASE-1)/.pt, colour=INK,
           label="Darker = higher response. The grid darkens\ngoing UP (ICI) more than going RIGHT (GEP).")

save_fig(p, "ICI_Fig_addedvalue_grid", width_mm=174, height_mm=104, dir="../results/figures/composite")
cat(sprintf("done — 2x2 grid. cells: %s (base %.0f%%)\n",
    paste(sprintf("%s/%s=%.0f%%", cell$GEP, cell$ICI, cell$pct), collapse=", "), base))
