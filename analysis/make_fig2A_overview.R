# =============================================================================
# Figure 2A - how the three scores are DERIVED and APPLIED (end to end).
#   Top    : score each tumor on its own against the locked reference medians
#   Middle : the three scores, each on its own range, with the meaning of extremes
#   Bottom : application -> (1) PREDICT ICB response  (2) LOCALIZE the bottleneck
# House style: Arial, warm layer palette, geometric square connectors.
# RULE: all text BLACK (INK), no gray; no em-dashes.
# -----------------------------------------------------------------------------
# TUNABLE CONFIG
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R")
.pt <- 2.845276
col <- pal_layer

EX <- list(S = 0.75, I = 0.34, E = 0.74)              # example tumor scores

G <- list(
  list(sym="S", name="Sensory",     col=unname(col["Sense"]),
       lo=0,  hi=1,  val=EX$S, mid=NA, valfmt="%.2f",
       lo_lab="0", hi_lab="1", mid_lab=NA,
       lo_mean="not detected", mid_mean=NA,         hi_mean="fully detected"),
  list(sym="I", name="Integration", col=unname(col["Integration"]),
       lo=-1, hi=1,  val=EX$I, mid=0,  valfmt="%+.2f",
       lo_lab="−1", hi_lab="+1", mid_lab="0",
       lo_mean="suppressed",   mid_mean="balanced", hi_mean="open access"),
  list(sym="E", name="Execution",   col=unname(col["Execution"]),
       lo=-1, hi=1,  val=EX$E, mid=0,  valfmt="%+.2f",
       lo_lab="−1", hi_lab="+1", mid_lab="0",
       lo_mean="exhausted",    mid_mean="balanced", hi_mean="cytotoxic"))
GY  <- c(91, 80, 69); TX0 <- 48; TX1 <- 176

# Localize rows: the lowest score names the failing layer (which layer, not the therapy).
LOC <- list(
  list(sym="S", col=unname(col["Sense"]),       why="tumor not detected"),
  list(sym="I", col=unname(col["Integration"]), why="immune access blocked"),
  list(sym="E", col=unname(col["Execution"]),   why="no tumor-cell killing"))

# Archetypes: canonical S/I/E high(1)/low(0) patterns -> named phenotype (name only; the
# ICB expectation per archetype is a result and lives in a separate figure).
# Four data-grounded archetypes (IMvigor210 check, 2026-10-31): Sensing Failed (0,1,1)
# is only 2.9% and all boundary cases -> dropped. Excluded = tumour detected (S+) but
# immune access blocked (I-), so effectors never reach the tumour and cannot kill (E-):
# S+ I- E- (1,0,0). Ordered STEP-WISE down the cascade (S->I->E): the failure point moves
# one layer earlier each row, so the + cells form a descending staircase.
ARCH <- list(
  list(pat=c(1,1,1), name="Inflamed"),          # cascade completes
  list(pat=c(1,1,0), name="Execution Failed"),  # fails at E
  list(pat=c(1,0,0), name="Excluded"),           # fails at I
  list(pat=c(0,0,0), name="Desert"))             # fails at S
AY  <- c(30, 23.5, 17, 10.5)          # archetype row y-centers (inside card 3)
ASX <- c(10, 16, 22)                  # x of the S/I/E signature-code columns (archetype card = leftmost)
ASCOL <- unname(col[c("Sense","Integration","Execution")])

# =============================================================================
FSCALE <- 1.16   # global font bump (legibility prioritized)
txt <- function(x,y,lab,size=9,face=1,colour=INK,hjust=0.5,parse=FALSE)
  annotate("text", x=x, y=y, label=lab, family=PUB_FONT, size=size*FSCALE/.pt, fontface=face, colour=colour, hjust=hjust, parse=parse)
rect <- function(xmin,xmax,ymin,ymax,fill,colour,lwd=0.6)
  annotate("rect", xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=fill,colour=colour,linewidth=lwd)
harrow <- function(x,xend,y) annotate("segment", x=x,xend=xend,y=y,yend=y,
            colour=INK, linewidth=0.7, arrow=arrow(length=unit(1.8,"mm"), type="closed"))
xmap <- function(v,g) TX0 + (v-g$lo)/(g$hi-g$lo)*(TX1-TX0)
# simple patient pictogram (head + torso) centered at (cx,cy)
person <- function(cx, cy, fill) list(
  geom_polygon(data=data.frame(x=c(cx-2.6,cx-2.1,cx-1.0,cx+1.0,cx+2.1,cx+2.6),
                               y=c(cy-3.2,cy+0.6,cy+1.6,cy+1.6,cy+0.6,cy-3.2)),
               aes(x,y), fill=fill, colour=fill, linewidth=0),
  annotate("point", x=cx, y=cy+3.0, size=3.2, colour=fill))

gauge_layers <- list()
for (i in seq_along(G)) { g <- G[[i]]; y <- GY[i]
  gauge_layers <- c(gauge_layers, list(
    rect(8, 14, y-2.6, y+2.6, g$col, g$col),
    txt(11, y, g$sym, size=9, face=2, colour="white"),
    txt(17, y, g$name, size=8.8, face=2, hjust=0, colour=g$col),
    annotate("segment", x=TX0, xend=TX1, y=y, yend=y, colour=g$col, linewidth=0.8),
    annotate("segment", x=TX0, xend=TX0, y=y-1.5, yend=y+1.5, colour=g$col, linewidth=0.8),
    annotate("segment", x=TX1, xend=TX1, y=y-1.5, yend=y+1.5, colour=g$col, linewidth=0.8),
    txt(TX0, y-3.3, g$lo_lab, size=7.2, colour=g$col),
    txt(TX1, y-3.3, g$hi_lab, size=7.2, colour=g$col),
    txt(TX0, y-6.5, g$lo_mean, size=7.4, face=3, hjust=0, colour=INK),
    txt(TX1, y-6.5, g$hi_mean, size=7.4, face=3, hjust=1, colour=INK)))
  if (!is.na(g$mid)) gauge_layers <- c(gauge_layers, list(
    annotate("segment", x=xmap(g$mid,g), xend=xmap(g$mid,g), y=y-1.5, yend=y+1.5, colour=g$col, linewidth=0.4),
    txt(xmap(g$mid,g), y-3.3, g$mid_lab, size=7.4, colour=g$col),
    txt(xmap(g$mid,g), y-6.5, g$mid_mean, size=7.4, face=3, colour=INK))) }

# localize rows (card 3, now rightmost): layer tag + low + which layer fails (one column)
loc_layers <- list(); ly <- c(30, 23, 16)
for (i in seq_along(LOC)) { l <- LOC[[i]]; y <- ly[i]
  loc_layers <- c(loc_layers, list(
    rect(131, 135.5, y-2.4, y+2.4, pal_layer_fill[c("Sense","Integration","Execution")[i]], l$col, lwd=0.5),
    txt(133.25, y, l$sym, size=8.3, face=2, colour=l$col),
    txt(137, y, "low", size=7.8, hjust=0, colour=l$col),
    txt(145, y, l$why, size=7.5, hjust=0, colour=INK))) }

# archetype rows (card 1, leftmost): layer box (filled=high / open=low) with its +/- sign inside + name
arch_layers <- list()
for (i in seq_along(ARCH)) { a <- ARCH[[i]]; y <- AY[i]
  for (j in 1:3) { hi <- a$pat[j]==1
    arch_layers <- c(arch_layers, list(
      rect(ASX[j]-1.9, ASX[j]+1.9, y-1.9, y+1.9, if (hi) ASCOL[j] else "white", ASCOL[j], lwd=0.6),
      txt(ASX[j], y, if (hi) "+" else "−", size=6.5, face=2, colour=if (hi) "white" else ASCOL[j]))) }
  arch_layers <- c(arch_layers, list(
    txt(30, y, a$name, size=8, face=2, hjust=0, colour=INK))) }

pA <- ggplot() +
  annotate("text", x=1, y=129, label="A", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +

  # ===================== DERIVE =====================
  txt(93, 126, "Score each patient at each layer", size=10.5, face=2) +
  rect(5, 37, 108, 122, CREAM, MUTED) +
  txt(21, 117.5, "patient cohort", size=8.5, face=2) + txt(21, 113, "transcriptomes", size=8) +
  rect(50, 120, 107, 123, CREAM, MUTED) +
  txt(85, 119, "binarize each gene at its", size=8) +
  txt(85, 114.8, "cohort median", size=8, face=2) +
  txt(85, 110.5, "b[g]==1~~'if'~~x[g] >= m[g]", size=8.2, parse=TRUE) +
  rect(133, 181, 108, 122, CREAM, MUTED) +
  txt(157, 118, "combine per layer", size=8) + txt(157, 113.5, "(÷ gene count)", size=8) +
  harrow(37, 50, 115) + harrow(120, 133, 115) +
  annotate("segment", x=157, xend=157, y=108, yend=102, colour=INK, linewidth=0.7,
           arrow=arrow(length=unit(1.8,"mm"), type="closed")) +

  # ===================== SCORES + MEANING =====================
  txt(8, 101, "Three scores per patient, each on its own range:", size=9.5, face=2, hjust=0) +
  gauge_layers +

  # ===================== APPLICATION (three cards on one row) =====================
  txt(93, 58, "Application of the three scores", size=10.5, face=2) +

  # ---- Card 1 : ARCHETYPE (classify the tumor from the S/I/E pattern; step 1) ----
  rect(5, 63, 4, 52, "white", INK, lwd=0.9) +
  annotate("point", x=11, y=48, size=5.2, colour=INK) + txt(11, 48, "1", size=8, face=2, colour="white") +
  txt(16.5, 48, "Archetypes", size=9, face=2, hjust=0, colour=INK) +
  txt(8, 42, "Each +/− pattern is a type:", size=7.3, hjust=0, colour=INK) +
  txt(ASX[1], 36, "S", size=7.4, face=2, colour=ASCOL[1]) +
  txt(ASX[2], 36, "I", size=7.4, face=2, colour=ASCOL[2]) +
  txt(ASX[3], 36, "E", size=7.4, face=2, colour=ASCOL[3]) +
  txt(30, 36, "Archetype", size=7.6, face=2, hjust=0, colour=INK) +
  annotate("segment", x=8, xend=58, y=34, yend=34, colour=INK, linewidth=0.3) +
  arch_layers +

  # ---- Card 2 : STRATIFY benefit (layer-agnostic; the test decides which layer; step 2) ----
  rect(66, 124, 4, 52, "white", INK, lwd=0.9) +
  annotate("point", x=72, y=48, size=5.2, colour=INK) + txt(72, 48, "2", size=8, face=2, colour="white") +
  txt(77.5, 48, "Stratify ICB benefit", size=9, face=2, hjust=0, colour=INK) +
  txt(69, 42, "from the layer scores:", size=7.6, hjust=0, colour=INK) +
  person(70.5, 30, MUTED) + person(76.5, 30, MUTED) +
  txt(73.5, 23.5, "patients", size=7, hjust=0.5, colour=INK) +
  annotate("segment", x=80, xend=83.5, y=30, yend=30, colour=INK, linewidth=0.6, arrow=arrow(length=unit(1.3,"mm"), type="closed")) +
  annotate("segment", x=88, xend=88, y=15, yend=38, colour=INK, linewidth=0.6) +
  annotate("segment", x=88, xend=119, y=15, yend=15, colour=INK, linewidth=0.6) +
  annotate("segment", x=90, xend=117, y=18, yend=36, colour=INK, linewidth=1.1) +
  annotate("segment", x=90, xend=117, y=18, yend=18, colour=NEUTRAL, linewidth=1.0, linetype="22") +
  txt(110, 38, "ICB", size=6.3, face=2, colour=INK) +
  txt(111, 20, "control", size=6.3, hjust=0.5, colour=INK) +
  annotate("text", x=85.3, y=26.5, label="response", angle=90,
           family=PUB_FONT, size=6.8*FSCALE/.pt, colour=INK) +
  txt(90, 12.5, "low", size=6.8, colour=INK) + txt(117, 12.5, "high", size=6.8, colour=INK) +
  txt(103.5, 9.3, "layer score", size=7.4, face=2, hjust=0.5, colour=INK) +
  txt(69, 6.3, "layer × treatment model.", size=6.8, face=3, hjust=0, colour=INK) +

  # ---- Card 3 : LOCALIZE the bottleneck (for non-responders: which layer limits; step 3) ----
  rect(127, 182, 4, 52, "white", INK, lwd=0.9) +
  annotate("point", x=133, y=48, size=5.2, colour=INK) + txt(133, 48, "3", size=8, face=2, colour="white") +
  txt(138.5, 48, "Localize the bottleneck", size=8.3, face=2, hjust=0, colour=INK) +
  txt(130, 42, "Lowest score names the layer:", size=7.6, hjust=0, colour=INK) +
  loc_layers +

  coord_fixed(ratio=1, xlim=c(0,186), ylim=c(1,131), expand=FALSE, clip="off") +
  theme_void() + theme(plot.background=element_rect(fill="white", colour=NA),
                       plot.margin=margin(1,1,1,1,"mm"))

if (!exists("FIG2_NOSAVE")) save_fig(pA, "ICI_Fig2A_overview", width_mm=174, height_mm=174*130/186)
cat("done\n")
