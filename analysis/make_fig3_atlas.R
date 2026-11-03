# =============================================================================
# Figure 3 (anatomy pass): two SMART human figures (male + female) with the
# organs that can develop cancer placed at their anatomical sites, each tagged
# with a number.  Shared organs carry the SAME number on both figures; sex-
# specific organs have their own numbers.  A key decodes every number.
# Male urogenital + female reproductive systems shown INTERNALLY (SMART
# apparatus art); breast shown forward-facing on the anterior chest.
# No scores yet (added later); organs scored in this study are highlighted.
# All art: Servier Medical Art (smart.servier.com), CC BY 4.0.  House style ici_theme.R.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(grid); library(png) })
source("ici_theme.R"); .pt <- 2.845276
ICON_DIR <- "assets/icons"
ICO   <- function(f) rasterGrob(readPNG(file.path(ICON_DIR, f)), interpolate=TRUE)
.dimc <- new.env()
PDIM  <- function(f){ if(is.null(.dimc[[f]])){ a<-readPNG(file.path(ICON_DIR,f)); .dimc[[f]]<-c(w=dim(a)[2],h=dim(a)[1]) }; .dimc[[f]] }
imgwh <- function(f, cx, cy, w, h) annotation_custom(ICO(f), xmin=cx-w/2, xmax=cx+w/2, ymin=cy-h/2, ymax=cy+h/2)
imgw  <- function(f, cx, cy, w){ d<-PDIM(f); h<-w*d["h"]/d["w"]; imgwh(f,cx,cy,w,unname(h)) }
disc  <- function(cx, cy, r, n=56){ a<-seq(0,2*pi,length.out=n); data.frame(x=cx+r*cos(a), y=cy+r*sin(a)) }

## ---- figure geometry (both silhouettes share the vertical frame) ----
YTOP <- 158; HH <- 128
PY   <- function(f) YTOP - f*HH
MALE <- list(file="male_silhouette.png",   cx=50,  w=HH*673/1207)
FEM  <- list(file="female_silhouette.png", cx=150, w=HH*603/1208)
BADGE_DATA <- PLUM; BADGE_OTH <- "#6E6257"

## master legend (num, name, scope, scored-here) -------------------------------
# status per site: "val" = validated (ICB response predicted), "null" = tested but
# null (esophageal/pancreatic/renal), "none" = not studied here.
LEG <- list(
  c(1,"Brain","both","none"),    c(2,"Thyroid","both","none"), c(3,"Lung","both","val"),
  c(4,"Esophagus","both","null"),c(5,"Gastric","both","val"),  c(6,"Liver","both","none"),
  c(7,"Pancreas","both","null"), c(8,"Kidney","both","null"),  c(9,"Colon","both","none"),
  c(10,"Bladder","both","val"),  c(11,"Skin (melanoma)","both","val"), c(12,"Breast","female","val"),
  c(13,"Ovary","female","none"), c(14,"Uterus","female","none"), c(15,"Prostate","male","none"), c(16,"Testis","male","none"))
STAT <- setNames(vapply(LEG, `[`, "", 4), vapply(LEG, `[`, "", 1))   # num -> status

## collect badges, draw on top at the end -------------------------------------
BADGES <- list()
addbadge <- function(x,y,num,data) BADGES[[length(BADGES)+1]] <<- list(x=x,y=y,num=num,data=data)

## ---- shared viscera (placed identically on both figures), back -> front ----
O <- function(num,file,f,dx,w,data=FALSE,bdx=0,bdy=0) list(num=num,file=file,f=f,dx=dx,w=w,data=data,bdx=bdx,bdy=bdy)
shared <- list(
  O(4 ,"esophagus_smart.png",0.42, 5.5, 2.4,FALSE, bdx=3.4,bdy=5),
  O(7 ,"pancreas_smart.png" ,0.615,-1.0,15.0,FALSE, bdx=-2, bdy=-4),
  O(3 ,"lungs_smart.png"    ,0.41, 0.0,18.0,TRUE , bdx=-5, bdy=1),
  O(6 ,"liver_smart.png"    ,0.55,-8.0,16.0,FALSE, bdx=-5, bdy=2),
  O(5 ,"stomach_smart.png"  ,0.55, 7.0,13.0,FALSE, bdx=4,  bdy=3),
  O(9 ,"colon_smart.png"    ,0.72, 0.0,18.0,FALSE, bdx=-7, bdy=1),
  O(1 ,"brain_smart.png"    ,0.10,-1.0,11.0,FALSE, bdx=5,  bdy=-2),
  O(2 ,"thyroid_smart.png"  ,0.255,0.0, 6.0,FALSE, bdx=0,  bdy=0),
  O(11,"melanoma_lesion.png",0.35,16.0, 8.0,TRUE , bdx=4,  bdy=5))

p <- ggplot() +
  imgwh(MALE$file, MALE$cx, YTOP-HH/2, MALE$w, HH) +
  imgwh(FEM$file,  FEM$cx,  YTOP-HH/2, FEM$w,  HH)

for (fig in list(MALE, FEM)) for (g in shared){
  p <- p + imgw(g$file, fig$cx+g$dx, PY(g$f), g$w)
  addbadge(fig$cx+g$dx+g$bdx, PY(g$f)+g$bdy, g$num, g$data)
}

## ---- female breast (forward-facing): external areola (right) + gland ductal tree (left) ----
BR_SKIN<-"#F1CDC6"; BR_LINE<-"#C98F93"; DUCT<-"#BE7490"; NIP<-"#9B5E66"; AREO<-"#D7A2A6"
LOB <- cbind(c(0.6,-0.5,0.15), c(0.15,0.55,-0.5))         # fixed lobule offsets (deterministic)
draw_breast <- function(bx,by,r,gland){
  d<-disc(bx,by,r)
  L<-list(geom_polygon(data=d, aes(x,y), fill=BR_SKIN, colour=BR_LINE, linewidth=0.3, alpha=0.9))
  if(gland){
    for(a in seq(0,2*pi,length.out=9)[-9]){ ex<-bx+r*0.68*cos(a); ey<-by+r*0.68*sin(a)
      L<-c(L, list(annotate("segment", x=bx,y=by, xend=ex,yend=ey, linewidth=0.4, colour=DUCT)))
      for(j in 1:3) L<-c(L, list(annotate("point", x=ex+LOB[j,1], y=ey+LOB[j,2], size=0.7, colour=DUCT))) }
    L<-c(L, list(annotate("point", x=bx,y=by, size=1.5, colour=NIP)))
  } else {
    L<-c(L, list(annotate("point", x=bx,y=by, size=3.3, colour=AREO),
                 annotate("point", x=bx,y=by, size=1.2, colour=NIP)))
  }
  L
}
p <- p + draw_breast(FEM$cx+9, PY(0.47), 7, FALSE) +     # right: external (forward-facing)
         draw_breast(FEM$cx-9, PY(0.47), 7, TRUE)        # left: gland (forward-facing ducts)
addbadge(FEM$cx+13, PY(0.41), 12, TRUE)

## ---- male urinary+genital: kidneys-ureters-bladder-prostate (one tract) + testis ----
p <- p + imgw("urinary_male_np.png", MALE$cx, 63.5, 26) +
         imgw("testis_smart.png",    MALE$cx+1, 39, 7)
addbadge(MALE$cx-12, 80.0, 8,  FALSE)  # kidney
addbadge(MALE$cx,    52.0, 10, TRUE)   # bladder
addbadge(MALE$cx+1,  45.0, 15, FALSE)  # prostate
addbadge(MALE$cx+1,  39.0, 16, FALSE)  # testis

## ---- female urinary (kidneys-ureters-bladder) + reproductive (uterus/ovaries) ----
p <- p + imgw("urinary_female.png", FEM$cx, 64.0, 24) +
         imgw("femalerepro.png",    FEM$cx, 57.0, 15)
addbadge(FEM$cx-11, 80.0, 8,  FALSE)   # kidney
addbadge(FEM$cx,    52.0, 10, TRUE)    # bladder
addbadge(FEM$cx+6,  61.0, 13, FALSE)   # ovary
addbadge(FEM$cx-1,  59.0, 14, FALSE)   # uterus

## ---- draw all badges on top ----
badge <- function(x,y,num){ s <- STAT[[as.character(num)]]
  if (identical(s,"null"))        # hollow plum = tested but null
    list(annotate("point", x=x, y=y, size=5.0, shape=21, fill="white", colour=PLUM, stroke=0.9),
         annotate("text",  x=x, y=y, label=num, family=PUB_FONT, size=6/.pt, fontface=2, colour=PLUM))
  else                            # filled: plum = validated, grey = not studied
    list(annotate("point", x=x, y=y, size=5.0, colour=if(identical(s,"val")) BADGE_DATA else BADGE_OTH),
         annotate("point", x=x, y=y, size=5.0, shape=1, colour="white", stroke=0.5),
         annotate("text",  x=x, y=y, label=num, family=PUB_FONT, size=6/.pt, fontface=2, colour="white")) }
for (b in BADGES) p <- p + badge(b$x, b$y, b$num)

## ---- titles, sex labels ----
p <- p +
  annotate("text", x=2, y=163, label="A", family=PUB_FONT, size=PT_PANEL_TAG/.pt, fontface=2, hjust=0) +
  annotate("text", x=MALE$cx, y=23, label="Male",   family=PUB_FONT, size=8/.pt, fontface=2, hjust=0.5, colour=INK) +
  annotate("text", x=FEM$cx,  y=23, label="Female", family=PUB_FONT, size=8/.pt, fontface=2, hjust=0.5, colour=INK)

## ---- key / legend (numbered), right column ----
lx<-190; ly0<-152; dy<-7.2
p <- p + annotate("text", x=lx-2, y=ly0+9, label="Key", family=PUB_FONT, size=11.5/.pt, fontface=2, hjust=0, colour=INK)
for (i in seq_along(LEG)){ e<-LEG[[i]]; num<-as.integer(e[1]); nm<-e[2]; scope<-e[3]; st<-e[4]; y<-ly0-(i-1)*dy
  suffix <- if(scope=="female") "  (female)" else if(scope=="male") "  (male)" else ""
  p <- p + badge(lx, y, num) +
    annotate("text", x=lx+6, y=y, label=paste0(nm, suffix),
             family=PUB_FONT, size=9.5/.pt, fontface=if(st=="val") 2 else 1, hjust=0, colour=INK) }

## ---- status legend (below the key, in existing whitespace; no layout shift) ----
lgy <- ly0 - length(LEG)*dy - 1.5
p <- p +
  annotate("point", x=lx, y=lgy,      size=5.0, colour=PLUM) +
  annotate("point", x=lx, y=lgy,      size=5.0, shape=1, colour="white", stroke=0.5) +
  annotate("text",  x=lx+6, y=lgy,      label="validated (response predicted)", family=PUB_FONT, size=8/.pt, hjust=0, colour=INK) +
  annotate("point", x=lx, y=lgy-dy,   size=5.0, shape=21, fill="white", colour=PLUM, stroke=0.9) +
  annotate("text",  x=lx+6, y=lgy-dy,   label="tested, null (no prediction)",   family=PUB_FONT, size=8/.pt, hjust=0, colour=INK) +
  annotate("point", x=lx, y=lgy-2*dy, size=5.0, colour=BADGE_OTH) +
  annotate("point", x=lx, y=lgy-2*dy, size=5.0, shape=1, colour="white", stroke=0.5) +
  annotate("text",  x=lx+6, y=lgy-2*dy, label="not studied here",                family=PUB_FONT, size=8/.pt, hjust=0, colour=INK)

## ---- frame ----
p <- p +
  coord_fixed(ratio=1, xlim=c(0,250), ylim=c(0,168), expand=FALSE, clip="off") +
  theme_void() + theme(plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(1,1,1,1,"mm"))

save_fig(p, "ICI_Fig3_atlas", width_mm=174, height_mm=174*168/250, dir="../results/figures/composite")
cat("done\n")
