# =============================================================================
# Archetype (limiting-layer) distribution across cohorts + how archetypes map to
# response.  Panel A: per-cohort 100% stacked bars of the four archetypes.
# Panel B: alluvial of archetype -> responder/non-responder, pooled over the
# validated cohorts (where the cascade predicts response).
# Archetype = icb_eval per-patient limiting layer ($scores$limiting):
#   intact -> Inflamed; no-kill -> Execution-failed (E); blocked -> Excluded (I);
#   cold -> Desert/cold (S).  Colours follow the S/I/E palette of the failing layer.
# Needs per-sample scores+limiting in pancancer_validation.rds. Editable SVG.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(ggalluvial); library(patchwork); library(grid) })
source("ici_theme.R"); .pt <- 2.845276
R <- readRDS("pancancer_validation.rds")

WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
CTLA4 <- c("VanAllen2015"); NULLB <- c("GSE165252","PRINCE","Braun")
AUTHOR <- c(IMvigor210="IMvigor210", GSE78220="Hugo 2016", GSE91061="Riaz 2017",
          GSE126044="Cho 2020", GSE135222="Jung 2019", Gide2019="Gide 2019", Liu2019="Liu 2019",
          Puch="PUCH 2022", Kim2018gas="Kim 2018", Rose="Rose 2021", ISPY2="I-SPY2 trial",
          VanAllen2015="Van Allen 2015", GSE165252="PERFECT trial", PRINCE="PRINCE trial", Braun="Braun 2020")
CANC <- c(melanoma="Melanoma", urothelial="Urothelial", NSCLC="NSCLC", gastric="Gastric",
          breast="Breast", EAC="Esophageal", pancreatic="Pancreatic", ccRCC="Renal")
cdisp <- function(x){ v<-unname(CANC[x]); ifelse(is.na(v), x, v) }

## order cohorts by Integration AUC (desc), as in Fig 3b
WORKS <- WORKS[order(-sapply(WORKS, function(k) R[[k]]$int_ci_single[1]))]
ORDER <- c(WORKS, CTLA4, NULLB)
bucket <- setNames(c(rep("works",length(WORKS)),rep("ctla4",length(CTLA4)),rep("null",length(NULLB))), ORDER)

## ---- archetype levels, labels, colours (failing layer's S/I/E colour) ------
ARCH <- c("intact","no-kill","blocked","cold")                    # stack bottom -> top (cascade order)
ALAB <- c(intact="Inflamed", "no-kill"="Execution-failed", blocked="Excluded", cold="Desert / cold")
ACOL <- c(intact="#5E9C86", "no-kill"=ROSE, blocked=PLUM, cold=SAND_DEEP)
afac <- function(v) factor(ALAB[v], levels=unname(ALAB[ARCH]))

## ===================== Panel A: per-cohort stacked bars =====================
xpos <- numeric(0); xc <- 1; prev <- "works"
for(k in ORDER){ if(bucket[k]!=prev){ xc <- xc+1.0; prev <- bucket[k] }; xpos[k]<-xc; xc<-xc+1 }

## cumulative proportions -> rectangles
segs <- do.call(rbind, lapply(ORDER, function(k){ s<-R[[k]]$scores
  tb <- table(factor(s$limiting, ARCH)); pr <- as.numeric(tb)/sum(tb)
  hi <- cumsum(pr); lo <- hi-pr
  data.frame(cohort=k, x=xpos[k], arch=ARCH, lo=lo, hi=hi, bucket=bucket[k]) }))
segs$archf <- afac(segs$arch)
meta <- data.frame(k=ORDER, x=xpos[ORDER],
  label=sapply(ORDER, function(k) sprintf("%s · %s (n=%d)", cdisp(R[[k]]$cancer), AUTHOR[[k]], R[[k]]$n)))

bw <- 0.42; yb <- c(-0.02, 1.02)
bandA <- function(ks, fill){ xs<-xpos[ks]; annotate("rect", xmin=min(xs)-0.55, xmax=max(xs)+0.55,
                 ymin=-0.45, ymax=1.18, fill=fill, colour=NA) }
pA <- ggplot() +
  bandA(WORKS,"#F7F2F5") + bandA(CTLA4,"#F8F3E6") + bandA(NULLB,"#F1EEE9") +
  annotate("text", x=mean(range(xpos[WORKS])), y=1.14, label="Validated (anti-PD-1/PD-L1)",
           family=PUB_FONT, size=(PT_BASE-1)/.pt, fontface=2, colour=INK) +
  annotate("text", x=xpos[CTLA4], y=1.14, label="CTLA-4", family=PUB_FONT, size=(PT_BASE-1)/.pt, fontface=2, colour=INK) +
  annotate("text", x=mean(range(xpos[NULLB])), y=1.14, label="Null", family=PUB_FONT, size=(PT_BASE-1)/.pt, fontface=2, colour=INK) +
  geom_rect(data=segs, aes(xmin=x-bw, xmax=x+bw, ymin=lo, ymax=hi, fill=archf), colour="white", linewidth=0.25) +
  scale_fill_manual(values=setNames(unname(ACOL[ARCH]), unname(ALAB[ARCH])), name=NULL) +
  geom_text(data=meta, aes(x=x, y=-0.03, label=label), angle=90, hjust=1, vjust=0.5,
            family=PUB_FONT, size=(PT_BASE-2)/.pt, colour=INK) +
  annotate("segment", x=-0.0, xend=0.0, y=0, yend=1, linewidth=0.3, colour=RULE) +
  lapply(c(0,.25,.5,.75,1), function(t) list(
    annotate("segment", x=-0.08, xend=0, y=t, yend=t, linewidth=0.3, colour=RULE),
    annotate("text", x=-0.2, y=t, label=paste0(t*100), family=PUB_FONT, size=(PT_BASE-2.5)/.pt, colour=INK, hjust=1))) +
  annotate("text", x=-1.0, y=0.5, label="Patients (%)", angle=90, family=PUB_FONT, size=PT_AXIS_TITLE/.pt, colour=INK) +
  coord_cartesian(xlim=c(-1.1, max(xpos)+0.6), ylim=c(-3.2, 1.4), clip="off") +
  theme_void() + theme(legend.position="none",
        plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(2,2,2,2,"mm"))

## ===================== Panel B: archetype -> response alluvial ==============
al <- do.call(rbind, lapply(WORKS, function(k){ s<-R[[k]]$scores; s$cohort<-k; s }))
al$Archetype <- afac(al$limiting)
al$Response  <- factor(ifelse(al$y==1,"Responder","Non-responder"), c("Responder","Non-responder"))
alc <- as.data.frame(table(Archetype=al$Archetype, Response=al$Response)); names(alc)[3]<-"n"
alc <- alc[alc$n>0, ]
RESCOL <- c(Responder=PLUM_DEEP, "Non-responder"="#C9BEB0")

## axis1=Response, axis2=Archetype so after coord_flip the ARCHETYPE row is on top
pB <- ggplot(alc, aes(y=n, axis1=Response, axis2=Archetype)) +
  geom_alluvium(aes(fill=Archetype), width=0.26, alpha=0.75, colour=NA) +
  geom_stratum(width=0.26, colour="grey45", fill="grey96", linewidth=0.3) +
  geom_text(stat="stratum", aes(label=after_stat(stratum)), family=PUB_FONT, size=(PT_BASE-3)/.pt) +
  scale_fill_manual(values=setNames(unname(ACOL[ARCH]), unname(ALAB[ARCH])), guide="none") +
  scale_x_continuous(breaks=1:2, labels=c("Response","Archetype"), expand=expansion(add=c(0.22,0.22))) +
  coord_flip(clip="off") +
  labs(subtitle="Archetype → response, validated cohorts (pooled)") +
  theme_void(base_family=PUB_FONT) +
  theme(plot.subtitle=element_text(size=PT_ANNOT, colour=MUTED, hjust=0.5),
        plot.background=element_rect(fill="white", colour=NA), plot.margin=margin(1,3,2,3,"mm"))

## shared archetype legend (built from pA's scale), placed under the combo
leg <- cowplot::get_legend(pA + theme(legend.position="bottom") +
         guides(fill=guide_legend(nrow=1, override.aes=list(colour="white"))))

combo <- (pA / pB) + plot_layout(heights=c(1.5, 1)) +
  plot_annotation(tag_levels="A") &
  theme(plot.tag=element_text(family=PUB_FONT, size=PT_PANEL_TAG, face=2, colour=INK))
final <- cowplot::plot_grid(combo, leg, ncol=1, rel_heights=c(1, 0.05))

save_fig2 <- function(){ dir.create("../results/figures/composite", showWarnings=FALSE, recursive=TRUE)
  w<-FIG_W_DOUBLE*MM_TO_IN; h<-158*MM_TO_IN
  ggsave("../results/figures/composite/ICI_Fig_archetypes.png", final, width=w, height=h, dpi=600,
         device=if(requireNamespace("ragg",quietly=TRUE)) ragg::agg_png else "png")
  ggsave("../results/figures/composite/ICI_Fig_archetypes.svg", final, width=w, height=h,
         device=if(requireNamespace("svglite",quietly=TRUE)) svglite::svglite else "svg")
  if(file.exists("svg_join_text.py")) system("python3 svg_join_text.py ../results/figures/composite/ICI_Fig_archetypes.svg", ignore.stdout=TRUE, ignore.stderr=TRUE)
}
save_fig2()

## ---- supplementary table: archetype composition per cohort -----------------
comp <- do.call(rbind, lapply(ORDER, function(k){ s<-R[[k]]$scores
  tb <- table(factor(s$limiting, ARCH))
  data.frame(cohort=R[[k]]$cohort, cancer=R[[k]]$cancer, bucket=bucket[k], n=sum(tb),
             Inflamed=as.integer(tb["intact"]), Execution_failed=as.integer(tb["no-kill"]),
             Excluded=as.integer(tb["blocked"]), Desert=as.integer(tb["cold"]),
             pct_Inflamed=round(100*tb["intact"]/sum(tb),1), row.names=NULL) }))
write.csv(comp, "../results/tables/ICI_SuppTable_archetype_composition.csv", row.names=FALSE)
cat("done — ICI_Fig_archetypes + ICI_SuppTable_archetype_composition.csv\n")
