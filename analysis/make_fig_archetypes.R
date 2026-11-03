# =============================================================================
# Three-stage alluvial: CANCER TYPE -> ARCHETYPE -> RESPONSE (validated cohorts
# pooled by cancer type). One read shows which cancers feed which archetypes and
# where those archetypes end up. Ribbons coloured by archetype.
# Archetype classified STRICTLY from per-patient layer scores (classify()):
#   Inflamed = S+ I+ E+ (all layers at/above the cohort median, per Fig 2A); else the
#   lowest sub-median layer names the failure (cold=S, blocked=I, no-kill=E).
# NB these are ABSOLUTE patient counts (base ICB response rate ~29%): read ribbon
# widths as counts, not rates -- per-archetype response RATE is in the supp table.
# Needs per-sample layer scores in pancancer_validation.rds. Editable SVG.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(ggalluvial) })
source("ici_theme.R"); .pt <- 2.845276
R <- readRDS("pancancer_validation.rds")

WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
CTLA4 <- c("VanAllen2015"); NULLB <- c("GSE165252","PRINCE","Braun")
ORDER <- c(WORKS, CTLA4, NULLB)
bucket <- setNames(c(rep("works",length(WORKS)),rep("ctla4",length(CTLA4)),rep("null",length(NULLB))), ORDER)
CANC <- c(melanoma="Melanoma", urothelial="Urothelial", NSCLC="NSCLC", gastric="Gastric",
          breast="Breast", EAC="Esophageal", pancreatic="Pancreatic", ccRCC="Renal")
cdisp <- function(x){ v<-unname(CANC[x]); ifelse(is.na(v), x, v) }

## ---- archetype levels, labels, colours (failing layer's S/I/E colour) ------
ARCH <- c("intact","no-kill","blocked","cold")
ALAB <- c(intact="Inflamed", "no-kill"="Execution-failed", blocked="Excluded", cold="Desert / cold")
ACOL <- c(intact="#5E9C86", "no-kill"=ROSE, blocked=PLUM, cold=SAND_DEEP)
afac <- function(v) factor(ALAB[v], levels=unname(ALAB[ARCH]))
classify <- function(s){ z<-function(v)(v-mean(v))/sd(v)
  M <- cbind(cold=z(s$sens_cohort), blocked=z(s$integ_cohort), "no-kill"=z(s$exec_cohort))
  ifelse(M[,1]>=0 & M[,2]>=0 & M[,3]>=0, "intact", colnames(M)[max.col(-M, ties.method="first")]) }

## ---- per-patient Cancer / Archetype / Response (validated cohorts) ----------
dd <- do.call(rbind, lapply(WORKS, function(k){ s<-R[[k]]$scores
  data.frame(Cancer=cdisp(R[[k]]$cancer), arch=classify(s), y=s$y, stringsAsFactors=FALSE) }))
base <- mean(dd$y)
ord_c <- names(sort(tapply(dd$y, dd$Cancer, mean)))                 # ascending response rate
dd$Cancer    <- factor(dd$Cancer, levels=ord_c)                    # best-responding cancer ends on top
dd$Archetype <- afac(dd$arch)
dd$Response  <- factor(ifelse(dd$y==1,"Responder","Non-responder"), c("Responder","Non-responder"))
alc <- as.data.frame(table(Cancer=dd$Cancer, Archetype=dd$Archetype, Response=dd$Response)); names(alc)[4]<-"n"
alc <- alc[alc$n>0, ]

## ---- three-stage alluvial --------------------------------------------------
pal <- setNames(unname(ACOL[ARCH]), unname(ALAB[ARCH]))
p <- ggplot(alc, aes(y=n, axis1=Cancer, axis2=Archetype, axis3=Response)) +
  geom_alluvium(aes(fill=Archetype), width=0.24, alpha=0.78, colour=NA) +
  geom_stratum(width=0.24, fill="grey96", colour="grey45", linewidth=0.3) +
  ## vertical labels where the stratum is tall enough (archetype/response, and the
  ## larger cancer strata); small cancer strata (NSCLC, Gastric) stay horizontal,
  ## centred in their box. All labels centred (hjust 0.5).
  geom_text(stat="stratum", aes(label=after_stat(stratum),
                                angle=after_stat(ifelse(x > 1 | count > 55, 90, 0))),
            family=PUB_FONT, size=(PT_BASE-0.5)/.pt, colour=INK) +
  scale_fill_manual(values=pal, name=NULL) +
  scale_x_continuous(breaks=1:3, labels=c("Cancer type","Archetype","ICB response"),
                     expand=expansion(add=c(0.24,0.24)), position="top") +
  coord_cartesian(clip="off") +
  guides(fill=guide_legend(nrow=1, override.aes=list(colour="white"))) +
  theme_void(base_family=PUB_FONT) +
  theme(axis.text.x.top=element_text(size=PT_AXIS_TITLE, colour=INK, face=2, margin=margin(b=2,t=1)),
        legend.position="bottom", legend.text=element_text(size=PT_LEGEND, colour=INK),
        legend.margin=margin(t=0), plot.margin=margin(2,5,2,5,"mm"),
        plot.background=element_rect(fill="white", colour=NA))

dir.create("../results/figures/composite", showWarnings=FALSE, recursive=TRUE)
w <- FIG_W_DOUBLE*MM_TO_IN; h <- 158*MM_TO_IN
ggsave("../results/figures/composite/ICI_Fig_archetypes.png", p, width=w, height=h, dpi=600,
       device=if(requireNamespace("ragg",quietly=TRUE)) ragg::agg_png else "png")
ggsave("../results/figures/composite/ICI_Fig_archetypes.svg", p, width=w, height=h,
       device=if(requireNamespace("svglite",quietly=TRUE)) svglite::svglite else "svg")
if(file.exists("svg_join_text.py")) system("python3 svg_join_text.py ../results/figures/composite/ICI_Fig_archetypes.svg", ignore.stdout=TRUE, ignore.stderr=TRUE)

## ---- supplementary table: archetype composition + response rate per cohort --
comp <- do.call(rbind, lapply(ORDER, function(k){ s<-R[[k]]$scores; a<-classify(s)
  tb <- table(factor(a, ARCH))
  data.frame(cohort=R[[k]]$cohort, cancer=R[[k]]$cancer, bucket=bucket[k], n=sum(tb),
             Inflamed=as.integer(tb["intact"]), Execution_failed=as.integer(tb["no-kill"]),
             Excluded=as.integer(tb["blocked"]), Desert=as.integer(tb["cold"]),
             pct_Inflamed=round(100*tb["intact"]/sum(tb),1),
             pctR_Inflamed=round(100*mean(s$y[a=="intact"]),1), row.names=NULL) }))
write.csv(comp, "../results/tables/ICI_SuppTable_archetype_composition.csv", row.names=FALSE)
cat(sprintf("done — 3-stage alluvial (base ICB response %.0f%%) + supp table\n", 100*base))
