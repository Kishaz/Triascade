# =============================================================================
# The ICI is a dissociable, forward-coherent immune cascade (pooled 11 validated
# anti-PD-(L)1 cohorts, single-patient layer scores z-scored within cohort).
#   (a) DISSOCIABILITY: the 3 layers are near-independent axes (mean pairwise
#       R^2 = 9%  ->  ~91% of each layer's variance is independent).
#   (b) CASCADE STAIRCASE: forward-concordant S/I/E states are enriched and
#       forward-breaking (discordant) states depleted vs independence (chi-sq).
#   (c) CASCADE COMPLETION: ICB response rate climbs with the number of intact
#       layers (logistic trend) -- a graded dose-response of cascade integrity.
# Writes results/tables/ICI_cascade_{dissociability,patterns,completion}.csv.
# House style; single canvas; editable SVG.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("ici_theme.R"); .pt <- 2.845276
R <- readRDS("pancancer_validation.rds")
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
OUT <- "../results/tables"; dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
z <- function(v){ s<-sd(v,na.rm=TRUE); if(is.na(s)||s==0) v*0 else (v-mean(v,na.rm=TRUE))/s }
D <- do.call(rbind, lapply(WORKS, function(k){ s<-R[[k]]$scores
  data.frame(coh=k, y=s$y, S=z(s$sens_single), I=z(s$integ_single), X=z(s$exec_single)) }))
hib <- function(v) as.integer(v > median(v))
D$Sb<-ave(D$S,D$coh,FUN=hib); D$Ib<-ave(D$I,D$coh,FUN=hib); D$Xb<-ave(D$X,D$coh,FUN=hib)

## ---- (a) dissociability ----------------------------------------------------
rr <- c(`Sensory vs Integration`=cor(D$S,D$I), `Sensory vs Execution`=cor(D$S,D$X),
        `Integration vs Execution`=cor(D$I,D$X))
indep <- 100*(1 - mean(rr^2))
write.csv(data.frame(pair=names(rr), r=round(rr,3)), file.path(OUT,"ICI_cascade_dissociability.csv"), row.names=FALSE)
PA <- rbind(data.frame(pair=names(rr)[1], x=D$S, y=D$I),
            data.frame(pair=names(rr)[2], x=D$S, y=D$X),
            data.frame(pair=names(rr)[3], x=D$I, y=D$X))
PA$pair <- factor(PA$pair, levels=names(rr))
rlab <- data.frame(pair=factor(names(rr),levels=names(rr)), r=sprintf("r = %.2f", rr))
pa <- ggplot(PA, aes(x,y)) +
  geom_point(colour=PLUM, alpha=0.12, size=0.5, stroke=0) +
  geom_smooth(method="lm", se=FALSE, colour=SAND_DEEP, linewidth=0.6, formula=y~x) +
  geom_text(data=rlab, aes(x=-2.7, y=2.9, label=r), hjust=0, vjust=1, inherit.aes=FALSE,
            family=PUB_FONT, size=(PT_BASE-1)/.pt, fontface=2, colour=INK) +
  facet_wrap(~pair, nrow=1) +
  coord_cartesian(xlim=c(-3,3), ylim=c(-3,3)) +
  scale_x_continuous("layer score (z)", breaks=seq(-2,2,2)) +
  scale_y_continuous("layer score (z)", breaks=seq(-2,2,2)) +
  theme_publication(grid="both") +
  theme(axis.text=element_text(size=PT_BASE-1, colour=INK),
        axis.title=element_text(size=PT_AXIS_TITLE-0.5, colour=INK),
        strip.text=element_text(size=PT_BASE-0.5, face=2, colour=INK),
        panel.spacing=unit(3,"mm"), plot.margin=margin(2,3,1,2,"mm"))

## ---- (b) cascade staircase -------------------------------------------------
valid <- c("000","100","110","111")
pat <- paste0(D$Sb,D$Ib,D$Xb)
lv  <- c("000","100","110","111","010","001","011","101")
tab <- table(factor(pat, levels=lv))
exp_cnt <- setNames(rep(0,8), lv)
for(k in WORKS){ d<-D[D$coh==k,]; ps<-mean(d$Sb); pi<-mean(d$Ib); px<-mean(d$Xb); n<-nrow(d)
  pr<-c("000"=(1-ps)*(1-pi)*(1-px),"100"=ps*(1-pi)*(1-px),"110"=ps*pi*(1-px),"111"=ps*pi*px,
        "010"=(1-ps)*pi*(1-px),"001"=(1-ps)*(1-pi)*px,"011"=(1-ps)*pi*px,"101"=ps*(1-pi)*px)
  exp_cnt <- exp_cnt + n*pr[lv] }
PB <- data.frame(pattern=lv, type=ifelse(lv%in%valid,"forward-concordant","discordant"),
                 observed=as.integer(tab), expected=as.numeric(exp_cnt))
PB$type <- factor(PB$type, c("forward-concordant","discordant"))
PB$pattern <- factor(PB$pattern, levels=lv)
write.csv(transform(PB, expected=round(expected,1)), file.path(OUT,"ICI_cascade_patterns.csv"), row.names=FALSE)
chip <- chisq.test(PB$observed, p=PB$expected/sum(PB$expected))$p.value
pb <- ggplot(PB, aes(pattern, observed, fill=type)) +
  geom_col(width=0.72, colour="white", linewidth=0.3) +
  geom_errorbar(aes(ymin=expected, ymax=expected), width=0.52, colour=INK, linewidth=0.55) +
  scale_fill_manual(values=c(`forward-concordant`=PLUM, discordant=SAND_DEEP), name=NULL) +
  scale_y_continuous("tumours", expand=expansion(mult=c(0,0.08))) +
  scale_x_discrete("S / I / E high-low pattern") +
  theme_publication(grid="y") +
  theme(axis.text.x=element_text(size=PT_BASE-1, colour=INK, family="mono"),
        axis.text.y=element_text(size=PT_BASE-1, colour=INK),
        axis.title=element_text(size=PT_AXIS_TITLE-0.5, colour=INK),
        legend.position="top", legend.justification="right", legend.box.spacing=unit(0.5,"mm"),
        legend.text=element_text(size=PT_BASE-1.5, colour=INK), plot.margin=margin(2,3,2,2,"mm"))

## ---- (c) cascade completion dose-response ----------------------------------
wilson <- function(k,n){ zc<-1.96; p<-k/n; d<-1+zc^2/n
  c(100*((p+zc^2/(2*n))/d - zc*sqrt(p*(1-p)/n+zc^2/(4*n^2))/d),
    100*((p+zc^2/(2*n))/d + zc*sqrt(p*(1-p)/n+zc^2/(4*n^2))/d)) }
D$nhigh <- D$Sb + D$Ib + D$Xb
PC <- do.call(rbind, lapply(0:3, function(k){ v<-D$y[D$nhigh==k]; ci<-wilson(sum(v),length(v))
  data.frame(nhigh=k, k=sum(v), n=length(v), pct=100*mean(v), lo=ci[1], hi=ci[2]) }))
write.csv(PC, file.path(OUT,"ICI_cascade_completion.csv"), row.names=FALSE)
g <- glm(y~nhigh, binomial, D); orT <- exp(coef(g)[2]); pT <- summary(g)$coef[2,4]
pc <- ggplot(PC, aes(factor(nhigh), pct, fill=nhigh)) +
  geom_col(width=0.68, colour="white", linewidth=0.3) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.18, linewidth=0.45, colour=INK) +
  geom_text(aes(label=sprintf("%.0f%%", pct), y=hi+2.5), family=PUB_FONT, size=(PT_BASE-0.5)/.pt, fontface=2, colour=INK) +
  scale_fill_gradient(low="#E7D7E2", high=PLUM, guide="none") +
  scale_y_continuous("ICB response rate (%)", limits=c(0,NA), expand=expansion(mult=c(0,0.12))) +
  scale_x_discrete("number of intact (high) layers") +
  theme_publication(grid="y") +
  theme(axis.text=element_text(size=PT_BASE-1, colour=INK),
        axis.title=element_text(size=PT_AXIS_TITLE-0.5, colour=INK), plot.margin=margin(2,2,2,2,"mm"))

tag_th <- theme(plot.tag=element_text(size=PT_PANEL_TAG, face="bold", family=PUB_FONT, colour=INK, hjust=0, vjust=1))
pa<-pa+labs(tag="a")+tag_th; pb<-pb+labs(tag="b")+tag_th; pc<-pc+labs(tag="c")+tag_th
fig <- wrap_elements(full=pa) / (wrap_elements(full=pb) + wrap_elements(full=pc) + plot_layout(widths=c(1.25,1))) +
  plot_layout(heights=c(1, 1.1)) & theme(plot.margin=margin(2,2,2,2,"mm"))
save_fig(fig, "ICI_Fig_cascade_biology", width_mm=174, height_mm=135, dir="../results/figures/composite")
cat(sprintf("done. independent=%.0f%%; staircase chi-sq p=%.1e; completion OR=%.2f p=%.1e\n", indep, chip, orT, pT))
