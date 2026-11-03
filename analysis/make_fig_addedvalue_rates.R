# =============================================================================
# Added value in PLAIN RESPONSE RATES (intuitive companion to the OR dumbbell).
# Pooled validated patients, each biomarker AND the ICI Integration score split at
# the within-cohort median. Within every biomarker stratum, the ICB response rate
# for ICI-low vs ICI-high: the ICI split separates responders even among patients
# the established marker calls "high" -> that is the added value, in percentages.
# Wilson 95% CIs; dashed = overall response rate. Reads pancancer_validation.rds.
# House style; single canvas; editable SVG (no textLength).
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R"); .pt <- 2.845276
R <- readRDS("pancancer_validation.rds")
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")
MK <- c(GEP="T-cell-inflamed GEP", IFNG="IFN-γ signature",
        CYT="Cytolytic activity", CD274="PD-L1 mRNA (CD274)")

hi <- function(v) ifelse(v > median(v), "high", "low")      # within-cohort median split
D <- do.call(rbind, lapply(WORKS, function(k){ s<-R[[k]]$scores
  data.frame(y=s$y, ICI=hi(s$integ_single),
             GEP=hi(s$GEP), IFNG=hi(s$IFNG), CYT=hi(s$CYT), CD274=hi(s$CD274)) }))
base <- 100*mean(D$y)

wilson <- function(k,n){ z<-1.96; p<-k/n; d<-1+z^2/n
  c(100*(p+z^2/(2*n))/d - 100*z*sqrt(p*(1-p)/n+z^2/(4*n^2))/d,
    100*(p+z^2/(2*n))/d + 100*z*sqrt(p*(1-p)/n+z^2/(4*n^2))/d) }

## cross-classified response rate per (marker stratum) x (ICI stratum) ----------
grid <- expand.grid(m=names(MK), ms=c("low","high"), is=c("low","high"), stringsAsFactors=FALSE)
AG <- do.call(rbind, Map(function(m, ms, is){
  v <- D$y[D[[m]]==ms & D$ICI==is]; ci <- wilson(sum(v), length(v))
  data.frame(marker=MK[m], mstrat=ms, ICI=is, n=length(v),
             pct=100*mean(v), lo=ci[1], hi=ci[2], row.names=NULL)
}, grid$m, grid$ms, grid$is))
AG$marker <- factor(AG$marker, levels=unname(MK))
AG$mstrat <- factor(AG$mstrat, levels=c("low","high"), labels=c("Low","High"))
AG$ICI    <- factor(AG$ICI,    levels=c("low","high"), labels=c("ICI low","ICI high"))
write.csv(AG, "../results/tables/ICI_added_value_rates.csv", row.names=FALSE)

p <- ggplot(AG, aes(mstrat, pct, fill=ICI)) +
  geom_hline(yintercept=base, linetype="22", colour=RULE, linewidth=0.35) +
  geom_col(position=position_dodge(width=0.74), width=0.66, colour="white", linewidth=0.3) +
  geom_errorbar(aes(ymin=lo, ymax=hi), position=position_dodge(width=0.74), width=0.18,
                linewidth=0.45, colour=INK) +
  geom_text(aes(label=sprintf("%.0f%%", pct), y=hi+3), position=position_dodge(width=0.74),
            family=PUB_FONT, size=(PT_BASE-0.5)/.pt, fontface=2, colour=INK) +
  geom_text(aes(label=n, y=3.2), position=position_dodge(width=0.74),
            family=PUB_FONT, size=(PT_BASE-3)/.pt, colour="white") +
  facet_wrap(~marker, nrow=1) +
  scale_fill_manual(values=c(`ICI low`=SAND, `ICI high`=PLUM), name=NULL) +
  scale_y_continuous("ICB response rate (%)", limits=c(0,62), breaks=seq(0,60,20),
                     expand=expansion(mult=c(0,0.02))) +
  scale_x_discrete("Established-biomarker stratum  (within-cohort median split)") +
  annotate("text", x=0.52, y=base+2.6, hjust=0, family=PUB_FONT, size=(PT_BASE-2.5)/.pt,
           colour=MUTED, label=sprintf("all %.0f%%", base)) +
  theme_publication(grid="y") +
  theme(axis.title.x=element_text(size=PT_AXIS_TITLE, colour=INK),
        axis.text=element_text(size=PT_BASE, colour=INK),
        strip.text=element_text(size=PT_BASE, face=2, colour=INK),
        legend.position="top", legend.justification="right",
        legend.text=element_text(size=PT_BASE+0.5, colour=INK), legend.box.spacing=unit(1,"mm"))

save_fig(p, "ICI_Fig_addedvalue_rates", width_mm=174, height_mm=78, dir="../results/figures/composite")
cat(sprintf("done — within-stratum response rates (base %.0f%%)\n", base))
