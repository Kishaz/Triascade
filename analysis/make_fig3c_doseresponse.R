# =============================================================================
# Figure 3 performance panel (dose-response): does a higher Integration score mean
# a higher ICB response rate?  Validated trials pooled; patients binned into
# within-cohort Integration tertiles (low/medium/high); ICB response rate per bin
# with Wilson 95% CIs, vs the overall base rate.  Reads per-sample scores from
# pancancer_validation.rds. House style; editable SVG.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
source("ici_theme.R"); .pt <- 2.845276
R <- readRDS("pancancer_validation.rds")
WORKS <- c("IMvigor210","GSE78220","GSE91061","GSE126044","GSE135222",
           "Gide2019","Liu2019","Puch","Kim2018gas","Rose","ISPY2")

## pooled validated patients: Integration score z-scored within its cohort
pv <- do.call(rbind, lapply(WORKS, function(k){ s<-R[[k]]$scores
  data.frame(z=as.numeric(scale(s$integ_single)), y=s$y) }))
base <- mean(pv$y)
qs <- quantile(pv$z, c(0, 1/3, 2/3, 1))
pv$bin <- cut(pv$z, breaks=qs, include.lowest=TRUE, labels=c("Low","Medium","High"))

## response rate + Wilson 95% CI per tertile
wilson <- function(k, n){ z<-1.96; p<-k/n; d<-1+z^2/n
  c(100*(p+z^2/(2*n))/d - 100*z*sqrt(p*(1-p)/n + z^2/(4*n^2))/d,
    100*(p+z^2/(2*n))/d + 100*z*sqrt(p*(1-p)/n + z^2/(4*n^2))/d) }
agg <- do.call(rbind, lapply(levels(pv$bin), function(b){ v<-pv$y[pv$bin==b]; ci<-wilson(sum(v),length(v))
  data.frame(bin=b, n=length(v), pct=100*mean(v), lo=ci[1], hi=ci[2]) }))
agg$bin <- factor(agg$bin, levels=c("Low","Medium","High"))
fillcols <- c(Low="#D9B8CC", Medium="#A86A8E", High=PLUM)

p <- ggplot(agg, aes(bin, pct, fill=bin)) +
  geom_hline(yintercept=100*base, linetype="22", colour=RULE, linewidth=0.35) +
  annotate("text", x=0.58, y=100*base+2.4, hjust=0, family=PUB_FONT, size=(PT_BASE-2.5)/.pt, colour=MUTED,
           label=sprintf("all patients  %.0f%%", 100*base)) +
  geom_col(width=0.66, colour="white", linewidth=0.3) +
  geom_errorbar(aes(ymin=lo, ymax=hi), width=0.16, linewidth=0.45, colour=INK) +
  geom_text(aes(label=sprintf("%.0f%%", pct), y=hi+3), family=PUB_FONT, size=(PT_BASE-0.5)/.pt, fontface=2, colour=INK) +
  geom_text(aes(label=sprintf("n=%d", n), y=2.5), family=PUB_FONT, size=(PT_BASE-2.5)/.pt, colour="white") +
  scale_fill_manual(values=fillcols, guide="none") +
  scale_y_continuous("ICB response rate (%)", limits=c(0,100), breaks=seq(0,100,20), expand=expansion(mult=c(0,0.01))) +
  scale_x_discrete("Integration score (within-trial tertile)") +
  theme_publication(grid="y") +
  theme(plot.tag=element_text(size=PT_PANEL_TAG, face="bold", family=PUB_FONT, colour=INK)) +
  labs(tag="d")

save_fig(p, "ICI_Fig3c_doseresponse", width_mm=85, height_mm=84, dir="../results/figures/composite")
cat(sprintf("done — response rate by Integration tertile: %s (base %.0f%%)\n",
            paste(sprintf("%s %.0f%%(n=%d)", agg$bin, agg$pct, agg$n), collapse=", "), 100*base))
