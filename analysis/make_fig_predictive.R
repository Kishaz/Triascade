# =============================================================================
# Fig 5 — ICI is PREDICTIVE of checkpoint-blockade benefit, not merely prognostic.
#  (a) Randomized interaction: pCR rate for Integration Low vs High within each
#      arm of two randomized trials (per-trial median-split interaction p).
#  (b) Pooled continuous confirmation: P(pCR) vs Integration score diverges
#      between the ICB and control arms across both checkpoint trials (n=348,
#      trial-stratified Integration x arm p=0.035, OR/SD 1.82).
#  (c) Cross-context specificity: Integration single-patient AUC by treatment
#      context across 16 cohorts. Works for checkpoint blockade, partial for the
#      related CTLA-4 axis, and falls to/below chance for non-checkpoint
#      immunotherapy and the ccRCC negative control.
# Reads the run_predictive_prognostic.R tables. House style; editable SVG.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
source("ici_theme.R"); .pt <- 2.845276
T  <- "../results/tables"
RT <- read.csv(file.path(T,"ICI_predictive_interaction_rates.csv"), check.names=FALSE)
CK <- read.csv(file.path(T,"ICI_predictive_checkpoint_patients.csv"), check.names=FALSE)
CN <- read.csv(file.path(T,"ICI_predictive_interaction_continuous.csv"), check.names=FALSE)
IA <- read.csv(file.path(T,"ICI_predictive_interaction.csv"), check.names=FALSE)
CC <- read.csv(file.path(T,"ICI_predictive_crosscontext.csv"), check.names=FALSE)

## ---- panel a: randomized interaction, rate view ----------------------------
trial_lab <- c("I-SPY2 (pembro vs pac)"="I-SPY2", "GeparNuevo (durva vs placebo)"="GeparNuevo")
cmp_lab   <- c("I-SPY2 (pembro vs pac)"="pembrolizumab vs paclitaxel",
               "GeparNuevo (durva vs placebo)"="durvalumab vs placebo")
RT$trialf <- factor(trial_lab[RT$trial], levels=unname(trial_lab))
RT$arm    <- factor(RT$arm, levels=c("Control arm","ICB arm"))
RT$integ  <- factor(RT$integ, levels=c("Low","High"))
pint <- IA[IA$layer=="Integration" & IA$scoring=="locked single-patient", c("trial","interaction_p")]
pint$trialf <- factor(trial_lab[pint$trial], levels=unname(trial_lab))
pint$lab <- sprintf("%s\ninteraction p = %.3f", cmp_lab[pint$trial], pint$interaction_p)

pa <- ggplot(RT, aes(arm, pct, fill=integ)) +
  geom_col(position=position_dodge(width=0.72), width=0.62, colour="white", linewidth=0.3) +
  geom_errorbar(aes(ymin=lo, ymax=hi), position=position_dodge(width=0.72),
                width=0.16, linewidth=0.45, colour=INK) +
  geom_text(aes(label=sprintf("%.0f%%", pct), y=hi+6), position=position_dodge(width=0.72),
            family=PUB_FONT, size=(PT_BASE-0.5)/.pt, fontface=2, colour=INK) +
  geom_text(aes(label=sprintf("%d/%d", k, n), y=4), position=position_dodge(width=0.72),
            family=PUB_FONT, size=(PT_BASE-3.5)/.pt, colour="white") +
  geom_text(data=pint, aes(x=0.5, y=97, label=lab), inherit.aes=FALSE, hjust=0, vjust=1,
            family=PUB_FONT, size=(PT_BASE-1.5)/.pt, fontface=3, colour=INK, lineheight=0.95) +
  facet_wrap(~trialf, nrow=1) +
  scale_fill_manual(values=c(Low=SAND, High=PLUM), name="Integration score",
                    labels=c("Low (≤ median)","High (> median)")) +
  scale_y_continuous("pCR rate (%)", limits=c(0,100), breaks=seq(0,100,25),
                     expand=expansion(mult=c(0,0.02))) +
  scale_x_discrete(NULL) +
  theme_publication(grid="y") +
  theme(axis.text=element_text(size=PT_BASE, colour=INK),
        axis.title.y=element_text(size=PT_AXIS_TITLE, colour=INK, margin=margin(r=1)),
        strip.text=element_text(size=PT_BASE+0.5, face=2, colour=INK),
        legend.position="top", legend.justification="right", legend.box.spacing=unit(1,"mm"),
        legend.title=element_text(size=PT_BASE, colour=INK),
        panel.spacing=unit(6,"mm"), plot.margin=margin(2,4,2,3,"mm"))

## ---- panel b: cross-context specificity gradient ---------------------------
ctx_lv  <- c("anti-PD-(L)1","anti-CTLA4","chemo+ICB combo",
             "negative control (ccRCC)","specificity control (bispecific)")
ctx_sh  <- c("anti-PD-(L)1"="anti-PD-(L)1", "anti-CTLA4"="anti-CTLA-4",
             "chemo+ICB combo"="chemo+ICB", "negative control (ccRCC)"="neg-ctrl (ccRCC)",
             "specificity control (bispecific)"="spec-ctrl (bispecific)")
CC$context <- factor(ctx_sh[CC$context], levels=unname(ctx_sh[ctx_lv]))
CC <- CC[order(CC$context, CC$AUC), ]
CC$cohort <- factor(CC$cohort, levels=CC$cohort)        # sorted within context
CC$side   <- ifelse(CC$AUC >= 0.5, "at/above chance", "below chance")

pooled <- "Integration pooled AUC 0.72 [0.66–0.78]"   # from Fig 3 random-effects meta
pb <- ggplot(CC, aes(AUC, cohort, colour=side)) +
  geom_vline(xintercept=0.5, linetype="22", colour=RULE, linewidth=0.4) +
  geom_errorbarh(aes(xmin=lo, xmax=hi), height=0.55, linewidth=0.7) +
  geom_point(size=2.7) +
  facet_grid(context~., scales="free_y", space="free_y", switch="y") +
  scale_colour_manual(values=c(`at/above chance`=PLUM, `below chance`=SAND_DEEP), guide="none") +
  scale_x_continuous("Integration single-patient AUC (ICB response)",
                     limits=c(0,1.02), breaks=seq(0,1,0.25), expand=expansion(mult=0)) +
  scale_y_discrete(NULL) +
  theme_publication(grid="x") +
  theme(axis.text.y=element_text(size=PT_BASE-0.5, colour=INK),
        axis.text.x=element_text(size=PT_BASE, colour=INK),
        strip.placement="outside",
        strip.text.y.left=element_text(size=PT_BASE-1, face=2, colour=INK, angle=0, hjust=1),
        panel.spacing=unit(1.6,"mm"), plot.margin=margin(2,5,2,2,"mm"))

## ---- panel b (new): pooled continuous interaction, logistic curves ----------
CK$arm <- factor(CK$arm, c("ctrl","icb"))
mC <- glm(y ~ Integration*arm + trial, binomial, CK)
gx <- seq(min(CK$Integration), max(CK$Integration), length.out=80)
gr <- expand.grid(Integration=gx, arm=factor(c("ctrl","icb"), c("ctrl","icb")), trial="I-SPY2")
pr <- predict(mC, gr, type="link", se.fit=TRUE)
gr$p  <- plogis(pr$fit); gr$lo <- plogis(pr$fit-1.96*pr$se.fit); gr$hi <- plogis(pr$fit+1.96*pr$se.fit)
pint2 <- sprintf("Integration × arm\np = %.3f",
                 CN$interaction_p[CN$layer=="Integration"])
pc <- ggplot(gr, aes(Integration, p, colour=arm, fill=arm)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.14, colour=NA) +
  geom_line(linewidth=1.1) +
  geom_point(data=transform(CK, yv=ifelse(y==1, 0.97, 0.03)), aes(Integration, yv),
             position=position_jitter(height=0.015, width=0), size=0.7, alpha=0.5, show.legend=FALSE) +
  annotate("text", x=min(CK$Integration), y=0.97, hjust=0, vjust=1, lineheight=0.95,
           family=PUB_FONT, size=(PT_BASE-1.5)/.pt, fontface=3, colour=INK, label=pint2) +
  scale_colour_manual(values=c(ctrl=SAND_DEEP, icb=PLUM), name=NULL,
                      labels=c("Control arm","ICB arm")) +
  scale_fill_manual(values=c(ctrl=SAND_DEEP, icb=PLUM), guide="none") +
  scale_y_continuous("P(pCR)", limits=c(0,1), breaks=seq(0,1,0.25), expand=expansion(mult=0)) +
  scale_x_continuous("Integration score (z)", breaks=seq(-2,2,1), expand=expansion(mult=0.02)) +
  theme_publication(grid="y") +
  theme(axis.title=element_text(size=PT_AXIS_TITLE, colour=INK),
        axis.text=element_text(size=PT_BASE, colour=INK),
        legend.position="top", legend.justification="right", legend.box.spacing=unit(1,"mm"),
        legend.text=element_text(size=PT_BASE-0.5, colour=INK),
        plot.margin=margin(2,5,2,2,"mm"))

tag_th <- theme(plot.tag=element_text(size=PT_PANEL_TAG, face="bold", family=PUB_FONT,
                                      colour=INK, hjust=0, vjust=1))
pa <- pa + labs(tag="a") + tag_th
pc <- pc + labs(tag="b") + tag_th
pb <- pb + labs(tag="c") + tag_th
# wrap the top-row panels so their y-axis columns are not stretched to match
# panel c's long cohort labels (keeps each y-title hugging its axis)
top <- wrap_elements(full=pa) + wrap_elements(full=pc) + plot_layout(widths=c(0.95, 1.18))
fig <- top / pb + plot_layout(heights=c(0.95, 1.6)) & theme(plot.margin=margin(2,2,2,2,"mm"))

save_fig(fig, "ICI_Fig_predictive", width_mm=174, height_mm=192, dir="../results/figures/composite")
cat("done — predictive-vs-prognostic figure;", pooled, "\n")
