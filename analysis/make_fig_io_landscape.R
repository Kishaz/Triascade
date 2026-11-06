# =============================================================================
# IO therapy landscape vs the ICI atlas. Counts DISTINCT approved checkpoint-
# inhibitor REGIMENS (drug x combination-partner class, collapsing line/setting)
# per cancer site, FDA and EMA, keyed to the 16 atlas organs.
# Source: primary FDA prescribing information via Drugs@FDA / openFDA + FDA OCE
# approval notifications, cross-checked vs the Cancer Research Institute IO
# approval landscape; EMA approvals from EPARs. Compiled as-of 2026-11-05.
# Backing per-indication table (with per-row sources): ICI_io_approval_landscape.csv
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2) })
setwd("/Users/smwambu1/Desktop/Projects/ICI/analysis")
source("ici_theme.R"); .pt <- 2.845276
D <- read.csv("../results/tables/ICI_io_approval_landscape.csv", stringsAsFactors=FALSE)

ATLAS <- data.frame(
  organ = c("Brain","Thyroid","Esophagus","Lung","Gastric","Liver","Pancreas","Kidney",
            "Colon","Bladder","Skin","Breast","Ovary","Uterus","Prostate","Testis"),
  ici_status = c("none","none","null","val","val","none","null","null",
                 "none","val","val","val","none","none","none","none"),
  stringsAsFactors=FALSE)

active_fda <- function(x) !(x %in% c("no","withdrawn"))
active_ema <- function(x) x != "no"

# ---- collapse each indication to a distinct REGIMEN: drug + combination-partner class ----
canon <- function(reg){
  r <- tolower(reg); p <- character(0)
  if (grepl("ipilimumab|tremelimumab", r))            p <- c(p,"+CTLA4")
  if (grepl("relatlimab", r))                          p <- c(p,"+LAG3")
  if (grepl("axitinib|lenvatinib|cabozantinib|bevacizumab|\\bbev\\b|enfortumab|sacituzumab|cobimetinib|vemurafenib|belzutifan|brentuximab|trastuzumab|olaparib|lurbinectedin|\\bbcg\\b", r))
    p <- c(p,"+targeted")
  if (grepl("chemo|platinum|pemetrexed|paclitaxel|carbo|nab-|fluoropyrimidine|\\bfu\\b|flot|gemcitabine|cisplatin|etoposide|folfox|\\bavd\\b|chemoradiotherapy|\\bcrt\\b", r))
    p <- c(p,"+chemo")
  if (length(p)==0) p <- "mono"
  paste(sort(unique(p)), collapse="")
}
D$regimen_key <- paste(D$drug, vapply(D$regimen, canon, ""))

reg_counts <- function(df){
  f <- unique(df$regimen_key[active_fda(df$fda)])
  e <- unique(df$regimen_key[active_ema(df$ema)])
  data.frame(n_fda=length(f), n_ema=length(e),
             n_drugs=length(unique(df$drug)),
             n_withdrawn=sum(df$fda=="withdrawn"))
}
S <- do.call(rbind, by(D, D$organ, reg_counts)); S$organ <- rownames(S); rownames(S) <- NULL

atlasS <- merge(ATLAS, S, by="organ", all.x=TRUE)
for (cc in c("n_fda","n_ema","n_drugs","n_withdrawn")) atlasS[[cc]][is.na(atlasS[[cc]])] <- 0
atlasS <- atlasS[order(-atlasS$n_fda, -atlasS$n_ema), ]
extras <- S[!S$organ %in% ATLAS$organ, ]; extras <- extras[order(-extras$n_fda),]

write.csv(rbind(
  cbind(on_atlas="yes", atlasS[,c("organ","ici_status","n_fda","n_ema","n_drugs","n_withdrawn")]),
  cbind(on_atlas="no",  data.frame(organ=extras$organ, ici_status=NA,
        n_fda=extras$n_fda, n_ema=extras$n_ema, n_drugs=extras$n_drugs, n_withdrawn=extras$n_withdrawn))),
  "../results/tables/ICI_io_regimen_counts_by_cancer.csv", row.names=FALSE)

cat("=== DISTINCT approved checkpoint-inhibitor REGIMENS per atlas organ (as-of 2026-11-05) ===\n")
print(atlasS[,c("organ","ici_status","n_fda","n_ema","n_drugs","n_withdrawn")], row.names=FALSE)
cat("\n=== off-atlas sites (shown in footnote) ===\n"); print(extras[,c("organ","n_fda","n_ema","n_drugs")], row.names=FALSE)

# ---- bar: distinct FDA + EMA regimen tally per organ ----
P <- atlasS; P$organ <- factor(P$organ, levels=rev(P$organ))
L <- rbind(data.frame(organ=P$organ, reg="FDA", n=P$n_fda),
           data.frame(organ=P$organ, reg="EMA", n=P$n_ema))
L$reg <- factor(L$reg, levels=c("FDA","EMA"))
xmax <- max(L$n)

p <- ggplot() +
  geom_col(data=L, aes(x=n, y=organ, fill=reg), position=position_dodge(width=0.74), width=0.72) +
  geom_text(data=L, aes(x=n, y=organ, label=n, group=reg), position=position_dodge(width=0.74),
            hjust=-0.3, size=(PT_BASE-1.5)/.pt, family=PUB_FONT, colour=INK) +
  scale_fill_manual(values=c("FDA"="#1F6F8B","EMA"="#E07B39"), name=NULL, breaks=c("FDA","EMA")) +
  guides(fill=guide_legend(order=1)) +
  scale_x_continuous("Distinct approved checkpoint-inhibitor regimens (drug + combination partner)",
                     limits=c(0, xmax*1.1), breaks=seq(0, 20, 5), expand=expansion(mult=c(0,0))) +
  scale_y_discrete("Cancer type") +
  theme_publication(grid="x") +
  theme(axis.text.y=element_text(size=PT_BASE, colour=INK, face=2),
        axis.text.x=element_text(size=PT_BASE-1, colour=INK),
        axis.title.x=element_text(size=PT_AXIS_TITLE-0.5, colour=INK, margin=margin(t=2)),
        axis.title.y=element_text(size=PT_AXIS_TITLE-0.5, colour=INK, face=2, margin=margin(r=2)),
        legend.position="top", legend.justification="left", legend.text=element_text(size=PT_BASE-1.5, colour=INK),
        legend.box.spacing=unit(1,"mm"),
        plot.margin=margin(3,6,3,3,"mm"))

save_fig(p, "ICI_Fig_io_landscape", width_mm=174, height_mm=120, dir="../results/figures/composite")
cat("\ndone - wrote ICI_io_regimen_counts_by_cancer.csv + ICI_Fig_io_landscape.{png,svg}\n")
