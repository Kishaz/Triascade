# =====================================================================
# FIGURE 2 — Integration predicts benefit; Sensory and Execution are
# prognostic. I-SPY2, the only randomized cohort.
#
#   a  layer AUC in each arm, joined, with the interaction p beneath
#   b  odds ratio per SD by arm, with CIs
#   c  receptor subgroup forest (TNBC vs HR+/HER2-)
#   d  pCR rate by Integration tertile, by arm
#
# Everything is recomputed here from the score file and GSE194040 so
# the figure cannot drift from the text.
#
# Output: Aug_1_ICI_FINAL_PAPER/Figures/Fig2_ispy2_predictive.pdf/.png
#         + Fig2_source_values.csv
# =====================================================================
SCRIPT_VERSION <- "fig2_ispy2_predictive v1.0  2026-08-23"
cat(SCRIPT_VERSION, "\n\n")
suppressPackageStartupMessages({ library(ggplot2); library(pROC) })

MAN  <- "/Users/dbedi/Desktop/R data FINAL/ICI_Manuscript"
OUT  <- file.path(MAN,"Aug_1_ICI_FINAL_PAPER","Figures")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
SCOR <- file.path(MAN,"MASTER_GRID","ICI_binary_ispy2_scores.csv")
SER  <- file.path(MAN,"ISPY2", c("GSE194040-GPL20078_series_matrix.txt.gz",
                                 "GSE194040-GPL30493_series_matrix.txt.gz"))
stopifnot(file.exists(SCOR))

d <- read.csv(SCOR, stringsAsFactors=FALSE)
d$pid_str <- as.character(d$pid_str)
d$y   <- as.integer(d$pcr_binary)
d$arm <- factor(d$arm, levels=c("control","pembro"))
d <- d[!is.na(d$y) & !is.na(d$arm), ]
LAY <- c("Sensory","Integration","Execution")
for (L in LAY) d[[paste0(L,"_z")]] <- as.numeric(scale(d[[L]]))
cat(sprintf("n=%d (control %d, pembro %d)\n", nrow(d),
            sum(d$arm=="control"), sum(d$arm=="pembro")))

# receptor subtype from the series matrices
grab <- function(f){ L <- readLines(gzfile(f), warn=FALSE)
  ch <- grep("^!Sample_characteristics_ch1", L, value=TRUE)
  m <- list(pid = gsub('"','',strsplit(grep("^!Sample_title",L,value=TRUE)[1],"\t")[[1]][-1]))
  for (l in ch){ v <- gsub('"','',strsplit(l,"\t")[[1]][-1])
    m[[tolower(gsub("[^a-z0-9]+","_",tolower(sub(":.*$","",v[1]))))]] <-
      trimws(sub("^[^:]*:\\s*","",v)) }
  as.data.frame(m, stringsAsFactors=FALSE) }
have_sub <- all(file.exists(SER))
if (have_sub) {
  md <- do.call(rbind, lapply(SER, grab))
  key <- if ("patient_id" %in% names(md)) "patient_id" else "pid"
  md$pid2 <- as.character(md[[key]])
  md$hr <- suppressWarnings(as.integer(md$hr)); md$her2 <- suppressWarnings(as.integer(md$her2))
  md$subtype <- ifelse(md$her2==1,"HER2+", ifelse(md$hr==1,"HR+/HER2-","TNBC"))
  d <- merge(d, md[,c("pid2","subtype")], by.x="pid_str", by.y="pid2", all.x=TRUE)
  cat("subtype merged:", sum(!is.na(d$subtype)), "\n")
}

auc1 <- function(y,v) as.numeric(auc(roc(y,v,quiet=TRUE,direction="<")))
inter <- function(dat, col) {
  m1 <- glm(as.formula(sprintf("y ~ %s * arm", col)), dat, family=binomial)
  m0 <- glm(as.formula(sprintf("y ~ %s + arm", col)), dat, family=binomial)
  a  <- anova(m0,m1,test="Chisq"); cf <- coef(summary(m1)); k <- grep(":",rownames(cf))
  ci <- suppressMessages(confint(m1))
  list(p=a$`Pr(>Chi)`[2],
       or_c=exp(cf[2,1]), or_c_lo=exp(ci[2,1]), or_c_hi=exp(ci[2,2]),
       or_p=exp(cf[2,1]+cf[k,1]),
       or_p_lo=exp(cf[2,1]+cf[k,1]-1.96*sqrt(cf[2,2]^2+cf[k,2]^2)),
       or_p_hi=exp(cf[2,1]+cf[k,1]+1.96*sqrt(cf[2,2]^2+cf[k,2]^2))) }

PAL <- c(control="#5F5E5A", pembro="#B03A2E")
src <- list()

# ---- panel a --------------------------------------------------------
A <- do.call(rbind, lapply(LAY, function(L){
  r <- inter(d, paste0(L,"_z"))
  data.frame(layer=L,
    arm=c("control","pembro"),
    auc=c(auc1(d$y[d$arm=="control"], d[[paste0(L,"_z")]][d$arm=="control"]),
          auc1(d$y[d$arm=="pembro"],  d[[paste0(L,"_z")]][d$arm=="pembro"])),
    p=r$p, stringsAsFactors=FALSE) }))
A$layer <- factor(A$layer, levels=LAY)
lab <- unique(A[,c("layer","p")])
lab$txt <- sprintf("interaction\np = %.3f", lab$p)
src$panel_a <- A

pa <- ggplot(A, aes(arm, auc, group=layer)) +
  geom_hline(yintercept=0.5, linetype=2, linewidth=0.3, colour="grey70") +
  geom_line(colour="grey55", linewidth=0.6) +
  geom_point(aes(colour=arm), size=3) +
  geom_text(data=lab, aes(x=1.5, y=0.44, label=txt), inherit.aes=FALSE,
            size=2.7, colour="grey25", lineheight=0.95) +
  scale_colour_manual(values=PAL, guide="none") +
  scale_y_continuous(limits=c(0.40,0.85), breaks=seq(0.4,0.8,0.1)) +
  facet_wrap(~layer, nrow=1) +
  labs(x=NULL, y="AUC for pCR", subtitle="a   Discrimination by arm") +
  theme_minimal(base_size=9) +
  theme(panel.grid.minor=element_blank(), strip.text=element_text(face="bold", size=9),
        plot.subtitle=element_text(face="bold", size=9))

# ---- panel b --------------------------------------------------------
B <- do.call(rbind, lapply(LAY, function(L){
  r <- inter(d, paste0(L,"_z"))
  data.frame(layer=L, arm=c("control","pembro"),
    or=c(r$or_c,r$or_p), lo=c(r$or_c_lo,r$or_p_lo), hi=c(r$or_c_hi,r$or_p_hi),
    stringsAsFactors=FALSE) }))
B$layer <- factor(B$layer, levels=rev(LAY))
src$panel_b <- B

pb <- ggplot(B, aes(or, layer, colour=arm)) +
  geom_vline(xintercept=1, linetype=2, linewidth=0.3, colour="grey70") +
  geom_errorbarh(aes(xmin=lo,xmax=hi), height=0, linewidth=0.5,
                 position=position_dodge(width=0.55)) +
  geom_point(size=2.4, position=position_dodge(width=0.55)) +
  scale_colour_manual(values=PAL, name=NULL) +
  scale_x_log10() +
  labs(x="odds ratio for pCR per SD", y=NULL,
       subtitle="b   Effect size by arm") +
  theme_minimal(base_size=9) +
  theme(legend.position="top", panel.grid.minor=element_blank(),
        plot.subtitle=element_text(face="bold", size=9))

# ---- panel c: receptor subgroups ------------------------------------
pc <- NULL
if (have_sub && sum(!is.na(d$subtype))>100) {
  C <- list()
  for (st in c("TNBC","HR+/HER2-")) {
    s <- d[!is.na(d$subtype) & d$subtype==st, ]
    if (nrow(s) < 40) next
    s$Iz <- as.numeric(scale(s$Integration))
    r <- inter(s, "Iz")
    C[[st]] <- data.frame(subtype=st, arm=c("control","pembro"),
      or=c(r$or_c,r$or_p), lo=c(r$or_c_lo,r$or_p_lo), hi=c(r$or_c_hi,r$or_p_hi),
      p=r$p, n=nrow(s), stringsAsFactors=FALSE) }
  if (length(C)) {
    C <- do.call(rbind, C); src$panel_c <- C
    C$lab <- sprintf("%s (n=%d)\ninteraction p = %.3f", C$subtype, C$n, C$p)
    pc <- ggplot(C, aes(or, lab, colour=arm)) +
      geom_vline(xintercept=1, linetype=2, linewidth=0.3, colour="grey70") +
      geom_errorbarh(aes(xmin=lo,xmax=hi), height=0, linewidth=0.5,
                     position=position_dodge(width=0.55)) +
      geom_point(size=2.4, position=position_dodge(width=0.55)) +
      scale_colour_manual(values=PAL, guide="none") + scale_x_log10() +
      labs(x="odds ratio for pCR per SD of Integration", y=NULL,
           subtitle="c   Integration by receptor subtype") +
      theme_minimal(base_size=9) +
      theme(panel.grid.minor=element_blank(),
            plot.subtitle=element_text(face="bold", size=9)) } }

# ---- panel d: tertiles ----------------------------------------------
q <- quantile(d$Integration, c(1/3,2/3))
d$tert <- cut(d$Integration, c(-Inf,q,Inf), labels=c("low","mid","high"))
D <- aggregate(y ~ tert + arm, d, function(v) c(rate=mean(v), n=length(v)))
D <- data.frame(tert=D$tert, arm=D$arm, rate=D$y[,"rate"], n=D$y[,"n"])
mt <- glm(y ~ tert*arm, d, family=binomial); m0 <- glm(y ~ tert+arm, d, family=binomial)
pt <- anova(m0,mt,test="Chisq")$`Pr(>Chi)`[2]
src$panel_d <- D; src$tertile_p <- pt
cat(sprintf("\ntertile x arm interaction p = %.4f\n", pt))

pd <- ggplot(D, aes(tert, rate, fill=arm)) +
  geom_col(position=position_dodge(width=0.75), width=0.65) +
  geom_text(aes(label=sprintf("%.0f%%\n(n=%d)", 100*rate, n)),
            position=position_dodge(width=0.75), vjust=-0.25, size=2.4, colour="grey25") +
  scale_fill_manual(values=PAL, name=NULL) +
  scale_y_continuous(labels=function(v) paste0(100*v,"%"), limits=c(0,0.72), expand=c(0,0)) +
  labs(x="Integration tertile", y="pCR rate",
       subtitle=sprintf("d   pCR by Integration tertile (interaction p = %.4f)", pt)) +
  theme_minimal(base_size=9) +
  theme(legend.position="top", panel.grid.major.x=element_blank(),
        panel.grid.minor=element_blank(),
        plot.subtitle=element_text(face="bold", size=9))

# ---- assemble --------------------------------------------------------
fp <- file.path(OUT,"Fig2_ispy2_predictive")
combine <- function(a,b,c,e){
  if (requireNamespace("patchwork", quietly=TRUE)) { library(patchwork)
    if (is.null(c)) (a)/(b|e) else (a)/(b|c)/(e) + plot_layout(heights=c(1,1,1))
  } else if (requireNamespace("gridExtra", quietly=TRUE)) {
    gridExtra::arrangeGrob(a, gridExtra::arrangeGrob(b, if(is.null(c)) e else c, ncol=2),
                           if (is.null(c)) NULL else e, ncol=1)
  } else NULL }
g <- combine(pa,pb,pc,pd)
if (is.null(g)) {
  cat("\npatchwork/gridExtra missing - writing panels separately\n")
  ggsave(paste0(fp,"_a.pdf"), pa, width=7.5, height=2.8, device=cairo_pdf)
  ggsave(paste0(fp,"_b.pdf"), pb, width=4.6, height=2.8, device=cairo_pdf)
  if (!is.null(pc)) ggsave(paste0(fp,"_c.pdf"), pc, width=4.6, height=2.8, device=cairo_pdf)
  ggsave(paste0(fp,"_d.pdf"), pd, width=5.0, height=3.0, device=cairo_pdf)
} else {
  ggsave(paste0(fp,".pdf"), g, width=9.0, height=9.0, device=cairo_pdf)
  ggsave(paste0(fp,".png"), g, width=9.0, height=9.0, dpi=400)
  cat("\nWrote:\n ", paste0(fp,".pdf"), "\n ", paste0(fp,".png"), "\n") }

for (k in names(src)[grepl("panel", names(src))]) {
  write.csv(src[[k]], file.path(OUT, sprintf("Fig2_source_%s.csv", k)), row.names=FALSE) }
write.csv(data.frame(item="tertile_interaction_p", value=src$tertile_p),
          file.path(OUT,"Fig2_source_tertile_p.csv"), row.names=FALSE)
cat(" ", file.path(OUT,"Fig2_source_values.csv"), "\n")
cat("\nEvery number in the figure is in that CSV - use it to write the legend.\n")
