# =============================================================================
# Archetype-distribution check on IMvigor210 (Mariathasan 2018).
# Question: under the cascade (Sensory -> Integration -> Execution), which S/I/E
# high/low patterns actually occur? Do the "discordant" patterns that break the
# forward dependency (downstream high while upstream low) -- above all
# Sensing Failed (0,1,1) -- show up, or is the data a clean cascade staircase?
#
# Scoring = LOCKED v1 panel + Fig 2 formula:
#   binarize each gene b_g = 1 if expr >= cohort median (IMvigor210 cohort),
#   S = mean(b over 24 Sensory),  I = mean(b over 26 activating) - mean(b over 31 suppressive),
#   E = mean(b over 27 Execution).
# Layer called HIGH at the principled midpoint: S,E >= 0.5 ; I >= 0.
# =============================================================================
suppressPackageStartupMessages({ library(data.table) })

obj <- readRDS("data/imvigor210/imvigor210_clean.rds")
E  <- obj$E                                   # genes (symbols) x samples
pd <- obj$pd
ip <- factor(pd[["Immune phenotype"]], c("desert","excluded","inflamed"))

# ---- LOCKED v1 panel --------------------------------------------------------
p <- fread("ICI_panel_LOCKED_v1.csv", data.table = FALSE)
g_up <- toupper(p$gene); rn_up <- toupper(rownames(E))
pick <- function(mask) { g <- g_up[mask]; rownames(E)[match(intersect(g, rn_up), rn_up)] }
PAN <- list(
  Sensory      = pick(p$layer=="Sense"       & p$direction=="positive"),
  Int_activate = pick(p$layer=="Integration" & p$direction=="positive"),
  Int_suppress = pick(p$layer=="Integration" & p$direction=="negative"),
  Execution    = pick(p$layer=="Execution"   & p$direction=="positive"))
cat("panel genes measured in IMvigor210 (of locked count):\n")
print(data.frame(layer=names(PAN), measured=sapply(PAN,length),
                 locked=c(24,26,31,27)), row.names=FALSE)

# ---- binarize genes at cohort median, then layer scores ---------------------
prop_high <- function(g){ g<-g[!is.na(g)]; if(!length(g)) return(rep(NA_real_, ncol(E)))
  b <- t(apply(E[g,,drop=FALSE],1,function(v) as.integer(v >= median(v,na.rm=TRUE))))
  colMeans(b, na.rm=TRUE) }
S <- prop_high(PAN$Sensory)
I <- prop_high(PAN$Int_activate) - prop_high(PAN$Int_suppress)
Ex<- prop_high(PAN$Execution)

# ---- layer high/low (principled midpoints) ----------------------------------
Sh <- as.integer(S  >= 0.5)
Ih <- as.integer(I  >= 0)
Eh <- as.integer(Ex >= 0.5)
pat <- paste0(Sh, Ih, Eh)

NAME <- c("111"="Inflamed","110"="Execution Failed","101"="Excluded (1,0,1)",
          "100"="Excluded (1,0,0)","011"="Sensing Failed","010"="I-only (0,1,0)",
          "001"="E-only (0,0,1)","000"="Desert")
staircase <- c("111","110","100","000")        # cascade-consistent (ones form a prefix)

n <- length(pat)
tab <- as.data.frame(table(factor(pat, names(NAME))), responseName="n")
names(tab)[1] <- "pattern"
tab$archetype <- NAME[as.character(tab$pattern)]
tab$pct <- round(100*tab$n/n, 1)
tab$cascade_ok <- ifelse(as.character(tab$pattern) %in% staircase, "staircase", "DISCORDANT")
tab <- tab[order(-tab$n), c("pattern","archetype","n","pct","cascade_ok")]

cat(sprintf("\nN tumors scored: %d\n", n))
cat("\n==== S/I/E high-low pattern frequencies ====\n")
print(tab, row.names = FALSE)

disc <- sum(tab$n[tab$cascade_ok=="DISCORDANT"])
cat(sprintf("\nStaircase (cascade-consistent): %d (%.1f%%)   Discordant: %d (%.1f%%)\n",
            n-disc, 100*(n-disc)/n, disc, 100*disc/n))
cat(sprintf("Sensing Failed (0,1,1) specifically: %d (%.1f%%)\n",
            sum(pat=="011"), 100*mean(pat=="011")))

# ---- forward-dependency probabilities (does downstream-high need upstream-high?) ----
cat("\n==== forward-dependency check  P(downstream high | upstream low) ====\n")
cat(sprintf("P(I high | S low) = %.2f    P(I high | S high) = %.2f\n",
            mean(Ih[Sh==0]), mean(Ih[Sh==1])))
cat(sprintf("P(E high | I low) = %.2f    P(E high | I high) = %.2f\n",
            mean(Eh[Ih==0]), mean(Eh[Ih==1])))
cat(sprintf("P(E high | S low) = %.2f    P(E high | S high) = %.2f\n",
            mean(Eh[Sh==0]), mean(Eh[Sh==1])))
cat("\nSpearman correlations among continuous layer scores:\n")
cc <- cor(cbind(S=S,I=I,E=Ex), method="spearman", use="complete.obs")
print(round(cc,2))

# ---- cross-tab with Mariathasan immune phenotype ----------------------------
cat("\n==== pattern x Mariathasan immune phenotype (counts) ====\n")
print(addmargins(table(archetype=NAME[pat], phenotype=ip)))
cat("\nmean layer score by immune phenotype (sanity: S,E rise desert->inflamed):\n")
print(round(rbind(S=tapply(S,ip,mean,na.rm=TRUE),
                  I=tapply(I,ip,mean,na.rm=TRUE),
                  E=tapply(Ex,ip,mean,na.rm=TRUE)),3))

saveRDS(data.frame(sample=colnames(E), S=S, I=I, E=Ex, pattern=pat,
                   phenotype=as.character(ip)), "archetype_scores_IMvigor210.rds")
cat("\ndone\n")
