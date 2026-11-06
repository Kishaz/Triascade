# =============================================================================
# SUPPLEMENTARY: co-expression structure of the locked 98-gene ICI panel.
# Empirical (data-driven) co-expression of the panel genes; reference = IMvigor210
# (n=348). Spearman gene-gene correlation; positive co-expression network (igraph)
# + correlation heatmap. Produced in TWO colourings:
#   *_coexpression          : nodes coloured by signed direction (activating / suppressive)
#   *_coexpression_modules   : nodes coloured by the 6 designed modules (M1-M6)
# Expression recovers ~2 coarse programs (immune-activating manifold vs stromal /
# suppressive); the 6 modules are a finer, literature-defined partition (provenance
# to be re-derived separately). House style; bottom-centred legend; editable SVG.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(igraph); library(ggraph); library(patchwork) })
source("ici_theme.R"); source("ici_locked_reference.R"); .pt <- 2.845276
PANEL <- load_panel("ICI_98gene_panel.csv")
asg <- read.csv("ICI_98_module_assignment.csv", stringsAsFactors=FALSE); asg$gene <- toupper(asg$gene)
asg <- asg[!duplicated(asg$gene), ]
asg$Layer <- c(Sense="Sensory", Integration="Integration", Execution="Execution")[asg$layer]
asg$dir <- ifelse(asg$direction=="positive", "activating (+)", "suppressive (−)")
obj <- readRDS("data/imvigor210/imvigor210_clean.rds"); E <- obj$E; rownames(E) <- toupper(rownames(E))
g98 <- intersect(asg$gene, rownames(E)); asg <- asg[match(g98, asg$gene), ]
M <- E[g98, ]; M <- M[apply(M,1,sd)>0, ]; g98 <- rownames(M); asg <- asg[match(g98, asg$gene), ]
C <- cor(t(M), method="spearman")

DIRCOL <- c(`activating (+)`=PLUM, `suppressive (−)`=SAND_DEEP)
MODCOL <- c(M1="#C9A66B", M2="#6A2C5A", M3="#8A6D30", M4="#E79FC4", M5="#C44E6E", M6="#2B1C33")
MODSH  <- c(M1="M1 antigen", M2="M2 trafficking", M3="M3 stroma",
            M4="M4 angiogenic", M5="M5 effector", M6="M6 exhaustion")

## shared network graph + layout (colouring-independent) ----------------------
thr <- 0.40; A <- C; A[A < thr] <- 0; diag(A) <- 0
gr <- graph_from_adjacency_matrix(A, mode="undirected", weighted=TRUE)
V(gr)$name <- asg$gene; V(gr)$dir <- asg$dir; V(gr)$module <- asg$module; V(gr)$deg <- degree(gr)
set.seed(11); lay <- create_layout(gr, layout="fr")
MARK <- c("HLA-A","B2M","STAT1","TAP1","CXCL9","CXCL10","CXCR3","CCL5","CD8A","GZMB","PRF1","IFNG",
          "PDCD1","LAG3","CTLA4","TOX","IDO1","TGFB1","FAP","COL1A1","VEGFA","MMP9")
labdat <- lay[lay$name %in% MARK, ]

## centred bottom legend (hand-placed) ---------------------------------------
make_leg <- function(items, cols, ncol){
  n <- length(items); nrow <- ceiling(n/ncol)
  xs <- seq(0.5 - (ncol-1)*0.16/1, 0.5 + (ncol-1)*0.16/1, length.out=ncol)  # column centres (dot)
  if(ncol==2) xs <- c(0.40, 0.56); if(ncol==3) xs <- c(0.17, 0.44, 0.71)
  ys <- if(nrow==1) 0.5 else seq(0.72, 0.28, length.out=nrow)
  df <- data.frame(lab=names(items), fill=unname(items),
                   col=xs[((seq_len(n)-1) %% ncol)+1], row=ys[((seq_len(n)-1) %/% ncol)+1])
  ggplot(df) + xlim(0,1) + ylim(0,1) + theme_void() +
    geom_point(aes(col, row), shape=21, size=3.2, fill=df$fill, colour="white", stroke=0.3) +
    geom_text(aes(col+0.03, row, label=lab), hjust=0, family=PUB_FONT, size=(PT_BASE-1.5)/.pt, colour=INK) +
    theme(plot.margin=margin(0,0,0,0,"mm"))
}

build <- function(scheme){
  if(scheme=="dir"){ fillvar<-asg$dir; cmap<-DIRCOL; items<-DIRCOL; lcol<-2
    ordf <- order(factor(asg$dir,names(DIRCOL)), factor(asg$Layer,c("Sensory","Integration","Execution")), asg$module, asg$gene)
    barcol <- function(go) DIRCOL[asg$dir[match(go,asg$gene)]]; nm<-"ICI_SuppFig_coexpression"
  } else { fillvar<-asg$module; cmap<-MODCOL; items<-setNames(MODCOL, MODSH[names(MODCOL)]); lcol<-3
    ordf <- order(factor(asg$module,paste0("M",1:6)), asg$gene)
    barcol <- function(go) MODCOL[asg$module[match(go,asg$gene)]]; nm<-"ICI_SuppFig_coexpression_modules" }
  lay$fillv <- fillvar
  pa <- ggraph(lay) +
    geom_edge_link(aes(alpha=weight), edge_colour=MUTED, edge_width=0.25, show.legend=FALSE) +
    geom_node_point(aes(fill=I(cmap[fillvar]), size=deg), shape=21, colour="white", stroke=0.2) +
    ggrepel::geom_text_repel(data=labdat, aes(x=x, y=y, label=name), size=4.2/.pt, family=PUB_FONT,
              colour=INK, fontface=2, max.overlaps=Inf, segment.size=0.12, segment.colour=MUTED,
              box.padding=0.3, point.padding=0.1, min.segment.length=0) +
    scale_edge_alpha(range=c(0.04,0.4)) + scale_size(range=c(0.7,2.6), guide="none") +
    coord_equal(clip="off") + theme_void() + theme(legend.position="none", plot.margin=margin(1,2,1,2,"mm"))
  go <- g98[ordf]; HC <- C[go,go]
  HD <- data.frame(gi=factor(rep(go,each=length(go)),levels=go), gj=factor(rep(go,times=length(go)),levels=rev(go)), r=as.vector(HC))
  bar <- data.frame(g=factor(go,levels=go))
  pb <- ggplot(HD, aes(gi,gj,fill=r)) + geom_tile() +
    scale_fill_gradient2(low=SAND_DEEP, mid="white", high=PLUM, midpoint=0, limits=c(-1,1),
      name="Spearman r", guide=guide_colourbar(barheight=unit(20,"mm"), barwidth=unit(2.4,"mm"))) +
    geom_tile(data=bar, aes(x=g,y=-1.5,fill=NULL), fill=barcol(go), height=3, inherit.aes=FALSE) +
    coord_equal(clip="off") + scale_x_discrete(NULL,breaks=NULL) + scale_y_discrete(NULL,breaks=NULL) +
    theme_void() + theme(legend.title=element_text(size=PT_BASE-1,colour=INK,family=PUB_FONT),
      legend.text=element_text(size=PT_BASE-2,colour=INK,family=PUB_FONT), plot.margin=margin(1,1,1,4,"mm"))
  tag_th <- theme(plot.tag=element_text(size=PT_PANEL_TAG,face="bold",family=PUB_FONT,colour=INK,hjust=0,vjust=1))
  pa<-pa+labs(tag="a")+tag_th; pb<-pb+labs(tag="b")+tag_th
  leg <- make_leg(items, cmap, lcol)
  top <- wrap_elements(full=pa) + wrap_elements(full=pb) + plot_layout(widths=c(1.25,1))
  fig <- top / leg + plot_layout(heights=c(1, if(lcol==3) 0.12 else 0.08)) & theme(plot.margin=margin(2,2,2,2,"mm"))
  save_fig(fig, nm, width_mm=174, height_mm=if(lcol==3) 116 else 112, dir="../results/figures/composite")
}
build("dir"); build("module")
# Supplementary module-membership table = the FULL authoritative 98-gene assignment
# (all 100 layer assignments, incl. the IFNG/TGFB1 dual-layer rows), independent of
# which genes survive into the co-expression network above (that set is filtered to
# genes expressed/variable in the IMvigor210 reference, so it is legitimately smaller).
full <- read.csv("ICI_98_module_assignment.csv", stringsAsFactors=FALSE); full$gene <- toupper(full$gene)
write.csv(data.frame(gene=full$gene, module=full$module,
                     layer=c(Sense="Sensory", Integration="Integration", Execution="Execution")[full$layer],
                     direction=ifelse(full$direction=="positive", "activating (+)", "suppressive (−)")),
          "../results/tables/ICI_coexpression_modules.csv", row.names=FALSE)
cat("done - co-expression figures (direction + 6-module), ref IMvigor210 n=", ncol(E), "\n")
