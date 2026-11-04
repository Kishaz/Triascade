# =============================================================================
# SUPPLEMENTARY: co-expression structure of the locked 98-gene ICI panel.
# Empirical (data-driven) check that the panel genes co-express in blocks that
# recover the designed 6 modules / 3-layer cascade. Reference = IMvigor210
# (n=348, uniform platform). Spearman gene-gene correlation; a positive
# co-expression network (igraph) with Louvain community detection; concordance
# of data-driven communities vs the designed modules.
#   (a) co-expression network, nodes = genes coloured by designed module
#   (b) gene x gene correlation heatmap, ordered by module (block structure)
# NOTE: literature/PRISMA provenance for the 98 will be re-derived separately;
# this figure is expression-based only. House style; editable SVG.
# =============================================================================
suppressPackageStartupMessages({ library(ggplot2); library(igraph); library(ggraph); library(patchwork) })
source("ici_theme.R"); source("ici_locked_reference.R"); .pt <- 2.845276
PANEL <- load_panel("ICI_98gene_panel.csv")
asg <- read.csv("ICI_98_module_assignment.csv", stringsAsFactors=FALSE); asg$gene <- toupper(asg$gene)
asg <- asg[!duplicated(asg$gene), ]                       # one module per gene for layout
asg$Layer <- c(Sense="Sensory", Integration="Integration", Execution="Execution")[asg$layer]
obj <- readRDS("data/imvigor210/imvigor210_clean.rds"); E <- obj$E; rownames(E) <- toupper(rownames(E))
g98 <- intersect(asg$gene, rownames(E)); asg <- asg[match(g98, asg$gene), ]
M <- E[g98, ]; M <- M[apply(M,1,sd)>0, ]; g98 <- rownames(M); asg <- asg[match(g98, asg$gene), ]
C <- cor(t(M), method="spearman")
asg$dir <- ifelse(asg$direction=="positive", "activating (+)", "suppressive (−)")

DIRCOL <- c(`activating (+)`=PLUM, `suppressive (−)`=SAND_DEEP)

## ---- network: positive co-expression edges + Louvain communities -----------
thr <- 0.40
A <- C; A[A < thr] <- 0; diag(A) <- 0
gr <- graph_from_adjacency_matrix(A, mode="undirected", weighted=TRUE)
V(gr)$dir <- asg$dir; V(gr)$layer <- asg$Layer; V(gr)$name <- asg$gene
set.seed(7); comm <- cluster_louvain(gr, weights=E(gr)$weight)
ct  <- table(community=membership(comm), module=asg$module)
ari <- igraph::compare(membership(comm), as.integer(factor(asg$module)), method="adjusted.rand")
# degree for node sizing / hub labels
V(gr)$deg <- degree(gr)
set.seed(11); lay <- create_layout(gr, layout="fr")
MARK <- c("HLA-A","B2M","STAT1","TAP1","CXCL9","CXCL10","CXCR3","CCL5","CD8A","GZMB","PRF1","IFNG",
          "PDCD1","LAG3","CTLA4","TOX","IDO1","TGFB1","FAP","COL1A1","VEGFA","MMP9")
labdat <- lay[lay$name %in% MARK, ]

pa <- ggraph(lay) +
  geom_edge_link(aes(alpha=weight), edge_colour=MUTED, edge_width=0.25, show.legend=FALSE) +
  geom_node_point(aes(fill=dir, size=deg), shape=21, colour="white", stroke=0.2) +
  ggrepel::geom_text_repel(data=labdat, aes(x=x, y=y, label=name), size=4.2/.pt, family=PUB_FONT,
            colour=INK, fontface=2, max.overlaps=Inf, segment.size=0.12, segment.colour=MUTED,
            box.padding=0.3, point.padding=0.1, min.segment.length=0) +
  scale_fill_manual(values=DIRCOL, guide="none") +
  scale_edge_alpha(range=c(0.04,0.4)) + scale_size(range=c(0.7,2.6), guide="none") +
  coord_equal(clip="off") + theme_void() +
  theme(legend.position="none", plot.margin=margin(1,2,1,2,"mm"))

## ---- heatmap: correlation ordered by module --------------------------------
ord <- order(factor(asg$dir, c("activating (+)","suppressive (−)")),
             factor(asg$Layer, c("Sensory","Integration","Execution")), asg$module, asg$gene)
go <- g98[ord]; HC <- C[go, go]
HD <- data.frame(gi=factor(rep(go, each=length(go)), levels=go),
                 gj=factor(rep(go, times=length(go)), levels=rev(go)),
                 r =as.vector(HC))
dirbar <- data.frame(g=factor(go, levels=go), dir=asg$dir[match(go,asg$gene)])
pb <- ggplot(HD, aes(gi, gj, fill=r)) +
  geom_tile() +
  scale_fill_gradient2(low=SAND_DEEP, mid="white", high=PLUM, midpoint=0, limits=c(-1,1),
                       name="Spearman r", guide=guide_colourbar(barheight=unit(20,"mm"), barwidth=unit(2.4,"mm"))) +
  geom_tile(data=dirbar, aes(x=g, y=-1.5, fill=NULL), fill=DIRCOL[dirbar$dir], height=3, inherit.aes=FALSE) +
  coord_equal(clip="off") +
  scale_x_discrete(NULL, breaks=NULL) + scale_y_discrete(NULL, breaks=NULL) +
  theme_void() +
  theme(legend.title=element_text(size=PT_BASE-1, colour=INK, family=PUB_FONT),
        legend.text=element_text(size=PT_BASE-2, colour=INK, family=PUB_FONT),
        plot.margin=margin(1,1,1,4,"mm"))

# centered two-row legend (activating / suppressive) for the bottom
leg <- ggplot() + xlim(0,1) + ylim(0,1) + theme_void() +
  annotate("point", x=0.47, y=0.72, shape=21, size=3.2, fill=PLUM, colour="white", stroke=0.3) +
  annotate("text",  x=0.49, y=0.72, label="activating (+)", hjust=0, family=PUB_FONT, size=(PT_BASE-1)/.pt, colour=INK) +
  annotate("point", x=0.47, y=0.28, shape=21, size=3.2, fill=SAND_DEEP, colour="white", stroke=0.3) +
  annotate("text",  x=0.49, y=0.28, label="suppressive (−)", hjust=0, family=PUB_FONT, size=(PT_BASE-1)/.pt, colour=INK) +
  theme(plot.margin=margin(0,0,0,0,"mm"))

tag_th <- theme(plot.tag=element_text(size=PT_PANEL_TAG, face="bold", family=PUB_FONT, colour=INK, hjust=0, vjust=1))
pa <- pa + labs(tag="a") + tag_th; pb <- pb + labs(tag="b") + tag_th
top <- wrap_elements(full=pa) + wrap_elements(full=pb) + plot_layout(widths=c(1.25,1))
fig <- top / leg + plot_layout(heights=c(1, 0.08)) & theme(plot.margin=margin(2,2,2,2,"mm"))
save_fig(fig, "ICI_SuppFig_coexpression", width_mm=174, height_mm=112, dir="../results/figures/composite")
write.csv(data.frame(gene=asg$gene, module=asg$module, layer=asg$Layer,
          coexpr_community=as.integer(membership(comm))), "../results/tables/ICI_coexpression_modules.csv", row.names=FALSE)
cat(sprintf("done. genes=%d, ref=IMvigor210 n=%d; Louvain communities=%d; ARI(community vs designed module)=%.2f\n",
            length(g98), ncol(E), length(comm), ari))
print(ct)
