# =============================================================================
# 03_phylogenetics.R
# Tree construction, annotation, and plotting
# =============================================================================

suppressPackageStartupMessages({
  library(ape)
  library(phangorn)
  library(ggtree)
  library(treeio)
  library(tidyverse)
})

snp_distance_matrix <- function(gt_matrix) {
  message("Computing SNP distance matrix ...")
  dist(t(gt_matrix), method="manhattan") / 2
}

gt_to_phydat <- function(gt_matrix) {
  message("Converting genotype matrix to phyDat ...")
  t_mat <- t(gt_matrix); t_mat[is.na(t_mat)] <- "?"
  phyDat(t_mat, type="USER", levels=c("0","1","2","?"))
}

build_nj_tree <- function(dist_mat, outgroup=NULL) {
  message("Building NJ tree ...")
  tree <- nj(dist_mat)
  if (!is.null(outgroup)) { tree <- root(tree, outgroup=outgroup, resolve.root=TRUE); message("  Rooted on: ", outgroup) }
  tree
}

build_upgma_tree <- function(dist_mat) {
  message("Building UPGMA tree ..."); upgma(dist_mat)
}

bootstrap_nj <- function(phydat, tree, nboot=100) {
  message(sprintf("Running %d bootstrap replicates ...", nboot))
  bs <- bootstrap.phyDat(phydat, FUN=function(x) nj(dist.hamming(x)), bs=nboot)
  tree_bs <- plotBS(tree, bs, type="none", p=0)
  message("  Bootstrap support added."); tree_bs
}

read_iqtree  <- function(path) { message("Reading IQ-TREE: ", path); treeio::read.iqtree(path) }
read_beast   <- function(path) { message("Reading BEAST: ",   path); treeio::read.beast(path) }
read_newick  <- function(path) { message("Reading Newick: ",  path); ape::read.tree(path) }
read_nexus   <- function(path) { message("Reading Nexus: ",   path); ape::read.nexus(path) }

annotate_tree <- function(tree, meta) {
  message("Annotating tree ...")
  if (inherits(tree, "phylo")) tree <- as.treedata(tree)
  tree %<+% meta
}

plot_tree_basic <- function(tree, color_by="population") {
  ggtree(tree, layout="rectangular") + aes(color=.data[[color_by]]) +
    geom_tiplab(size=2.5, align=TRUE) + theme_tree2() + labs(color=color_by)
}

plot_tree_circular <- function(tree, meta, heat_cols) {
  p <- ggtree(tree, layout="circular") + geom_tiplab(size=2, align=TRUE) + theme_tree2()
  for (col in heat_cols) p <- gheatmap(p, meta[,col,drop=FALSE], width=0.1, colnames_angle=90, font.size=2)
  p
}

plot_tree_support <- function(tree, bs_threshold=70) {
  ggtree(tree) + geom_tiplab(size=2.5) +
    geom_label2(aes(subset=!isTip & !is.na(as.numeric(label)) & as.numeric(label)>=bs_threshold, label=label),
                size=2.5, fill="white", label.padding=unit(0.1,"lines")) + theme_tree2()
}

save_tree <- function(tree, prefix) {
  write.tree(tree,  file=paste0(prefix,".nwk"))
  write.nexus(tree, file=paste0(prefix,".nexus"))
  message("Tree saved: ", prefix, ".nwk / .nexus")
}
