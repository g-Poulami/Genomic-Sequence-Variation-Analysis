# =============================================================================
# 03_phylogenetics.R
# Phylogenetic Tree Construction & Visualisation
# Domains: Bacterial Genomics | Ancient DNA
# =============================================================================

library(vcfR)
library(ape)       # read/write/manipulate trees
library(phangorn)  # distance matrices, NJ, parsimony
library(ggtree)    # ggplot2-based tree visualisation
library(treeio)    # read annotated tree formats (Newick, Nexus, BEAST)
library(tidyverse)


# --- 1. VCF → SNP distance matrix ---------------------------------------------

#' Build a pairwise SNP distance matrix from a genotype dosage matrix.
#' Counts sites differing between each pair of samples.
snp_distance_matrix <- function(gt_matrix) {
  # gt_matrix: samples as columns, sites as rows (from get_genotype_matrix())
  gt_t <- t(gt_matrix)               # samples × SNPs
  n    <- nrow(gt_t)
  ids  <- rownames(gt_t)
  mat  <- matrix(0, n, n, dimnames = list(ids, ids))

  for (i in seq_len(n - 1)) {
    for (j in (i + 1):n) {
      both <- !is.na(gt_t[i, ]) & !is.na(gt_t[j, ])
      diff <- sum(gt_t[i, both] != gt_t[j, both])
      mat[i, j] <- mat[j, i] <- diff
    }
  }
  as.dist(mat)
}


# --- 2. VCF → phyDat (for ML / parsimony) ------------------------------------

#' Convert genotype dosage matrix to phangorn::phyDat (binary encoding)
gt_to_phydat <- function(gt_matrix) {
  gt_t <- t(gt_matrix)
  # Encode: 0 → "0", 2 → "1", NA → "?"
  gt_bin <- apply(gt_t, 2, function(x) {
    x_chr           <- as.character(x / 2)
    x_chr[is.na(x)] <- "?"
    x_chr
  })
  phangorn::phyDat(
    data   = as.data.frame(gt_bin),
    type   = "USER",
    levels = c("0", "1")
  )
}


# --- 3. Neighbour-joining tree ------------------------------------------------

build_nj_tree <- function(dist_mat, outgroup = NULL) {
  tree <- ape::nj(dist_mat)
  tree <- ladderize(tree)
  if (!is.null(outgroup)) {
    tree <- root(tree, outgroup = outgroup, resolve.root = TRUE)
    tree <- ladderize(tree)
  }
  tree
}


# --- 4. UPGMA tree (useful for bacterial clonals) -----------------------------

build_upgma_tree <- function(dist_mat) {
  tree <- phangorn::upgma(dist_mat)
  ladderize(tree)
}


# --- 5. Bootstrap support -----------------------------------------------------

#' Add bootstrap values to an NJ tree
#' @param phydat   phyDat object (from gt_to_phydat)
#' @param tree     Starting tree (ape phylo)
#' @param nboot    Number of bootstrap replicates
bootstrap_nj <- function(phydat, tree, nboot = 100) {
  message("Bootstrapping (", nboot, " replicates) ...")
  boot <- phangorn::bootstrap.phyDat(
    phydat,
    FUN    = function(x) ape::nj(phangorn::dist.hamming(x)),
    bs     = nboot
  )
  tree_bs <- phangorn::plotBS(
    tree, boot, type = "phylogram", p = 0, digits = 0
  )
  tree_bs
}


# --- 6. Read external trees (IQ-TREE, BEAST, FastTree) -----------------------

read_iqtree  <- function(path) treeio::read.iqtree(path)
read_beast   <- function(path) treeio::read.beast(path)
read_newick  <- function(path) ape::read.tree(path)
read_nexus   <- function(path) ape::read.nexus(path)


# --- 7. Annotate tree with metadata ------------------------------------------

#' Join sample metadata onto a treedata / phylo object
#' @param tree  treedata or phylo object
#' @param meta  data.frame with column ind_id plus any annotation columns
annotate_tree <- function(tree, meta) {
  # treeio full_join requires treedata; convert phylo if needed
  if (inherits(tree, "phylo")) {
    tree <- treeio::as.treedata(tree)
  }
  treeio::full_join(tree, meta, by = c("label" = "ind_id"))
}


# --- 8. Plot: basic tree ------------------------------------------------------

plot_tree_basic <- function(tree, layout = "rectangular",
                            tip_size = 3, tip_col = "population") {
  p <- ggtree(tree, layout = layout) +
    geom_tiplab(size = tip_size, align = TRUE, linesize = 0.2) +
    theme_tree2()

  if (!is.null(tip_col) && inherits(tree, "treedata")) {
    p <- ggtree(tree, layout = layout) %<+%
      as.data.frame(tree@data) +
      aes(color = .data[[tip_col]]) +
      geom_tiplab(size = tip_size, align = TRUE, linesize = 0.2) +
      geom_tippoint(size = 2) +
      scale_color_brewer(palette = "Set2", name = tip_col) +
      theme_tree2() +
      theme(legend.position = "right")
  }
  p
}


# --- 9. Plot: circular tree with metadata heatmap ----------------------------

#' Circular tree annotated with a trait heatmap
#' @param tree        treedata object (annotated)
#' @param meta        data.frame (ind_id + columns to plot as heatmap)
#' @param heat_cols   Character vector of column names for heatmap strips
plot_tree_circular <- function(tree, meta, heat_cols) {
  p <- ggtree(tree, layout = "circular", branch.length = "none") %<+% meta +
    geom_tippoint(aes(color = population), size = 2) +
    geom_tiplab(aes(label = label), size = 2.5, align = TRUE,
                offset = 0.5, linesize = 0.2) +
    scale_color_brewer(palette = "Set1", name = "Population") +
    theme_tree2(legend.position = "right")

  for (col in heat_cols) {
    p <- gheatmap(
      p, meta[, col, drop = FALSE],
      offset = 1, width = 0.08,
      colnames_angle = 90, colnames_offset_y = 0.25,
      font.size = 3
    ) +
      scale_fill_viridis_c(name = col) +
      new_scale_fill()
  }
  p
}


# --- 10. Plot: bootstrap / posterior support labels ---------------------------

plot_tree_support <- function(tree, bs_threshold = 70) {
  ggtree(tree) +
    geom_tiplab(size = 3) +
    geom_nodelab(
      aes(
        label  = ifelse(
          !is.na(as.numeric(label)) & as.numeric(label) >= bs_threshold,
          label, ""
        ),
        subset = !isTip
      ),
      size  = 2.5,
      nudge_x = -0.002,
      color = "firebrick"
    ) +
    theme_tree2() +
    labs(caption = paste0("Values shown: bootstrap ≥ ", bs_threshold))
}


# --- 11. Save tree outputs ----------------------------------------------------

save_tree <- function(tree, out_prefix) {
  nwk_path <- paste0(out_prefix, ".nwk")
  nex_path <- paste0(out_prefix, ".nexus")

  if (inherits(tree, "treedata")) tree_phylo <- tree@phylo
  else tree_phylo <- tree

  ape::write.tree(tree_phylo, file = nwk_path)
  ape::write.nexus(tree_phylo, file = nex_path)
  message("Saved: ", nwk_path, " | ", nex_path)
}


# =============================================================================
# USAGE EXAMPLE
# =============================================================================
# source("01_vcf_snp_analysis.R")
#
# vcf    <- load_vcf("data/filtered_bacterial.vcf.gz")
# gt_mat <- get_genotype_matrix(vcf)
#
# ## Distance-based NJ tree
# dist_mat <- snp_distance_matrix(gt_mat)
# nj_tree  <- build_nj_tree(dist_mat, outgroup = "SampleA")
# save_tree(nj_tree, "output/bacterial_nj")
#
# ## Bootstrap
# phydat   <- gt_to_phydat(gt_mat)
# nj_bs    <- bootstrap_nj(phydat, nj_tree, nboot = 100)
#
# ## Annotate and plot
# meta      <- read.csv("data/sample_metadata.csv")  # ind_id, population, ...
# tree_ann  <- annotate_tree(nj_bs, meta)
# p_basic   <- plot_tree_basic(tree_ann, tip_col = "population")
# p_support <- plot_tree_support(nj_bs, bs_threshold = 70)
# p_circ    <- plot_tree_circular(tree_ann, meta, heat_cols = c("year", "country"))
#
# ## IQ-TREE output
# iqtree <- read_iqtree("output/iqtree.treefile")
