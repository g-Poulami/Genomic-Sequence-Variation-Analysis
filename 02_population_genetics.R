# =============================================================================
# 02_population_genetics.R
# PCA, FST, ADMIXTURE parsing
# =============================================================================

suppressPackageStartupMessages({
  library(vcfR)
  library(adegenet)
  library(hierfstat)
  library(tidyverse)
})

vcf_to_genlight <- function(vcf) {
  message("Converting VCF to genlight ...")
  gl <- vcfR2genlight(vcf)
  message(sprintf("  genlight: %d individuals x %d SNPs", nInd(gl), nLoc(gl)))
  gl
}

ld_prune <- function(gl, window=50, r2_thresh=0.2) {
  message(sprintf("LD pruning (window=%d, r2=%.2f) ...", window, r2_thresh))
  mat <- as.matrix(gl); n_loci <- ncol(mat); keep <- rep(TRUE, n_loci)
  for (i in seq(1, n_loci-1)) {
    if (!keep[i]) next
    for (j in seq(i+1, min(i+window-1, n_loci))) {
      if (!keep[j]) next
      r2 <- suppressWarnings(cor(mat[,i], mat[,j], use="pairwise.complete.obs")^2)
      if (!is.na(r2) && r2 > r2_thresh) keep[j] <- FALSE
    }
  }
  message(sprintf("  Retained %d / %d SNPs", sum(keep), n_loci))
  gl[, keep]
}

run_pca <- function(gl, ncomp=10) {
  message("Running PCA ...")
  pca <- glPca(gl, nf=ncomp, parallel=FALSE)
  var_exp <- (pca$eig / sum(pca$eig)) * 100
  for (i in seq_len(min(5, length(var_exp)))) message(sprintf("  PC%d: %.2f%%", i, var_exp[i]))
  pca
}

pca_scores_df <- function(pca, meta) {
  as.data.frame(pca$scores) %>% rownames_to_column("ind_id") %>% left_join(meta, by="ind_id")
}

compute_fst <- function(vcf, pop_map) {
  message("Computing pairwise FST ...")
  gl <- vcfR2genlight(vcf); pop(gl) <- pop_map[indNames(gl)]
  pairwise.WCfst(genlight2hierfstat(gl))
}

overall_fst <- function(vcf, pop_map) {
  gl <- vcfR2genlight(vcf); pop(gl) <- pop_map[indNames(gl)]
  fst_val <- basic.stats(genlight2hierfstat(gl))$overall["Fst"]
  message(sprintf("Overall FST = %.4f", fst_val)); invisible(fst_val)
}

read_admixture <- function(q_files, ind_ids) {
  message("Reading ADMIXTURE Q files ...")
  purrr::map_dfr(q_files, function(f) {
    k_val <- as.integer(gsub(".*\\.([0-9]+)\\.Q$", "\\1", basename(f)))
    q_mat <- read.table(f, header=FALSE)
    colnames(q_mat) <- paste0("Cluster", seq_len(ncol(q_mat)))
    q_mat %>% mutate(ind_id=ind_ids, K_value=k_val) %>%
      pivot_longer(starts_with("Cluster"), names_to="cluster", values_to="proportion")
  })
}

admixture_cv <- function(log_files) {
  purrr::map_dfr(log_files, function(f) {
    lines <- readLines(f); cv_line <- grep("CV error", lines, value=TRUE)
    if (!length(cv_line)) return(NULL)
    data.frame(K=as.integer(gsub(".*K=(\\d+).*","\\1",cv_line[1])),
               cv_error=as.numeric(gsub(".*:\\s*([0-9.]+).*","\\1",cv_line[1])))
  }) %>% arrange(K)
}

pseudo_haploidise <- function(gt_matrix) {
  mat <- gt_matrix; het <- which(mat==1, arr.ind=TRUE)
  if (nrow(het) > 0) mat[het] <- sample(c(0L,2L), nrow(het), replace=TRUE)
  message(sprintf("  Pseudo-haploidised %d heterozygous calls", nrow(het)))
  mat
}
