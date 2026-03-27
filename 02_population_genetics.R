# =============================================================================
# 02_population_genetics.R
# PCA, FST & ADMIXTURE from VCF
# Domains: Bacterial Genomics | Ancient DNA
# =============================================================================

library(vcfR)
library(adegenet)      # genind / genlight objects
library(hierfstat)     # pairwise FST
library(tidyverse)
library(ggplot2)


# --- 1. VCF → genlight (memory-efficient SNP matrix) -------------------------

vcf_to_genlight <- function(vcf) {
  gl <- vcfR2genlight(vcf)
  ploidy(gl) <- 2        # diploid; set to 1 for aDNA pseudo-haploid
  gl
}


# --- 2. PCA -------------------------------------------------------------------

#' Run PCA on a genlight object
#' @param gl        genlight object
#' @param ncomp     Number of principal components to retain (default 10)
#' @param n_cores   Threads for glPca (default 1)
run_pca <- function(gl, ncomp = 10, n_cores = 1) {
  message("Running PCA on ", nInd(gl), " individuals x ", nLoc(gl), " SNPs ...")
  pca <- glPca(gl, nf = ncomp, n.cores = n_cores)
  message("Done. Variance explained:")
  var_pct <- round(100 * pca$eig / sum(pca$eig), 2)
  print(head(var_pct, ncomp))
  pca
}


#' Extract tidy PCA scores data frame
#' @param pca    glPca object
#' @param meta   data.frame with at minimum columns: ind_id, population
pca_scores_df <- function(pca, meta = NULL) {
  scores <- as.data.frame(pca$scores)
  scores$ind_id <- rownames(scores)

  if (!is.null(meta)) {
    scores <- left_join(scores, meta, by = "ind_id")
  }
  scores
}


# --- 3. FST -------------------------------------------------------------------

#' Compute pairwise FST (Weir & Cockerham) using hierfstat
#' @param vcf   Filtered vcfR object
#' @param pop_map  data.frame with columns: sample, population
compute_fst <- function(vcf, pop_map) {
  # Build genind
  gi <- vcfR2genind(vcf)

  # Assign populations
  samp_order <- indNames(gi)
  pop_map_ord <- pop_map[match(samp_order, pop_map$sample), ]
  pop(gi) <- pop_map_ord$population

  message("Computing pairwise FST across ",
          length(unique(pop_map$population)), " populations ...")

  # Convert to hierfstat format
  hf <- genind2hierfstat(gi)
  fst_mat <- pairwise.WCfst(hf)

  message("Pairwise FST matrix:")
  print(round(fst_mat, 4))
  fst_mat
}


#' Overall Fst summary (basic)
overall_fst <- function(vcf, pop_map) {
  gi  <- vcfR2genind(vcf)
  samp_order <- indNames(gi)
  pop_map_ord <- pop_map[match(samp_order, pop_map$sample), ]
  pop(gi)  <- pop_map_ord$population
  hf  <- genind2hierfstat(gi)
  bs  <- basic.stats(hf)
  message(sprintf("Overall Fst = %.4f", bs$overall["Fst"]))
  invisible(bs)
}


# --- 4. ADMIXTURE (Q-matrix) parsing -----------------------------------------

#' Read ADMIXTURE .Q files and tidy into long format
#' ADMIXTURE must be run externally (CLI); this function imports results.
#' @param q_files   Named vector: names = K values, values = .Q file paths
#' @param ind_ids   Character vector of individual IDs (same order as .Q rows)
read_admixture <- function(q_files, ind_ids) {
  purrr::map_dfr(names(q_files), function(k) {
    q  <- read.table(q_files[[k]], header = FALSE)
    colnames(q) <- paste0("K", k, "_", seq_len(ncol(q)))
    q$ind_id <- ind_ids
    q$K      <- as.integer(k)
    q %>%
      pivot_longer(
        cols      = starts_with("K"),
        names_to  = "cluster",
        values_to = "proportion"
      )
  })
}


#' Cross-validation error from ADMIXTURE log files
#' @param log_files  Named vector: names = K, values = log file paths
admixture_cv <- function(log_files) {
  purrr::map_dfr(names(log_files), function(k) {
    lines <- readLines(log_files[[k]])
    cv_line <- grep("CV error", lines, value = TRUE)
    cv_val  <- as.numeric(sub(".*CV error \\(K=(\\d+)\\): ([0-9.]+).*",
                              "\\2", cv_line))
    data.frame(K = as.integer(k), CV_error = cv_val)
  })
}


# --- 5. aDNA: pseudo-haploid conversion ---------------------------------------

#' For aDNA: randomly sample one allele per heterozygous site (pseudo-haploid)
#' This avoids biases introduced by heterozygous calls on low-coverage data.
pseudo_haploidise <- function(gt_matrix) {
  apply(gt_matrix, c(1, 2), function(x) {
    if (is.na(x)) return(NA_real_)
    if (x == 1) return(sample(c(0, 2), 1))   # het → random homozygote
    x
  })
}


# --- 6. LD pruning (before PCA / ADMIXTURE) -----------------------------------

#' Simple window-based LD pruning on a genlight object
#' Uses correlation (r^2) between SNP pairs; keeps one SNP per correlated pair.
#' @param gl       genlight object
#' @param window   Window size in SNPs
#' @param r2_thresh r^2 threshold for pruning
ld_prune <- function(gl, window = 50, r2_thresh = 0.2) {
  message("LD pruning: ", nLoc(gl), " input SNPs ...")
  snp_names <- locNames(gl)

  keep <- rep(TRUE, nLoc(gl))
  mat  <- as.matrix(gl)

  for (i in seq_len(nLoc(gl) - 1)) {
    if (!keep[i]) next
    end <- min(i + window - 1, nLoc(gl))
    for (j in (i + 1):end) {
      if (!keep[j]) next
      pair <- cor(mat[, i], mat[, j], use = "pairwise.complete.obs")^2
      if (!is.na(pair) && pair > r2_thresh) keep[j] <- FALSE
    }
  }

  message("  Retained: ", sum(keep), " SNPs after pruning.")
  gl[, keep]
}


# =============================================================================
# USAGE EXAMPLE
# =============================================================================
# source("01_vcf_snp_analysis.R")
#
# vcf  <- load_vcf("data/filtered_bacterial.vcf.gz")
# gl   <- vcf_to_genlight(vcf)
# gl_p <- ld_prune(gl, window = 50, r2_thresh = 0.2)
#
# ## PCA
# pca     <- run_pca(gl_p, ncomp = 10)
# meta    <- read.csv("data/sample_metadata.csv")   # cols: ind_id, population
# scores  <- pca_scores_df(pca, meta)
#
# ## FST
# pop_map <- meta %>% rename(sample = ind_id)
# fst_mat <- compute_fst(vcf, pop_map)
#
# ## ADMIXTURE (after running externally)
# q_files <- c("2" = "admixture/out.2.Q",
#              "3" = "admixture/out.3.Q",
#              "4" = "admixture/out.4.Q")
# q_long  <- read_admixture(q_files, ind_ids = meta$ind_id)
#
# ## aDNA pseudo-haploid
# gt_mat   <- get_genotype_matrix(vcf)          # from 01_vcf_snp_analysis.R
# gt_haplo <- pseudo_haploidise(gt_mat)
