# =============================================================================
# 01_vcf_snp_analysis.R
# VCF Loading, SNP Filtering & Variant Analysis
# Domains: Bacterial Genomics | Ancient DNA
# =============================================================================

# --- Dependencies -------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation", "GenomicRanges"), ask = FALSE)
install.packages(c("vcfR", "tidyverse", "data.table"), repos = "https://cloud.r-project.org")

library(vcfR)
library(tidyverse)
library(data.table)


# --- 1. Load VCF --------------------------------------------------------------

load_vcf <- function(vcf_path, verbose = TRUE) {
  if (verbose) message("Loading VCF: ", vcf_path)
  vcf <- read.vcfR(vcf_path, verbose = verbose)
  if (verbose) {
    message("  Variants : ", nrow(vcf@fix))
    message("  Samples  : ", ncol(vcf@gt) - 1)
  }
  vcf
}


# --- 2. VCF to tidy data frame ------------------------------------------------

vcf_to_tidy <- function(vcf) {
  fix  <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  info <- INFO2df(vcf)
  gt   <- extract.gt(vcf, element = "GT")

  fix <- fix %>%
    mutate(
      POS    = as.integer(POS),
      QUAL   = as.numeric(QUAL),
      is_snp = nchar(REF) == 1 & nchar(ALT) == 1
    )

  list(fix = fix, info = info, gt = gt)
}


# --- 3. Standard SNP quality filters ------------------------------------------

#' @param vcf       vcfR object
#' @param min_qual  Minimum QUAL score         (default 30)
#' @param min_dp    Minimum mean depth          (default 5)
#' @param max_miss  Maximum missingness / site  (default 0.2)
#' @param min_maf   Minimum minor allele freq   (default 0.01)
filter_snps <- function(vcf,
                        min_qual = 30,
                        min_dp   = 5,
                        max_miss = 0.2,
                        min_maf  = 0.01) {

  qual       <- as.numeric(vcf@fix[, "QUAL"])
  pass_qual  <- !is.na(qual) & qual >= min_qual

  dp         <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
  mean_dp    <- rowMeans(dp, na.rm = TRUE)
  pass_dp    <- !is.na(mean_dp) & mean_dp >= min_dp

  miss       <- apply(dp, 1, function(x) mean(is.na(x)))
  pass_miss  <- miss <= max_miss

  maf_vals   <- maf(vcf)
  pass_maf   <- !is.na(maf_vals[, "Frequency"]) &
                maf_vals[, "Frequency"] >= min_maf

  pass_snp   <- is.biallelic(vcf)

  keep <- pass_qual & pass_dp & pass_miss & pass_maf & pass_snp

  message(sprintf(
    "Filter summary:\n  Input      : %d\n  QUAL >= %d  : %d\n  DP >= %d    : %d\n  miss <= %.0f%% : %d\n  MAF >= %.2f  : %d\n  biallelic  : %d\n  Retained   : %d",
    nrow(vcf@fix), min_qual, sum(pass_qual),
    min_dp, sum(pass_dp),
    max_miss * 100, sum(pass_miss),
    min_maf, sum(pass_maf),
    sum(pass_snp), sum(keep)
  ))

  vcf[keep, ]
}


# --- 4. Ancient DNA-specific filters ------------------------------------------

#' Extra filters for aDNA — lower DP threshold; flag/remove C>T and G>A
#' transitions attributable to post-mortem deamination damage.
#' @param remove_ct_ga  If TRUE, exclude deamination-prone transitions
filter_adna <- function(vcf,
                        min_qual     = 30,
                        min_dp       = 3,
                        max_miss     = 0.5,
                        remove_ct_ga = FALSE) {

  # aDNA: no MAF cut (pseudo-haploid calling)
  vcf_f <- filter_snps(vcf,
                       min_qual = min_qual,
                       min_dp   = min_dp,
                       max_miss = max_miss,
                       min_maf  = 0)

  fix <- as.data.frame(vcf_f@fix, stringsAsFactors = FALSE)

  ct_ga <- with(fix,
    (REF == "C" & ALT == "T") | (REF == "G" & ALT == "A") |
    (REF == "T" & ALT == "C") | (REF == "A" & ALT == "G")
  )

  n_ct_ga <- sum(ct_ga, na.rm = TRUE)
  message(sprintf(
    "C>T / G>A transitions (deamination candidates): %d (%.1f%%)",
    n_ct_ga, 100 * n_ct_ga / nrow(fix)
  ))

  if (remove_ct_ga) {
    message("  Removing deamination-prone transitions.")
    vcf_f <- vcf_f[!ct_ga, ]
  } else {
    message("  Retaining (set remove_ct_ga = TRUE to exclude).")
  }

  vcf_f
}


# --- 5. Genotype dosage matrix (0 / 1 / 2 / NA) --------------------------------

get_genotype_matrix <- function(vcf) {
  gt <- extract.gt(vcf, element = "GT", convertNA = TRUE)
  gt_num <- apply(gt, 2, function(col) {
    col <- gsub("\\|", "/", col)
    sapply(strsplit(col, "/"), function(alleles) {
      if (any(is.na(alleles)) || "." %in% alleles) return(NA_real_)
      sum(as.integer(alleles))
    })
  })
  rownames(gt_num) <- paste0(vcf@fix[, "CHROM"], ":", vcf@fix[, "POS"])
  gt_num
}


# --- 6. Per-site summary statistics -------------------------------------------

site_stats <- function(vcf) {
  gt_num <- get_genotype_matrix(vcf)
  fix    <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)

  data.frame(
    CHROM     = fix$CHROM,
    POS       = as.integer(fix$POS),
    REF       = fix$REF,
    ALT       = fix$ALT,
    QUAL      = as.numeric(fix$QUAL),
    n_called  = rowSums(!is.na(gt_num)),
    prop_miss = rowMeans(is.na(gt_num)),
    mean_dp   = rowMeans(
      extract.gt(vcf, "DP", as.numeric = TRUE), na.rm = TRUE),
    alt_freq  = rowMeans(gt_num / 2, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}


# --- 7. Transition / Transversion ratio ----------------------------------------

ts_tv_ratio <- function(vcf) {
  fix <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  ref <- fix$REF; alt <- fix$ALT

  ts <- sum(
    (ref == "A" & alt == "G") | (ref == "G" & alt == "A") |
    (ref == "C" & alt == "T") | (ref == "T" & alt == "C"),
    na.rm = TRUE)
  tv <- sum(
    (ref %in% c("A","G") & alt %in% c("C","T")) |
    (ref %in% c("C","T") & alt %in% c("A","G")),
    na.rm = TRUE)

  ratio <- ts / tv
  message(sprintf("Ts: %d  |  Tv: %d  |  Ts/Tv: %.3f", ts, tv, ratio))
  invisible(list(Ts = ts, Tv = tv, ratio = ratio))
}


# --- 8. Write filtered VCF ----------------------------------------------------

write_filtered_vcf <- function(vcf, out_path) {
  write.vcf(vcf, file = out_path)
  message("Written: ", out_path)
}


# =============================================================================
# USAGE EXAMPLE
# =============================================================================
# vcf_raw  <- load_vcf("data/raw_variants.vcf.gz")
#
# ## Bacterial
# vcf_bact <- filter_snps(vcf_raw, min_qual = 30, min_dp = 5, max_miss = 0.2)
# stats_b  <- site_stats(vcf_bact)
# ts_tv_ratio(vcf_bact)
# write_filtered_vcf(vcf_bact, "data/filtered_bacterial.vcf.gz")
#
# ## Ancient DNA
# vcf_adna <- filter_adna(vcf_raw, min_dp = 3, remove_ct_ga = FALSE)
# stats_a  <- site_stats(vcf_adna)
# gt_mat   <- get_genotype_matrix(vcf_adna)
