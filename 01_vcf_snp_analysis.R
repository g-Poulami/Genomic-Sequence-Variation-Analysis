# =============================================================================
# 01_vcf_snp_analysis.R
# VCF loading, SNP filtering, and site statistics
# =============================================================================

suppressPackageStartupMessages({
  library(vcfR)
  library(tidyverse)
  library(data.table)
})

load_vcf <- function(path) {
  message("Loading VCF: ", path)
  vcf <- read.vcfR(path, verbose = FALSE)
  message(sprintf("  %d variants x %d samples", nrow(vcf@fix), ncol(vcf@gt) - 1))
  vcf
}

vcf_to_tidy <- function(vcf) {
  fix_df <- as.data.frame(vcf@fix, stringsAsFactors = FALSE) %>%
    mutate(POS = as.integer(POS), QUAL = as.numeric(QUAL))
  gt_df <- extract_gt_tidy(vcf, format_fields = c("GT", "DP", "GQ"), verbose = FALSE)
  list(fix = fix_df, gt = gt_df)
}

filter_snps <- function(vcf, min_qual=30, min_dp=5, max_miss=0.2, min_maf=0.01) {
  message("Filtering SNPs ...")
  n_before <- nrow(vcf@fix)
  qual_vals <- as.numeric(vcf@fix[, "QUAL"])
  vcf <- vcf[!is.na(qual_vals) & qual_vals >= min_qual, ]
  vcf <- vcf[is.biallelic(vcf), ]
  miss <- prop.miss(vcf)
  vcf  <- vcf[miss <= max_miss, ]
  dp_mat <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
  mean_dp <- rowMeans(dp_mat, na.rm = TRUE)
  vcf <- vcf[mean_dp >= min_dp, ]
  maf_v <- maf(vcf)
  vcf <- vcf[maf_v[, "Frequency"] >= min_maf, ]
  message(sprintf("  Retained %d / %d variants", nrow(vcf@fix), n_before))
  vcf
}

filter_adna <- function(vcf, min_qual=20, min_dp=3, max_miss=0.5, min_maf=0, remove_ct_ga=FALSE) {
  message("Applying aDNA filters ...")
  n_before <- nrow(vcf@fix)
  qual_vals <- as.numeric(vcf@fix[, "QUAL"])
  vcf <- vcf[!is.na(qual_vals) & qual_vals >= min_qual, ]
  vcf <- vcf[is.biallelic(vcf), ]
  miss <- prop.miss(vcf)
  vcf  <- vcf[miss <= max_miss, ]
  dp_mat <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
  mean_dp <- rowMeans(dp_mat, na.rm = TRUE)
  vcf <- vcf[mean_dp >= min_dp, ]
  if (min_maf > 0) { maf_vals <- maf(vcf); vcf <- vcf[maf_vals[, "Frequency"] >= min_maf, ] }
  ref <- vcf@fix[, "REF"]; alt <- vcf@fix[, "ALT"]
  ct_ga <- (ref == "C" & alt == "T") | (ref == "G" & alt == "A")
  message(sprintf("  C>T / G>A sites detected: %d", sum(ct_ga, na.rm = TRUE)))
  if (remove_ct_ga) { vcf <- vcf[!ct_ga, ]; message("  C>T / G>A sites removed.") }
  message(sprintf("  Retained %d / %d variants", nrow(vcf@fix), n_before))
  vcf
}

get_genotype_matrix <- function(vcf) {
  gt <- extract.gt(vcf, element = "GT", as.numeric = FALSE)
  dosage <- apply(gt, 2, function(x) {
    x <- gsub("\\|", "/", x)
    sapply(strsplit(x, "/"), function(alleles) {
      if (any(is.na(alleles)) || any(alleles == ".")) return(NA_integer_)
      sum(as.integer(alleles), na.rm = TRUE)
    })
  })
  storage.mode(dosage) <- "integer"
  rownames(dosage) <- paste0(vcf@fix[, "CHROM"], ":", vcf@fix[, "POS"])
  dosage
}

site_stats <- function(vcf) {
  fix_df <- as.data.frame(vcf@fix, stringsAsFactors = FALSE) %>%
    mutate(POS = as.integer(POS), QUAL = as.numeric(QUAL))
  dp_mat   <- extract.gt(vcf, element = "DP", as.numeric = TRUE)
  fix_df %>% select(CHROM, POS, QUAL) %>%
    mutate(mean_DP = rowMeans(dp_mat, na.rm=TRUE), missingness = prop.miss(vcf), AF = maf(vcf)[,"Frequency"])
}

ts_tv_ratio <- function(vcf) {
  ref <- vcf@fix[, "REF"]; alt <- vcf@fix[, "ALT"]
  purines <- c("A", "G")
  transitions <- (ref %in% purines & alt %in% purines) | (!ref %in% purines & !alt %in% purines)
  ts <- sum(transitions, na.rm=TRUE); tv <- sum(!transitions, na.rm=TRUE)
  message(sprintf("Ts = %d | Tv = %d | Ts/Tv = %.3f", ts, tv, ts/tv))
  invisible(ts/tv)
}

write_filtered_vcf <- function(vcf, path) {
  message("Writing filtered VCF to: ", path)
  write.vcf(vcf, file = path)
  message("  Done.")
}
