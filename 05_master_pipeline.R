# =============================================================================
# 05_master_pipeline.R
# End-to-end orchestration — start here
# =============================================================================

set.seed(42)
source("01_vcf_snp_analysis.R")
source("02_population_genetics.R")
source("03_phylogenetics.R")
source("04_visualisation.R")

CONFIG <- list(
  vcf_path     = "data/raw_variants.vcf.gz",
  metadata     = "data/sample_metadata.csv",
  domain       = "bacterial",
  outgroup     = NULL,
  min_qual     = 30, min_dp = 5, max_miss = 0.2, min_maf = 0.01,
  remove_ct_ga = FALSE,
  pca_n_comp   = 10, ld_window = 50, ld_r2 = 0.2,
  admixture_k       = 2:5,
  admixture_q_dir   = "admixture/",
  admixture_log_dir = "admixture/logs/",
  n_bootstrap  = 100, bs_threshold = 70
)

dir.create("figures", showWarnings=FALSE)
dir.create("output",  showWarnings=FALSE)

message("\n=== Loading metadata ===")
meta    <- read.csv(CONFIG$metadata, stringsAsFactors=FALSE)
pop_map <- setNames(meta$population, meta$ind_id)

message("\n=== Step 1: VCF Loading and SNP Filtering ===")
vcf_raw <- load_vcf(CONFIG$vcf_path)
vcf_filtered <- if (CONFIG$domain == "adna") {
  filter_adna(vcf_raw, CONFIG$min_qual, CONFIG$min_dp, CONFIG$max_miss, CONFIG$min_maf, CONFIG$remove_ct_ga)
} else {
  filter_snps(vcf_raw, CONFIG$min_qual, CONFIG$min_dp, CONFIG$max_miss, CONFIG$min_maf)
}
ts_tv_ratio(vcf_filtered)
stats_df <- site_stats(vcf_filtered)
write_filtered_vcf(vcf_filtered, "output/filtered.vcf.gz")
write.csv(stats_df, "output/site_stats.csv", row.names=FALSE)

message("\n=== Step 2: QC Figures ===")
save_figure(combine_panels(list(plot_qual_distribution(stats_df), plot_missingness(stats_df),
  plot_afs(stats_df), plot_mutation_spectrum(vcf_filtered)), ncol=2), "figures/01_qc_panel", 14, 10)
save_figure(plot_snp_density(stats_df), "figures/02_snp_density")

message("\n=== Step 3: PCA ===")
gt_matrix <- get_genotype_matrix(vcf_filtered)
if (CONFIG$domain == "adna") gt_matrix <- pseudo_haploidise(gt_matrix)
gl      <- vcf_to_genlight(vcf_filtered)
gl_ld   <- ld_prune(gl, CONFIG$ld_window, CONFIG$ld_r2)
pca_obj <- run_pca(gl_ld, CONFIG$pca_n_comp)
scores  <- pca_scores_df(pca_obj, meta)
write.csv(scores, "output/pca_scores.csv", row.names=FALSE)
save_figure(combine_panels(list(plot_pca(scores,pca_obj), plot_pca_scree(pca_obj)), ncol=2), "figures/03_pca")

message("\n=== Step 4: FST ===")
fst_mat <- compute_fst(vcf_filtered, pop_map)
overall_fst(vcf_filtered, pop_map)
write.csv(as.matrix(fst_mat), "output/fst_matrix.csv")
save_figure(plot_fst_heatmap(fst_mat), "figures/04_fst_heatmap")

message("\n=== Step 5: ADMIXTURE ===")
q_files   <- list.files(CONFIG$admixture_q_dir,   pattern="\\.Q$",   full.names=TRUE)
log_files <- list.files(CONFIG$admixture_log_dir,  pattern="\\.log$", full.names=TRUE)
if (length(q_files) > 0) {
  q_long <- read_admixture(q_files, meta$ind_id)
  save_figure(plot_admixture(q_long, meta, CONFIG$admixture_k), "figures/05_admixture_barplot", 14, 8)
} else message("  No ADMIXTURE .Q files found — run ADMIXTURE externally first.")
if (length(log_files) > 0) save_figure(plot_admixture_cv(admixture_cv(log_files)), "figures/05b_admixture_cv")

message("\n=== Step 6: Phylogenetics ===")
dist_mat <- snp_distance_matrix(gt_matrix)
nj_tree  <- build_nj_tree(dist_mat, CONFIG$outgroup)
phydat   <- gt_to_phydat(gt_matrix)
nj_bs    <- bootstrap_nj(phydat, nj_tree, CONFIG$n_bootstrap)
save_tree(nj_tree, "output/nj_tree")
save_figure(plot_tree_basic(annotate_tree(nj_tree, meta)), "figures/06_nj_tree")
save_figure(plot_tree_support(nj_bs, CONFIG$bs_threshold), "figures/06b_nj_tree_bootstrap")

if (CONFIG$domain == "adna" && file.exists("data/damage_profile.csv")) {
  message("\n=== Step 7: aDNA Damage ===")
  save_figure(plot_adna_damage(read.csv("data/damage_profile.csv")), "figures/07_adna_damage")
}

message("\n=== Pipeline complete ===")
