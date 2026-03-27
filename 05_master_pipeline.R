# =============================================================================
# 05_master_pipeline.R
# End-to-end orchestration: load → filter → analyse → visualise
# Domains: Bacterial Genomics | Ancient DNA
# =============================================================================

# --- 0. Source modules --------------------------------------------------------
source("01_vcf_snp_analysis.R")
source("02_population_genetics.R")
source("03_phylogenetics.R")
source("04_visualisation.R")

dir.create("figures", showWarnings = FALSE)
dir.create("output",  showWarnings = FALSE)


# =============================================================================
# USER CONFIG — edit these paths and parameters
# =============================================================================

CONFIG <- list(
  # Input
  vcf_path     = "data/raw_variants.vcf.gz",
  metadata     = "data/sample_metadata.csv",   # cols: ind_id, population, ...
  domain       = "bacterial",                  # "bacterial" | "adna"
  outgroup     = NULL,                          # sample ID for tree rooting

  # Filter thresholds
  min_qual     = 30,
  min_dp       = 5,       # use 3 for aDNA
  max_miss     = 0.2,     # use 0.5 for aDNA
  min_maf      = 0.01,    # use 0   for aDNA
  remove_ct_ga = FALSE,   # aDNA: TRUE to strip deamination-prone sites

  # PCA
  pca_n_comp   = 10,
  ld_window    = 50,
  ld_r2        = 0.2,

  # ADMIXTURE (set paths after running ADMIXTURE externally)
  admixture_k  = 2:5,
  admixture_q_dir  = "admixture/",
  admixture_log_dir = "admixture/logs/",

  # Phylogenetics
  n_bootstrap  = 100,
  bs_threshold = 70
)


# =============================================================================
# STEP 1 — Load VCF
# =============================================================================

message("\n=== STEP 1: Load VCF ===")
vcf_raw <- load_vcf(CONFIG$vcf_path)
meta    <- read.csv(CONFIG$metadata, stringsAsFactors = FALSE)


# =============================================================================
# STEP 2 — Filter SNPs
# =============================================================================

message("\n=== STEP 2: Filter SNPs ===")

if (CONFIG$domain == "adna") {
  vcf_filt <- filter_adna(
    vcf_raw,
    min_qual     = CONFIG$min_qual,
    min_dp       = CONFIG$min_dp,
    max_miss     = CONFIG$max_miss,
    remove_ct_ga = CONFIG$remove_ct_ga
  )
} else {
  vcf_filt <- filter_snps(
    vcf_raw,
    min_qual = CONFIG$min_qual,
    min_dp   = CONFIG$min_dp,
    max_miss = CONFIG$max_miss,
    min_maf  = CONFIG$min_maf
  )
}

write_filtered_vcf(vcf_filt, "output/filtered.vcf.gz")
stats_df <- site_stats(vcf_filt)
write.csv(stats_df, "output/site_stats.csv", row.names = FALSE)


# =============================================================================
# STEP 3 — QC visualisation
# =============================================================================

message("\n=== STEP 3: QC figures ===")

p_qual  <- plot_qual_distribution(stats_df, min_qual = CONFIG$min_qual)
p_miss  <- plot_missingness(stats_df,      max_miss = CONFIG$max_miss)
p_afs   <- plot_afs(stats_df)
p_spec  <- plot_mutation_spectrum(vcf_filt)
p_dens  <- plot_snp_density(stats_df, window_kb = 10)

p_qc <- combine_panels(list(p_qual, p_miss, p_afs, p_spec),
                        ncol = 2, title = "SNP QC overview")
save_figure(p_qc,  "figures/01_qc_panel",       width = 12, height = 10)
save_figure(p_dens, "figures/02_snp_density",   width = 14, height = 4)

ts_tv_ratio(vcf_filt)


# =============================================================================
# STEP 4 — Population genetics: PCA
# =============================================================================

message("\n=== STEP 4: PCA ===")

gl    <- vcf_to_genlight(vcf_filt)
gl_ld <- ld_prune(gl, window = CONFIG$ld_window, r2_thresh = CONFIG$ld_r2)
pca   <- run_pca(gl_ld, ncomp = CONFIG$pca_n_comp)

scores   <- pca_scores_df(pca, meta)
write.csv(scores, "output/pca_scores.csv", row.names = FALSE)

p_pca    <- plot_pca(scores, pca, colour_col = "population",
                     label_col = "ind_id")
p_scree  <- plot_pca_scree(pca)

p_pca_panel <- combine_panels(list(p_pca, p_scree),
                               ncol = 2, title = "PCA")
save_figure(p_pca_panel, "figures/03_pca", width = 14, height = 6)


# =============================================================================
# STEP 5 — Population genetics: FST
# =============================================================================

message("\n=== STEP 5: FST ===")

pop_map <- meta %>% rename(sample = ind_id)
fst_mat <- compute_fst(vcf_filt, pop_map)
write.csv(as.data.frame(as.matrix(fst_mat)), "output/fst_matrix.csv")

p_fst <- plot_fst_heatmap(fst_mat)
save_figure(p_fst, "figures/04_fst_heatmap", width = 7, height = 6)


# =============================================================================
# STEP 6 — ADMIXTURE (reads pre-computed .Q files)
# =============================================================================

message("\n=== STEP 6: ADMIXTURE ===")

k_vals   <- CONFIG$admixture_k
q_files  <- setNames(
  file.path(CONFIG$admixture_q_dir, paste0("out.", k_vals, ".Q")),
  k_vals
)
log_files <- setNames(
  file.path(CONFIG$admixture_log_dir, paste0("k", k_vals, ".log")),
  k_vals
)

existing_q   <- q_files[file.exists(q_files)]
existing_log <- log_files[file.exists(log_files)]

if (length(existing_q) > 0) {
  q_long <- read_admixture(existing_q, ind_ids = meta$ind_id)
  p_admx <- plot_admixture(q_long, meta)
  save_figure(p_admx, "figures/05_admixture_barplot", width = 14, height = 5)

  if (length(existing_log) > 0) {
    cv_df <- admixture_cv(existing_log)
    p_cv  <- plot_admixture_cv(cv_df)
    save_figure(p_cv, "figures/05b_admixture_cv", width = 6, height = 5)
  }
} else {
  message("No ADMIXTURE .Q files found — skipping.")
}


# =============================================================================
# STEP 7 — Phylogenetics
# =============================================================================

message("\n=== STEP 7: Phylogenetics ===")

gt_mat   <- get_genotype_matrix(vcf_filt)
dist_mat <- snp_distance_matrix(gt_mat)
nj_tree  <- build_nj_tree(dist_mat, outgroup = CONFIG$outgroup)
save_tree(nj_tree, "output/nj_tree")

# Bootstrap
phydat <- gt_to_phydat(gt_mat)
nj_bs  <- bootstrap_nj(phydat, nj_tree, nboot = CONFIG$n_bootstrap)

# Annotate and plot
tree_ann  <- annotate_tree(nj_bs, meta)
p_tree    <- plot_tree_basic(tree_ann, tip_col = "population")
p_tree_bs <- plot_tree_support(nj_bs, bs_threshold = CONFIG$bs_threshold)

save_figure(p_tree,    "figures/06_nj_tree",            width = 8, height = 10)
save_figure(p_tree_bs, "figures/06b_nj_tree_bootstrap", width = 8, height = 10)


# =============================================================================
# STEP 8 — aDNA-specific: damage profile
# =============================================================================

if (CONFIG$domain == "adna") {
  message("\n=== STEP 8: aDNA damage profile ===")

  dmg_path <- "data/damage_profile.csv"
  if (file.exists(dmg_path)) {
    damage_df <- read.csv(dmg_path)
    p_dmg     <- plot_adna_damage(damage_df)
    save_figure(p_dmg, "figures/07_adna_damage", width = 9, height = 5)
  } else {
    message("  damage_profile.csv not found — run mapDamage2 first.")
  }
}


# =============================================================================
# DONE
# =============================================================================

message("\n=== Pipeline complete ===")
message("Figures: ./figures/")
message("Data:    ./output/")
