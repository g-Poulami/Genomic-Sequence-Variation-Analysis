# =============================================================================
# 04_visualisation.R
# Publication-quality figures for Genomic Sequence Variation Analysis
# Domains: Bacterial Genomics | Ancient DNA
# =============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)     # compose multi-panel figures
library(viridis)
library(RColorBrewer)
library(ggrepel)       # non-overlapping labels
library(scales)

# Consistent theme across all figures
theme_genomics <- function(base_size = 12) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.background  = element_rect(fill = "grey92", colour = NA),
      legend.position   = "right",
      plot.title        = element_text(face = "bold", size = base_size + 1),
      axis.title        = element_text(size = base_size),
      legend.title      = element_text(size = base_size - 1),
      legend.text       = element_text(size = base_size - 2)
    )
}

PALETTE_POP <- brewer.pal(8, "Set2")   # population colours
PALETTE_DIV <- brewer.pal(9, "RdBu")  # diverging (FST, etc.)


# =============================================================================
# A. VARIANT-LEVEL PLOTS
# =============================================================================

# --- A1. SNP quality distribution (QUAL scores) --------------------------------

plot_qual_distribution <- function(stats_df, min_qual = 30) {
  ggplot(stats_df, aes(x = QUAL)) +
    geom_histogram(bins = 60, fill = "#4292c6", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = min_qual, colour = "firebrick",
               linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = min_qual + 2, y = Inf,
             label = paste0("Threshold: ", min_qual),
             hjust = 0, vjust = 1.5, colour = "firebrick", size = 3.5) +
    scale_x_continuous(limits = c(0, NA)) +
    labs(title = "SNP quality score distribution",
         x = "QUAL score", y = "Count") +
    theme_genomics()
}


# --- A2. Missingness per site --------------------------------------------------

plot_missingness <- function(stats_df, max_miss = 0.2) {
  ggplot(stats_df, aes(x = prop_miss)) +
    geom_histogram(bins = 50, fill = "#74c476", colour = "white", linewidth = 0.2) +
    geom_vline(xintercept = max_miss, colour = "firebrick",
               linetype = "dashed", linewidth = 0.8) +
    scale_x_continuous(labels = percent_format()) +
    labs(title = "Per-site missingness",
         x = "Proportion missing", y = "Count") +
    theme_genomics()
}


# --- A3. Depth distribution per sample ----------------------------------------

plot_depth_per_sample <- function(vcf) {
  dp <- extract.gt(vcf, element = "DP", as.numeric = TRUE) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("site") %>%
    pivot_longer(-site, names_to = "sample", values_to = "depth") %>%
    filter(!is.na(depth))

  ggplot(dp, aes(x = depth, y = sample, fill = sample)) +
    ggridges::geom_density_ridges(scale = 1.5, alpha = 0.75,
                                  show.legend = FALSE) +
    scale_x_continuous(limits = c(0, quantile(dp$depth, 0.99, na.rm = TRUE))) +
    scale_fill_viridis_d(option = "C") +
    labs(title = "Read depth distribution per sample",
         x = "Depth (DP)", y = NULL) +
    theme_genomics()
}


# --- A4. SNP density along genome (Manhattan-style) ---------------------------

plot_snp_density <- function(stats_df,
                              window_kb = 10,
                              highlight_chrom = NULL) {
  df <- stats_df %>%
    mutate(window = floor(POS / (window_kb * 1e3))) %>%
    count(CHROM, window) %>%
    mutate(pos_mid = (window + 0.5) * window_kb * 1e3 / 1e6)

  p <- ggplot(df, aes(x = pos_mid, y = n, colour = CHROM)) +
    geom_line(linewidth = 0.6, alpha = 0.8) +
    facet_wrap(~CHROM, scales = "free_x", ncol = 1) +
    scale_colour_viridis_d(guide = "none") +
    labs(title   = paste0("SNP density (", window_kb, " kb windows)"),
         x = "Position (Mb)", y = "SNP count") +
    theme_genomics() +
    theme(strip.text = element_text(face = "bold"))

  if (!is.null(highlight_chrom)) {
    p <- p + geom_rect(
      data = df %>% filter(CHROM == highlight_chrom),
      aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
      fill = "yellow", alpha = 0.05, inherit.aes = FALSE
    )
  }
  p
}


# --- A5. Allele frequency spectrum (AFS / SFS) --------------------------------

plot_afs <- function(stats_df, n_bins = 20) {
  ggplot(stats_df %>% filter(!is.na(alt_freq)),
         aes(x = alt_freq)) +
    geom_histogram(bins = n_bins, fill = "#9e9ac8",
                   colour = "white", linewidth = 0.2) +
    scale_x_continuous(labels = percent_format(), limits = c(0, 1)) +
    labs(title = "Allele frequency spectrum",
         x = "Alternative allele frequency", y = "SNP count") +
    theme_genomics()
}


# --- A6. Transition / Transversion composition (aDNA damage context) ----------

plot_mutation_spectrum <- function(vcf) {
  fix <- as.data.frame(vcf@fix, stringsAsFactors = FALSE)
  fix <- fix %>%
    filter(nchar(REF) == 1 & nchar(ALT) == 1) %>%
    mutate(change = paste0(REF, ">", ALT)) %>%
    count(change) %>%
    mutate(
      class = ifelse(
        change %in% c("C>T","G>A","A>G","T>C"), "Transition", "Transversion"
      ),
      deamination = change %in% c("C>T","G>A")
    )

  ggplot(fix, aes(x = reorder(change, -n), y = n,
                  fill = class, alpha = deamination)) +
    geom_col(colour = "white", linewidth = 0.3) +
    scale_fill_manual(values = c("Transition" = "#2171b5",
                                 "Transversion" = "#cb181d"),
                      name = "Type") +
    scale_alpha_manual(values = c("FALSE" = 0.6, "TRUE" = 1.0),
                       labels = c("No", "Yes"),
                       name   = "Deamination\n(C>T / G>A)") +
    labs(title = "Mutation spectrum",
         x = "Substitution type", y = "Count") +
    theme_genomics() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}


# =============================================================================
# B. POPULATION GENETICS PLOTS
# =============================================================================

# --- B1. PCA scatter ----------------------------------------------------------

#' @param scores_df  Output of pca_scores_df() — must contain PC1, PC2, population
#' @param pca        glPca object (for variance explained axis labels)
#' @param label_col  Column to label points with (NULL = no labels)
plot_pca <- function(scores_df, pca, colour_col = "population",
                     shape_col = NULL, label_col = NULL,
                     pc_x = 1, pc_y = 2) {

  var_pct <- round(100 * pca$eig / sum(pca$eig), 1)
  xlab    <- sprintf("PC%d (%.1f%%)", pc_x, var_pct[pc_x])
  ylab    <- sprintf("PC%d (%.1f%%)", pc_y, var_pct[pc_y])

  xcol <- paste0("PC", pc_x)
  ycol <- paste0("PC", pc_y)

  aes_base <- aes(x = .data[[xcol]], y = .data[[ycol]],
                  colour = .data[[colour_col]])
  if (!is.null(shape_col))
    aes_base <- modifyList(aes_base, aes(shape = .data[[shape_col]]))

  p <- ggplot(scores_df, aes_base) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
    geom_point(size = 3, alpha = 0.85) +
    scale_colour_brewer(palette = "Set2", name = colour_col) +
    labs(title = "Principal Component Analysis",
         x = xlab, y = ylab) +
    theme_genomics() +
    coord_equal()

  if (!is.null(label_col)) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = .data[[label_col]]),
      size = 2.8, max.overlaps = 20, colour = "grey30"
    )
  }
  p
}


# --- B2. PCA scree plot -------------------------------------------------------

plot_pca_scree <- function(pca, n_show = 15) {
  var_pct <- 100 * pca$eig / sum(pca$eig)
  df <- data.frame(PC = seq_along(var_pct), var = var_pct) %>%
    slice(1:min(n_show, n()))

  ggplot(df, aes(x = PC, y = var)) +
    geom_col(fill = "#6baed6", width = 0.7) +
    geom_line(colour = "grey40", linewidth = 0.6) +
    geom_point(size = 2, colour = "grey30") +
    scale_x_continuous(breaks = df$PC) +
    labs(title = "PCA scree plot",
         x = "Principal component", y = "Variance explained (%)") +
    theme_genomics()
}


# --- B3. Pairwise FST heatmap -------------------------------------------------

plot_fst_heatmap <- function(fst_mat) {
  fst_df <- as.data.frame(as.matrix(fst_mat)) %>%
    tibble::rownames_to_column("pop1") %>%
    pivot_longer(-pop1, names_to = "pop2", values_to = "fst")

  ggplot(fst_df, aes(x = pop1, y = pop2, fill = fst)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    geom_text(aes(label = round(fst, 3)), size = 3.2, colour = "grey20") +
    scale_fill_gradientn(
      colours  = rev(PALETTE_DIV),
      na.value = "grey90",
      limits   = c(0, max(fst_df$fst, na.rm = TRUE)),
      name     = expression(italic(F)[ST])
    ) +
    labs(title = expression(paste("Pairwise ", italic(F)[ST], " (Weir & Cockerham)")),
         x = NULL, y = NULL) +
    theme_genomics() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.grid  = element_blank())
}


# --- B4. ADMIXTURE barplot ----------------------------------------------------

#' @param q_long    Long data.frame from read_admixture()
#' @param meta      data.frame with ind_id and population (for sorting / faceting)
#' @param k_vals    Integer vector of K values to plot (default: all)
plot_admixture <- function(q_long, meta, k_vals = NULL) {
  if (!is.null(k_vals)) q_long <- q_long %>% filter(K %in% k_vals)

  q_long <- q_long %>%
    left_join(meta, by = "ind_id") %>%
    arrange(population, ind_id)

  q_long$ind_id <- factor(q_long$ind_id,
                           levels = unique(q_long$ind_id))

  # Cluster colour palette (max 12 clusters)
  n_clust <- q_long %>% group_by(K) %>%
    summarise(n = n_distinct(cluster)) %>% pull(n) %>% max()
  clust_pal <- colorRampPalette(brewer.pal(min(n_clust, 9), "Set1"))(n_clust)

  ggplot(q_long, aes(x = ind_id, y = proportion, fill = cluster)) +
    geom_col(width = 1, colour = NA) +
    facet_grid(K ~ population, scales = "free_x", space = "free_x",
               labeller = labeller(K = function(x) paste0("K = ", x))) +
    scale_fill_manual(values = clust_pal, guide = "none") +
    scale_y_continuous(expand = c(0, 0)) +
    labs(title = "ADMIXTURE ancestry proportions",
         x = NULL, y = "Ancestry proportion") +
    theme_genomics() +
    theme(
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      panel.spacing.x  = unit(0.5, "mm"),
      strip.text.y     = element_text(angle = 0),
      strip.text.x     = element_text(size = 8, face = "bold")
    )
}


# --- B5. ADMIXTURE cross-validation plot --------------------------------------

plot_admixture_cv <- function(cv_df) {
  best_k <- cv_df$K[which.min(cv_df$CV_error)]

  ggplot(cv_df, aes(x = K, y = CV_error)) +
    geom_line(colour = "grey50", linewidth = 0.8) +
    geom_point(aes(colour = K == best_k), size = 3.5) +
    geom_vline(xintercept = best_k, colour = "firebrick",
               linetype = "dashed", linewidth = 0.7) +
    annotate("text", x = best_k + 0.2, y = max(cv_df$CV_error),
             label = paste0("Best K = ", best_k),
             hjust = 0, colour = "firebrick", size = 3.5) +
    scale_colour_manual(values = c("FALSE" = "steelblue", "TRUE" = "firebrick"),
                        guide = "none") +
    scale_x_continuous(breaks = cv_df$K) +
    labs(title = "ADMIXTURE cross-validation error",
         x = "K (number of ancestral populations)",
         y = "CV error") +
    theme_genomics()
}


# =============================================================================
# C. ANCIENT DNA PLOTS
# =============================================================================

# --- C1. mapDamage-style damage plot -----------------------------------------

#' Plot C>T substitution frequency at read termini (deamination signature)
#' @param damage_df  data.frame with columns: position, ct_freq, ga_freq, end
#'   end: "5prime" | "3prime"
plot_adna_damage <- function(damage_df) {
  df_long <- damage_df %>%
    pivot_longer(c(ct_freq, ga_freq),
                 names_to  = "type",
                 values_to = "frequency") %>%
    mutate(type = recode(type,
                         "ct_freq" = "C>T (5′ end)",
                         "ga_freq" = "G>A (3′ end)"))

  ggplot(df_long, aes(x = position, y = frequency, colour = type)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~end, scales = "free_x",
               labeller = labeller(end = c("5prime" = "5′ end",
                                           "3prime" = "3′ end"))) +
    scale_colour_manual(
      values = c("C>T (5′ end)" = "#cb181d", "G>A (3′ end)" = "#2171b5"),
      name   = "Substitution"
    ) +
    scale_y_continuous(labels = percent_format(), limits = c(0, NA)) +
    labs(title = "Post-mortem DNA damage profile",
         subtitle = "C>T / G>A substitution frequency at read termini",
         x = "Position from read end (bp)",
         y = "Substitution frequency") +
    theme_genomics()
}


# --- C2. Fragment length distribution (aDNA characteristic short fragments) ---

#' @param frag_df  data.frame with columns: length, count (or sample)
plot_fragment_length <- function(frag_df, expected_peak = 60) {
  ggplot(frag_df, aes(x = length, y = count)) +
    geom_col(fill = "#74c476", colour = "white", linewidth = 0.1) +
    geom_vline(xintercept = expected_peak, colour = "firebrick",
               linetype = "dashed", linewidth = 0.8) +
    annotate("text", x = expected_peak + 2, y = Inf,
             label = paste0("Mode: ", expected_peak, " bp"),
             hjust = 0, vjust = 1.5, colour = "firebrick", size = 3.5) +
    scale_x_continuous(limits = c(0, 300)) +
    labs(title = "Fragment length distribution",
         subtitle = "Characteristic short fragments expected in aDNA",
         x = "Fragment length (bp)", y = "Count") +
    theme_genomics()
}


# =============================================================================
# D. SAVE UTILITIES
# =============================================================================

#' Save a ggplot to file (PNG + PDF)
save_figure <- function(p, prefix, width = 8, height = 6, dpi = 300) {
  ggsave(paste0(prefix, ".png"), p, width = width, height = height, dpi = dpi)
  ggsave(paste0(prefix, ".pdf"), p, width = width, height = height)
  message("Saved: ", prefix, ".png / .pdf")
}


#' Combine multiple ggplots into a labelled panel (patchwork)
combine_panels <- function(plots, labels = LETTERS, ncol = 2,
                           title = NULL) {
  pw <- patchwork::wrap_plots(plots, ncol = ncol) +
    patchwork::plot_annotation(
      title = title,
      tag_levels = list(labels[seq_along(plots)])
    )
  pw
}


# =============================================================================
# USAGE EXAMPLE
# =============================================================================
# source("01_vcf_snp_analysis.R")
# source("02_population_genetics.R")
#
# vcf      <- load_vcf("data/filtered_bacterial.vcf.gz")
# stats_df <- site_stats(vcf)
# meta     <- read.csv("data/sample_metadata.csv")
#
# ## Variant-level plots
# p1 <- plot_qual_distribution(stats_df)
# p2 <- plot_missingness(stats_df)
# p3 <- plot_afs(stats_df)
# p4 <- plot_mutation_spectrum(vcf)
# p_qc <- combine_panels(list(p1, p2, p3, p4), ncol = 2, title = "QC overview")
# save_figure(p_qc, "figures/qc_panel", width = 12, height = 10)
#
# ## PCA
# gl     <- vcf_to_genlight(vcf)
# pca    <- run_pca(gl, ncomp = 10)
# scores <- pca_scores_df(pca, meta)
# p_pca  <- plot_pca(scores, pca, colour_col = "population", label_col = "ind_id")
# save_figure(p_pca, "figures/pca", width = 7, height = 6)
#
# ## FST heatmap
# pop_map <- meta %>% rename(sample = ind_id)
# fst_mat <- compute_fst(vcf, pop_map)
# p_fst   <- plot_fst_heatmap(fst_mat)
# save_figure(p_fst, "figures/fst_heatmap")
#
# ## ADMIXTURE
# q_long <- read_admixture(q_files, ind_ids = meta$ind_id)
# cv_df  <- admixture_cv(log_files)
# p_cv   <- plot_admixture_cv(cv_df)
# p_admx <- plot_admixture(q_long, meta, k_vals = c(2, 3, 4))
# save_figure(p_admx, "figures/admixture_barplot", width = 14, height = 6)
#
# ## aDNA damage
# damage_df <- read.csv("data/damage_profile.csv")  # from mapDamage2 output
# p_dmg     <- plot_adna_damage(damage_df)
# save_figure(p_dmg, "figures/adna_damage")
