# =============================================================================
# 04_visualisation.R
# All ggplot2 figure functions
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(patchwork); library(viridis)
  library(RColorBrewer); library(ggrepel); library(ggridges)
  library(scales); library(tidyverse)
})

theme_genomics <- function(base_size=12) {
  theme_bw(base_size=base_size) +
    theme(panel.grid.minor=element_blank(),
          strip.background=element_rect(fill="grey92", colour="grey70"))
}

plot_qual_distribution <- function(stats_df) {
  ggplot(stats_df, aes(x=QUAL)) +
    geom_histogram(bins=60, fill="#4E79A7", colour="white", linewidth=0.2) +
    scale_x_continuous(labels=comma) +
    labs(title="QUAL Score Distribution", x="QUAL", y="Number of sites") + theme_genomics()
}

plot_missingness <- function(stats_df) {
  ggplot(stats_df, aes(x=missingness)) +
    geom_histogram(bins=50, fill="#F28E2B", colour="white", linewidth=0.2) +
    scale_x_continuous(labels=percent_format()) +
    labs(title="Per-site Missingness", x="Proportion missing", y="Number of sites") + theme_genomics()
}

plot_depth_per_sample <- function(vcf) {
  dp_mat <- vcfR::extract.gt(vcf, element="DP", as.numeric=TRUE)
  as.data.frame(dp_mat) %>% rownames_to_column("site") %>%
    pivot_longer(-site, names_to="sample", values_to="depth") %>% filter(!is.na(depth)) %>%
    ggplot(aes(x=depth, y=sample, fill=after_stat(x))) +
    geom_density_ridges_gradient(scale=1.5, rel_min_height=0.01) +
    scale_fill_viridis_c(option="plasma", name="Depth") +
    labs(title="Depth per Sample", x="Read depth", y=NULL) + theme_genomics()
}

plot_snp_density <- function(stats_df, window_kb=100) {
  stats_df %>% mutate(window=floor(POS/(window_kb*1000))*(window_kb*1000)) %>%
    count(CHROM, window, name="n_snps") %>%
    ggplot(aes(x=window/1e6, y=n_snps, fill=CHROM)) +
    geom_col(width=window_kb/1000, show.legend=FALSE) +
    facet_wrap(~CHROM, scales="free_x", ncol=1) +
    scale_fill_brewer(palette="Set2") +
    labs(title=sprintf("SNP Density (window=%d kb)", window_kb), x="Position (Mb)", y="SNPs per window") +
    theme_genomics()
}

plot_afs <- function(stats_df) {
  ggplot(stats_df, aes(x=AF)) +
    geom_histogram(bins=50, fill="#59A14F", colour="white", linewidth=0.2) +
    scale_x_continuous(labels=percent_format(), limits=c(0,0.5)) +
    labs(title="Allele Frequency Spectrum", x="Minor allele frequency", y="Number of sites") + theme_genomics()
}

plot_mutation_spectrum <- function(vcf) {
  data.frame(mutation=paste0(vcf@fix[,"REF"],">",vcf@fix[,"ALT"])) %>%
    filter(nchar(vcf@fix[,"REF"])==1, nchar(vcf@fix[,"ALT"])==1) %>%
    count(mutation) %>% arrange(desc(n)) %>%
    ggplot(aes(x=reorder(mutation,-n), y=n, fill=mutation)) +
    geom_col(show.legend=FALSE) + scale_fill_brewer(palette="Set3") +
    labs(title="Mutation Spectrum", x="Substitution type", y="Count") + theme_genomics()
}

plot_pca <- function(scores, pca, pc_x=1, pc_y=2, color_by="population") {
  var_exp <- (pca$eig/sum(pca$eig))*100
  ggplot(scores, aes(x=.data[[paste0("PC",pc_x)]], y=.data[[paste0("PC",pc_y)]],
                     colour=.data[[color_by]], label=ind_id)) +
    geom_point(size=3, alpha=0.8) + geom_text_repel(size=2.5, max.overlaps=20, show.legend=FALSE) +
    scale_colour_brewer(palette="Set1") +
    labs(title="PCA", x=sprintf("PC%d (%.1f%%)",pc_x,var_exp[pc_x]),
         y=sprintf("PC%d (%.1f%%)",pc_y,var_exp[pc_y]), colour=color_by) + theme_genomics()
}

plot_pca_scree <- function(pca) {
  data.frame(PC=seq_along(pca$eig), variance=(pca$eig/sum(pca$eig))*100) %>%
    ggplot(aes(x=PC, y=variance)) + geom_col(fill="#4E79A7") + geom_line(colour="grey40") +
    geom_point(colour="grey20") +
    labs(title="PCA Scree Plot", x="Principal component", y="Variance explained (%)") + theme_genomics()
}

plot_fst_heatmap <- function(fst_mat) {
  as.data.frame(as.table(as.matrix(fst_mat))) %>% rename(pop1=Var1, pop2=Var2, Fst=Freq) %>%
    ggplot(aes(x=pop1, y=pop2, fill=Fst)) +
    geom_tile(colour="white") + geom_text(aes(label=round(Fst,3)), size=3) +
    scale_fill_viridis_c(option="magma", na.value="grey90", limits=c(0,NA), name="FST") +
    labs(title="Pairwise FST Heatmap", x=NULL, y=NULL) + theme_genomics() +
    theme(axis.text.x=element_text(angle=45, hjust=1))
}

plot_admixture <- function(q_long, meta, k_vals=NULL) {
  df <- if (!is.null(k_vals)) filter(q_long, K_value %in% k_vals) else q_long
  df <- df %>% left_join(meta, by="ind_id") %>% arrange(population, ind_id)
  df$ind_id <- factor(df$ind_id, levels=unique(df$ind_id))
  ggplot(df, aes(x=ind_id, y=proportion, fill=cluster)) +
    geom_col(width=1, position="stack") +
    facet_grid(K_value~population, scales="free_x", space="free_x") +
    scale_fill_brewer(palette="Set2") +
    labs(title="ADMIXTURE Ancestry Proportions", x=NULL, y="Ancestry proportion", fill="Cluster") +
    theme_genomics() + theme(axis.text.x=element_text(angle=90, hjust=1, size=6))
}

plot_admixture_cv <- function(cv_df) {
  best_k <- cv_df$K[which.min(cv_df$cv_error)]
  ggplot(cv_df, aes(x=K, y=cv_error)) + geom_line(colour="#4E79A7") + geom_point(size=3, colour="#4E79A7") +
    geom_vline(xintercept=best_k, linetype="dashed", colour="firebrick") +
    scale_x_continuous(breaks=cv_df$K) +
    labs(title="ADMIXTURE CV Error", x="K", y="CV error") + theme_genomics()
}

plot_adna_damage <- function(damage_df) {
  damage_df %>% pivot_longer(c(CtoT,GtoA), names_to="substitution", values_to="frequency") %>%
    ggplot(aes(x=pos, y=frequency, colour=substitution)) + geom_line(linewidth=1) + geom_point(size=2) +
    scale_colour_manual(values=c(CtoT="#E15759", GtoA="#4E79A7")) +
    scale_y_continuous(labels=percent_format()) +
    labs(title="aDNA Damage Profile", x="Position from read end (bp)",
         y="Substitution frequency", colour="Substitution") + theme_genomics()
}

plot_fragment_length <- function(frag_df) {
  ggplot(frag_df, aes(x=length)) +
    geom_histogram(bins=80, fill="#B07AA1", colour="white", linewidth=0.2) +
    labs(title="Fragment Length Distribution", x="Fragment length (bp)", y="Count") + theme_genomics()
}

save_figure <- function(p, prefix, width=10, height=7) {
  dir.create(dirname(prefix), recursive=TRUE, showWarnings=FALSE)
  ggsave(paste0(prefix,".png"), plot=p, width=width, height=height, dpi=300)
  ggsave(paste0(prefix,".pdf"), plot=p, width=width, height=height)
  message("Saved: ", prefix, ".png / .pdf")
}

combine_panels <- function(plots, ncol=NULL, ...) {
  wrap_plots(plots, ncol=ncol, ...) + plot_annotation(tag_levels="A")
}
