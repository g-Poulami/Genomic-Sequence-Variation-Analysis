# Genomic Sequence Variation Analysis — R Pipeline

R scripts for analysing SNP variation in bacterial and ancient human genomic datasets from VCF input. The pipeline covers quality filtering, population genetics, phylogenetics, and publication-quality visualisation.

---

## Scripts

| File | Purpose |
|------|---------|
| `01_vcf_snp_analysis.R` | VCF loading, SNP filtering, site statistics |
| `02_population_genetics.R` | PCA, FST, ADMIXTURE parsing |
| `03_phylogenetics.R` | Tree construction, annotation, plotting |
| `04_visualisation.R` | All ggplot2 figure functions |
| `05_master_pipeline.R` | End-to-end orchestration — start here |

---

## Requirements

### R version
R ≥ 4.2.0

### CRAN packages
```r
install.packages(c(
  "vcfR", "tidyverse", "data.table", "ape", "phangorn",
  "adegenet", "hierfstat", "ggtree", "treeio", "patchwork",
  "viridis", "RColorBrewer", "ggrepel", "ggridges", "scales", "purrr"
))
```

### Bioconductor packages
```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("VariantAnnotation", "GenomicRanges"))
```

> **Note for ancient DNA:** [mapDamage2](https://ginolhac.github.io/mapDamage/) must be run externally before using `plot_adna_damage()`. [ADMIXTURE](https://dalexander.github.io/admixture/) must also be run externally before using the ADMIXTURE plotting and CV functions.

---

## Input data

| File | Description |
|------|-------------|
| `data/raw_variants.vcf.gz` | Indexed, bgzipped VCF (multi-sample) |
| `data/sample_metadata.csv` | Sample metadata — must include `ind_id` and `population` columns |
| `data/damage_profile.csv` | mapDamage2 output (aDNA only) |
| `admixture/*.Q` | ADMIXTURE Q-matrix files, one per K value |
| `admixture/logs/k*.log` | ADMIXTURE log files containing CV error lines |

### Metadata format (`sample_metadata.csv`)

```
ind_id,population,country,year
SampleA,PopA,Germany,1200
SampleB,PopA,Germany,1350
SampleC,PopB,France,800
```

---

## Quick start

Open `05_master_pipeline.R` and edit the `CONFIG` block at the top:

```r
CONFIG <- list(
  vcf_path     = "data/raw_variants.vcf.gz",
  metadata     = "data/sample_metadata.csv",
  domain       = "bacterial",   # or "adna"
  outgroup     = NULL,          # sample ID to root the tree, or NULL

  # SNP filters
  min_qual     = 30,
  min_dp       = 5,
  max_miss     = 0.2,
  min_maf      = 0.01,
  remove_ct_ga = FALSE,   # aDNA: set TRUE to remove deamination-prone sites

  # PCA
  pca_n_comp   = 10,
  ld_window    = 50,
  ld_r2        = 0.2,

  # ADMIXTURE
  admixture_k       = 2:5,
  admixture_q_dir   = "admixture/",
  admixture_log_dir = "admixture/logs/",

  # Phylogenetics
  n_bootstrap  = 100,
  bs_threshold = 70
)
```

Then source the script:

```r
source("05_master_pipeline.R")
```

---

## Domain-specific settings

### Bacterial genomics (default)

```r
domain       = "bacterial"
min_qual     = 30
min_dp       = 5
max_miss     = 0.2
min_maf      = 0.01
remove_ct_ga = FALSE
```

### Ancient DNA

```r
domain       = "adna"
min_dp       = 3       # lower — coverage is typically sparse
max_miss     = 0.5     # more permissive — many sites will be partially missing
min_maf      = 0       # no MAF filter — pseudo-haploid calling
remove_ct_ga = FALSE   # set TRUE to exclude deamination-prone C>T / G>A sites
```

When `domain = "adna"` the pipeline automatically applies `filter_adna()` instead of `filter_snps()`, reports C>T/G>A transition counts, applies pseudo-haploid conversion before PCA, and generates the aDNA damage profile figure if `data/damage_profile.csv` is present.

---

## Outputs

### Figures (`figures/`)

| File | Content |
|------|---------|
| `01_qc_panel` | QUAL distribution, missingness, AFS, mutation spectrum |
| `02_snp_density` | SNP density along genome in sliding windows |
| `03_pca` | PCA scatter coloured by population + scree plot |
| `04_fst_heatmap` | Pairwise FST matrix heatmap |
| `05_admixture_barplot` | Ancestry proportion barplot across K values |
| `05b_admixture_cv` | Cross-validation error by K |
| `06_nj_tree` | Neighbour-joining tree coloured by population |
| `06b_nj_tree_bootstrap` | NJ tree with bootstrap support labels |
| `07_adna_damage` | C>T / G>A substitution frequency at read termini (aDNA only) |

All figures are saved as both `.png` (300 dpi) and `.pdf`.

### Data (`output/`)

| File | Content |
|------|---------|
| `filtered.vcf.gz` | Quality-filtered VCF |
| `site_stats.csv` | Per-site CHROM, POS, QUAL, depth, missingness, allele frequency |
| `pca_scores.csv` | PC scores per individual |
| `fst_matrix.csv` | Pairwise FST values |
| `nj_tree.nwk` | Neighbour-joining tree (Newick) |
| `nj_tree.nexus` | Neighbour-joining tree (Nexus) |

---

## Module reference

### `01_vcf_snp_analysis.R`

| Function | Description |
|----------|-------------|
| `load_vcf(path)` | Read VCF into vcfR object |
| `vcf_to_tidy(vcf)` | Extract fix, info, and GT fields as data frames |
| `filter_snps(vcf, ...)` | Apply QUAL / DP / missingness / MAF / biallelic filters |
| `filter_adna(vcf, ...)` | aDNA filters + deamination reporting |
| `get_genotype_matrix(vcf)` | Dosage matrix (0/1/2/NA), samples as columns |
| `site_stats(vcf)` | Per-site summary data frame |
| `ts_tv_ratio(vcf)` | Compute and print Ts/Tv ratio |
| `write_filtered_vcf(vcf, path)` | Write filtered VCF to disk |

### `02_population_genetics.R`

| Function | Description |
|----------|-------------|
| `vcf_to_genlight(vcf)` | Convert to adegenet genlight object |
| `ld_prune(gl, window, r2_thresh)` | Window-based LD pruning |
| `run_pca(gl, ncomp)` | Run glPca, print variance explained |
| `pca_scores_df(pca, meta)` | Tidy PC scores joined to metadata |
| `compute_fst(vcf, pop_map)` | Pairwise Weir & Cockerham FST via hierfstat |
| `overall_fst(vcf, pop_map)` | Overall Fst summary |
| `read_admixture(q_files, ind_ids)` | Parse .Q files into long format |
| `admixture_cv(log_files)` | Extract CV error per K from log files |
| `pseudo_haploidise(gt_matrix)` | Randomly sample one allele at het sites |

### `03_phylogenetics.R`

| Function | Description |
|----------|-------------|
| `snp_distance_matrix(gt_matrix)` | Pairwise SNP distance matrix |
| `gt_to_phydat(gt_matrix)` | Convert to phangorn phyDat (binary) |
| `build_nj_tree(dist_mat, outgroup)` | Neighbour-joining tree |
| `build_upgma_tree(dist_mat)` | UPGMA tree |
| `bootstrap_nj(phydat, tree, nboot)` | Add bootstrap support to NJ tree |
| `read_iqtree / read_beast / read_newick / read_nexus` | Import external trees |
| `annotate_tree(tree, meta)` | Join metadata onto treedata object |
| `plot_tree_basic(tree, ...)` | Rectangular tree coloured by population |
| `plot_tree_circular(tree, meta, heat_cols)` | Circular tree with metadata heatmap |
| `plot_tree_support(tree, bs_threshold)` | Tree with bootstrap labels |
| `save_tree(tree, prefix)` | Write Newick + Nexus |

### `04_visualisation.R`

| Function | Description |
|----------|-------------|
| `plot_qual_distribution(stats_df)` | QUAL score histogram |
| `plot_missingness(stats_df)` | Per-site missingness histogram |
| `plot_depth_per_sample(vcf)` | Depth ridgeline plot per sample |
| `plot_snp_density(stats_df, window_kb)` | Sliding-window SNP density |
| `plot_afs(stats_df)` | Allele frequency spectrum |
| `plot_mutation_spectrum(vcf)` | Substitution type bar chart |
| `plot_pca(scores, pca, ...)` | PCA scatter with variance labels |
| `plot_pca_scree(pca)` | Scree plot |
| `plot_fst_heatmap(fst_mat)` | Annotated FST heatmap |
| `plot_admixture(q_long, meta, k_vals)` | Stacked ancestry barplot |
| `plot_admixture_cv(cv_df)` | CV error line plot |
| `plot_adna_damage(damage_df)` | C>T / G>A damage profile |
| `plot_fragment_length(frag_df)` | Fragment length distribution |
| `save_figure(p, prefix, ...)` | Save PNG + PDF |
| `combine_panels(plots, ...)` | Multi-panel patchwork layout |

---

## Directory layout

```
project/
├── 01_vcf_snp_analysis.R
├── 02_population_genetics.R
├── 03_phylogenetics.R
├── 04_visualisation.R
├── 05_master_pipeline.R
├── data/
│   ├── raw_variants.vcf.gz
│   ├── raw_variants.vcf.gz.tbi
│   ├── sample_metadata.csv
│   └── damage_profile.csv        # aDNA only
├── admixture/
│   ├── out.2.Q
│   ├── out.3.Q
│   └── logs/
│       ├── k2.log
│       └── k3.log
├── figures/                      # created automatically
└── output/                       # created automatically
```

---

## Notes

**LD pruning before PCA and ADMIXTURE** is applied by default (`ld_window = 50`, `ld_r2 = 0.2`). Disable by setting `ld_r2 = 1` in CONFIG.

**ADMIXTURE** must be run externally on the filtered VCF before steps 6 in the pipeline can execute. A typical command line would be:
```bash
for K in 2 3 4 5; do
  admixture --cv filtered.vcf.gz $K | tee logs/k${K}.log
done
```

**External phylogenetic software** (IQ-TREE2, FastTree, BEAST2) can be integrated by passing their output tree files to `read_iqtree()`, `read_newick()`, or `read_beast()` and then calling the annotation and plotting functions from `03_phylogenetics.R` as normal.

**Reproducibility:** Set a random seed at the top of your session before running PCA or bootstrap steps:
```r
set.seed(42)
```
