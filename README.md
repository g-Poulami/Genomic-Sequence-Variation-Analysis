# 🧬 Genomic Sequence Variation Analysis

---

## Overview

This repository contains tools, pipelines, and analyses developed for investigating genomic sequence variation. The work spans three major domains:

- **Bacterial genomics** — population structure, SNP calling, comparative genomics
- **Viral genomics** — within-host diversity, phylodynamics, recombination detection
- **Ancient human genomics (aDNA)** — damage assessment, population history, ancestry inference

---

## Repository Structure

```
genomics-sequence-analysis/
├── pipelines/
│   ├── bacterial/          # Bacterial variant calling & pan-genome workflows
│   ├── viral/              # Viral diversity & phylogenetic pipelines
│   └── ancient_dna/        # aDNA damage, filtering & genotyping workflows
├── scripts/
│   ├── preprocessing/      # QC, trimming, alignment wrappers
│   ├── variant_calling/    # SNP/indel calling utilities
│   ├── population_genetics/# FST, PCA, admixture analyses
│   └── visualisation/      # Plotting and figure generation
├── notebooks/              # Jupyter notebooks for exploratory analyses
├── envs/                   # Conda environment YAML files
├── config/                 # Snakemake/Nextflow config templates
├── tests/                  # Unit tests for core scripts
└── docs/                   # Extended documentation and methods notes
```

---

## Key Features

### Bacterial Genomics
- Whole-genome SNP calling against reference assemblies (BWA-MEM2 + GATK/Snippy)
- Pan-genome construction and accessory genome analysis (Roary / Panaroo)
- Recombination filtering (Gubbins / ClonalFrameML)
- Phylogenetic reconstruction with bootstrapping (IQ-TREE2 / FastTree)

### Viral Genomics
- Deep-sequencing variant calling with strand-bias and frequency filters
- Within-host diversity metrics (π, Tajima's D, dN/dS)
- Recombination screening (RDP4 / 3SEQ)
- Time-resolved phylodynamics (BEAST2 / TreeTime)

### Ancient DNA
- Damage pattern quantification (mapDamage2 / DamageProfiler)
- Low-coverage genotyping and pseudo-haploid calling
- Population structure: PCA, ADMIXTURE, f-statistics (EIGENSOFT / AdmixTools2)
- Contamination estimation (ANGSD / Schmutzi)

---

## Getting Started

### Prerequisites

- Python ≥ 3.9
- Conda / Mamba (recommended for environment management)
- Snakemake ≥ 7.0 (for automated pipelines)

### Installation

```bash
git clone https://github.com/<your-username>/genomics-sequence-analysis.git
cd genomics-sequence-analysis

# Create and activate the base environment
mamba env create -f envs/base.yaml
conda activate genomics-env
```

### Running a Pipeline

```bash
# Bacterial SNP calling pipeline
snakemake --snakefile pipelines/bacterial/Snakefile \
          --configfile config/bacterial_config.yaml \
          --cores 16

# Ancient DNA damage & genotyping
snakemake --snakefile pipelines/ancient_dna/Snakefile \
          --configfile config/adna_config.yaml \
          --cores 8
```

---

## Dependencies

Key software used across pipelines (see `envs/` for pinned versions):

| Tool | Version | Purpose |
|------|---------|---------|
| BWA-MEM2 | 2.2.1 | Short-read alignment |
| GATK | 4.x | Variant calling |
| Snippy | 4.6 | Bacterial SNP calling |
| IQ-TREE2 | 2.2 | Phylogenetic inference |
| mapDamage2 | 2.2 | aDNA damage quantification |
| ANGSD | 0.940 | Population genomics (low-coverage) |
| BEAST2 | 2.7 | Bayesian phylodynamics |
| Pandas / NumPy | latest | Data manipulation |
| Matplotlib / Seaborn | latest | Visualisation |

---

## Data

Raw sequencing data associated with published analyses are deposited at:
- **NCBI SRA**: accession numbers listed per study in `docs/data_accessions.md`
- **European Nucleotide Archive (ENA)**: linked from relevant notebooks

Example/test datasets are provided in `tests/data/` for pipeline validation.

---

## Documentation

Extended methods notes and parameter justifications are in `docs/`. Jupyter notebooks in `notebooks/` walk through key analytical decisions with inline commentary.

---

## Citation

If you use code or methods from this repository, please cite:

```
[Author(s)]. Genomic Sequence Variation Analysis. GitHub (year).
https://github.com/<your-username>/genomics-sequence-analysis
```

For specific pipelines, see citation guidance in the relevant subdirectory `README`.

---

## Contributing

Pull requests and issues are welcome. Please open an issue first to discuss proposed changes. See `CONTRIBUTING.md` for style guidelines and testing requirements.

---

## Licence

MIT — see `LICENSE` for details.
