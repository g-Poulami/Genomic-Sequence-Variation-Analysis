# Snakemake-Genomic-Pipeline

[![Snakemake](https://img.shields.io/badge/Snakemake-%E2%89%A57.0-brightgreen?style=flat-square)](https://snakemake.readthedocs.io)
[![License](https://img.shields.io/badge/License-MIT-blue?style=flat-square)](LICENSE)


A reproducible, scalable genomic data processing pipeline built with **Snakemake**, demonstrating platform engineering best practices for large-scale short-read sequencing analysis. Covers the full journey from raw FASTQs to quality-controlled, aligned, and variant-called outputs.

---

## Overview

This pipeline showcases Snakemake as a workflow management system for genomics — using its DAG-based dependency resolution, conda environment integration, and cluster submission capabilities to handle multi-sample genomic datasets reproducibly. It is designed to be readable, modular, and straightforward to extend with additional rules.

---

## Pipeline DAG

```
Raw paired FASTQs
       |
       v
  FastQC (raw)           -- per-read quality metrics
       |
       v
  Trimmomatic            -- adapter removal, quality trimming
       |
       v
  FastQC (trimmed)       -- post-trim QC confirmation
       |
       v
  BWA index (once)
       |
  BWA MEM                -- read alignment, @RG tag embedding
       |
       v
  SAMtools sort          -- SAM to coordinate-sorted BAM
  SAMtools index         -- BAI index file
  SAMtools flagstat      -- per-sample alignment statistics
       |
       v
  GATK HaplotypeCaller   -- per-sample gVCF variant calls
       |
       v
  MultiQC                -- aggregated HTML QC report
```

---

## Key Features

- **Snakemake DAG**: automatic dependency resolution — only re-runs rules whose inputs have changed
- **Conda integration**: each rule specifies its own `conda:` environment, eliminating version conflicts between tools
- **Cluster-ready**: submit to SLURM or SGE by passing `--cluster` flags without modifying the Snakefile
- **Dry-run support**: `snakemake -n` previews execution plan before committing compute
- **Scalable**: processes multiple samples in parallel automatically
- **Reproducible**: `--use-conda` and pinned environment files ensure identical results across machines

---

## Requirements

```bash
# Install Snakemake (via mamba recommended)
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

All per-rule tool environments are managed automatically via `envs/*.yaml`.

---

## Quick Start

### Clone and configure

```bash
git clone https://github.com/g-Poulami/Snakemake-Genomic-Pipeline.git
cd Snakemake-Genomic-Pipeline
```

Edit `config/config.yaml`:

```yaml
samples:
  - SAMPLE_001
  - SAMPLE_002

reads_dir: data/raw/
genome: ref/hg38.fa
outdir: results/

trimming:
  leading: 3
  trailing: 3
  sliding_window: "4:15"
  min_len: 36
```

### Dry run

```bash
snakemake --use-conda -n
```

### Run locally

```bash
snakemake --use-conda --cores 8
```

### Run on a SLURM cluster

```bash
snakemake --use-conda \
  --cluster "sbatch --mem={resources.mem_mb}M --cpus-per-task={threads}" \
  --jobs 50
```

---

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `samples` | required | List of sample IDs |
| `reads_dir` | `data/raw/` | Directory containing paired FASTQs |
| `genome` | required | Path to reference genome FASTA |
| `outdir` | `results/` | Output directory |
| `trimming.min_len` | `36` | Minimum read length post-trimming |

---

## Outputs

| Directory | Contents |
|-----------|----------|
| `results/fastqc/` | FastQC HTML reports (raw and trimmed) |
| `results/trimmed/` | Trimmed FASTQ files |
| `results/bam/` | Sorted, indexed BAM files |
| `results/flagstat/` | Per-sample alignment rate summaries |
| `results/gvcf/` | Per-sample GATK gVCF files |
| `results/multiqc/` | `multiqc_report.html` |

---

## Project Structure

```
Snakemake-Genomic-Pipeline/
├── Snakefile                  # Main workflow definition
├── config/
│   └── config.yaml            # User-editable parameters
├── envs/
│   ├── fastqc.yaml
│   ├── trimmomatic.yaml
│   ├── bwa.yaml
│   ├── samtools.yaml
│   ├── gatk.yaml
│   └── multiqc.yaml
├── rules/
│   ├── qc.smk
│   ├── trim.smk
│   ├── align.smk
│   └── variant_call.smk
├── data/
│   └── raw/                   # Input FASTQs (not tracked)
├── ref/                       # Reference genome (not tracked)
└── .github/
    └── workflows/
        └── ci.yml
```

---

## Design Notes

**Why Snakemake over shell scripts?** Shell scripts become unmanageable with multiple samples and tools. Snakemake's DAG tracks which outputs exist and which rules need re-running, prevents redundant computation, and provides a readable, auditable record of the full analysis.

**Why per-rule conda environments?** Tools in genomics pipelines often have conflicting dependencies (e.g. different Python versions). Isolating each rule in its own conda environment avoids this entirely and makes the pipeline portable across machines without any manual installation.

---

## License

MIT

---

## Author

Poulami Ghosh — [LinkedIn](https://linkedin.com/in/poulami-ghosh-879439304) | [Google Scholar](https://scholar.google.com)
