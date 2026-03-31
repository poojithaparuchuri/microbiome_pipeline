# Microbiome Pipeline — Celiac Disease Biomarker Discovery

A cloud-based 16S rRNA amplicon sequencing analysis pipeline for microbiome 
biomarker discovery in Celiac disease, built on Google Cloud Platform.

## Overview

This pipeline analyses gut microbiome data from infants at familial risk of 
Celiac disease (PRJEB23313, Olivares et al. 2018, *Microbiome*). It identifies 
microbial biomarkers and their associations with immune dysfunction.

## Platform

- **Infrastructure:** Google Cloud Platform (e2-standard-4, Ubuntu 22.04)
- **Pipeline:** QIIME2 v2024.10, DADA2, Silva 138 classifier
- **Analysis:** Python (pandas, matplotlib, seaborn, scipy)
- **Version control:** Git/GitHub

## Dataset

- **Study:** PRJEB23313 - Gut microbiota trajectory in early life predicts Celiac disease
- **Samples:** 20 (10 Celiac cases, 10 healthy controls)
- **Type:** 16S rRNA V1-V2, Illumina MiSeq 2x250bp
- **Tissue:** Stool, 6 months of age

## Key Findings

| Finding | Significance |
|---------|-------------|
| Christensenellaceae R-7 group depleted in Celiac cases | q < 0.05, LFC = 4.3 |
| Acinetobacter correlates with IL-6 elevation | Spearman r = 0.8, p < 0.05 |
| Streptococcus negatively correlates with sIgA | p < 0.05 |
| Bray-Curtis PCoA shows case/control separation | Visible clustering |

## Figures

| Figure | Description |
|--------|-------------|
| figure1_shannon_diversity.png | Alpha diversity - Shannon index |
| figure2_bray_curtis_pcoa.png | Beta diversity - Bray-Curtis PCoA |
| figure3_taxonomy_barplot.png | Taxonomic composition top 10 genera |
| figure4_differential_abundance.png | ANCOM-BC biomarker discovery |
| figure5_immune_correlation.png | Microbiome - immune marker correlations |

## Methods

Raw FASTQ files downloaded from ENA (PRJEB23313). Quality control and 
denoising with DADA2 (trunc-len-f 200, trunc-len-r 180). Taxonomic 
classification against Silva 138 99% database. Alpha diversity (Shannon), 
beta diversity (Bray-Curtis PCoA), differential abundance (ANCOM-BC, q < 0.05), 
and Spearman correlations with immune markers (sIgA, TNF-α, IL-6).

## Repository Structure
```
microbiome-pipeline/
├── data/
│   └── metadata/          # Sample metadata and immune markers
├── scripts/
│   ├── download_samples.sh
│   ├── generate_figures.py
│   ├── figure4_differential_abundance.py
│   └── figure5_immune_correlation.py
└── results/
    └── figures/           # 5 publication-quality figures
```

## Roadmap

- [ ] Multi-study integration (466 samples, 6 published studies)
- [ ] Batch correction across sequencing platforms
- [ ] ML biomarker discovery (random forest, LASSO)
- [ ] PICRUSt2 functional pathway prediction
- [ ] AI-powered clinical interpretation layer
- [ ] R Shiny interactive dashboard

## Reference

Olivares et al. (2018) Gut microbiota trajectory in early life may predict 
development of celiac disease. *Microbiome* 6:36.
