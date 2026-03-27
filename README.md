# Microbiome Pipeline - Celiac Disease Analysis

A 16S rRNA amplicon sequencing analysis pipeline for microbiome biomarker discovery in Celiac disease.

## Overview

This pipeline analyses gut microbiome data from infants at familial risk of Celiac disease (PRJEB23313, Olivares et al. 2018, *Microbiome*). It reproduces and extends the published analysis using QIIME2 and Python.

## Dataset

- **Study:** PRJEB23313 — Gut microbiota trajectory in early life may predict development of celiac disease
- **Samples:** 20 (10 Celiac cases, 10 healthy controls)
- **Type:** 16S rRNA V1-V2 amplicon sequencing (Illumina MiSeq, 2x250bp)
- **Tissue:** Stool, collected at 6 months of age

## Key Findings

- Bray-Curtis PCoA shows partial separation between Celiac cases and controls
- Shannon diversity similar between groups, consistent with published results
- Bifidobacterium dominant genus in both groups, with compositional differences

## Methods

Raw FASTQ files were downloaded from the European Nucleotide Archive. Quality control and denoising were performed using DADA2 (v1.30) within QIIME2 (v2024.10), with forward reads truncated at 200bp and reverse reads at 180bp. Taxonomic classification was performed using a Naive Bayes classifier trained on the Silva 138 99% OTU reference database. Alpha diversity was assessed using the Shannon diversity index. Beta diversity was assessed using Bray-Curtis dissimilarity with PCoA ordination. Figures were generated using Python (matplotlib, seaborn).

## Results

| Figure | Description |
|--------|-------------|
| figure1_shannon_diversity.png | Alpha diversity comparison |
| figure2_bray_curtis_pcoa.png | Beta diversity PCoA ordination |
| figure3_taxonomy_barplot.png | Taxonomic composition top 10 genera |

## Repository Structure
```
microbiome-pipeline/
├── data/
│   └── metadata/          # Sample metadata and manifest
├── scripts/
│   ├── download_samples.sh
│   └── generate_figures.py
└── results/
    └── figures/           # Publication-quality figures
```

## Reference

Olivares et al. (2018) Gut microbiota trajectory in early life may predict development of celiac disease. *Microbiome* 6:36.
