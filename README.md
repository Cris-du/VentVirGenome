# 🧬 Supporting Materials for Manuscript: **VentVirGenomes: A Global Ocean Hydrothermal Vent Virus Genome Dataset for Exploring Diversity, Function, and Host Interactions**

---

## 📖 Overview

This repository contains scripts, configurations, and auxiliary resources supporting the analyses presented in the manuscript:

> **VentVirGenomes: A Global Ocean Hydrothermal Vent Virus Genome Dataset for Exploring Diversity, Function, and Host Interactions**

The workflows described here cover viral and microbial genome assembly, quality control, viral prediction, host association, phylogenetic analysis and functional annotation.

---

## 📂 Data Availability

- **VentVirGenomes — Global Ocean Hydrothermal Vent Virus Genome Dataset**  
  [VentVirGenomes](https://doi.org/10.6084/m9.figshare.32593182)

- **VentProkGenomes — Global Ocean Hydrothermal Vent Prokaryotes Genome Dataset**  
  [VentProkGenomes](https://doi.org/10.6084/m9.figshare.32593197)

---

## 💻 Code Availability and Workflow Summary

### 1. **[Metagenomic assembly and binning](./Metagenomic%20assembly%20and%20binning/README.md)**
Metagenomic reads were assembled into contigs and contig ≥ 1 kb were binned into MAG (metagenomic assembly genome).  
Microbial prediction, quality control and microbial taxonomic classification were performed for all bins.

### 2. **[Viral sequence identification and vOTU clustering](./Viral%20sequence%20identification%20and%20clustering%20clustering/README.md)**
Viral sequences prediction from contig ≥3 kb, and quality-checked following the pipeline described in the manuscript.  
Representative viral genomes were clustered at the species level (**vOTUs**) based on **Average Nucleotide Identity (ANI)**.

### 3. **[Viral protein prediction and functional annotation](./Viral%20protein%20prediction%20and%20functional%20annotation/README.md)**
Viral coding sequences were predicted and clustered into protein families, with functional annotation performed against known protein datasets for downstream comparative and functional analyses.

### 4. **[Viral taxonomic assignment](./Viral%20taxonomic%20assignment/README.md)**
Viral genomes with taxonomic assignment according to the **International Committee on Taxonomy of Viruses (ICTV)** framework.

### 5. **[Host prediction](./Host%20prediction/README.md)**
Viral–host infective relationship were inferred through:
- **CRISPR-Spacer sequence matches**
- **Whole-genome sequence matches** between VentVirGenomes and VentProkGenomes genomes.



