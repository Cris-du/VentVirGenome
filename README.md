# 🧬 Supporting Materials for Manuscript: **VentVirGenome: A Global Ocean Hydrothermal Vent Virus Genome Dataset for Exploring Diversity, Function, and Host Interactions**

---

## 📖 Overview

This repository contains scripts, configurations, and auxiliary resources supporting the analyses presented in the manuscript:

> **Global diversity mapping of hydrothermal vent viruses: Virus–host interactions pervade deep-sea "oases of life"**

The workflows described here cover viral and microbial genome assembly, quality control, viral prediction, host association, phylogenetic analysis, diversity assessment, and functional annotation.

---

## 📂 Data Availability

- **VentVirGenome — Global Ocean Hydrothermal Vent Virus Genome Dataset**  
  [OEZ00021625 (BioSino)](https://www.biosino.org/node/analysis/detail/OEZ00021625)

- **VentProkGenome — Global Ocean Hydrothermal Vent Prokaryotes Genome Dataset**  
  [OEZ00021644 (BioSino)](https://www.biosino.org/node/analysis/detail/OEZ00021644)

---

## 💻 Code Availability and Workflow Summary

### 1. **[Contig assembly and binning](https://github.com/Cris-du/RDPS/blob/main/Contig_assembly_and_binning/README.md)**
Metagenomic reads were assembled into contigs and contig ≥ 1 kb were binned into MAG (metagenomic assembly genome).  
Microbial prediction, quality control and microbial taxonomic classification were performed for all bins.

### 2. **[Viral prediction and vOTU clustering](https://github.com/Cris-du/RDPS/blob/main/Viral_prediction_and_vOTU_clustering/README.md)**
Viral sequences prediction from contig ≥3 kb, and quality-checked following the pipeline described in the manuscript.  
Representative viral genomes were clustered at the species level (**vOTUs**) based on **Average Nucleotide Identity (ANI)**.

### 3. **[Viral ORF prediction and protein clustering](https://github.com/Cris-du/RDPS/blob/main/Viral_ORF_prediction_and_protein_clustering/README.md)**
Viral coding sequences were predicted and clustered into protein families for downstream comparative and functional analyses.

### 4. **[Viral taxonomic assignment](https://github.com/Cris-du/RDPS/blob/main/Viral_taxonomic_assignment/README.md)**
Viral genomes were grouped into genera-level and family-level based on **Average Amino Acid Identity (AAI)**, with taxonomic assignment according to the **International Committee on Taxonomy of Viruses (ICTV)** framework.

### 5. **[Uniqueness and cross-Dataset Comparison of GOHVGD](https://github.com/Cris-du/RDPS/blob/main/Uniqueness_and_cross-Dataset_Comparison_of_GOHVGD/README.md)**
The uniqueness and sharedness of GOHVGD at **vOTU**, **genus**, **family**, and **protein-clusters (PCs)** levels was evaluated against surface-ocean virus datasets: **Global Ocean Virome 2 (GOV2.0)** and previously published hydrothermal virus datasets (Cheng et al,. 2022, Langwig et al,. 2025). 

### 6. **[Virus–Host infective relationship prediction](https://github.com/Cris-du/RDPS/blob/main/Virus%E2%80%93Host_infective_relationship_prediction/README.md)**
Viral–host infective relationship were inferred through:
- **CRISPR-Spacer sequence matches**
- **Whole-genome sequence matches**  
between GOHVGD and GOHMGD genomes.

### 7. **[Phylogenetic Analysis](https://github.com/Cris-du/RDPS/blob/main/Phylogenetic_Analysis/README.md)**
Perform phylogenetic analysis on the hallmark proteins of **Caudoviricetes** of GOHVGD and genome of GOHMGD.

---

## 🧩 Repository Structure

