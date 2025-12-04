# Genomic Analysis of 423,881 UK Biobank Participants Identifies Telomere Maintenance Dysfunction as a Key Driver of Multiple Primary Cancer Susceptibility

## Repository Overview
This repository contains all analysis code and scripts for the study investigating multiple primary cancer (MPC) susceptibility using UK Biobank data. The research focuses on identifying genetic factors, particularly telomere maintenance dysfunction, contributing to MPC risk.

## Project Structure

### 1. Data Processing Pipeline

#### 1.1 Population-Specific Processing
- `1.1.dataProcess-whitePeople.ipynb` - White population data processing
- `1.2.dataProcess-MPC.ipynb` - Main MPC data processing pipeline
- `1.2.1.dataProcess-otherMPC.Rmd` - R Markdown for non-primary MPC data processing
- `1.2.2.dataProcess-singleCancer.ipynb` - Single cancer case processing
- `1.2.3.dataProcess-3cancersMPC.ipynb` - 3-cancer MPC cases processing
- `1.2.4.dataProcess-4cancersMPC.ipynb` - 4-cancer MPC cases processing

#### 1.2 Genotype Data Extraction
- `1.3.1.dataProcess-3cancersMPCGenodataExtract-CommonSNP.sh` - Common SNP extraction for 3-cancer MPC
- `1.3.1.dataProcess-3cancersMPCGenodataExtract.sh` - Full genotype extraction for 3-cancer MPC
- `1.3.2.dataProcess-4cancersMPCGenodataExtract-CommonSNP.sh` - Common SNP extraction for 4-cancer MPC
- `1.3.dataProcess-mpcGenodataExtract.sh` - Main MPC genotype extraction pipeline

#### 1.3 GRM and Covariate Processing
- `1.4.dataProcess-mpcGRMcopy.sh` - Copy GRM files
- `1.4.dataProcess-mpcGRM.sh` - Generate GRM for MPC analysis
- `1.5.1.dataProcess-mpcCovar.ipynb` - Covariate processing for MPC
- `1.5.2.dataProcess-mpcQcovar.ipynb` - Quantitative covariate processing
- `1.5.3.dataProcess-mpcPheno.ipynb` - Phenotype data processing
- `1.5.4.dataProcess-otherMPCPheno.sh` - Phenotype processing for other MPC

### 2. GWAS Analysis

#### 2.1 Disease-Specific GWAS
- `2.1.GWAS-MPC.sh` - Main MPC GWAS pipeline
- `2.1.1.GWAS-LungDigestive.sh` - Lung and digestive system cancers GWAS
- `2.1.2.GWAS-otherMPC.sh` - Other MPC types GWAS
- `2.1.3.GWAS-3cancersMPC.sh` - 3-cancer MPC GWAS
- `2.1.4.GWAS-4cancersMPC.sh` - 4-cancer MPC GWAS
- `2.1.5.GWAS-5cancersMPC.sh` - 5-cancer MPC GWAS

### 3. Results Processing and Visualization

#### 3.1 GWAS Results Processing
- `3.1.dataProcess-GWAResults.ipynb` - Process GWAS results
- `3.2.Rplot-manhattaan.Rmd` - Manhattan plot generation (R Markdown)
- `3.3.MPC3CommomLoci5e_6_Manhattaan.R` - MPC3 common loci Manhattan plots (p<5e-6)
- `3.4.MPC3AllSites_MPC4AllSites_Manhattaan.R` - All sites Manhattan plots for MPC3 and MPC4
- `3.5.SpecificMPC_5e_8_Manhattan.R` - Specific MPC Manhattan plots (p<5e-8)
- `3.6.TWAS_result_plot.R` - TWAS results visualization

#### 3.2 Pattern Analysis
- `4.MPC3_MPC4_find_regular.ipynb` - Identify patterns in MPC3 and MPC4
- `MPC3_MPC4.ipynb` - Additional analysis for MPC3 and MPC4

### 4. Statistical and Epidemiological Analysis

- `5.First_second_MPC_heatmap_dotplot_and_diagnose_time_interval_dotplot.R` - Temporal analysis of MPC diagnosis intervals
- `Epidemiological-analysis.R` - Epidemiological analysis
- `MPC-recurrence-predict.R` - MPC recurrence prediction models
- `MPC-Stage-sample-size-statistic.R` - Sample size statistics by cancer stage
- `Three-line-table-computate.R` - Generate three-line tables for publication

### 5. Specialized Analysis Folders

#### CaseCaseGWAS/
Contains scripts for case-case GWAS analysis comparing different cancer groups.

#### Gene-basedGWAS/
Contains gene-based association analysis scripts and results.

#### GloblePattern/
Global pattern analysis across different populations or cancer types.

#### GWASSurvival/
Integration of GWAS results with survival analysis.

#### SingleCancerGWAS/
GWAS analysis for single primary cancers as controls.

#### TWAS/
Transcriptome-wide association study analysis pipeline and results.

## Prerequisites and Requirements

### Software Dependencies
- PLINK 1.9/2.0 for GWAS analysis
- R (≥ 4.0.0) with tidyverse, ggplot2, data.table packages
- Python (≥ 3.8) with pandas, numpy, matplotlib, jupyter
- BASH shell for pipeline execution

### Data Requirements
- UK Biobank genotype data (requires proper access authorization)
- Phenotype data from UK Biobank
- Reference genomes and annotation files
