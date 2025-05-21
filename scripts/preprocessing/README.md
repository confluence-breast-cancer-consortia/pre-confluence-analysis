# Pre-processing Scripts

This folder contains scripts to clean and preprocess summary statistics used in all analyses.

## Contents

- `AA_Overall_2024.R`  
  Cleans and harmonizes African American (AA) GWAS summary statistics, matches SNPs to HapMap3 and MEGA reference panels, applies allele frequency and effective sample size filters, and outputs cleaned summary statistics for downstream analysis.

- `AA_Overall_2024.sh`  
  Shell script to run `AA_Overall_2024.R` with the required environment and input files.

- `HM3_Sumstats_EAS_Latina.R`  
  Cleans and harmonizes East Asian (EAS) and Latina GWAS summary statistics, matches SNPs to HapMap3 reference, applies allele frequency and effective sample size filters, and outputs cleaned summary statistics for EAS and Latina cohorts.

- `HM3_Sumstats_EAS_Latina.sh`  
  Shell script to run `HM3_Sumstats_EAS_Latina.R` with the required environment and input files.

- `HM3_Sumstats_EUR_EAS.R`  
  Cleans and harmonizes European (EUR) and East Asian (EAS) GWAS summary statistics, matches SNPs to HapMap3 reference, applies allele frequency and effective sample size filters, and outputs cleaned summary statistics for EUR and EAS cohorts.

- `HM3_Sumstats_EUR_EAS.sh`  
  Shell script to run `HM3_Sumstats_EUR_EAS.R` with the required environment and input files.

- `HM3_MEGA_Sumstats_EAS_Latina.R`  
  Cleans and harmonizes East Asian (EAS) and Latina GWAS summary statistics, matches SNPs to the MEGA reference panel, applies allele frequency and effective sample size filters, and outputs cleaned summary statistics for EAS and Latina cohorts.

- `HM3_MEGA_Sumstats_EAS_Latina.sh`  
  Shell script to run `HM3_MEGA_Sumstats_EAS_Latina.R` with the required environment and input files.

- `HM3_MEGA_Sumstats_EUR_EAS.R`  
  Cleans and harmonizes European (EUR) and East Asian (EAS) GWAS summary statistics, matches SNPs to the MEGA reference panel, applies allele frequency and effective sample size filters, and outputs cleaned summary statistics for EUR and EAS cohorts.

- `HM3_MEGA_Sumstats_EUR_EAS.sh`  
  Shell script to run `HM3_MEGA_Sumstats_EUR_EAS.R` with the required environment and input files.

## Usage

The analysis was run on DNAnexus, to run the code locally, please follow these steps:

1. Ensure you have R installed, along with any required R packages.
2. To run a script, use the corresponding shell script. For example:
   ```bash
   bash HM3_Sumstats_EAS_Latina.sh
   ```
   Or run the R script directly:
   ```bash
   Rscript HM3_Sumstats_EAS_Latina.R
   ```

## Input/Output

- **Input:**
  - GWAS summary statistics files for various populations (AA, EAS, EUR, Latina)
  - Reference SNP lists: `w_hm3_updated.snplist`, `snpinfo_hm3_mega.csv`
- **Output:**
  - Cleaned and filtered summary statistics files for downstream analysis

## Dependencies

- R (version 4.2.0 or higher)
- R packages: `data.table`, `stringr`, `dplyr`
- bash
