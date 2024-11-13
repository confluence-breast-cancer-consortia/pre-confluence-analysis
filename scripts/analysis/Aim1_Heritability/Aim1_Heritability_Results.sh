#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=20G

# module purge
module load ldsc
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/Aim1_Heritability_Results.R > Aim1_Heritability_Results.Rout