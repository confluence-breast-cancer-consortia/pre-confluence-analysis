#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=96:00:00
#SBATCH --cpus-per-task=30
#SBATCH --mem-per-cpu=2G


# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/Aim2_Genesis/GENESIS_Code_AFR_AABCG.R > GENESIS_Code_AFR_AABCG.Rout