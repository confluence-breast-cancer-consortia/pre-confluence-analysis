#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=5G


# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/Aim2_Genesis/GENESIS_Code_EUR_Full1KG.R > GENESIS_Code_EUR_Full1KG.Rout