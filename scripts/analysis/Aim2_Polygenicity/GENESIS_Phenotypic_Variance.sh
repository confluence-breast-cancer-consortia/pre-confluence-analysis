#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=4:00:00
#SBATCH --array=1-8
#SBATCH --mem-per-cpu=10G

# module purge
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/Aim2_Genesis/GENESIS_Phenotypic_Variance.R ${SLURM_ARRAY_TASK_ID} > out"${SLURM_ARRAY_TASK_ID}".Rout