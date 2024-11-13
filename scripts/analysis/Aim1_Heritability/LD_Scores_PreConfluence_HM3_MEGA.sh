#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --array=1-110
#SBATCH --mem-per-cpu=100G

# module purge
module load ldsc
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/LD_Scores_PreConfluence_HM3_MEGA.R ${SLURM_ARRAY_TASK_ID} > LD_Scores_PreConfluence_HM3_MEGA"${SLURM_ARRAY_TASK_ID}".Rout