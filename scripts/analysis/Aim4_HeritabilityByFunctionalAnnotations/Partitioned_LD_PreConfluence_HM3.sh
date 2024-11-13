#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=48:00:00
#SBATCH --array=1-110
#SBATCH --mem-per-cpu=40G

# module purge
module load ldsc
module load R/4.3.2

Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/Partitioned_LD_PreConfluence_HM3.R ${SLURM_ARRAY_TASK_ID} > Partitioned_LD_PreConfluence_HM3"${SLURM_ARRAY_TASK_ID}".Rout