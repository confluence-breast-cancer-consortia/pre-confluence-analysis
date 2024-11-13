#!/bin/bash --login
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --time=24:00:00
#SBATCH --array=1-9
#SBATCH --mem-per-cpu=30G

# module purge
module load R/4.3.2

anc=${1} # EUR, AFR, AMR, SAS, EAS
extract='PreExtract' # PostExtract or PreExtract
MAF='05' # 01 or 05
hm3='HM3' # HM3 or HM3_MEGA

## step1:GWAS Summary Statistics from Train Data
Rscript --slave --no-restore --no-save /spin1/home/linux/williamsjacr/Aim2_Genesis/Genesis_LD/CompileFinal_Jake.R ${SLURM_ARRAY_TASK_ID} ${anc} ${extract} ${MAF} ${hm3} > CompileFinal_Jake"${SLURM_ARRAY_TASK_ID}".Rout
