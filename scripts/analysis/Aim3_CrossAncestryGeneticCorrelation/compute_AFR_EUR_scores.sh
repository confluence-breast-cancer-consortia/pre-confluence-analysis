#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mail-user=jahagirdarob@nih.gov
#SBATCH --mail-type=ALL
#SBATCH --mem=16g

# loads in python and installs popcorn
module load python/3.9
cd /data/Dutta_lab/Popcorn
# pip install .

# populations are numbered alphabetically
POP1=AFR
POP2=EUR

# precomputes scores using popcorn
POP1_REF=/data/Dutta_lab/REF/$POP1
POP2_REF=/data/Dutta_lab/REF/$POP2
OUT=/data/Dutta_lab/Popcorn/${POP1}_${POP2}_scores.txt

popcorn compute -v 1 --SNPs_to_store 100000 --bfile1 $POP1_REF --bfile2 $POP2_REF $OUT
