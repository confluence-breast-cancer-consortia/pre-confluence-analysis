#!/bin/bash
#SBATCH --time=4:00:00
#SBATCH --mem=16g

module load python/3.9		# loads in python3

# installs Popcorn (https://github.com/brielin/Popcorn)
git clone https://github.com/brielin/Popcorn.git
cd Popcorn
pip install .

POP1_REF=$1			# bed/bim/fam file prefix (before extension) for 1st population
POP2_REF=$2			# bed/bim/fam file prefix (before extension) for 2nd population
OUT=$3				# file path to write cross-ancestry LD scores

# runs Popcorn to compute the cross-ancestry scores
popcorn compute -v 1 --SNPs_to_store 100000 --bfile1 $POP1_REF --bfile2 $POP2_REF $OUT
