#!/bin/bash
#SBATCH --job-name=Run_subpop
#SBATCH --error=Run_subpop.err
#SBATCH --output=Run_subpop.out
#SBATCH --mem=60gb
#SBATCH --time=5:00:00
#SBATCH --cpus-per-task=1


Species=fprou #102293
echo "Extracting Species from data"
#bash Extract_species.sh $Species
echo "Running R script"
ml R/4.1.0-foss-2021a
Rscript Find_subpopulations.R $Species



