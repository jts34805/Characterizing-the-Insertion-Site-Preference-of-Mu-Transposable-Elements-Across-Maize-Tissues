#!/bin/bash
#SBATCH --job-name=HSnull
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --array=1-1000
#SBATCH --cpus-per-task=1
#SBATCH --time=72:00:00
#SBATCH --mem=64gb
#SBATCH --output=logs/HSnull_%A_%a.out
#SBATCH --error=logs/HSnull_%A_%a.err

module load R

Rscript run_hotspots_ARRAY.R
