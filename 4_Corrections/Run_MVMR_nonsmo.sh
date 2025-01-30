#!/bin/bash
#SBATCH --job-name=nonsmokers
#SBATCH --output=nonsmokers.txt
#SBATCH --nodes=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=100000M
#SBATCH --time=2:00:00
#SBATCH --account=sscm013902
#SBATCH --partition=mrcieu


module load languages/R/4.3.3

cd /mnt/storage/scratch/bu19525/BristolPhD/MVMR_follow_up/Code/MVMRalcsmo/4_Corrections/
Rscript "/mnt/storage/scratch/bu19525/BristolPhD/MVMR_follow_up/Code/MVMRalcsmo/4_Corrections/MR_for_adjusted_SNPs_nonsmo.R" nonsmokers.txt