#!/bin/bash
#SBATCH --job-name=MVMR_corrections
#SBATCH --output=MVMR_corrections.txt
#SBATCH --nodes=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=100000M
#SBATCH --time=20:00:00
#SBATCH --account=sscm013902
#SBATCH --partition=mrcieu


module load languages/R/4.3.3

cd /user/work/bu19525/BristolPhD/MVMR_follow_up/Code/MVMRalcsmo/3_MVMR/
Rscript "/user/work/bu19525/BristolPhD/MVMR_follow_up/Code/MVMRalcsmo/3_MVMR/MVMR_corrections_LSI.R" MVMR_corrections.txt