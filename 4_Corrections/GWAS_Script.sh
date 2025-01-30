#!/bin/bash
#SBATCH --job-name=GWAS
#SBATCH --output=GWAS.txt
#SBATCH --nodes=2
#SBATCH --cpus-per-task=1
#SBATCH --mem=100000M
#SBATCH --time=2:00:00
#SBATCH --account=sscm013902
#SBATCH --partition=mrcieu




plink --file yes_smo_file --linear --freq --ci 0.95 --covar Cov_smokers.txt --out FEV1_SNPs_GWAS_smokers