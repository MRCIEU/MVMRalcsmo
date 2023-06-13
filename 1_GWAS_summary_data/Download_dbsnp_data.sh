#!/bin/bash
#SBATCH --job-name=dbsnp
#SBATCH --output=dbsnp.txt
#SBATCH --nodes=4
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem=250000M
#SBATCH --time=8:00:00
#SBATCH --account=sscm013902
#SBATCH --partition=mrcieu



rsync -aP hgdownload.soe.ucsc.edu::gbdb/hg19/snp/dbSnp153.bb ./
