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





#chose only essential cols
 awk '{ print $1,  $2, $3, $4 }' dbsnp.txt >chr_snp_pos.txt

awk '{ print $1, $3 }' chr_snp_pos.txt >RSID_SNP.txt

awk '{ print $3, $4 }' chr_snp_pos.txt >RSID_SNP_names.txt



#truncate chr:bp

#get rid of chr
awk 'sub("[^[:digit:]]+", "", $0)' RSID_SNP.txt> RSID_SNP_noltr.txt

awk 'BEGIN{OFS=":"} {print $1,$2}' RSID_SNP_noltr.txt>RSID_SNP_Chr_ID.txt
paste -d" " RSID_SNP_names.txt RSID_SNP_Chr_ID.txt > RSID_final.txt

# cut only neeeded cols

awk '{ print $2, $3 }' RSID_final.txt >RSID_merge.txt
