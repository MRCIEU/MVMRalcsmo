#libraries and local clumpting paths

library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)

library(readr)
library(data.table)

#Smoking data read and clump
SmoInit<-read.table(gzfile("/mnt/storage/scratch/bu19525/BristolPhD/MVMR_follow_up/Clumping/SmokingInitiation.WithoutUKB.txt.gz"),sep="\t",header=TRUE)
head(SmoInit)
summary(SmoInit$EFFECTIVE_N)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
FORMATSmoInit<-format_data(SmoInit, type = "exposure", snp_col = "RSID",
                             beta_col = "BETA",
                             se_col = "SE",
                             eaf_col = "AF",
                             effect_allele_col = "ALT",
                             other_allele_col = "REF",
                             pval_col = "PVALUE",
                             samplesize_col = "EFFECTIVE_N",
                             chr_col = "CHROM",
                             pos_col = "POS"
)

  Pvallimit<-0.00000005
  SmoInit=FORMATSmoInit[FORMATSmoInit$pval.exposure<Pvallimit ,]
SmoClumped<-clump_data(
  SmoInit,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)
SmoClumped$Dir<-NULL
SmoClumped$Dir[ SmoClumped$pval.exposure > 0.00000005]<-0
SmoClumped$Dir[ SmoClumped$pval.exposure < 0.00000005]<-1
SmoClumped$Total<-1
(sum(SmoClumped$Dir)/sum(SmoClumped$Total))*100
table(SmoClumped$Dir) 

#cig per week 

CigPDay<-read.table(gzfile("/mnt/storage/scratch/bu19525/BristolPhD/MVMR_follow_up/Clumping/CigarettesPerDay.WithoutUKB.txt.gz"),sep="\t",header=TRUE)
head(CigPDay)
summary(CigPDay$EFFECTIVE_N)
summary(CigPDay$EFFECTIVE_N)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
FORMATCigPDay<-format_data(CigPDay, type = "exposure", snp_col = "RSID",
                             beta_col = "BETA",
                             se_col = "SE",
                             eaf_col = "AF",
                             effect_allele_col = "ALT",
                             other_allele_col = "REF",
                             pval_col = "PVALUE",
                             samplesize_col = "EFFECTIVE_N",
                             chr_col = "CHROM",
                             pos_col = "POS"
)
  Pvallimit<-0.00000005

  CigPDay=FORMATCigPDay[FORMATCigPDay$pval.exposure<Pvallimit ,]
CigClumped<-clump_data(
  CigPDay,
  clump_kb = 10000,
  clump_r2 = 0.01,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)
CigClumped$Dir<-NULL
CigClumped$Dir[ CigClumped$pval.exposure > 0.00000005]<-0
CigClumped$Dir[ CigClumped$pval.exposure < 0.00000005]<-1
CigClumped$Total<-1
(sum(CigClumped$Dir)/sum(CigClumped$Total))*100
table(CigClumped$Dir) 
# Alcohol data read and clump (1smr)

AUDITCSNPS<-fread(
  "/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/AUDIT2SMRdata.txt"
)
AUDSNPS<-fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/AUDTwoSampleMR.txt")
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
AUDITClumped<-clump_data(
  AUDITCSNPS,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)
AUDClumped<-clump_data(
  AUDSNPS,
  clump_kb = 10000,
  clump_r2 = 0.001,
  clump_p1 = 1,
  clump_p2 = 1,
  pop = "EUR"
)
# 1smr analyses alcohol
AUDITClumped<-data.frame(AUDITClumped)
AUDClumped<-data.frame(AUDClumped)
FORMATAUDIT<-format_data(AUDITClumped, type = "exposure", snp_col = "SNP",
                             beta_col = "beta",
                             se_col = "se",
                             eaf_col = "eaf",
                             effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele",
                             pval_col = "pval"
)
outcome_data <- extract_outcome_data(
  snps = FORMATAUDIT$SNP, outcomes = "ebi-a-GCST007432")
H_data <- harmonise_data(
  exposure_dat = FORMATAUDIT, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
mr_results$id.exposure<-"AUDIT-C"
mr_results
mr_AUDIT<-mr_results
FullResults<-mr_results
Presso<-run_mr_presso(H_data)
Presso
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))
#mr_res
# AUD
FORMATAUD<-format_data(AUDClumped, type = "exposure", snp_col = "SNP",
                             beta_col = "beta",
                             se_col = "se",
                             eaf_col = "eaf",
                             effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele",
                             pval_col = "pval"
)
Outcome<-fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Data/FEV1_SNPs_GWAS_smokers.txt")
Outcome<-subset(Outcome,TEST=="ADD")
Outcome$phen<-"FEV1"
Outcome<-data.frame(Outcome)
Outcome$samplesize<-180000
Outcome_maf<-fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Data/FEV1_MAF_smokers.txt")
Outcome<-merge(Outcome,Outcome_maf,by="SNP")
outcome_data<-format_data(Outcome,type="outcome",
 phenotype_col = "phen",
 snp_col = "SNP",
  beta_col = "BETA",
  effect_allele_col = "A1.x",
  other_allele_col = "A2",
  eaf_col = "MAF",
  pval_col = "P",
  chr_col = "CHR.x",
  pos_col = "BP",
  samplesize_col = "samplesize",
  se_col = "SE",
  log_pval = FALSE
)

H_data <- harmonise_data(
  exposure_dat = FORMATAUD, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
mr_results$id.exposure<-"AUD"
mr_results
mr_AUD<-mr_results
FullResults<-merge(FullResults,mr_results)
Presso<-run_mr_presso(H_data)
Presso
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

#2SMR analyses smoking

#smoking -ukb 
Outcome<-fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Data/FEV1_SNPs_GWAS_smokers.txt")
Outcome<-subset(Outcome,TEST=="ADD")
Outcome$phen<-"FEV1"
Outcome<-data.frame(Outcome)
Outcome$samplesize<-180000
Outcome_maf<-fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Data/FEV1_MAF_smokers.txt")
Outcome<-merge(Outcome,Outcome_maf,by="SNP")
outcome_data<-format_data(Outcome,type="outcome",
 phenotype_col = "phen",
 snp_col = "SNP",
  beta_col = "BETA",
  effect_allele_col = "A1.x",
  other_allele_col = "A2",
  eaf_col = "MAF",
  pval_col = "P",
  chr_col = "CHR.x",
  pos_col = "BP",
  samplesize_col = "samplesize",
  se_col = "SE",
  log_pval = FALSE
)
outcome_data<-data.frame(outcome_data)
H_data <- harmonise_data(
  exposure_dat = SmoClumped, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
mr_results$id.exposure<-"Smoking_init"
mr_results
mr_smo<-mr_results
FullResults<-merge(FullResults,mr_results)
Presso<-run_mr_presso(H_data)
Presso
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

#Cigs per week - subset

#outcome_data <- extract_outcome_data(
 # snps = CigClumped$SNP, outcomes = "ebi-a-GCST007432")
H_data <- harmonise_data(
  exposure_dat = CigClumped, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_weighted_median","mr_weighted_mode"))
mr_results$id.exposure<-"Cigs"
mr_results
mr_cig<-mr_results
FullResults<-merge(FullResults,mr_results)
Presso<-run_mr_presso(H_data)
Presso
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))
 #concacenate
MR_merge<-NULL
MR_merge<-rbind(mr_smo,mr_AUD)
MR_merge<-rbind(MR_merge,mr_AUDIT)
MR_merge<-rbind(MR_merge,mr_cig)

setwd("/mnt/storage/scratch/bu19525/BristolPhD/MVMR_follow_up/Output/Smokers_corrections/")
write.csv(MR_merge,"Full2SMRrescorrect_smokers.csv")
write.table(MR_merge,"Full2SMRrescorrect_smokers.txt")

