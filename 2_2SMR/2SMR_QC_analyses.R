#libraries and local clumpting paths

library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(archive)
library(readr)
library(data.table)

#Smoking data read and clump
SmoInit<-read.table(gzfile("SmokingInitiation.WithoutUKB.txt.gz"),sep="\t",header=TRUE)
head(SmoInit)
summary(SmoInit$EFFECTIVE_N)
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

CigPDay<-read.table(gzfile("CigarettesPerDay.WithoutUKB.txt.gz"),sep="\t",header=TRUE)
head(CigPDay)
summary(CigPDay$EFFECTIVE_N)
summary(CigPDay$EFFECTIVE_N)
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
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))
mr_results
Presso<-run_mr_presso(H_data)
Presso
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

# AUD
FORMATAUD<-format_data(AUDClumped, type = "exposure", snp_col = "SNP",
                             beta_col = "beta",
                             se_col = "se",
                             eaf_col = "eaf",
                             effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele",
                             pval_col = "pval"
)
outcome_data <- extract_outcome_data(
  snps = FORMATAUD$SNP, outcomes = "ebi-a-GCST007432")
H_data <- harmonise_data(
  exposure_dat = FORMATAUD, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))
mr_results
Presso<-run_mr_presso(H_data)
Presso
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

#2SMR analyses smoking

#smoking -ukb 
outcome_data <- extract_outcome_data(
  snps = SmoClumped$SNP, outcomes = "ebi-a-GCST007432")
H_data <- harmonise_data(
  exposure_dat = SmoClumped, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))
mr_results
Presso<-run_mr_presso(H_data)
Presso
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))

#Cigs per week - subset

outcome_data <- extract_outcome_data(
  snps = CigClumped$SNP, outcomes = "ebi-a-GCST007432")
H_data <- harmonise_data(
  exposure_dat = CigClumped, 
  outcome_dat = outcome_data
)
mr_results<-mr(H_data, method_list=c("mr_wald_ratio", "mr_egger_regression", "mr_ivw", "mr_two_sample_ml","mr_simple_median","mr_weighted_median","mr_weighted_mode"))
mr_results
Presso<-run_mr_presso(H_data)
Presso
mr_pleiotropy_test(H_data)
mr_heterogeneity(H_data, method_list=c("mr_egger_regression", "mr_ivw"))



