#libraries and local clumpting paths

library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(readr)
library(data.table)
library(remotes)
library(MVMR)
options(ieugwasr_api = 'gwas-api.mrcieu.ac.uk/')
#Read_data_in
setwd("/user/work/bu19525/BristolPhD/MVMR_follow_up/Output/Corrections/")
AUDITC<-data.frame(fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/AUDIT_hg19_final.txt"))
AUD<-data.frame(fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/AUD_MVP1_MVP2_PGC_Dec11.txt"))
FORMATAUD<-format_data(AUD, type = "exposure", snp_col = "SNP_ID",
                             beta_col = "Effect",
                             se_col = "SE",
                             eaf_col = "EAF",
                             effect_allele_col = "Allele1",
                             other_allele_col = "Allele2",
                             pval_col = "PValue",
                             chr_col = "Chrosome",
                             pos_col = "Position"
)
FORMATAUDIT<-format_data(AUDITC, type = "exposure", snp_col = "RSID",
                             beta_col = "Beta",
                             se_col = "SE",
                             eaf_col = "EAF",
                             effect_allele_col = "A1",
                             other_allele_col = "A2",
                             pval_col = "PValue",
                             chr_col = "Chr",
                             pos_col = "Pos"
)





#Smoking data read and clump
SmoInit<-data.frame(read.table("/user/work/bu19525/BristolPhD/MVMR_follow_up/Data/Lifetime_smoking.txt",sep=" ",header=TRUE,row.names = NULL))
head(SmoInit)
summary(SmoInit$EFFECTIVE_N)
FORMATSmoInit<-format_data(SmoInit, type = "exposure", snp_col = "SNP",
                             beta_col = "BETA",
                             se_col = "SE",
                             eaf_col = "EAF",
                             effect_allele_col = "EFFECT_ALLELE",
                             other_allele_col = "OTHER_ALLELE",
                             pval_col = "P",
                             chr_col = "CHR",
                             pos_col = "BP"
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

CigPDay<-data.frame(read.table(gzfile("/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/CigarettesPerDay.WithoutUKB.txt.gz"),sep="\t",header=TRUE))
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
#SNP list
AUDITCSNPS<-fread(
  "/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/AUDIT2SMRdata.txt"
)
AUDSNPS<-fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/AUDTwoSampleMR.txt")
AuditCSNPs<-AUDITCSNPS$SNP
AudSNPs<-AUDSNPS$SNP
CigSNPs<-CigClumped$SNP
InitSNPs<-SmoClumped$SNP

AUD_Cig<-c(AudSNPs,CigSNPs)
AUD_Init<-c(AudSNPs,InitSNPs)
AUDIT_Cig<-c(AuditCSNPs,CigSNPs)
AUDIT_Init<-c(AuditCSNPs,InitSNPs)



#MVMR 
#AUD Cig analysis 1
FORMATAUDcigs <- FORMATAUD[FORMATAUD$SNP %in% AUD_Cig,]
FORMATAUDcigs$id.exposure<-"AUD"
FORMATCigsAUD <- FORMATCigPDay[FORMATCigPDay$SNP %in% AUD_Cig, ]
FORMATCigsAUD$id.exposure<-"CigsPerDay"
AUD_CIG_MVMR<-rbind(FORMATAUDcigs,FORMATCigsAUD)
write.table(AUD_CIG_MVMR,"AUD_CIG_MVMR.txt")
outcome_data <- extract_outcome_data(
  snps = FORMATAUDcigs$SNP, outcomes = "ebi-a-GCST007432")
  
H_data<-mv_harmonise_data(
  exposure_dat = AUD_CIG_MVMR, 
  outcome_dat = outcome_data
)
H_data_snps<-data.frame(H_data[1]$exposure_beta[,c(1,2)])
H_data_snps$SNP <- row.names(H_data_snps)  
#hdata<-data.frame(H_data)
#colnames(hdata)[1] <- "SNP"
MVMRdata<-format_mvmr(BXGs = H_data[1]$exposure_beta[,c(1,2)],BYG=H_data[4]$outcome_beta,seBXGs = H_data[3]$exposure_se[,c(1,2)],H_data[6]$outcome_se,RSID=H_data_snps$SNP)
sres<-strength_mvmr(r_input=MVMRdata,gencov=0)
sres<-strength_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov=0)
res<-ivw_mvmr(r_input=MVMRdata)
res<-data.frame(res)
row.names(res) <- c("AUD","Cigs P Day")
res
for(i in 1:2){
res$upper[i]<-(res$Estimate[i]+1.96*res$Std..Error[i])
res$upper[i]<-(res$Estimate[i]+1.96*res$Std..Error[i])
}
for(i in 1:2){
res$lower[i]<-(res$Estimate[i]-1.96*res$Std..Error[i])
res$lower[i]<-(res$Estimate[i]-1.96*res$Std..Error[i])
}
write.csv(res,"AUD_CIG_results_MVMR.csv")
Merge_MVMR<-res
#MRegger

#bx<-H_data[1]$exposure_beta[,c(1,2)]
#bxse<-H_data[3]$exposure_se[,c(1,2)]
#by<-H_data[4]$outcome_beta
#byse<-H_data[6]$outcome_se
#MVMR_egger_input<-mr_mvinput(
#  bx = bx,
#  bxse = bxse,
#  by = by,
#  byse = byse,
#  #correlation = 0,
#  exposure = c("AUD","Cig per day"),
#  outcome = "FEV1",
#  snps = H_data_snps$SNP,
#  effect_allele = NA,
#  other_allele = NA,
#  eaf = NA
#)
#AUD_CIG_IVW<-mr_mvivw(
#  MVMR_egger_input,
#  model = "default",
#  robust = FALSE,
#  correl = FALSE,
#  distribution = "normal",
#  alpha = 0.05
#)
#AUD_CIG_IVW
#AUD_CIG_EGGER<-mr_mvegger(
#  MVMR_egger_input,
#  orientate = 1,
#  correl = FALSE,
#  distribution = "normal",
#  alpha = 0.05
#)
#AUD_CIG_EGGER
#
#MVMR_results_IVW<-NULL
#MVMR_results_IVW$Exposure<-AUD_CIG_IVW$Exposure 
#MVMR_results_IVW$Beta<-AUD_CIG_IVW$Estimate 
#MVMR_results_IVW$SE<-AUD_CIG_IVW$StdError 
#MVMR_results_IVW$Lower<-AUD_CIG_IVW$CILower
#MVMR_results_IVW$Upper<-AUD_CIG_IVW$CIUpper
#MVMR_results_IVW$Pvalue<-AUD_CIG_IVW$Pvalue
#MVMR_results_IVW$Analysis<-"IVW"
#MVMR_results_IVW<-as.data.frame(MVMR_results_IVW)
#
#MVMR_results_Egger<-NULL
#MVMR_results_Egger$Exposure<-AUD_CIG_EGGER$Exposure 
#MVMR_results_Egger$Beta<-AUD_CIG_EGGER$Estimate 
#MVMR_results_Egger$SE<-AUD_CIG_EGGER$StdError.Est 
#MVMR_results_Egger$Lower<-AUD_CIG_EGGER$CILower.Est
#MVMR_results_Egger$Upper<-AUD_CIG_EGGER$CIUpper.Est
#MVMR_results_Egger$Pvalue<-AUD_CIG_EGGER$Pvalue.Est
#MVMR_results_Egger$Analysis<-"MR-Egger"
#MVMR_results_Egger<-as.data.frame(MVMR_results_Egger)
# Merge
#Merge_MVMR_IVW<-MVMR_results_IVW
#Merge_MVMR_egger<-MVMR_results_Egger

#AUDIT Cigs analysis 2
FORMATAUDITcigs <- FORMATAUDIT[FORMATAUDIT$SNP %in% AUDIT_Cig, ]
FORMATAUDITcigs$id.exposure<-"AUDIT"
FORMATCigsAUDIT <- FORMATCigPDay[FORMATCigPDay$SNP %in% AUDIT_Cig, ]
FORMATCigsAUDIT$id.exposure<-"CigsPerDay"
AUDIT_CIG_MVMR<-rbind(FORMATAUDITcigs,FORMATCigsAUDIT)
write.table(AUDIT_CIG_MVMR,"AUDIT_CIG_MVMR.txt")
write.csv(AUDIT_CIG_MVMR,"AUDIT_CIG_MVMR.csv")
outcome_data <- extract_outcome_data(
  snps = FORMATAUDITcigs$SNP, outcomes = "ebi-a-GCST007432")

H_data<-mv_harmonise_data(
  exposure_dat = AUDIT_CIG_MVMR, 
  outcome_dat = outcome_data
)
H_data_snps<-data.frame(H_data[1]$exposure_beta[,c(1,2)])
H_data_snps$SNP <- row.names(H_data_snps)  
#hdata<-data.frame(H_data)
#colnames(hdata)[1] <- "SNP"
MVMRdata<-format_mvmr(BXGs = H_data[1]$exposure_beta[,c(1,2)],BYG=H_data[4]$outcome_beta,seBXGs = H_data[3]$exposure_se[,c(1,2)],H_data[6]$outcome_se,RSID=H_data_snps$SNP)
sres<-strength_mvmr(r_input=MVMRdata,gencov=0)
sres<-strength_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov=0)
res<-ivw_mvmr(r_input=MVMRdata)
res<-data.frame(res)
row.names(res) <- c("AUDIT","Cigs P Day")
res
for(i in 1:2){
res$upper[i]<-(res$Estimate[i]+1.96*res$Std..Error[i])
res$upper[i]<-(res$Estimate[i]+1.96*res$Std..Error[i])
}
for(i in 1:2){
res$lower[i]<-(res$Estimate[i]-1.96*res$Std..Error[i])
res$lower[i]<-(res$Estimate[i]-1.96*res$Std..Error[i])
}
write.csv(res,"AUDIT_CIG_results_MVMR.csv")
Merge_MVMR<-rbind(Merge_MVMR,res)
#MVMRres<-mr_mvegger(MVMRdata)
#MRegger
#bx<-H_data[1]$exposure_beta[,c(1,2)]
#bxse<-H_data[3]$exposure_se[,c(1,2)]
#by<-H_data[4]$outcome_beta
#byse<-H_data[6]$outcome_se
#MVMR_egger_input<-mr_mvinput(
#  bx = bx,
#  bxse = bxse,
#  by = by,
#  byse = byse,
#  #correlation = 0,
#  exposure = c("AUDIT-C","Cig per day"),
#  outcome = "FEV1",
#  snps = H_data_snps$SNP,
#  effect_allele = NA,
#  other_allele = NA,
#  eaf = NA
#)
#AUDIT_CIG_IVW<-mr_mvivw(
#  MVMR_egger_input,
#  model = "default",
#  robust = FALSE,
#  correl = FALSE,
#  distribution = "normal",
#  alpha = 0.05
#)
#AUDIT_CIG_IVW
#AUDIT_CIG_EGGER<-mr_mvegger(
#  MVMR_egger_input,
#  orientate = 1,
#  correl = FALSE,
#  distribution = "normal",
#  alpha = 0.05
#)
#AUDIT_CIG_EGGER#
#
#
#MVMR_results_IVW<-NULL
#MVMR_results_IVW$Exposure<-AUDIT_CIG_IVW$Exposure 
#MVMR_results_IVW$Beta<-AUDIT_CIG_IVW$Estimate 
#MVMR_results_IVW$SE<-AUDIT_CIG_IVW$StdError 
#MVMR_results_IVW$Lower<-AUDIT_CIG_IVW$CILower
#MVMR_results_IVW$Upper<-AUDIT_CIG_IVW$CIUpper
#MVMR_results_IVW$Pvalue<-AUDIT_CIG_IVW$Pvalue
#MVMR_results_IVW$Analysis<-"IVW"
#MVMR_results_IVW<-as.data.frame(MVMR_results_IVW)
##
#MVMR_results_Egger<-NULL
#MVMR_results_Egger$Exposure<-AUDIT_CIG_EGGER$Exposure 
#MVMR_results_Egger$Beta<-AUDIT_CIG_EGGER$Estimate 
#MVMR_results_Egger$SE<-AUDIT_CIG_EGGER$StdError.Est 
#MVMR_results_Egger$Lower<-AUDIT_CIG_EGGER$CILower.Est
#MVMR_results_Egger$Upper<-AUDIT_CIG_EGGER$CIUpper.Est
#MVMR_results_Egger$Pvalue<-AUDIT_CIG_EGGER$Pvalue.Est
#MVMR_results_Egger$Analysis<-"MR-Egger"
#MVMR_results_Egger<-as.data.frame(MVMR_results_Egger)
## Merge
#Merge_MVMR_IVW<-rbind(Merge_MVMR_IVW,MVMR_results_IVW)
#Merge_MVMR_egger<-rbind(Merge_MVMR_egger,MVMR_results_Egger)#


# AUD init analysis 3

FORMATAUDinit <- FORMATAUD[FORMATAUD$SNP %in% AUD_Init, ]
FORMATAUDinit$id.exposure<-"AUD"
FORMATInitAUD <- FORMATSmoInit[FORMATSmoInit$SNP %in% AUD_Init, ]
FORMATInitAUD$id.exposure<-"Smo_score"
AUD_INIT_MVMR<-rbind(FORMATAUDinit,FORMATInitAUD)
write.table(AUD_INIT_MVMR,"AUD_INIT_MVMR.txt")
outcome_data <- extract_outcome_data(
  snps = FORMATAUDinit$SNP, outcomes = "ebi-a-GCST007432")
  
H_data<-mv_harmonise_data(
  exposure_dat = AUD_INIT_MVMR, 
  outcome_dat = outcome_data
)
H_data_snps<-data.frame(H_data[1]$exposure_beta[,c(1,2)])
H_data_snps$SNP <- row.names(H_data_snps)  
#hdata<-data.frame(H_data)
#colnames(hdata)[1] <- "SNP"
MVMRdata<-format_mvmr(BXGs = H_data[1]$exposure_beta[,c(1,2)],BYG=H_data[4]$outcome_beta,seBXGs = H_data[3]$exposure_se[,c(1,2)],H_data[6]$outcome_se,RSID=H_data_snps$SNP)
sres<-strength_mvmr(r_input=MVMRdata,gencov=0)
sres<-strength_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov=0)
res<-ivw_mvmr(r_input=MVMRdata)
res<-data.frame(res)
row.names(res) <- c("AUD","Initiation")
res
for(i in 1:2){
res$upper[i]<-(res$Estimate[i]+1.96*res$Std..Error[i])
res$upper[i]<-(res$Estimate[i]+1.96*res$Std..Error[i])
}
for(i in 1:2){
res$lower[i]<-(res$Estimate[i]-1.96*res$Std..Error[i])
res$lower[i]<-(res$Estimate[i]-1.96*res$Std..Error[i])
}
write.csv(res,"AUD_init_results_MVMR.csv")
Merge_MVMR<-rbind(Merge_MVMR,res)
#MRegger


#bx<-H_data[1]$exposure_beta[,c(1,2)]
#bxse<-H_data[3]$exposure_se[,c(1,2)]
#by<-H_data[4]$outcome_beta
#byse<-H_data[6]$outcome_se
#MVMR_egger_input<-mr_mvinput(
#  bx = bx,
#  bxse = bxse,
#  by = by,
#  byse = byse,
#  #correlation = 0,
#  exposure = c("AUD","Init"),
#  outcome = "FEV1",
#  snps = H_data_snps$SNP,
#  effect_allele = NA,
#  other_allele = NA,
#  eaf = NA
#)
#AUD_init_IVW<-mr_mvivw(
#  MVMR_egger_input,
#  model = "default",
#  robust = FALSE,
#  correl = FALSE,
#  distribution = "normal",
#  alpha = 0.05
#)
#AUD_init_IVW
#AUD_init_EGGER<-mr_mvegger(
#  MVMR_egger_input,
#  orientate = 1,
#  correl = FALSE,
#  distribution = "normal",
#  alpha = 0.05
#)
#AUD_init_EGGER
##MVMR combine
#MVMR_results_IVW<-NULL
#MVMR_results_IVW$Exposure<-AUD_init_IVW$Exposure 
#MVMR_results_IVW$Beta<-AUD_init_IVW$Estimate 
#MVMR_results_IVW$SE<-AUD_init_IVW$StdError 
#MVMR_results_IVW$Lower<-AUD_init_IVW$CILower
#MVMR_results_IVW$Upper<-AUD_init_IVW$CIUpper
#MVMR_results_IVW$Pvalue<-AUD_init_IVW$Pvalue
#MVMR_results_IVW$Analysis<-"IVW"
#MVMR_results_IVW<-as.data.frame(MVMR_results_IVW)
##
#MVMR_results_Egger<-NULL
#MVMR_results_Egger$Exposure<-AUD_init_EGGER$Exposure 
#MVMR_results_Egger$Beta<-AUD_init_EGGER$Estimate 
#MVMR_results_Egger$SE<-AUD_init_EGGER$StdError.Est 
#MVMR_results_Egger$Lower<-AUD_init_EGGER$CILower.Est
#MVMR_results_Egger$Upper<-AUD_init_EGGER$CIUpper.Est
#MVMR_results_Egger$Pvalue<-AUD_init_EGGER$Pvalue.Est
#MVMR_results_Egger$Analysis<-"MR-Egger"
#MVMR_results_Egger<-as.data.frame(MVMR_results_Egger)
# Merge
#Merge_MVMR_IVW<-rbind(Merge_MVMR_IVW,MVMR_results_IVW)
#Merge_MVMR_egger<-rbind(Merge_MVMR_egger,MVMR_results_Egger)

#
#AUDIT init analysis 4
FORMATAUDITinit <- FORMATAUDIT[FORMATAUDIT$SNP %in% AUDIT_Init, ]
FORMATAUDITinit$id.exposure<-"AUDIT-C"
FORMATInitAUDIT <- FORMATSmoInit[FORMATSmoInit$SNP %in% AUDIT_Init, ]
FORMATInitAUDIT$id.exposure<-"Smo Init"
AUDIT_INIT_MVMR<-rbind(FORMATAUDITinit,FORMATInitAUDIT)
write.table(AUDIT_INIT_MVMR,"AUDIT_INIT_MVMR.txt")

outcome_data <- extract_outcome_data(
  snps = FORMATAUDITinit$SNP, outcomes = "ebi-a-GCST007432")
  
H_data<-mv_harmonise_data(
  exposure_dat = AUDIT_INIT_MVMR, 
  outcome_dat = outcome_data
)
H_data_snps<-data.frame(H_data[1]$exposure_beta[,c(1,2)])
H_data_snps$SNP <- row.names(H_data_snps)  
#hdata<-data.frame(H_data)
#colnames(hdata)[1] <- "SNP"
MVMRdata<-format_mvmr(BXGs = H_data[1]$exposure_beta[,c(1,2)],BYG=H_data[4]$outcome_beta,seBXGs = H_data[3]$exposure_se[,c(1,2)],H_data[6]$outcome_se,RSID=H_data_snps$SNP)
sres<-strength_mvmr(r_input=MVMRdata,gencov=0)
sres<-strength_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov="cor01")
pres<-pleiotropy_mvmr(r_input=MVMRdata,gencov=0)
res<-ivw_mvmr(r_input=MVMRdata)
res<-data.frame(res)
row.names(res) <- c("AUDIT","Smo Init")
res
for(i in 1:2){
res$upper[i]<-(res$Estimate[i]+1.96*res$Std..Error[i])
res$upper[i]<-(res$Estimate[i]+1.96*res$Std..Error[i])
}
for(i in 1:2){
res$lower[i]<-(res$Estimate[i]-1.96*res$Std..Error[i])
res$lower[i]<-(res$Estimate[i]-1.96*res$Std..Error[i])
}
write.csv(res,"AUDIT_init_results_MVMR_LSI.csv")
Merge_MVMR<-rbind(Merge_MVMR,res)
write.table(Merge_MVMR,"MVMR_results_LSI.txt")
#MR CLust
#MRegger
#bx<-H_data[1]$exposure_beta[,c(1,2)]
#bxse<-H_data[3]$exposure_se[,c(1,2)]
#by<-H_data[4]$outcome_beta
#byse<-H_data[6]$outcome_se
#MVMR_egger_input<-mr_mvinput(
#  bx = bx,
#  bxse = bxse,
#  by = by,
#  byse = byse,
#  #correlation = 0,
#  exposure = c("AUDIT-C","Smoking Init"),
#  outcome = "FEV1",
#  snps = H_data_snps$SNP,
#  effect_allele = NA,
#  other_allele = NA,
#  eaf = NA
#)
#AUDIT_init_IVW<-mr_mvivw(
#  MVMR_egger_input,
#  model = "default",
#  robust = FALSE,
#  correl = FALSE,
#  distribution = "normal",
#  alpha = 0.05
#)
#AUDIT_init_IVW
#AUDIT_init_EGGER<-mr_mvegger(
#  MVMR_egger_input,
#  orientate = 1,
#  correl = FALSE,
#  distribution = "normal",
#  alpha = 0.05
#)
#AUDIT_init_EGGER
##MVMR
#MVMR_results_IVW<-NULL
#MVMR_results_IVW$Exposure<-AUDIT_init_IVW$Exposure 
#MVMR_results_IVW$Beta<-AUDIT_init_IVW$Estimate 
#MVMR_results_IVW$SE<-AUDIT_init_IVW$StdError 
#MVMR_results_IVW$Lower<-AUDIT_init_IVW$CILower
#MVMR_results_IVW$Upper<-AUDIT_init_IVW$CIUpper
#MVMR_results_IVW$Pvalue<-AUDIT_init_IVW$Pvalue
#MVMR_results_IVW$Analysis<-"IVW"
#MVMR_results_IVW<-as.data.frame(MVMR_results_IVW)
##
#MVMR_results_Egger<-NULL
#MVMR_results_Egger$Exposure<-AUDIT_init_EGGER$Exposure 
#MVMR_results_Egger$Beta<-AUDIT_init_EGGER$Estimate 
#MVMR_results_Egger$SE<-AUDIT_init_EGGER$StdError.Est 
#MVMR_results_Egger$Lower<-AUDIT_init_EGGER$CILower.Est
#MVMR_results_Egger$Upper<-AUDIT_init_EGGER$CIUpper.Est
#MVMR_results_Egger$Pvalue<-AUDIT_init_EGGER$Pvalue.Est
#MVMR_results_Egger$Analysis<-"MR-Egger"
#MVMR_results_Egger<-as.data.frame(MVMR_results_Egger)
# Merge
#Merge_MVMR_IVW<-rbind(Merge_MVMR_IVW,MVMR_results_IVW)
#Merge_MVMR_egger<-rbind(Merge_MVMR_egger,MVMR_results_Egger)
#Final_combine<-rbind(Merge_MVMR_IVW,Merge_MVMR_egger)
#write.table(Final_combine,"MVMR_Egger_results.txt")
