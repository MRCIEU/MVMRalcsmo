#libraries and local clumpting paths
library(MendelianRandomization)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(archive)
library(readr)
library(data.table)
library(remotes)
library(MVMR)
#Read_data_in

AUDITC<-fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/AUDIT_hg19_final.txt")
AUD<-fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/AUD_MVP1_MVP2_PGC_Dec11.txt")
FORMATAUD<-format_data(AUD, type = "exposure", snp_col = "SNP_ID",
                             beta_col = "Effect",
                             se_col = "SE",
                             eaf_col = "EAF",
                             effect_allele_col = "Allele1",
                             other_allele_col = "Allele2",
                             pval_col = "PValue",
                             samplesize_col = "SampleSize",
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
                             samplesize_col = "Sample_N",
                             chr_col = "Chr",
                             pos_col = "Pos"
)





SmoInit<-read.table(gzfile("/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/SmokingInitiation.WithoutUKB.txt.gz"),sep="\t",header=TRUE)
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
CigPDay<-read.table(gzfile("/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/CigarettesPerDay.WithoutUKB.txt.gz"),sep="\t",header=TRUE)
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
#AUD Cig
FORMATAUDcigs <- FORMATAUD[FORMATAUD$SNP %in% AUD_Cig, ]
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


#AUDIT Cigs
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
#MVMRres<-mr_mvegger(MVMRdata)

# AUD init

FORMATAUDinit <- FORMATAUD[FORMATAUD$SNP %in% AUD_Init, ]
FORMATAUDinit$id.exposure<-"AUD"
FORMATInitAUD <- FORMATSmoInit[FORMATSmoInit$SNP %in% AUD_Init, ]
FORMATInitAUD$id.exposure<-"Init"
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

#AUDIT init
FORMATAUDITinit <- FORMATAUDIT[FORMATAUDIT$SNP %in% AUDIT_Init, ]
FORMATAUDITinit$id.exposure<-"AUDIT"
FORMATInitAUDIT <- FORMATSmoInit[FORMATSmoInit$SNP %in% AUDIT_Init, ]
FORMATInitAUDIT$id.exposure<-"Init"
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
row.names(res) <- value(c("AUDIT","Initiation"))
res
for(i in 1:2){
res$upper[i]<-(res$Estimate[i]+1.96*res$Std..Error[i])
res$upper[i]<-(res$Estimate[i]+1.96*res$Std..Error[i])
}
for(i in 1:2){
res$lower[i]<-(res$Estimate[i]-1.96*res$Std..Error[i])
res$lower[i]<-(res$Estimate[i]-1.96*res$Std..Error[i])
}
write.csv(res,"AUDIT_init_results_MVMR.csv")

#MR CLust
