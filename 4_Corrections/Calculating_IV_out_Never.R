library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(readr)
library(data.table)
library(remotes)
library(MVMR)
library(plyr)

#Read in SNPs
setwd("/mnt/storage/scratch/bu19525/BristolPhD/MVMR_follow_up/Code/MVMRalcsmo/4_Corrections/")
AUD_Cig<-read.csv("/mnt/storage/scratch/bu19525/BristolPhD/MVMR_follow_up/Code/MVMRalcsmo/3_MVMR/AUD_CIgs_per_SNPS.csv")
AUD_INIT<-read.csv("/mnt/storage/scratch/bu19525/BristolPhD/MVMR_follow_up/Code/MVMRalcsmo/3_MVMR/AUD_INIT_SNPS.csv")
AUDIT_CIG<-read.csv("/mnt/storage/scratch/bu19525/BristolPhD/MVMR_follow_up/Code/MVMRalcsmo/3_MVMR/AUDIT_Cigs_per_SNPS.csv")
AUDIT_INIT<-read.csv("/mnt/storage/scratch/bu19525/BristolPhD/MVMR_follow_up/Code/MVMRalcsmo/3_MVMR/AUDIT_INIT_SNPS.csv")

#Combine into SNP list

SNP_list<-as.vector(c(AUD_Cig$SNP,AUD_INIT$SNP,AUDIT_CIG$SNP,AUDIT_INIT$SNP))
write.table(SNP_list,"SNP_list.txt")

Pheno<-read.csv("/user/work/bu19525/BristolPhD/MVMR_follow_up/Data/yes_smokers_data.csv")
Pheno_only<-Pheno[,c(2,5,6)]
#Read Genetics

MAPfile<-fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Data/SNPs_for_MR.ped",header=F)
colnames(MAPfile)[1] <- "IID"

MAPfile = merge(MAPfile, Pheno_only, by='IID', all.x=TRUE, all.y=FALSE, sort=FALSE)
MAPfile$V6<-MAPfile$FEV1
#Remove non smokers
MAPfile<-subset(MAPfile,Smo_status==1)
MAPfile<-MAPfile[,c(1:86)]
setwd("/user/work/bu19525/BristolPhD/MVMR_follow_up/Data/")
write.table(MAPfile,'yes_smo_file.ped',col.names = FALSE, row.names = FALSE,quote = FALSE)

#Covariates!!!


Cov<-Pheno[,c(2,2,3,4)]
colnames(Cov)[1] <- "FID"
colnames(Cov)[2] <- "IID"

Cov<-match_df(Cov,MAPfile,on="IID")
write.table(Cov,'Cov_smokers.txt',row.names=FALSE,quote = FALSE)