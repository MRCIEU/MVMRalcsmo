library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(readr)
library(data.table)
library(remotes)
library(MVMR)

setwd("/user/work/bu19525/BristolPhD/MVMR_follow_up/Data/")

Linker<-read.csv("/user/work/bu19525/BristolPhD/MVMR_follow_up/Data/linker.81499.csv")
Pheno<-fread("/user/work/bu19525/SRA_files/Cancer_epi_project/Data/UKB/data.51913.csv",check.names=TRUE)

#Extract only the ones i want!!
#ID
which(colnames(Pheno) == "eid")
# 1
#SEX
which(colnames(Pheno) == "X31.0.0")
#23
#AGE
which(colnames(Pheno) == "X21022.0.0")
# 9838
#Smo status
which(colnames(Pheno) == "X20160.0.0")
# 9034
#FEV1
which(colnames(Pheno) == "X20256.0.0")

#9034
Pheno<-Pheno[,c(1,23,9838,9034,9603)]
#Name Linker
colnames(Linker)[2] <- "eid"
colnames(Linker)[1] <- "IID"

#Name Pheno
colnames(Pheno)[1] <- "eid"
colnames(Pheno)[2] <- "Sex"
colnames(Pheno)[3] <- "Age"
colnames(Pheno)[4] <- "Smo_status"
colnames(Pheno)[5] <- "FEV1"

#Merge
phenos = merge(Linker, Pheno, by='eid', all.x=TRUE, all.y=FALSE, sort=FALSE)

Non_smokers<-subset(phenos,Smo_status==1)
#Non_smokers

write.table(Non_smokers,'yes_smokers_data.csv', sep=',', row.names=FALSE)