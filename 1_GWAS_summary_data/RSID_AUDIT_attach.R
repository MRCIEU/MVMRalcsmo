library(biomaRt)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(archive)
library(readr)
library(data.table)

AUDITCdata<-fread("/user/work/bu19525/BristolPhD/MVMR_follow_up/Clumping/EUR_rawadjavg_merged.txt")

RSID_audit<-AUDITCdata[,c(1,2)]
RSID_AUDIT_Final<-NULL
RSID_audit$Chr_ID<- apply(RSID_audit, 1, paste, collapse = ":")
AUDIT_ChrPos<-merge(AUDITCdata,RSID_audit,by="V2")

RSID_merge<-fread("RSID_merge.txt")
colnames(RSID_merge)[1] <- "RSID"
colnames(RSID_merge)[2] <- "Chr_ID"
colnames(AUDIT_ChrPos)[1] <- "Pos"
colnames(AUDIT_ChrPos)[2] <- "Chr"
colnames(AUDIT_ChrPos)[3] <- "A1"
colnames(AUDIT_ChrPos)[4] <- "A2"
colnames(AUDIT_ChrPos)[5] <- "EAF"
colnames(AUDIT_ChrPos)[6] <- "Sample_N"
colnames(AUDIT_ChrPos)[7] <- "Beta"
colnames(AUDIT_ChrPos)[8] <- "SE"
colnames(AUDIT_ChrPos)[9] <- "PValue"

AUDIT_with_RSID<-merge(RSID_merge,AUDIT_ChrPos,by="Chr_ID")
write.table(AUDIT_with_RSID,"AUDIT_with_RSID_hg19.txt", col.names = TRUE)