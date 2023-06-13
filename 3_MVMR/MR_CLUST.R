library(ggplot2)
library(TwoSampleMR)
library(MRInstruments)
library(ieugwasr)
library(archive)
library(readr)
library(data.table)
library(remotes)
library(MVMR)
library(mrclust)
#Read in data
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

FORMATAUDIT<-format_data(AUDITClumped, type = "exposure", snp_col = "SNP",
                             beta_col = "beta",
                             se_col = "se",
                             eaf_col = "eaf",
                             effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele",
                             pval_col = "pval"
)

FORMATAUD<-format_data(AUDClumped, type = "exposure", snp_col = "SNP",
                             beta_col = "beta",
                             se_col = "se",
                             eaf_col = "eaf",
                             effect_allele_col = "effect_allele",
                             other_allele_col = "other_allele",
                             pval_col = "pval"
)
#Smoking 

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

#cig per week 

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
CigClumped$Dir<-NULL
CigClumped$Dir[ CigClumped$pval.exposure > 0.00000005]<-0
CigClumped$Dir[ CigClumped$pval.exposure < 0.00000005]<-1
CigClumped$Total<-1
(sum(CigClumped$Dir)/sum(CigClumped$Total))*100
table(CigClumped$Dir) 

#AUD#################################
outcome_data <- extract_outcome_data(
  snps = FORMATAUD$SNP, outcomes = "ebi-a-GCST007432")
H_data <- harmonise_data(
  exposure_dat = FORMATAUD, 
  outcome_dat = outcome_data
)
AUD_Clust<-H_data
#MR CLUST
set.seed(2000030885)
bx = AUD_Clust$beta.exposure
bxse = AUD_Clust$se.exposure
by = AUD_Clust$beta.outcome
byse = AUD_Clust$se.outcome
ratio_est = by/bx
ratio_est_se = byse/abs(bx)
snp_names = AUD_Clust$SNP
res_em_AUD = mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, bx = bx,
                     by = by, bxse = bxse, byse = byse, obs_names = snp_names)
# Graph AUD
AUD_Cluster_1<-c(res_em_AUD$cluster_membership$cluster_1[1:3])
AUD_Clus_1 <- Reduce(c,AUD_Cluster_1)
AUD_Cluster_2<-c(res_em_AUD$cluster_membership$cluster_2[1:3])
AUD_Clus_2 <- Reduce(c,AUD_Cluster_2)
AUD_Cluster_Junk<-c(res_em_AUD$cluster_membership$cluster_Junk[1:3])
AUD_Clus_Junk <- Reduce(c,AUD_Cluster_Junk)
AUD_Cluster_Null<-c(res_em_AUD$cluster_membership$cluster_Null[1:3])
AUD_Clus_Null <- Reduce(c,AUD_Cluster_Null)

AUD_table1<-as.data.frame(AUD_Clus_1)
colnames(AUD_table1)[1] <- "SNP"
AUD_table1$Cluster<-"1"


AUD_table2<-as.data.frame(AUD_Clus_2)
colnames(AUD_table2)[1] <- "SNP"
AUD_table2$Cluster<-"2"
AUD_tableJunk<-as.data.frame(AUD_Clus_Junk)
colnames(AUD_tableJunk)[1] <- "SNP"
AUD_tableJunk$Cluster<-"Junk"
AUD_tableNull<-as.data.frame(AUD_Clus_Null)
colnames(AUD_tableNull)[1] <- "SNP"
AUD_tableNull$Cluster<-"Null"
AUD_excel<-rbind(AUD_table1,AUD_table2)
AUD_excel<-rbind(AUD_excel,AUD_tableNull)
AUD_excel<-rbind(AUD_excel,AUD_tableJunk)
write.csv(AUD_excel,"AUD_MR_Clust_excel.csv")
AUD_SNP2GENE<-merge(AUD_excel,FORMATAUD)
write.csv(AUD_SNP2GENE,"AUD_SNP2GENE.csv") 
colnames(AUD_SNP2GENE)[1] <- "SNP"
colnames(AUD_SNP2GENE)[4] <- "A2"
colnames(AUD_SNP2GENE)[3] <- "A1"
colnames(AUD_SNP2GENE)[8] <- "pvalue"
colnames(AUD_SNP2GENE)[5] <- "Beta"
colnames(AUD_SNP2GENE)[6] <- "SE"
write.table(AUD_SNP2GENE, "AUD_SNP2GENE.txt", append = FALSE, sep = "\t",quote=FALSE, dec = ".",
            row.names = FALSE, col.names = TRUE)
plot_AUD_best = res_em_AUD$plots$two_stage +
  ggplot2::xlim(0, max(abs(bx) + 2*bxse)) +
  ggplot2::xlab("Genetic association with AUD") +
  ggplot2::ylab("Genetic association with FEV1") +
  ggplot2::ggtitle("")
plot_AUD_best
ggsave("plot_AUD_best.tiff", width = 28, height = 35, dpi=50)
ggsave("plot_AUD_best_smaller.tiff", width = 8, height = 10, dpi=50)
#AUDIT
outcome_data <- extract_outcome_data(
  snps = FORMATAUDIT$SNP, outcomes = "ebi-a-GCST007432")
H_data <- harmonise_data(
  exposure_dat = FORMATAUDIT, 
  outcome_dat = outcome_data
)
AUDIT_Clust<-H_data
set.seed(2000030885)
bx = AUDIT_Clust$beta.exposure
bxse = AUDIT_Clust$se.exposure
by = AUDIT_Clust$beta.outcome
byse = AUDIT_Clust$se.outcome
ratio_est = by/bx
ratio_est_se = byse/abs(bx)
snp_names = AUDIT_Clust$SNP
res_em_AUDIT = mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, bx = bx,
                     by = by, bxse = bxse, byse = byse, obs_names = snp_names)
AUDIT_Cluster_1<-c(res_em_AUDIT$cluster_membership$cluster_1[1:3])
AUDIT_Clus_1 <- Reduce(c,AUDIT_Cluster_1)
AUDIT_Cluster_2<-c(res_em_AUDIT$cluster_membership$cluster_2[1:3])
AUDIT_Clus_2 <- Reduce(c,AUDIT_Cluster_2)
AUDIT_Cluster_Junk<-c(res_em_AUDIT$cluster_membership$cluster_Junk[1:3])
AUDIT_Clus_Junk <- Reduce(c,AUDIT_Cluster_Junk)
AUDIT_Cluster_Null<-c(res_em_AUDIT$cluster_membership$cluster_Null[1:3])
AUDIT_Clus_Null <- Reduce(c,AUDIT_Cluster_Null)

AUDIT_table1<-as.data.frame(AUDIT_Clus_1)
colnames(AUDIT_table1)[1] <- "SNP"
AUDIT_table1$Cluster<-"1"


AUDIT_table2<-as.data.frame(AUDIT_Clus_2)
colnames(AUDIT_table2)[1] <- "SNP"
AUDIT_table2$Cluster<-"2"
AUDIT_tableJunk<-as.data.frame(AUDIT_Clus_Junk)
colnames(AUDIT_tableJunk)[1] <- "SNP"
AUDIT_tableJunk$Cluster<-"Junk"
AUDIT_tableNull<-as.data.frame(AUDIT_Clus_Null)
colnames(AUDIT_tableNull)[1] <- "SNP"
AUDIT_tableNull$Cluster<-"Null"
AUDIT_excel<-rbind(AUDIT_table1,AUDIT_table2)
AUDIT_excel<-rbind(AUDIT_excel,AUDIT_tableNull)
AUDIT_excel<-rbind(AUDIT_excel,AUDIT_tableJunk)
write.csv(AUDIT_excel,"AUDIT_MR_Clust_excel.csv")
AUDIT_SNP2GENE<-merge(AUDIT_excel,FORMATAUDIT)   
write.csv(AUDIT_SNP2GENE,"AUDIT_SNP2GENE.csv")   
colnames(AUDIT_SNP2GENE)[1] <- "SNP"
colnames(AUDIT_SNP2GENE)[4] <- "A2"
colnames(AUDIT_SNP2GENE)[3] <- "A1"
colnames(AUDIT_SNP2GENE)[7] <- "pvalue"
colnames(AUDIT_SNP2GENE)[5] <- "Beta"
colnames(AUDIT_SNP2GENE)[6] <- "SE"         
write.table(AUD_SNP2GENE, "AUDIT_SNP2GENE.txt", append = FALSE, sep = "\t",quote=FALSE, dec = ".",
            row.names = FALSE, col.names = TRUE)
                     plot_AUDIT_best = res_em_AUDIT$plots$two_stage +
  ggplot2::xlim(0, max(abs(bx) + 2*bxse)) +
  ggplot2::xlab("Genetic association with AUDIT-C") +
  ggplot2::ylab("Genetic association with FEV1") +
  ggplot2::ggtitle("")
plot_AUDIT_best
ggsave("plot_AUDIT_best.tiff", width = 28, height = 35, dpi=50)
ggsave("plot_AUDIT_smaller.tiff", width = 8, height = 10, dpi=50)

#Cig



outcome_data <- extract_outcome_data(
  snps = CigClumped$SNP, outcomes = "ebi-a-GCST007432")
H_data <- harmonise_data(
  exposure_dat = CigClumped, 
  outcome_dat = outcome_data
)
CigClumped_Clust<-H_data
#MR CLUST
set.seed(2000030885)
bx = CigClumped_Clust$beta.exposure
bxse = CigClumped_Clust$se.exposure
by = CigClumped_Clust$beta.outcome
byse = CigClumped_Clust$se.outcome
ratio_est = by/bx
ratio_est_se = byse/abs(bx)
snp_names = CigClumped_Clust$SNP
res_em_CIG = mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, bx = bx,
                     by = by, bxse = bxse, byse = byse, obs_names = snp_names)
                     names(res_em_CIG$cluster_membership)

Cig_Cluster_1<-c(res_em_CIG$cluster_membership$cluster_1[1:3])
Cig_Clus_1 <- Reduce(c,Cig_Cluster_1)
Cig_Cluster_2<-c(res_em_CIG$cluster_membership$cluster_2[1:3])
Cig_Clus_2 <- Reduce(c,Cig_Cluster_2)
Cig_Cluster_Junk<-c(res_em_CIG$cluster_membership$cluster_Junk[1:3])
Cig_Clus_Junk <- Reduce(c,Cig_Cluster_Junk)
Cig_Cluster_Null<-c(res_em_CIG$cluster_membership$cluster_Null[1:3])
Cig_Clus_Null <- Reduce(c,Cig_Cluster_Null)

Cig_table1<-as.data.frame(Cig_Clus_1)
colnames(Cig_table1)[1] <- "SNP"
Cig_table1$Cluster<-"1"


Cig_table2<-as.data.frame(Cig_Clus_2)
colnames(Cig_table2)[1] <- "SNP"
Cig_table2$Cluster<-"2"
Cig_tableJunk<-as.data.frame(Cig_Clus_Junk)
colnames(Cig_tableJunk)[1] <- "SNP"
Cig_tableJunk$Cluster<-"Junk"
Cig_tableNull<-as.data.frame(Cig_Clus_Null)
colnames(Cig_tableNull)[1] <- "SNP"
Cig_tableNull$Cluster<-"Null"
Cig_excel<-rbind(Cig_table1,Cig_table2)
Cig_excel<-rbind(Cig_excel,Cig_tableNull)
Cig_excel<-rbind(Cig_excel,Cig_tableJunk)
write.csv(Cig_excel,"Cig_MR_Clust_excel.csv")
Cig_SNP2GENE<-merge(Cig_excel,CigClumped)       
write.csv(  Cig_SNP2GENE, "Cig_SNP2GENE.csv") 
colnames(Cig_SNP2GENE)[1] <- "SNP"
colnames(Cig_SNP2GENE)[3] <- "CHR"
colnames(Cig_SNP2GENE)[4] <- "pos"
colnames(Cig_SNP2GENE)[5] <- "A2"
colnames(Cig_SNP2GENE)[6] <- "A1"
colnames(Cig_SNP2GENE)[8] <- "pvalue"
colnames(Cig_SNP2GENE)[9] <- "Beta"
colnames(Cig_SNP2GENE)[10] <- "SE"
write.table(Cig_SNP2GENE, "Cig_SNP2GENE.txt", append = FALSE, sep = "\t",quote=FALSE, dec = ".",
            row.names = FALSE, col.names = TRUE)
                     plot_CIG_best = res_em_CIG$plots$two_stage +
  ggplot2::xlim(0, max(abs(bx) + 2*bxse)) +
  ggplot2::xlab("Genetic association with AUDIT-C") +
  ggplot2::ylab("Genetic association with FEV1") +
  ggplot2::ggtitle("")
plot_CIG_best
ggsave("plot_CIG_best.tiff", width = 28, height = 35, dpi=50)     
ggsave("plot_CIG_smaller.tiff", width = 8, height = 10, dpi=50)
#Init 
outcome_data <- extract_outcome_data(
  snps = SmoClumped$SNP, outcomes = "ebi-a-GCST007432")
H_data <- harmonise_data(
  exposure_dat = SmoClumped, 
  outcome_dat = outcome_data
)
SmoClumped_Clust<-H_data
#MR CLUST
set.seed(2000030885)
bx = SmoClumped_Clust$beta.exposure
bxse = SmoClumped_Clust$se.exposure
by = SmoClumped_Clust$beta.outcome
byse = SmoClumped_Clust$se.outcome
ratio_est = by/bx
ratio_est_se = byse/abs(bx)
snp_names = SmoClumped_Clust$SNP
res_em_Init = mr_clust_em(theta = ratio_est, theta_se = ratio_est_se, bx = bx,
                     by = by, bxse = bxse, byse = byse, obs_names = snp_names)
Init_Cluster_1<-c(res_em_Init$cluster_membership$cluster_1[1:3])
Init_Clus_1 <- Reduce(c,Init_Cluster_1)
Init_Cluster_2<-c(res_em_Init$cluster_membership$cluster_2[1:3])
Init_Clus_2 <- Reduce(c,Init_Cluster_2)
Init_Cluster_Junk<-c(res_em_Init$cluster_membership$cluster_Junk[1:3])
Init_Clus_Junk <- Reduce(c,Init_Cluster_Junk)
Init_Cluster_Null<-c(res_em_Init$cluster_membership$cluster_Null[1:3])
Init_Clus_Null <- Reduce(c,Init_Cluster_Null)

Init_table1<-as.data.frame(Init_Clus_1)
colnames(Init_table1)[1] <- "SNP"
Init_table1$Cluster<-"1"


Init_table2<-as.data.frame(Init_Clus_2)
colnames(Init_table2)[1] <- "SNP"
Init_table2$Cluster<-"2"
Init_tableJunk<-as.data.frame(Init_Clus_Junk)
colnames(Init_tableJunk)[1] <- "SNP"
Init_tableJunk$Cluster<-"Junk"
Init_tableNull<-as.data.frame(Init_Clus_Null)
colnames(Init_tableNull)[1] <- "SNP"
Init_tableNull$Cluster<-"Null"
Init_excel<-rbind(Init_table1,Init_table2)
Init_excel<-rbind(Init_excel,Init_tableNull)
Init_excel<-rbind(Init_excel,Init_tableJunk)
write.csv(Init_excel,"Init_MR_Clust_excel.csv")
Init_SNP2GENE<-merge(Init_excel,SmoClumped)  
write.csv(  Init_SNP2GENE, "Init_SNP2GENE.csv")
colnames(Init_SNP2GENE)[1] <- "SNP"
colnames(Init_SNP2GENE)[3] <- "CHR"
colnames(Init_SNP2GENE)[4] <- "pos"
colnames(Init_SNP2GENE)[5] <- "A2"
colnames(Init_SNP2GENE)[6] <- "A1"
colnames(Init_SNP2GENE)[8] <- "pvalue"
colnames(Init_SNP2GENE)[9] <- "Beta"
colnames(Init_SNP2GENE)[10] <- "SE"
write.table(Init_SNP2GENE, "Init_SNP2GENE.txt", append = FALSE,  dec = ".",
            row.names = FALSE, col.names = TRUE,sep = "\t",quote=FALSE)
plot_CIG_best = res_em_CIG$plots$two_stage +
  ggplot2::xlim(0, max(abs(bx) + 2*bxse)) +
  ggplot2::xlab("Genetic association with AUDIT-C") +
  ggplot2::ylab("Genetic association with FEV1") +
  ggplot2::ggtitle("")
plot_CIG_best
ggsave("plot_init_best.tiff", width = 28, height = 35, dpi=50)     
ggsave("plot_init_smaller.tiff", width = 8, height = 10, dpi=50)


