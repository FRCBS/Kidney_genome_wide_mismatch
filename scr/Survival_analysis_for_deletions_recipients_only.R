###############################################################################
# Surival analysis for all deletion-tagging SNP:s
# for RECIPIENT-ONLY data of 1025 individuals
# 210427
###############################################################################

library(tidyverse)
library(survival)
library(survminer)
library(janitor)

###############################################################################
### Creating a file for survival analysis of 40 deletion-tagging variants
# CONTINUING WITH THE COVARIATE FILE FROM PREVIOUS MM ANALYSIS OF ALL MISSENSE SNPs

# Import Covariate file for recipients
R_covariates_mm_all <- read_table2("~/Kidney_analyses/Results_new/Mm_and_deletion_analyses/R_covariates_mm_all.txt")

# Import dosage file for deletion-tagged SNPs (including 40 deletions)
KIDNEY_DELS_dosage <- read_table2("~/Kidney_analyses/Results_new/KIDNEY_DELS_dosage.raw")

# Join dosage file with covariate file
R_dos_pheno_dels <- inner_join(R_covariates_mm_all, KIDNEY_DELS_dosage, by = c("Family_ID" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)

# Recoding the SNPs accordingly:
# 1 = homozygous for major allele or heterozygous
# 2 = homozygous for minor allele

# In the dosage file, the homozygosity/heterozygosity is coded with 0, 1 or 2
# First, recoding 0 (homozygous for major allele) into value 1
# (the other values remain the same (1=1 and 2=2))
# for every SNP except for
# chr7_101357735_A_G_G (rs6943474) and chr11_48548630_A_G_A (rs4882017)
# because they have a different minor alle than in general population
R_dos_pheno_dels$rs10927864 <- recode(R_dos_pheno_dels$chr1_15839166_C_G_G , '0' = 1)
R_dos_pheno_dels$rs11209948 <- recode(R_dos_pheno_dels$chr1_72346221_G_T_G, '0' = 1)
R_dos_pheno_dels$rs11249248 <- recode(R_dos_pheno_dels$chr1_25426560_T_C_T, '0' = 1)
R_dos_pheno_dels$rs11587012 <- recode(R_dos_pheno_dels$chr1_152786728_A_C_C, '0' = 1)
R_dos_pheno_dels$rs158736 <- recode(R_dos_pheno_dels$chr1_222193326_G_C_C, '0' = 1)
R_dos_pheno_dels$rs6693105 <- recode(R_dos_pheno_dels$chr1_152618187_T_C_T, '0' = 1)
R_dos_pheno_dels$rs7542235 <- recode(R_dos_pheno_dels$chr1_196854483_A_G_G, '0' = 1) # CFH
R_dos_pheno_dels$rs7419565 <- recode(R_dos_pheno_dels$chr2_130195449_T_C_C, '0' = 1) 
R_dos_pheno_dels$rs893403 <- recode(R_dos_pheno_dels$chr2_108606280_G_A_G, '0' = 1) # LIMS1
R_dos_pheno_dels$rs10053292 <- recode(R_dos_pheno_dels$chr5_149896561_T_C_C, '0' = 1)
R_dos_pheno_dels$rs2387715 <- recode(R_dos_pheno_dels$chr5_180934266_A_T_T, '0' = 1)
R_dos_pheno_dels$rs7703761 <- recode(R_dos_pheno_dels$chr5_148177325_C_T_C, '0' = 1)
R_dos_pheno_dels$rs17654108 <- recode(R_dos_pheno_dels$chr6_121480387_A_T_T, '0' = 1)
R_dos_pheno_dels$rs2160195 <- recode(R_dos_pheno_dels$chr7_38365370_A_T_T, '0' = 1)
R_dos_pheno_dels$rs4621754 <- recode(R_dos_pheno_dels$chr7_158706303_A_G_G, '0' = 1)
R_dos_pheno_dels$rs4729606 <- recode(R_dos_pheno_dels$chr7_100724167_T_C_C, '0' = 1)
R_dos_pheno_dels$rs11985201 <- recode(R_dos_pheno_dels$chr8_39581402_G_A_A, '0' = 1)
R_dos_pheno_dels$rs4543566 <- recode(R_dos_pheno_dels$chr8_6968014_C_G_G, '0' = 1)
R_dos_pheno_dels$rs1523688 <- recode(R_dos_pheno_dels$chr9_104592374_T_G_G, '0' = 1)
R_dos_pheno_dels$rs2174926 <- recode(R_dos_pheno_dels$chr9_118763272_A_G_A, '0' = 1)
R_dos_pheno_dels$rs10885336 <- recode(R_dos_pheno_dels$chr10_112351444_G_A_A, '0' = 1)
R_dos_pheno_dels$rs2342606 <- recode(R_dos_pheno_dels$chr10_80034407_T_C_T, '0' = 1)
R_dos_pheno_dels$rs3793917 <- recode(R_dos_pheno_dels$chr10_122459759_C_G_G, '0' = 1)
R_dos_pheno_dels$rs11228868 <- recode(R_dos_pheno_dels$chr11_55252994_C_T_T, '0' = 1)
R_dos_pheno_dels$rs1944862 <- recode(R_dos_pheno_dels$chr11_55521460_G_A_A, '0' = 1)
R_dos_pheno_dels$rs1478309 <- recode(R_dos_pheno_dels$chr12_10389146_T_G_T, '0' = 1)
R_dos_pheno_dels$rs9318648 <- recode(R_dos_pheno_dels$chr13_24566498_A_G_A, '0' = 1)
R_dos_pheno_dels$rs11156875 <- recode(R_dos_pheno_dels$chr14_35150610_A_G_G, '0' = 1)
R_dos_pheno_dels$rs8007442 <- recode(R_dos_pheno_dels$chr14_21924807_T_C_T, '0' = 1)
R_dos_pheno_dels$rs8022070 <- recode(R_dos_pheno_dels$chr14_81411038_C_T_T, '0' = 1)
R_dos_pheno_dels$rs10521145 <- recode(R_dos_pheno_dels$chr16_28585563_G_A_A, '0' = 1)
R_dos_pheno_dels$rs2244613 <- recode(R_dos_pheno_dels$chr16_55810697_G_T_G, '0' = 1)
R_dos_pheno_dels$rs16966699 <- recode(R_dos_pheno_dels$chr17_41343914_C_G_G, '0' = 1)
R_dos_pheno_dels$rs8064493 <- recode(R_dos_pheno_dels$chr17_41237071_A_G_A, '0' = 1)
R_dos_pheno_dels$rs103294 <- recode(R_dos_pheno_dels$chr19_54293995_T_C_T, '0' = 1)
R_dos_pheno_dels$rs324121 <- recode(R_dos_pheno_dels$chr19_52391192_G_A_A, '0' = 1)
R_dos_pheno_dels$rs3810336 <- recode(R_dos_pheno_dels$chr19_56175525_G_A_A, '0' = 1)
R_dos_pheno_dels$rs4806152 <- recode(R_dos_pheno_dels$chr19_35395758_A_C_C, '0' = 1)

# THEN for the other two:
# In LIMS1 paper, the minor allele of chr7_101357735_A_G_G (rs6943474) is A
# and in our data the minor allele is G, which means that for this SNP the dosage 
# file is interpreted according to:
# 0 = homozygous for deletion-tag allele (AA)
# 1 = heterozygous (AG)
# 2 = homozygous for non-deletion tag allele (GG)
# And also, in LIMS1 paper the minor allele of chr11_48548630_A_G_A (rs4882017) is G
# and in our data the minor allele is A, which means that for this SNP the dosage 
# file is interpreted according to:
# 0 = homozygous for deletion-tag allele (GG)
# 1 = heterozygous (AG)
# 2 = homozygous for non-deletion tag allele (AA)
# Based on this knowledge, recode the SNPs accordingly:
# 1 = homozygous for major allele or heterozygous
# 2 = homozygous for minor allele

# So, recoding 2 (homozygous for non-deletion allele) into value 1
# and 0 (homozygous for deletion-tag allele) into value 2
# (the heterozygous value 1 remains the same 1)
# for SNPs
# chr7_101357735_A_G_G (rs6943474) and chr11_48548630_A_G_A (rs4882017)
R_dos_pheno_dels$rs6943474 <- recode(R_dos_pheno_dels$chr7_101357735_A_G_G, '2' = 1)
R_dos_pheno_dels$rs6943474 <- recode(R_dos_pheno_dels$rs6943474, '0' = 2)

R_dos_pheno_dels$rs4882017 <- recode(R_dos_pheno_dels$chr11_48548630_A_G_A, '2' = 1)
R_dos_pheno_dels$rs4882017 <- recode(R_dos_pheno_dels$rs4882017, '0' = 2)

# Write out the table
write.table(R_dos_pheno_dels, file = "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/R_covariates_deletions.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

###############################################################################
### Analysing the association of risk genotype vs nonrisk with genotype to
### acute rejection with Kaplan-meier plot

## rs10927864
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs10927864, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs10927864)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10927864, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10927864 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs11209948
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11209948, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs11209948)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11209948, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11209948 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
## rs11249248
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11249248, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs11249248)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11249248, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11249248 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
## rs11587012
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11587012, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs11587012)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11587012, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11587012 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
## rs158736
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs158736, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs158736)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs158736, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs158736 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
## rs6693105
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs6693105, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs6693105)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs6693105, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs6693105 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
## CFH
## rs7542235
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs7542235, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability",
  risk.table = TRUE) 

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs7542235)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7542235, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7542235 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs7419565
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs7419565, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs7419565)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7419565, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7419565 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# LIMS1 (Steers et al.)
## rs893403
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs893403, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs893403)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs893403, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs893403 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs10053292
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs10053292, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs10053292)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10053292, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10053292 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs2387715
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs2387715, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs2387715)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2387715, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2387715 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs7703761
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs7703761, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs7703761)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7703761, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7703761 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs17654108
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs17654108, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs17654108)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs17654108, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs17654108 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs2160195
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs2160195, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs2160195)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2160195, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2160195 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs4621754
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs4621754, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs4621754)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4621754, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4621754 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs4729606
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs4729606, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs4729606)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4729606, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4729606 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs11985201
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11985201, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs11985201)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11985201, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11985201 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs4543566
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs4543566, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs4543566)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4543566, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4543566 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 


###############################################################################
## rs1523688
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs1523688, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs1523688)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1523688, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1523688 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 


###############################################################################
## rs2174926
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs2174926, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs2174926)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2174926, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2174926 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs10885336
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs10885336, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs10885336)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10885336, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10885336 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs2342606
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs2342606, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs2342606)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2342606, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2342606 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs3793917
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs3793917, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs3793917)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs3793917, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs3793917 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs11228868
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11228868, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs11228868)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11228868, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11228868 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs1944862 
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs1944862, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs1944862)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1944862, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1944862 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs1478309 
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs1478309, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs1478309)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1478309, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1478309 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs9318648 
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs9318648, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs9318648)
# No rejection: non-risk 760, risk 66

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs9318648, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs9318648 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs11156875
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11156875, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs11156875)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11156875, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11156875 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs8007442
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs8007442, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs8007442)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8007442, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8007442 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs8022070 
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs8022070, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs8022070)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8022070, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8022070 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs10521145
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs10521145, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs10521145)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10521145, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10521145 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs2244613
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs2244613, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs2244613)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2244613, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2244613 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs16966699
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs16966699, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs16966699)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs16966699, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs16966699 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs8064493
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs8064493, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs8064493)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8064493, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8064493 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs103294
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs103294, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs103294)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs103294, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs103294 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs324121
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs324121, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs324121)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs324121, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs324121 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs3810336
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs3810336, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs3810336)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs3810336, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs3810336 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs4806152
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs4806152, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs4806152)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4806152, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
testi <- coxph(Surv(Time_to_event_months, Rejection) ~ rs4806152 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs6943474
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs6943474, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs6943474)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs6943474, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs6943474 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
## rs4882017
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs4882017, data = R_dos_pheno_dels), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Amount of rejections among non-risk and risk genotype
tabyl(R_dos_pheno_dels, Rejection, rs4882017)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4882017, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4882017 + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# Creating the exactly the same covariate file for donors also for recipient-donor
# analysis

# Import Covariate file for donors
D_covariates_mm_all <- read_table2("~/Kidney_analyses/Results_new/Mm_and_deletion_analyses/D_covariates_mm_all.txt")

# Import dosage file for deletion-tagged SNPs (including 40 deletions)
KIDNEY_DELS_dosage <- read_table2("~/Kidney_analyses/Results_new/KIDNEY_DELS_dosage.raw")

# Join dosage file with covariate file
D_dos_pheno_dels <- inner_join(D_covariates_mm_all, KIDNEY_DELS_dosage, by = c("Family_ID" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)

# Recoding the SNPs accordingly:
# 1 = homozygous for major allele or heterozygous
# 2 = homozygous for minor allele

# Recode 0 (homozygous for major allele) into value 1
# (the other values remain the same (1=1 and 2=2))
# for every SNP except for
# chr7_101357735_A_G_G (rs6943474) and chr11_48548630_A_G_A (rs4882017)
D_dos_pheno_dels$rs10927864 <- recode(D_dos_pheno_dels$chr1_15839166_C_G_G , '0' = 1)
D_dos_pheno_dels$rs11209948 <- recode(D_dos_pheno_dels$chr1_72346221_G_T_G, '0' = 1)
D_dos_pheno_dels$rs11249248 <- recode(D_dos_pheno_dels$chr1_25426560_T_C_T, '0' = 1)
D_dos_pheno_dels$rs11587012 <- recode(D_dos_pheno_dels$chr1_152786728_A_C_C, '0' = 1)
D_dos_pheno_dels$rs158736 <- recode(D_dos_pheno_dels$chr1_222193326_G_C_C, '0' = 1)
D_dos_pheno_dels$rs6693105 <- recode(D_dos_pheno_dels$chr1_152618187_T_C_T, '0' = 1)
D_dos_pheno_dels$rs7542235 <- recode(D_dos_pheno_dels$chr1_196854483_A_G_G, '0' = 1) # CFH
D_dos_pheno_dels$rs7419565 <- recode(D_dos_pheno_dels$chr2_130195449_T_C_C, '0' = 1) 
D_dos_pheno_dels$rs893403 <- recode(D_dos_pheno_dels$chr2_108606280_G_A_G, '0' = 1) # LIMS1
D_dos_pheno_dels$rs10053292 <- recode(D_dos_pheno_dels$chr5_149896561_T_C_C, '0' = 1)
D_dos_pheno_dels$rs2387715 <- recode(D_dos_pheno_dels$chr5_180934266_A_T_T, '0' = 1)
D_dos_pheno_dels$rs7703761 <- recode(D_dos_pheno_dels$chr5_148177325_C_T_C, '0' = 1)
D_dos_pheno_dels$rs17654108 <- recode(D_dos_pheno_dels$chr6_121480387_A_T_T, '0' = 1)
D_dos_pheno_dels$rs2160195 <- recode(D_dos_pheno_dels$chr7_38365370_A_T_T, '0' = 1)
D_dos_pheno_dels$rs4621754 <- recode(D_dos_pheno_dels$chr7_158706303_A_G_G, '0' = 1)
D_dos_pheno_dels$rs4729606 <- recode(D_dos_pheno_dels$chr7_100724167_T_C_C, '0' = 1)
D_dos_pheno_dels$rs11985201 <- recode(D_dos_pheno_dels$chr8_39581402_G_A_A, '0' = 1)
D_dos_pheno_dels$rs4543566 <- recode(D_dos_pheno_dels$chr8_6968014_C_G_G, '0' = 1)
D_dos_pheno_dels$rs1523688 <- recode(D_dos_pheno_dels$chr9_104592374_T_G_G, '0' = 1)
D_dos_pheno_dels$rs2174926 <- recode(D_dos_pheno_dels$chr9_118763272_A_G_A, '0' = 1)
D_dos_pheno_dels$rs10885336 <- recode(D_dos_pheno_dels$chr10_112351444_G_A_A, '0' = 1)
D_dos_pheno_dels$rs2342606 <- recode(D_dos_pheno_dels$chr10_80034407_T_C_T, '0' = 1)
D_dos_pheno_dels$rs3793917 <- recode(D_dos_pheno_dels$chr10_122459759_C_G_G, '0' = 1)
D_dos_pheno_dels$rs11228868 <- recode(D_dos_pheno_dels$chr11_55252994_C_T_T, '0' = 1)
D_dos_pheno_dels$rs1944862 <- recode(D_dos_pheno_dels$chr11_55521460_G_A_A, '0' = 1)
D_dos_pheno_dels$rs1478309 <- recode(D_dos_pheno_dels$chr12_10389146_T_G_T, '0' = 1)
D_dos_pheno_dels$rs9318648 <- recode(D_dos_pheno_dels$chr13_24566498_A_G_A, '0' = 1)
D_dos_pheno_dels$rs11156875 <- recode(D_dos_pheno_dels$chr14_35150610_A_G_G, '0' = 1)
D_dos_pheno_dels$rs8007442 <- recode(D_dos_pheno_dels$chr14_21924807_T_C_T, '0' = 1)
D_dos_pheno_dels$rs8022070 <- recode(D_dos_pheno_dels$chr14_81411038_C_T_T, '0' = 1)
D_dos_pheno_dels$rs10521145 <- recode(D_dos_pheno_dels$chr16_28585563_G_A_A, '0' = 1)
D_dos_pheno_dels$rs2244613 <- recode(D_dos_pheno_dels$chr16_55810697_G_T_G, '0' = 1)
D_dos_pheno_dels$rs16966699 <- recode(D_dos_pheno_dels$chr17_41343914_C_G_G, '0' = 1)
D_dos_pheno_dels$rs8064493 <- recode(D_dos_pheno_dels$chr17_41237071_A_G_A, '0' = 1)
D_dos_pheno_dels$rs103294 <- recode(D_dos_pheno_dels$chr19_54293995_T_C_T, '0' = 1)
D_dos_pheno_dels$rs324121 <- recode(D_dos_pheno_dels$chr19_52391192_G_A_A, '0' = 1)
D_dos_pheno_dels$rs3810336 <- recode(D_dos_pheno_dels$chr19_56175525_G_A_A, '0' = 1)
D_dos_pheno_dels$rs4806152 <- recode(D_dos_pheno_dels$chr19_35395758_A_C_C, '0' = 1)

# THEN, the same for donors also
# In LIMS1 paper, the minor allele of chr7_101357735_A_G_G (rs6943474) is A
# and in our data the minor allele is G, which means that for this SNP the dosage 
# file is interpreted according to:
# 0 = homozygous for deletion-tag allele (AA)
# 1 = heterozygous (AG)
# 2 = homozygous for non-deletion tag allele (GG)
# And also, in LIMS1 paper the minor allele of chr11_48548630_A_G_A (rs4882017) is G
# and in our data the minor allele is A, which means that for this SNP the dosage 
# file is interpreted according to:
# 0 = homozygous for deletion-tag allele (GG)
# 1 = heterozygous (AG)
# 2 = homozygous for non-deletion tag allele (AA)
# Based on this knowledge, recode the SNPs accordingly:
# 1 = homozygous for major allele or heterozygous
# 2 = homozygous for minor allele

# Recode 2 (homozygous for non-deletion allele) into value 1
# and 0 (homozygous for deletion-tag allele) into value 2
# (the heterozygous value 1 remains the same 1)
# for SNPs
# chr7_101357735_A_G_G (rs6943474) and chr11_48548630_A_G_A (rs4882017)
D_dos_pheno_dels$rs6943474 <- recode(D_dos_pheno_dels$chr7_101357735_A_G_G, '2' = 1)
D_dos_pheno_dels$rs6943474 <- recode(D_dos_pheno_dels$rs6943474, '0' = 2)

D_dos_pheno_dels$rs4882017 <- recode(D_dos_pheno_dels$chr11_48548630_A_G_A, '2' = 1)
D_dos_pheno_dels$rs4882017 <- recode(D_dos_pheno_dels$rs4882017, '0' = 2)

# Write out the table
write.table(D_dos_pheno_dels, file = "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/D_covariates_deletions.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)
