#!/bin/env Rscript --no-save
###############################################################################
# Mismatch analysis for (genomic collision) for all deletion-tagging SNPs
# between recipient and donor

# CONTINUING WITH THE COVARIATE FILES FROM PREVIOUS script for
# Create covariate file for all deletions before genomic collision
###############################################################################

library(tidyverse)
library(survival)
library(survminer)
library(janitor)

###############################################################################

# Importing the covariate files for both recipients and donors
R_covariates <- read_table2("~/Kidney_analyses/Results_new/Mm_and_deletion_analyses/R_covariates_deletions_0_1.txt")
D_covariates <- read_table2("~/Kidney_analyses/Results_new/Mm_and_deletion_analyses/D_covariates_deletions_0_1.txt")

# Remove extra columns from both covariate files (so that we only have the 'pairs' column and the deletion-columns with either 
# value 1 or 0)
R_dos_pheno_dels_collision <- select(R_covariates, "Pair", 167:206)
D_dos_pheno_dels_collision <- select(D_covariates, "Pair", 167:206)

# The files for both recipients and donors are coded accordingly:
# 0 = heterozygous/homozygous for non-risk allele
# 1 = homozygous for risk-allele

# The next thing to do is to compare the patient and donor genotype
# I want mismatch when patient is 1 (homozygous for risk-allele) and donor is 0 (either 
# heterozygous or homozygous for the non-risk allele)

# Identify R/D pairs
Collision_pairs <- inner_join(R_dos_pheno_dels_collision, D_dos_pheno_dels_collision, by = "Pair") %>% select(Pair)

# Select paired data
R_paired_dosage_collision <- inner_join(Collision_pairs, R_dos_pheno_dels_collision, by = "Pair")
missmap(R_paired_dosage_collision)

D_paired_dosage_collision <- inner_join(Collision_pairs, D_dos_pheno_dels_collision, by = "Pair")
missmap(D_paired_dosage_collision)

# Calculate mismatch sum
# The function to calculate the sum is:
Mismatch <- function(R,D) {
  ifelse((D==0 & R==1), 1, 0) 
}

Collision <- sapply(2:ncol(R_paired_dosage_collision), 
                           function(i) {Mismatch(R_paired_dosage_collision[,i], D_paired_dosage_collision[,i])})

Collision_df <- data.frame(Pair=R_paired_dosage_collision$Pair, Collision)

# Rename the columns in genomic collision-dataframe
Collision_df <- rename(Collision_df, rs6943474_col = X1)
Collision_df <- rename(Collision_df, rs4882017_col = X2)
Collision_df <- rename(Collision_df, rs10927864_col = X3)
Collision_df <- rename(Collision_df, rs11209948_col = X4)
Collision_df <- rename(Collision_df, rs11249248_col = X5)
Collision_df <- rename(Collision_df, rs11587012_col = X6)
Collision_df <- rename(Collision_df, rs158736_col = X7)
Collision_df <- rename(Collision_df, rs6693105_col = X8)
Collision_df <- rename(Collision_df, rs7542235_col = X9)
Collision_df <- rename(Collision_df, rs7419565_col = X10)
Collision_df <- rename(Collision_df, rs893403_col = X11)
Collision_df <- rename(Collision_df, rs10053292_col = X12)
Collision_df <- rename(Collision_df, rs2387715_col = X13)
Collision_df <- rename(Collision_df, rs7703761_col = X14)
Collision_df <- rename(Collision_df, rs17654108_col = X15)
Collision_df <- rename(Collision_df, rs2160195_col = X16)
Collision_df <- rename(Collision_df, rs4621754_col = X17)
Collision_df <- rename(Collision_df, rs4729606_col = X18)
Collision_df <- rename(Collision_df, rs11985201_col = X19)
Collision_df <- rename(Collision_df, rs4543566_col = X20)
Collision_df <- rename(Collision_df, rs1523688_col = X21)
Collision_df <- rename(Collision_df, rs2174926_col = X22)
Collision_df <- rename(Collision_df, rs10885336_col = X23)
Collision_df <- rename(Collision_df, rs2342606_col = X24)
Collision_df <- rename(Collision_df, rs3793917_col = X25)
Collision_df <- rename(Collision_df, rs11228868_col = X26)
Collision_df <- rename(Collision_df, rs1944862_col = X27)
Collision_df <- rename(Collision_df, rs1478309_col = X28)
Collision_df <- rename(Collision_df, rs9318648_col = X29)
Collision_df <- rename(Collision_df, rs11156875_col = X30)
Collision_df <- rename(Collision_df, rs8007442_col = X31)
Collision_df <- rename(Collision_df, rs8022070_col = X32)
Collision_df <- rename(Collision_df, rs10521145_col = X33)
Collision_df <- rename(Collision_df, rs2244613_col = X34)
Collision_df <- rename(Collision_df, rs16966699_col = X35)
Collision_df <- rename(Collision_df, rs8064493_col = X36)
Collision_df <- rename(Collision_df, rs103294_col = X37)
Collision_df <- rename(Collision_df, rs324121_col = X38)
Collision_df <- rename(Collision_df, rs3810336_col = X39)
Collision_df <- rename(Collision_df, rs4806152_col = X40)

# Match covariate files with genomic collision result (for both recipient and donor)
R_dos_pheno_dels_collision <- inner_join(R_covariates, Collision_df, by = "Pair")
D_dos_pheno_dels_collision <- inner_join(D_covariates, Collision_df, by = "Pair")

# Write out the tables
write.table(R_dos_pheno_dels_collision, 
            file = "results/Mm_and_deletion_analyses/R_covariates_deletions_collision.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(D_dos_pheno_dels_collision, 
            file = "results/Mm_and_deletion_analyses/D_covariates_deletions_collision.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

###############################################################################
### Analysing the association of mismatch vs nonmismatch to
### acute rejection with Kaplan-meier plot

# rs6943474
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs6943474_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability",
  risk.table = TRUE)

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs6943474_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs6943474_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs6943474_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs4882017
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs4882017_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs4882017_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4882017_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4882017_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs10927864
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs10927864_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs10927864_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10927864_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10927864_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs11209948
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11209948_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs11209948_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11209948_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11209948_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs11249248
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11249248_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs11249248_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11249248_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11249248_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs11587012
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11587012_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs11587012_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11587012_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11587012_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs158736
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs158736_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs158736_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs158736_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs158736_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs6693105
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs6693105_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs6693105_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs6693105_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs6693105_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# CFH
# rs7542235
rs7542235_del_col <- ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs7542235_col, data = R_dos_pheno_dels_collision), 
  xlab = "Time to event (months)", 
  ylab = "Probability of graft survival without rejection",
  font.x = c(face = "bold"),
  font.y = c(face = "bold"),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("no mismatch", "mismatch"),
  legend = c(0.20,0.25))

rs7542235_del_col$plot <- rs7542235_del_col$plot+ 
  ggplot2::annotate("text", 
                    x = 100, y = 0.25, # x and y coordinates of the text
                    label = "HR 3.10, 95% CI 1.53-6.29 \n P-value 0.002",
                    size = 5)
rs7542235_del_col

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs7542235_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7542235_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7542235_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs7419565
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs7419565_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs7419565_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7419565_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7419565_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs893403
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs893403_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs893403_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs893403_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs893403_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs10053292
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs10053292_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs10053292_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10053292_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10053292_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs2387715
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs2387715_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs2387715_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2387715_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2387715_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs7703761
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs7703761_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs7703761_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7703761_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs7703761_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs17654108
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs17654108_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs17654108_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs17654108_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs17654108_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs2160195
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs2160195_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs2160195_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2160195_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2160195_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs4621754
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs4621754_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs4621754_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4621754_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4621754_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs4729606
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs4729606_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs4729606_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4729606_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4729606_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs11985201
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11985201_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs11985201_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11985201_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11985201_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs4543566
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs4543566_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs4543566_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4543566_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4543566_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs1523688
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs1523688_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs1523688_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1523688_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1523688_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs2174926
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs2174926_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs2174926_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2174926_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2174926_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs10885336
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs10885336_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs10885336_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10885336_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10885336_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs2342606
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs2342606_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs2342606_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2342606_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2342606_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs3793917
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs3793917_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs3793917_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs3793917_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs3793917_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs11228868
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11228868_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs11228868_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11228868_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11228868_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
# rs1944862
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs1944862_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs1944862_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1944862_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1944862_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs1478309
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs1478309_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs1478309_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1478309_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs1478309_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs9318648
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs9318648_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs9318648_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs9318648_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs9318648_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)


###############################################################################
# rs11156875
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs11156875_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs11156875_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11156875_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs11156875_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs8007442
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs8007442_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs8007442_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8007442_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8007442_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs8022070
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs8022070_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs8022070_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8022070_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8022070_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs10521145
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs10521145_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs10521145_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10521145_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs10521145_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs2244613
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs2244613_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs2244613_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2244613_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs2244613_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs16966699
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs16966699_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs16966699_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs16966699_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs16966699_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs8064493
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs8064493_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs8064493_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8064493_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs8064493_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs103294
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs103294_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs103294_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs103294_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs103294_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs324121
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs324121_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs324121_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs324121_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs324121_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs3810336
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs3810336_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs3810336_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs3810336_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs3810336_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

###############################################################################
# rs4806152
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ rs4806152_col, data = R_dos_pheno_dels_collision), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Number of mismatches in rejection and no rejection group
tabyl(R_dos_pheno_dels_collision, Rejection, rs4806152_col)

# Hazard ratio, 95% CI and P-value with unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4806152_col, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# Hazard ratio, 95% CI and P-value with adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ rs4806152_col + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_dels_collision) %>% 
  gtsummary::tbl_regression(exp = TRUE)
