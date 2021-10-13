#!/bin/env Rscript --no-save

###############################################################################
### Recipient and donor genome-wide matching: imputed missense variants 
# kidney-expressed proteins
# 210426
###############################################################################
# Required packages

library(data.table)
library(tidyverse)
library(magrittr)
library(Amelia)
library(broom)

###############################################################################
### Performing the mismatch analyses for kidney-expressed proteins
### with final dataset of 1025 patients and their donors
# CONTINUING WITH THE COVARIATE FILES FROM PREVIOUS MM ANALYSIS OF TRANSMEMBRANE PROTEINS

# Importing the phenotype files for both recipients and donors
R_covariates_mm_transmembrane <- read_table2("~/Kidney_analyses/Results_new/Mm_and_deletion_analyses/R_covariates_mm_transmembrane.txt")
D_covariates_mm_transmembrane <- read_table2("~/Kidney_analyses/Results_new/Mm_and_deletion_analyses/D_covariates_mm_transmembrane.txt")

# Read dosage files containing imputed missense SNPs for kidney-expressed proteins
KIDNEY_missense_KIDNEY_dosage <- read_table2("~/Kidney_analyses/Results_new/KIDNEY_missense_KIDNEY_dosage.raw")

# Join dosage files with phenotype files (both recipient and donor files)
R_dos_pheno_kidney <- inner_join(R_covariates_mm_transmembrane, KIDNEY_missense_KIDNEY_dosage, by = c("Family_ID" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)

D_dos_pheno_kidney <- inner_join(D_covariates_mm_transmembrane, KIDNEY_missense_KIDNEY_dosage, by = c("Family_ID" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)

# Remove extra columns (leave only 'pair' and dosage-columns)
R_dosage_kidney <- select(R_dos_pheno_kidney, -Family_ID, -Individual_ID, -Rejection, -Kad_pseudo, -R_Gender, -D_Gender, -Cold_ischemia, -Crea, 
                          -R_Age, -PRAI, -PRAII, -R_CMV, -D_CMV, -D_Age, -A_mm_high, -A_mm_low, -B_mm_high, -B_mm_low, -C_mm_high, -C_mm_low, 
                          -DPB1_mm_high, -DPB1_mm_low, -DQA1_mm_high, -DQA1_mm_low, -DQB1_mm_high, -DQB1_mm_low, -DRB1_mm_high,
                          -DRB1_mm_low, -A_mm_high_0_1, -A_mm_low_0_1, -B_mm_high_0_1, -B_mm_low_0_1, -C_mm_high_0_1, -C_mm_low_0_1, -DPB1_mm_high_0_1, 
                          -DPB1_mm_low_0_1, -DQA1_mm_high_0_1, -DQA1_mm_low_0_1, -DQB1_mm_high_0_1, -DQB1_mm_low_0_1, -DRB1_mm_high_0_1, -DRB1_mm_low_0_1, 
                          -HLAI_mm_low_total, -HLAII_mm_low_total, -HLAI_mm_high_total, -HLAII_mm_high_total, -PC1, -PC2,
                          -PC3, -PC4, -PC5, -PC6, -PC7, -PC8, -PC9, -PC10, -PC11, -PC12, -PC13, -PC14, -PC15, -PC16, -PC17, -PC18, -PC19, -PC20, -R_D,
                          -Eplets_total_HLAI, -Eplets_AbVer_HLAI, -Eplets_total_HLAII, -Eplets_AbVer_HLAII, -Eplets_TOTAL, -Eplets_AbVer_TOTAL, 
                          -Mismatch, -Time_to_rejection_months, -Time_to_rejection_days, -Time_to_last_fup_months, -Time_to_last_fup_days, 
                          -Time_to_event_months, -Time_to_event_days, -Rejektiotyyppi1, -Mm_secr_transm, -Mm_transmembrane)

D_dosage_kidney <- select(D_dos_pheno_kidney, -Family_ID, -Individual_ID, -Rejection, -Potilas_pseudo, -R_Gender, -D_Gender, -Cold_ischemia, -Crea, 
                          -R_Age, -PRAI, -PRAII, -R_CMV, -D_CMV, -D_Age, -A_mm_high, -A_mm_low, -B_mm_high, -B_mm_low, -C_mm_high, -C_mm_low, 
                          -DPB1_mm_high, -DPB1_mm_low, -DQA1_mm_high, -DQA1_mm_low, -DQB1_mm_high, -DQB1_mm_low, -DRB1_mm_high,
                          -DRB1_mm_low, -A_mm_high_0_1, -A_mm_low_0_1, -B_mm_high_0_1, -B_mm_low_0_1, -C_mm_high_0_1, -C_mm_low_0_1, -DPB1_mm_high_0_1, 
                          -DPB1_mm_low_0_1, -DQA1_mm_high_0_1, -DQA1_mm_low_0_1, -DQB1_mm_high_0_1, -DQB1_mm_low_0_1, -DRB1_mm_high_0_1, -DRB1_mm_low_0_1, 
                          -HLAI_mm_low_total, -HLAII_mm_low_total, -HLAI_mm_high_total, -HLAII_mm_high_total, -PC1, -PC2,
                          -PC3, -PC4, -PC5, -PC6, -PC7, -PC8, -PC9, -PC10, -PC11, -PC12, -PC13, -PC14, -PC15, -PC16, -PC17, -PC18, -PC19, -PC20, -R_D,
                          -Eplets_total_HLAI, -Eplets_AbVer_HLAI, -Eplets_total_HLAII, -Eplets_AbVer_HLAII, -Eplets_TOTAL, -Eplets_AbVer_TOTAL, 
                          -Mismatch, -Time_to_rejection_months, -Time_to_rejection_days, -Time_to_last_fup_months, -Time_to_last_fup_days, 
                          -Time_to_event_months, -Time_to_event_days, -Rejektiotyyppi1, -Mm_secr_transm, -Mm_transmembrane)

# Identify R/D pairs
R_D_pairs_kidney <- inner_join(R_dosage_kidney, D_dosage_kidney, by = "Pair") %>% select(Pair)

# Select paired data
R_paired_dosage_kidney <- inner_join(R_D_pairs_kidney, R_dosage_kidney, by = "Pair")
missmap(R_paired_dosage_kidney)

D_paired_dosage_kidney <- inner_join(R_D_pairs_kidney, D_dosage_kidney, by = "Pair")
missmap(D_paired_dosage_kidney)

# Calculate mismatch sum
# The function to calculate the sum is:
Mismatch <- function(R,D) {
  ifelse((D>0 & R==0) | (D<2 & R==2), T, F) 
}

Mm_kidney_result <- sapply(2:ncol(R_paired_dosage_kidney), 
                          function(i) {Mismatch(R_paired_dosage_kidney[,i], D_paired_dosage_kidney[,i])})

Mm_kidney_result_df <- data.frame(Pair=R_paired_dosage_kidney$Pair, Mm_kidney=rowSums(Mm_kidney_result))


# Match phenotype files with adjusted R/D mismatches, recipients:
R_covariates_mm_kidney <- inner_join(R_covariates_mm_transmembrane, Mm_kidney_result_df, by = "Pair")

# Also, the same for donors:
D_covariates_mm_kidney <- inner_join(D_covariates_mm_transmembrane, Mm_kidney_result_df, by = "Pair")

# Writing out the covariate table including mm sum of kidney-expressed proteins, recipients:
write.table(R_covariates_mm_kidney, file = "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/R_covariates_mm_kidney.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# Writing out the covariate table including mm sum of transmembrane proteins, donors:
write.table(D_covariates_mm_kidney, file = "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/D_covariates_mm_kidney.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)


###############################################################################
### Logistic regression analyses for mismatch sum of kidney-expressed
### proteins and acute rejection endpoint

# Rejection vs R/D mismatches: logistic regression model with all covariates
Glm_kidney <- glm(Rejection ~ Mm_kidney + R_Gender + D_Gender + R_Age + D_Age + Cold_ischemia + PRAI + PRAII + 
                           Eplets_total_HLAI + Eplets_total_HLAII,
                         family = binomial(link = "logit"), data = R_covariates_mm_kidney)
summary(Glm_kidney)
write.table(tidy(Glm_kidney), 
            "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/Glm_kidney",
            sep = "\t", quote = F, row.names = F)

# Odds ratio and 95% CI for adjusted logistic regression model
exp(cbind(OR = coef(Glm_kidney), confint(Glm_kidney)))

# Rejection vs R/D mismatches: logistic regression model without covariates
Glm_kidney_only_sum <- glm(Rejection ~ Mm_kidney,
                             family = binomial(link = "logit"), data = R_covariates_mm_kidney)
summary(Glm_kidney_only_sum)
write.table(tidy(Glm_kidney_only_sum), 
            "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/Glm_kidney_only_sum",
            sep = "\t", quote = F, row.names = F)

# Saving the R data
save.image("~/Kidney_analyses/Kidney_genetics_analyses/src_for_mm_and_deletion_analyses/FINAL_GW_mm_KIDNEY.RData")

###############################################################################
### The survival analysis: the mismatch sum association to acute rejection event

# Cox proportional hazards model for adjusted data
cox_kidney <- coxph(Surv(Time_to_event_months, Rejection) ~ Mm_kidney + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + PRAI + PRAII +
                              Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_kidney)
summary(cox_kidney)
write.table(tidy(cox_kidney), 
            "/home/markkinens/Kidney_analyses/Results_new/Mm_and_deletion_analyses/Cox_kidney",
            sep = "\t", quote = F, row.names = F)

# Hazard ratio and CI 95% for adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ Mm_kidney + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + PRAI + PRAII +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_kidney) %>% 
  gtsummary::tbl_regression(exp = TRUE) 


###############################################################################
### Analysing the quartiles of the mismatch data in kidney-expressed proteins

# Dividing the mismatch sum into quartiles
quantile(R_covariates_mm_kidney$Mm_kidney)

# Creating a new column from quartiles (1 = the lowest quartile, 2 = the second, 3 = the third, 4 = the highest quartile)
R_covariates_mm_kidney <- within(R_covariates_mm_kidney, quartile <- as.integer(cut(Mm_kidney, 
                                                                                    quantile(Mm_kidney, probs = 0:4/4), include.lowest = TRUE)))

# The number of pairs per each quartile
survdiff(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_kidney)

# Kapplan-Meier survival analysis of quartiles of mismatches in transmembrane and secretory proteins
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_kidney), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# HR and CI 95% for adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ quartile + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + PRAI + PRAII +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_kidney) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# HR and CI 95% for unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_kidney) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# Cox proportional hazards model for adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ quartile + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + PRAI + PRAII +
                                Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_kidney)

###############################################################################
# Redoing the Kapplan-Meier plots for manuscript

# Kapplan-Meier plot with risk table
kidney_quartiles_plot <- ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_kidney), 
  xlab = "Time to event (months)", 
  ylab = "Probability of graft survival without rejection",
  font.x = c(face = "bold"),
  font.y = c(face = "bold"),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Q1: MM sum 509-584", "Q2: MM sum 585-604", "Q3: MM sum 605-626", "Q4: MM sum 627-757"),
  legend = c(0.20,0.25))
kidney_quartiles_plot

# The same plot with additional details written in the picture (HR and 95% CI)
kidney_quartiles_plot$plot <- kidney_quartiles_plot$plot+ 
  ggplot2::annotate("text", 
                    x = 100, y = 0.25, # x and y coordinates of the text
                    label = "HR 1.15, 95% CI 1.01-1.30 \n P-value 0.029",
                    size = 5)
kidney_quartiles_plot
