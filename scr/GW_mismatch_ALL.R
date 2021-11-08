#!/bin/env Rscript --no-save
###############################################################################
### Recipient and donor genome-wide matching: imputed missense variants
# all missense variants
###############################################################################

library(data.table)
library(tidyverse)
library(magrittr)
library(Amelia)
library(broom)

###############################################################################
### Performing the mismatch analyses for all genome-wide missense variants
### with final dataset of 1025 patients and their donors
# CONTINUING WITH THE COVARIATE FILES FROM PREVIOUS MM ANALYSIS OF KIDNEY-RELATED PROTEINS

# Importing the phenotype files for both recipients and donors
R_covariates <- read_table2("results/Mm_and_deletion_analyses/R_covariates_mm_kidney.txt")
D_covariates <- read_table2("results/Mm_and_deletion_analyses/D_covariates_mm_kidney.txt")

# Read dosage files containing imputed missense SNPs only
KIDNEY_missense_dosage <- read_table2("results/KIDNEY_missense_dosage_without_X_MHC.raw")


# Join dosage files with phenotype files (both recipient and donor files)
R_dos_pheno_all <- inner_join(R_covariates, KIDNEY_missense_dosage, by = c("Family_ID" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)

D_dos_pheno_all <- inner_join(D_covariates, KIDNEY_missense_dosage, by = c("Family_ID" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)

# Remove extra columns (leave only 'pair' and dosage-columns)
R_dosage_all <- select(R_dos_pheno_all, -Family_ID, -Individual_ID, -Rejection, -Kad_pseudo, -R_Gender, -D_Gender, -Cold_ischemia, -Crea, 
                       -R_Age, -PRAI, -PRAII, -R_CMV, -D_CMV, -D_Age, -A_mm_high, -A_mm_low, -B_mm_high, -B_mm_low, -C_mm_high, -C_mm_low, 
                       -DPB1_mm_high, -DPB1_mm_low, -DQA1_mm_high, -DQA1_mm_low, -DQB1_mm_high, -DQB1_mm_low, -DRB1_mm_high,
                       -DRB1_mm_low, -A_mm_high_0_1, -A_mm_low_0_1, -B_mm_high_0_1, -B_mm_low_0_1, -C_mm_high_0_1, -C_mm_low_0_1, -DPB1_mm_high_0_1, 
                       -DPB1_mm_low_0_1, -DQA1_mm_high_0_1, -DQA1_mm_low_0_1, -DQB1_mm_high_0_1, -DQB1_mm_low_0_1, -DRB1_mm_high_0_1, -DRB1_mm_low_0_1, 
                       -HLAI_mm_low_total, -HLAII_mm_low_total, -HLAI_mm_high_total, -HLAII_mm_high_total, -PC1, -PC2,
                       -PC3, -PC4, -PC5, -PC6, -PC7, -PC8, -PC9, -PC10, -PC11, -PC12, -PC13, -PC14, -PC15, -PC16, -PC17, -PC18, -PC19, -PC20, -R_D,
                       -Eplets_total_HLAI, -Eplets_AbVer_HLAI, -Eplets_total_HLAII, -Eplets_AbVer_HLAII, -Eplets_TOTAL, -Eplets_AbVer_TOTAL, 
                       -Mismatch, -Time_to_rejection_months, -Time_to_rejection_days, -Time_to_last_fup_months, -Time_to_last_fup_days, 
                       -Time_to_event_months, -Time_to_event_days, -Rejektiotyyppi1, -Mm_secr_transm, -Mm_transmembrane, -Mm_kidney)

D_dosage_all <- select(D_dos_pheno_all, -Family_ID, -Individual_ID, -Rejection, -Potilas_pseudo, -R_Gender, -D_Gender, -Cold_ischemia, -Crea, 
                       -R_Age, -PRAI, -PRAII, -R_CMV, -D_CMV, -D_Age, -A_mm_high, -A_mm_low, -B_mm_high, -B_mm_low, -C_mm_high, -C_mm_low, 
                       -DPB1_mm_high, -DPB1_mm_low, -DQA1_mm_high, -DQA1_mm_low, -DQB1_mm_high, -DQB1_mm_low, -DRB1_mm_high,
                       -DRB1_mm_low, -A_mm_high_0_1, -A_mm_low_0_1, -B_mm_high_0_1, -B_mm_low_0_1, -C_mm_high_0_1, -C_mm_low_0_1, -DPB1_mm_high_0_1, 
                       -DPB1_mm_low_0_1, -DQA1_mm_high_0_1, -DQA1_mm_low_0_1, -DQB1_mm_high_0_1, -DQB1_mm_low_0_1, -DRB1_mm_high_0_1, -DRB1_mm_low_0_1, 
                       -HLAI_mm_low_total, -HLAII_mm_low_total, -HLAI_mm_high_total, -HLAII_mm_high_total, -PC1, -PC2,
                       -PC3, -PC4, -PC5, -PC6, -PC7, -PC8, -PC9, -PC10, -PC11, -PC12, -PC13, -PC14, -PC15, -PC16, -PC17, -PC18, -PC19, -PC20, -R_D,
                       -Eplets_total_HLAI, -Eplets_AbVer_HLAI, -Eplets_total_HLAII, -Eplets_AbVer_HLAII, -Eplets_TOTAL, -Eplets_AbVer_TOTAL, 
                       -Mismatch, -Time_to_rejection_months, -Time_to_rejection_days, -Time_to_last_fup_months, -Time_to_last_fup_days, 
                       -Time_to_event_months, -Time_to_event_days, -Rejektiotyyppi1, -Mm_secr_transm, -Mm_transmembrane, -Mm_kidney)

# Identify R/D pairs
R_D_pairs_all <- inner_join(R_dosage_all, D_dosage_all, by = "Pair") %>% select(Pair)

# Select paired data
R_paired_dosage_all <- inner_join(R_D_pairs_all, R_dosage_all, by = "Pair")
missmap(R_paired_dosage_all)

D_paired_dosage_all <- inner_join(R_D_pairs_all, D_dosage_all, by = "Pair")
missmap(D_paired_dosage_all)

# Calculate mismatch sum
# The function to calculate the sum is:
Mismatch <- function(R,D) {
  ifelse((D>0 & R==0) | (D<2 & R==2), T, F) 
}

Mm_all_result <- sapply(2:ncol(R_paired_dosage_all), 
                          function(i) {Mismatch(R_paired_dosage_all[,i], D_paired_dosage_all[,i])})

Mm_all_result_df <- data.frame(Pair=R_paired_dosage_all$Pair, Mm_all=rowSums(Mm_all_result))


# Match phenotype files with adjusted R/D mismatches, recipients:
R_covariates_mm_all <- inner_join(R_covariates, Mm_all_result_df, by = "Pair")

# Also, the same for donors:
D_covariates_mm_all <- inner_join(D_covariates, Mm_all_result_df, by = "Pair")

# Writing out the covariate table including mm sum of all missense variants, recipients:
write.table(R_covariates_mm_all, file = "results/Mm_and_deletion_analyses/R_covariates_mm_all.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# Writing out the covariate table including mm sum of all missense variants, donors:
write.table(D_covariates_mm_all, file = "results/Mm_and_deletion_analyses/D_covariates_mm_all.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)


###############################################################################
### Logistic regression analyses for mismatch sum of all missense variants
### and acute rejection endpoint

# Rejection vs R/D mismatches: logistic regression model with all covariates
Glm_all <- glm(Rejection ~ Mm_all + R_Gender + D_Gender + R_Age + D_Age + Cold_ischemia + PRAI + PRAII + 
                    Eplets_total_HLAI + Eplets_total_HLAII,
                  family = binomial(link = "logit"), data = R_covariates_mm_all)
summary(Glm_all)
write.table(tidy(Glm_all), 
            "results/Mm_and_deletion_analyses/Glm_all",
            sep = "\t", quote = F, row.names = F)

# Odds ratio and 95% CI for adjusted logistic regression model
exp(cbind(OR = coef(Glm_all), confint(Glm_all)))

# Rejection vs R/D mismatches: logistic regression model without covariates
Glm_all_only_sum <- glm(Rejection ~ Mm_all,
                           family = binomial(link = "logit"), data = R_covariates_mm_all)
summary(Glm_all_only_sum)
write.table(tidy(Glm_all_only_sum), 
            "results/Mm_and_deletion_analyses/Glm_all_only_sum",
            sep = "\t", quote = F, row.names = F)

###############################################################################
### The survival analysis: the mismatch sum association to acute rejection event

# Cox proportional hazards model for adjusted data
cox_all <- coxph(Surv(Time_to_event_months, Rejection) ~ Mm_all + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + PRAI + PRAII +
                      Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_all)
summary(cox_all)
write.table(tidy(cox_all), 
            "results/Mm_and_deletion_analyses/Cox_all",
            sep = "\t", quote = F, row.names = F)

# Hazard ratio and CI 95% for adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ Mm_all + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + PRAI + PRAII +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_all) %>% 
  gtsummary::tbl_regression(exp = TRUE) 


###############################################################################
### Analysing the quartiles of the mismatch data in all missense variants

# Dividing the mismatch sum into quartiles
quantile(R_covariates_mm_all$Mm_all)

# Creating a new column from quartiles (1 = the lowest quartile, 2 = the second, 3 = the third, 4 = the highest quartile)
R_covariates_mm_all <- within(R_covariates_mm_all, quartile <- as.integer(cut(Mm_all, 
                                                                              quantile(Mm_all, probs = 0:4/4), include.lowest = TRUE)))

# The number of pairs per each quartile
survdiff(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_all)

# Kapplan-Meier survival analysis of quartiles of mismatches in all missense variants
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_all), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# HR and CI 95% for adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ quartile + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + PRAI + PRAII +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_all) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# HR and CI 95% for unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_all) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# Cox proportional hazards model for adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ quartile + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + PRAI + PRAII +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_all)

###############################################################################
# Redoing the Kapplan-Meier plots for manuscript

# Kapplan-Meier plot with risk table
all_quartiles_plot <- ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_all), 
  xlab = "Time to event (months)", 
  ylab = "Probability of graft survival without rejection",
  font.x = c(face = "bold"),
  font.y = c(face = "bold"),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Q1: MM sum 4528-4860", "Q2: MM sum 4861-4934", "Q3: MM sum 4935-5011", "Q4: MM sum 5012-5958"),
  legend = c(0.20,0.25))
all_quartiles_plot

# The same plot with additional details written in the picture (HR and 95% CI)
all_quartiles_plot$plot <- all_quartiles_plot$plot+ 
  ggplot2::annotate("text", 
                    x = 100, y = 0.25, # x and y coordinates of the text
                    label = "HR 1.08, 95% CI 0.95-1.22 \n P-value 0.2",
                    size = 5)
all_quartiles_plot
