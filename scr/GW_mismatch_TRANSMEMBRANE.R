#!/bin/env Rscript --no-save
###############################################################################
### Recipient and donor genome-wide matching: imputed missense variants 
# transmembrane proteins
###############################################################################

library(data.table)
library(tidyverse)
library(magrittr)
library(Amelia)
library(broom)

###############################################################################
### Performing the mismatch analyses for transmembrane proteins
### with final dataset of 1025 patients and their donors
### CONTINUING WITH THE COVARIATE FILES FROM PREVIOUS MM ANALYSIS OF SECR AND TRANSM PROTEINS

# Importing the phenotype files for both recipients and donors
R_covariates <- read_table2("results/Mm_and_deletion_analyses/R_covariates_mm_secr_transm.txt")
D_covariates <- read_table2("results/Mm_and_deletion_analyses/D_covariates_mm_secr_transm.txt")

# Read dosage files containing imputed missense SNPs for transmembrane proteins
KIDNEY_missense_transmembrane_dosage <- read_table2("results/KIDNEY_missense_transmembrane_dosage.raw")

# Join dosage files with phenotype files (both recipient and donor files)
R_dos_pheno_transmembrane <- inner_join(R_covariates, KIDNEY_missense_transmembrane_dosage, by = c("Family_ID" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)

D_dos_pheno_transmembrane <- inner_join(D_covariates, KIDNEY_missense_transmembrane_dosage, by = c("Family_ID" = "IID")) %>% 
  select(-FID, -PAT, -MAT, -SEX, -PHENOTYPE)

# Remove extra columns (leave only 'pair' and dosage-columns)
R_dosage_transmembrane <- select(R_dos_pheno_transmembrane, -Family_ID, -Individual_ID, -Rejection, -Kad_pseudo, -R_Gender, -D_Gender, -Cold_ischemia, -Crea, 
                                 -R_Age, -PRAI, -PRAII, -R_CMV, -D_CMV, -D_Age, -A_mm_high, -A_mm_low, -B_mm_high, -B_mm_low, -C_mm_high, -C_mm_low, 
                                 -DPB1_mm_high, -DPB1_mm_low, -DQA1_mm_high, -DQA1_mm_low, -DQB1_mm_high, -DQB1_mm_low, -DRB1_mm_high,
                                 -DRB1_mm_low, -A_mm_high_0_1, -A_mm_low_0_1, -B_mm_high_0_1, -B_mm_low_0_1, -C_mm_high_0_1, -C_mm_low_0_1, -DPB1_mm_high_0_1, 
                                 -DPB1_mm_low_0_1, -DQA1_mm_high_0_1, -DQA1_mm_low_0_1, -DQB1_mm_high_0_1, -DQB1_mm_low_0_1, -DRB1_mm_high_0_1, -DRB1_mm_low_0_1, 
                                 -HLAI_mm_low_total, -HLAII_mm_low_total, -HLAI_mm_high_total, -HLAII_mm_high_total, -PC1, -PC2,
                                 -PC3, -PC4, -PC5, -PC6, -PC7, -PC8, -PC9, -PC10, -PC11, -PC12, -PC13, -PC14, -PC15, -PC16, -PC17, -PC18, -PC19, -PC20, -R_D,
                                 -Eplets_total_HLAI, -Eplets_AbVer_HLAI, -Eplets_total_HLAII, -Eplets_AbVer_HLAII, -Eplets_TOTAL, -Eplets_AbVer_TOTAL, 
                                 -Mismatch, -Time_to_rejection_months, -Time_to_rejection_days, -Time_to_last_fup_months, -Time_to_last_fup_days, 
                                 -Time_to_event_months, -Time_to_event_days, -Rejektiotyyppi1, -Mm_secr_transm)

D_dosage_transmembrane <- select(D_dos_pheno_transmembrane, -Family_ID, -Individual_ID, -Rejection, -Potilas_pseudo, -R_Gender, -D_Gender, -Cold_ischemia, -Crea, 
                                 -R_Age, -PRAI, -PRAII, -R_CMV, -D_CMV, -D_Age, -A_mm_high, -A_mm_low, -B_mm_high, -B_mm_low, -C_mm_high, -C_mm_low, 
                                 -DPB1_mm_high, -DPB1_mm_low, -DQA1_mm_high, -DQA1_mm_low, -DQB1_mm_high, -DQB1_mm_low, -DRB1_mm_high,
                                 -DRB1_mm_low, -A_mm_high_0_1, -A_mm_low_0_1, -B_mm_high_0_1, -B_mm_low_0_1, -C_mm_high_0_1, -C_mm_low_0_1, -DPB1_mm_high_0_1, 
                                 -DPB1_mm_low_0_1, -DQA1_mm_high_0_1, -DQA1_mm_low_0_1, -DQB1_mm_high_0_1, -DQB1_mm_low_0_1, -DRB1_mm_high_0_1, -DRB1_mm_low_0_1, 
                                 -HLAI_mm_low_total, -HLAII_mm_low_total, -HLAI_mm_high_total, -HLAII_mm_high_total, -PC1, -PC2,
                                 -PC3, -PC4, -PC5, -PC6, -PC7, -PC8, -PC9, -PC10, -PC11, -PC12, -PC13, -PC14, -PC15, -PC16, -PC17, -PC18, -PC19, -PC20, -R_D,
                                 -Eplets_total_HLAI, -Eplets_AbVer_HLAI, -Eplets_total_HLAII, -Eplets_AbVer_HLAII, -Eplets_TOTAL, -Eplets_AbVer_TOTAL, 
                                 -Mismatch, -Time_to_rejection_months, -Time_to_rejection_days, -Time_to_last_fup_months, -Time_to_last_fup_days, 
                                 -Time_to_event_months, -Time_to_event_days, -Rejektiotyyppi1, -Mm_secr_transm)

# Identify R/D pairs
R_D_pairs_transmembrane <- inner_join(R_dosage_transmembrane, D_dosage_transmembrane, by = "Pair") %>% select(Pair)

# Select paired data
R_paired_dosage_transmembrane <- inner_join(R_D_pairs_transmembrane, R_dosage_transmembrane, by = "Pair")
missmap(R_paired_dosage_transmembrane)

D_paired_dosage_transmembrane <- inner_join(R_D_pairs_transmembrane, D_dosage_transmembrane, by = "Pair")
missmap(D_paired_dosage_transmembrane)

# Calculate mismatch sum
# The function to calculate the sum is:
Mismatch <- function(R,D) {
  ifelse((D>0 & R==0) | (D<2 & R==2), T, F) 
}

Mm_transm_result <- sapply(2:ncol(R_paired_dosage_transmembrane), 
                          function(i) {Mismatch(R_paired_dosage_transmembrane[,i], D_paired_dosage_transmembrane[,i])})

Mm_transm_result_df <- data.frame(Pair=R_paired_dosage_transmembrane$Pair, Mm_transmembrane=rowSums(Mm_transm_result))


# Match phenotype files with adjusted R/D mismatches, recipients:
R_covariates_mm_transmembrane <- inner_join(R_covariates, Mm_transm_result_df, by = "Pair")

# Also, the same for donors:
D_covariates_mm_transmembrane <- inner_join(D_covariates, Mm_transm_result_df, by = "Pair")

# Writing out the covariate table including mm sum of transmembrane proteins, recipients:
write.table(R_covariates_mm_transmembrane, file = "results/Mm_and_deletion_analyses/R_covariates_mm_transmembrane.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

# Writing out the covariate table including mm sum of transmembrane proteins, donors:
write.table(D_covariates_mm_transmembrane, file = "results/Mm_and_deletion_analyses/D_covariates_mm_transmembrane.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)


###############################################################################
### Logistic regression analyses for mismatch sum of transmembrane
### proteins and acute rejection endpoint

# Rejection vs R/D mismatches: logistic regression model with all covariates
Glm_transmembrane <- glm(Rejection ~ Mm_transmembrane + R_Gender + D_Gender + R_Age + D_Age + Cold_ischemia + PRAI + PRAII + 
                           Eplets_total_HLAI + Eplets_total_HLAII,
                         family = binomial(link = "logit"), data = R_covariates_mm_transmembrane)
summary(Glm_transmembrane)
write.table(tidy(Glm_transmembrane), 
            "results/Mm_and_deletion_analyses/Glm_transmembrane",
            sep = "\t", quote = F, row.names = F)

# Odds ratio and 95% CI for adjusted logistic regression model
exp(cbind(OR = coef(Glm_transmembrane), confint(Glm_transmembrane)))

# Rejection vs R/D mismatches: logistic regression model without covariates
Glm_transmembrane_only_sum <- glm(Rejection ~ Mm_transmembrane,
                             family = binomial(link = "logit"), data = R_covariates_mm_transmembrane)
summary(Glm_transmembrane_only_sum)
write.table(tidy(Glm_transmembrane_only_sum), 
            "results/Mm_and_deletion_analyses/testi",
            sep = "\t", quote = F, row.names = F)

###############################################################################
### The survival analysis: the mismatch sum association to acute rejection event

# Cox proportional hazards model for adjusted data
cox_transmembrane <- coxph(Surv(Time_to_event_months, Rejection) ~ Mm_transmembrane + R_Age + D_Age + R_Gender + D_Gender + 
                             Cold_ischemia + PRAI + PRAII +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_transmembrane)

summary(cox_transmembrane)
write.table(tidy(cox_transmembrane), 
            "results/Mm_and_deletion_analyses/Cox_transmembrane",
            sep = "\t", quote = F, row.names = F)

# Hazard ratio and CI 95% for adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ Mm_transmembrane + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + PRAI + PRAII +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_transmembrane) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

###############################################################################
### Analysing the quartiles of the mismatch data in transmembrane proteins

# Dividing the mismatch sum into quartiles
quantile(R_covariates_mm_transmembrane$Mm_transmembrane)

# Creating a new column from quartiles (1 = the lowest quartile, 2 = the second, 3 = the third, 4 = the highest quartile)
R_covariates_mm_transmembrane <- within(R_covariates_mm_transmembrane, quartile 
                                        <- as.integer(cut(Mm_transmembrane, quantile(Mm_transmembrane, probs = 0:4/4), include.lowest = TRUE)))

# The number of pairs per each quartile
survdiff(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_transmembrane)

# Kapplan-Meier survival analysis of quartiles of mismatches in transmembrane and secretory proteins
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_transmembrane), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# HR and CI 95% for adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ quartile + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + PRAI + PRAII +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_transmembrane) %>% 
  gtsummary::tbl_regression(exp = TRUE)

# HR and CI 95% for unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_transmembrane) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# Cox proportional hazards model for adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ quartile + R_Age + D_Age + R_Gender + D_Gender + PRAI + PRAII +
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_covariates_mm_transmembrane)

###############################################################################
# Redoing the Kapplan-Meier plots for manuscript

# Kapplan-Meier plot with risk table
transm_quartiles_plot <- ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ quartile, data = R_covariates_mm_transmembrane), 
  xlab = "Time to event (months)", 
  ylab = "Probability of graft survival without rejection",
  font.x = c(face = "bold"),
  font.y = c(face = "bold"),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Q1: MM sum 1191-1291", "Q2: MM sum 1292-1333", "Q3: MM sum 1334-1373", "Q4: MM sum 1374-1671"),
  legend = c(0.20,0.25))
transm_quartiles_plot

# The same plot with additional details written in the picture (HR and 95% CI)
transm_quartiles_plot$plot <- transm_quartiles_plot$plot+ 
  ggplot2::annotate("text", 
                    x = 100, y = 0.25, # x and y coordinates of the text
                    label = "HR 1.02, 95% CI 0.9-1.16 \n P-value 0.7",
                    size = 5)
transm_quartiles_plot
