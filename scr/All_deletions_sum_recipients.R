#!/bin/env Rscript --no-save
###############################################################################
# Deletion sum for patients and donors

# CONTINUING WITH THE COVARIATE FILES FROM PREVIOUS script in which we created the files
# for genomic collision
###############################################################################
# Calculating the sum of deletions 

# Importing the covariate files 
R_covariates <- read_table2("results/Mm_and_deletion_analyses/R_covariates_deletions_0_1.txt")
D_covariates <- read_table2("results/Mm_and_deletion_analyses/D_covariates_deletions_0_1.txt")

# In these files, the SNPs are coded accordingly:
# 0 = homozygous for major allele or heterozygous (no deletions)
# 1 = homozygous for minor allele (deletion)

# Make a copy from the data
R_dos_pheno_sum <- R_dos_pheno_dels
D_dos_pheno_sum <- D_dos_pheno_dels

# Calculate the sum of deletions (40 SNPs) per individual
R_dos_pheno_sum$Del_sum <- c(Del_sum = rowSums(R_dos_pheno_dels[167:206]))
D_dos_pheno_sum$Del_sum <- c(Del_sum = rowSums(D_dos_pheno_dels[167:206]))

# Write out the tables
write.table(R_dos_pheno_sum, 
            file = "results/Mm_and_deletion_analyses/R_covariates_deletions_sum.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(D_dos_pheno_sum, 
            file = "results/Mm_and_deletion_analyses/D_covariates_deletions_sum.txt", 
            quote = FALSE, row.names = FALSE, col.names = TRUE)

###############################################################################
### Logistic regression analyses for deletion sum
### and acute rejection endpoint among recipients

# Logistic regression for deletion sum, unadjusted
glm_deletionsum <- glm(Rejection ~ Del_sum,
                       family = binomial(link = "logit"), data = R_dos_pheno_sum)
summary(glm_deletionsum)

# Logistic regression for deletion sum, adjusted
glm_delsum_adjusted <- glm(Rejection ~ Del_sum + R_Gender + D_Gender + R_Age + D_Age + Cold_ischemia + 
                             Eplets_total_HLAI + Eplets_total_HLAII,
                           family = binomial(link = "logit"), data = R_dos_pheno_sum)
summary(glm_delsum_adjusted)

write.table(tidy(glm_delsum_adjusted), 
            "results/Mm_and_deletion_analyses/glm_delsum_adjusted",
            sep = "\t", quote = F, row.names = F)

# OR and 95% CI
exp(cbind(OR = coef(glm_delsum_adjusted), confint(glm_delsum_adjusted)))

###############################################################################
### The survival analysis: the association of deletion sum to acute rejection 
### event

# Kaplan meier for deletion sum 
ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ Del_sum, data = R_dos_pheno_sum), 
  xlab = "Months", 
  ylab = "Overall survival probability")

# Kapplan-Meier plot with risk table for the manuscript
del_sum <- ggsurvplot(
  fit = survfit(Surv(Time_to_event_months, Rejection) ~ Del_sum, data = R_dos_pheno_sum), 
  xlab = "Time to event (months)", 
  ylab = "Probability of graft survival without rejection",
  font.x = c(face = "bold"),
  font.y = c(face = "bold"),
  risk.table = TRUE,
  risk.table.y.text = FALSE,
  legend.title = "",
  legend.labs = c("Del sum 0", "Del sum 1", "Del sum 2", "Del sum 3", "Del sum 4", "Del sum 5", 
                  "Del sum 6", "Del sum 7", "Del sum 8", "Del sum 9", "Del sum 10"),
  legend = c(0.25,0.25))

del_sum$plot <- del_sum$plot+ 
  ggplot2::annotate("text", 
                    x = 100, y = 0.25, # x and y coordinates of the text
                    label = "HR 1.02, 95% CI 0.94-1.10 \n P-value 0.6",
                    size = 5)
del_sum

# HR and CI 95% for unadjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ Del_sum, data = R_dos_pheno_sum) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# HR and CI 95% for adjusted data
coxph(Surv(Time_to_event_months, Rejection) ~ Del_sum + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia + 
        Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_sum) %>% 
  gtsummary::tbl_regression(exp = TRUE) 

# Cox proportional hazards model for adjusted data
cox_delsum_patients <- coxph(Surv(Time_to_event_months, Rejection) ~ Del_sum + R_Age + D_Age + R_Gender + D_Gender + Cold_ischemia +
                               Eplets_total_HLAI + Eplets_total_HLAII, data = R_dos_pheno_sum)
summary(cox_delsum_patients)
